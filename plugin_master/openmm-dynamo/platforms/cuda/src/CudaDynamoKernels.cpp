/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors: V. Ovchinnikov                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "CudaDynamoKernels.h"
#include "CudaDynamoKernelSources.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include <cstring>
#include <map>

using namespace DynamoPlugin;
using namespace OpenMM;
using namespace std;

class CudaCalcDynamoForceKernel::StartCalculationPreComputation : public CudaContext::ForcePreComputation {
public:
    StartCalculationPreComputation(CudaCalcDynamoForceKernel& owner) : owner(owner) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        owner.beginComputation(includeForces, includeEnergy, groups);
    }
    CudaCalcDynamoForceKernel& owner;
};

class CudaCalcDynamoForceKernel::ExecuteTask : public CudaContext::WorkTask {
public:
    ExecuteTask(CudaCalcDynamoForceKernel& owner) : owner(owner) {
    }
    void execute() {
        owner.executeOnWorkerThread();
    }
    CudaCalcDynamoForceKernel& owner;
};

class CudaCalcDynamoForceKernel::AddForcesPostComputation : public CudaContext::ForcePostComputation {
public:
    AddForcesPostComputation(CudaCalcDynamoForceKernel& owner) : owner(owner) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        return owner.addForces(includeForces, includeEnergy, groups);
    }
    CudaCalcDynamoForceKernel& owner;
};

CudaCalcDynamoForceKernel::~CudaCalcDynamoForceKernel() {
    if (hasInitialized) {
        master_done_plugin();
        free(r);
        free(fr);
    }
    cu.setAsCurrent();
    if (dynamoForces != NULL)
        delete dynamoForces;
    cuStreamDestroy(stream);
    cuEventDestroy(syncEvent);
}

void CudaCalcDynamoForceKernel::initialize(const System& system, const DynamoForce& force) {
    natoms = system.getNumParticles();
    //
    cu.setAsCurrent();
    cuStreamCreate(&stream, CU_STREAM_NON_BLOCKING);
    cuEventCreate(&syncEvent, CU_EVENT_DISABLE_TIMING);
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    dynamoForces = new CudaArray(cu, 3*natoms, elementSize, "dynamoForces");
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaDynamoKernelSources::dynamoForce, defines);
    addForcesKernel = cu.getKernel(module, "addForces");
    forceGroupFlag = (1<<force.getForceGroup());
    cu.addPreComputation(new StartCalculationPreComputation(*this));
    cu.addPostComputation(new AddForcesPostComputation(*this));

    //script name
    string inputfile__=force.getScript();
    int ilen=inputfile__.length();
    std::vector<char> inputfile_(ilen+1);
    std::strcpy(&inputfile_[0], inputfile__.c_str());
    char* inputfile = &inputfile_[0];
    //log name
    string logfile__ = force.getLog();
    int llen=logfile__.length();
    std::vector<char> logfile_(llen+1);
    std::strcpy(&logfile_[0], logfile__.c_str());
    char* logfile = &logfile_[0];
    // allocate position and force arrays
    r=(double*) calloc(3 * natoms, sizeof(double));
    fr=(double*) calloc(3 * natoms, sizeof(double));
    //PBC
    usesPeriodic = system.usesPeriodicBoundaryConditions();
    if (usesPeriodic) {
     system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
     box[0]=boxVectors[0][0]*nm2A;
     box[1]=boxVectors[0][1]*nm2A;
     box[2]=boxVectors[0][2]*nm2A;
     box[3]=boxVectors[1][0]*nm2A;
     box[4]=boxVectors[1][1]*nm2A;
     box[5]=boxVectors[1][2]*nm2A;
     box[6]=boxVectors[2][0]*nm2A;
     box[7]=boxVectors[2][1]*nm2A;
     box[8]=boxVectors[2][2]*nm2A;
    } else {
     for ( int i=0 ; i < 9 ; i++ ) { box[i]=0.0 ; } // initialize "by hand" for compatibility with older compilers
    }
    // Get particle masses and charges (if available)
    double *m=NULL; //mass
    double *q=NULL; //charge
    m = (double*) malloc(natoms * sizeof(double)); // allocate memory
    for (int i = 0; i < natoms; i++)
        m[i] = system.getParticleMass(i);

    q = (double*) calloc(natoms, sizeof(double));
    // If there's a NonbondedForce, get charges from it (otherwise, they will remain zero)

    for (int j = 0; j < system.getNumForces(); j++) {
        const NonbondedForce* nonbonded = dynamic_cast<const NonbondedForce*>(&system.getForce(j));
        if (nonbonded != NULL) {
            double sigma, epsilon;
            for (int i = 0; i < natoms; i++)
                nonbonded->getParticleParameters(i, q[i], sigma, epsilon);
        }
    }
    // initialize dynamo
    int ierr=master_init_plugin(natoms, m, q, inputfile, ilen, logfile, llen, &atomlist, usesPeriodic, box);
    free(m);
    free(q);
    pos.resize(natoms);
    hasInitialized = true;
    if (ierr) throw OpenMMException("Could not initialize DYNAMO plugin");
} // initialize

double CudaCalcDynamoForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    // This method does nothing.  The actual calculation is started by the pre-computation, continued on
    // the worker thread, and finished by the post-computation.
    return 0;
}

void CudaCalcDynamoForceKernel::beginComputation(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return;
    contextImpl.getPositions(pos);
    // The actual force computation will be done on a different thread.
    cu.getWorkThread().addTask(new ExecuteTask(*this));
}

void CudaCalcDynamoForceKernel::executeOnWorkerThread() {
//NOTE : 1/2019 : changes interface to match that of OPENCL (serialized) ; not clear that the removal of extra parallel copy 
// is always beneficial
    long int iteration = cu.getStepCount();
    double* rptr; // pointer to coordinate array
    int* aptr; // pointer to atom index array
    int i, j, ierr;
    // buffer for uploading forces to the device:
    bool qdble=cu.getUseDoublePrecision();
    // update periodic vectors
    if (usesPeriodic) {
     contextImpl.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
     box[0]=boxVectors[0][0]*nm2A;
     box[1]=boxVectors[0][1]*nm2A;
     box[2]=boxVectors[0][2]*nm2A;
     box[3]=boxVectors[1][0]*nm2A;
     box[4]=boxVectors[1][1]*nm2A;
     box[5]=boxVectors[1][2]*nm2A;
     box[6]=boxVectors[2][0]*nm2A;
     box[7]=boxVectors[2][1]*nm2A;
     box[8]=boxVectors[2][2]*nm2A;
    }
    // copy coordinates :
    if (atomlist==NULL) { // atomlist is not defined; therefore, provide all coords
//    if (1) {
     for (i=0, rptr=r ; i < natoms ; i++) {
      *(rptr++) = pos[i][0]*nm2A; //units
      *(rptr++) = pos[i][1]*nm2A;
      *(rptr++) = pos[i][2]*nm2A;
     }
    // compute plugin forces and energy
// cout << "CALLING DYNAMO"<<endl;
     ierr=master_dyna_plugin(iteration, r, fr, NULL, 0, &master_energy, &atomlist, usesPeriodic, box); // might return valid atomlist
    // copy plugin forces
//=============
     if (qdble) { // double precision version
      double *frc = (double*) cu.getPinnedBuffer();
      memset(frc,0,3*natoms*sizeof(double));
      if (atomlist!=NULL) { // atom indices provided; use them for adding forces
       for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
        i=*aptr - 1; // for zero offset (e.g. first coordinate lives in r[0]
        j=3*i;
        frc[j]= fr[j]*str2omm_f; j++; //units
// cout << "x-force on atom: "<<j<<"="<<frc[j]<<endl;
        frc[j]= fr[j]*str2omm_f; j++;
        frc[j]= fr[j]*str2omm_f;
       }
      } else { // no atomlist provided; loop over all atoms
       for (j=0 ; j < 3*natoms ; j++) {
        frc[j]= fr[j]*str2omm_f; //units
       }
      } // atomlist
//=============
     } else { // single precision, identical code
//      cout << "CUDA DOUBLE PREC: "<<qdble<<endl;
      float *frc = (float*) cu.getPinnedBuffer();
      memset(frc,0,3*natoms*sizeof(float));
      if (atomlist!=NULL) { // atom indices provided; use them for adding forces
//      if (0) { // atom indices provided; use them for adding forces
       for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
        i=*aptr - 1; // for zero offset (e.g. first coordinate lives in r[0]
        j=3*i;
        frc[j]= (fr[j]*str2omm_f);j++; //units
        frc[j]= (fr[j]*str2omm_f);j++;
        frc[j]= (fr[j]*str2omm_f);
       }
      } else { // no atomlist provided; loop over all atoms
       for (j=0 ; j < 3*natoms ; j++) {
        frc[j]= fr[j]*str2omm_f; //units
 //cout << "force on atom: "<<j<<"="<<frc[j]<<endl;
       }
      } // atomlist
     } //qdble
    } else { // atomlist not null : loop over only the desired indices
     for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
      j=*aptr - 1;
      rptr=r + 3*j ;
      *(rptr++) = pos[j][0]*nm2A; //units
      *(rptr++) = pos[j][1]*nm2A;
      *(rptr)   = pos[j][2]*nm2A;
     }
//
     ierr=master_dyna_plugin(iteration, r, fr, NULL, 0, &master_energy, &atomlist, usesPeriodic, box); // atomlist should not be modified in this call
//
     if (qdble) { // double
      double *frc = (double*) cu.getPinnedBuffer();
      memset(frc,0,3*natoms*sizeof(double));
      for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
       i=*aptr - 1; // zero offset (see above)
       j=3*i ;
       frc[j]= fr[j]*str2omm_f;j++; //units
       frc[j]= fr[j]*str2omm_f;j++;
       frc[j]= fr[j]*str2omm_f;
      }
     } else { // single
      float *frc = (float*) cu.getPinnedBuffer();
      memset(frc,0,3*natoms*sizeof(float));
      for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
       i=*aptr - 1; // zero offset (see above)
       j=3*i ;
       frc[j]= fr[j]*str2omm_f;j++; //units
       frc[j]= fr[j]*str2omm_f;j++;
       frc[j]= fr[j]*str2omm_f;
      }
     } // qdble
    } // atomlist == NULL
    // upload forces to device
//    cu.setAsCurrent();
    cuMemcpyHtoDAsync(dynamoForces->getDevicePointer(), cu.getPinnedBuffer(), dynamoForces->getSize()*dynamoForces->getElementSize(), stream);
    cuEventRecord(syncEvent, stream);
}

double CudaCalcDynamoForceKernel::addForces(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return 0;
    // Wait until executeOnWorkerThread() is finished.
    cu.getWorkThread().flush();
    cuStreamWaitEvent(cu.getCurrentStream(), syncEvent, 0);
    // Add in the forces.
    if (includeForces) {
        void* args[] = {&dynamoForces->getDevicePointer(), &cu.getForce().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer()};
        cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
    }
    // Return plugin energy.
    master_energy*=str2omm_e;
    return master_energy;
}
