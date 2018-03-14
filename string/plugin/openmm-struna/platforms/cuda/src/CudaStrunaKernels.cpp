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

#include "CudaStrunaKernels.h"
#include "CudaStrunaKernelSources.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include <cstring>
#include <map>

using namespace StrunaPlugin;
using namespace OpenMM;
using namespace std;

class CudaCalcStrunaForceKernel::StartCalculationPreComputation : public CudaContext::ForcePreComputation {
public:
    StartCalculationPreComputation(CudaCalcStrunaForceKernel& owner) : owner(owner) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        owner.beginComputation(includeForces, includeEnergy, groups);
    }
    CudaCalcStrunaForceKernel& owner;
};

class CudaCalcStrunaForceKernel::ExecuteTask : public CudaContext::WorkTask {
public:
    ExecuteTask(CudaCalcStrunaForceKernel& owner) : owner(owner) {
    }
    void execute() {
        owner.executeOnWorkerThread();
    }
    CudaCalcStrunaForceKernel& owner;
};

class CudaCalcStrunaForceKernel::CopyForcesTask : public ThreadPool::Task {
public:
    CopyForcesTask(CudaContext& cu, vector<Vec3>& forces) : cu(cu), forces(forces) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        // Copy the forces applied by STRUNA to a buffer for uploading.  This is done in parallel for speed.

        int numParticles = cu.getNumAtoms();
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        if (cu.getUseDoublePrecision()) {
            double* buffer = (double*) cu.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                const Vec3& p = forces[i];
                buffer[3*i] = p[0];
                buffer[3*i+1] = p[1];
                buffer[3*i+2] = p[2];
            }
        }
        else {
            float* buffer = (float*) cu.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                const Vec3& p = forces[i];
                buffer[3*i] = (float) p[0];
                buffer[3*i+1] = (float) p[1];
                buffer[3*i+2] = (float) p[2];
            }
        }
    }
    CudaContext& cu;
    vector<Vec3>& forces;
};

class CudaCalcStrunaForceKernel::AddForcesPostComputation : public CudaContext::ForcePostComputation {
public:
    AddForcesPostComputation(CudaCalcStrunaForceKernel& owner) : owner(owner) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        return owner.addForces(includeForces, includeEnergy, groups);
    }
    CudaCalcStrunaForceKernel& owner;
};

CudaCalcStrunaForceKernel::~CudaCalcStrunaForceKernel() {
    if (hasInitialized) {
        sm_done_plugin();
        free(r);
        free(fr);
    }
    cu.setAsCurrent();
    if (strunaForces != NULL)
        delete strunaForces;
    cuStreamDestroy(stream);
    cuEventDestroy(syncEvent);
}

void CudaCalcStrunaForceKernel::initialize(const System& system, const StrunaForce& force) {
    natoms = system.getNumParticles();
    //
    cu.setAsCurrent();
    cuStreamCreate(&stream, CU_STREAM_NON_BLOCKING);
    cuEventCreate(&syncEvent, CU_EVENT_DISABLE_TIMING);
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    strunaForces = new CudaArray(cu, 3*natoms, elementSize, "strunaForces");
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaStrunaKernelSources::strunaForce, defines);
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
    // initialize struna
    int ierr=sm_init_plugin(natoms, m, q, inputfile, ilen, logfile, llen, &atomlist, usesPeriodic, box);
    free(m);
    free(q);
    pos.resize(natoms);
    frc.resize(natoms);
    hasInitialized = true;
} // initialize

double CudaCalcStrunaForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    // This method does nothing.  The actual calculation is started by the pre-computation, continued on
    // the worker thread, and finished by the post-computation.
    return 0;
}

void CudaCalcStrunaForceKernel::beginComputation(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return;
    contextImpl.getPositions(pos);
    // The actual force computation will be done on a different thread.
    cu.getWorkThread().addTask(new ExecuteTask(*this));
}

void CudaCalcStrunaForceKernel::executeOnWorkerThread() {
    int iteration = cu.getStepCount();
    double *rptr; // pointer to positions array
    double *fptr; // pointer to force array
    int *aptr; // pointer to atom index array
    int i, j, ierr;
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
     for (i=0, rptr=r ; i < natoms ; i++) {
      *(rptr++) = pos[i][0]*nm2A;
      *(rptr++) = pos[i][1]*nm2A;
      *(rptr++) = pos[i][2]*nm2A;
     }
    // compute plugin forces and energy
     ierr=sm_dyna_plugin(iteration, r, fr, NULL, 0, &sm_energy, &atomlist, usesPeriodic, box); // might return valid atomlist
    // copy plugin forces
     if (atomlist!=NULL) { // atom indices provided; use them for adding forces
      for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
       j=*aptr - 1; // for zero offset (e.g. first coordinate lives in r[0]
       fptr=fr + 3*j ;
       frc[j][0]= (*(fptr++))*str2omm_f;
       frc[j][1]= (*(fptr++))*str2omm_f;
       frc[j][2]= (*(fptr))*str2omm_f;
      }
     } else { // no atomlist provided; loop over all atoms
      for (i=0, fptr=fr ; i < natoms ; i++) {
       frc[i][0]= (*(fptr++))*str2omm_f;
       frc[i][1]= (*(fptr++))*str2omm_f;
       frc[i][2]= (*(fptr++))*str2omm_f;
      }
     } // atomlist
    } else { // atomlist not null : loop over only the desired indices
     for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
      j=*aptr - 1;
      rptr=r + 3*j ;
      *(rptr++) = pos[j][0]*nm2A; // units
      *(rptr++) = pos[j][1]*nm2A;
      *(rptr)   = pos[j][2]*nm2A;
     }
//
     ierr=sm_dyna_plugin(iteration, r, fr, NULL, 0, &sm_energy, &atomlist, usesPeriodic, box); // atomlist should not be modified in this call
//
     for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
      j=*aptr - 1;
      fptr=fr + 3*j ;
      frc[j][0]= (*(fptr++))*str2omm_f; // convert units
      frc[j][1]= (*(fptr++))*str2omm_f;
      frc[j][2]= (*(fptr))*str2omm_f;
     }
    } // atomlist == NULL
    //
    // Upload the forces to the device.
    CopyForcesTask task(cu, frc);
    cu.getPlatformData().threads.execute(task);
    cu.getPlatformData().threads.waitForThreads();
    cu.setAsCurrent();
    cuMemcpyHtoDAsync(strunaForces->getDevicePointer(), cu.getPinnedBuffer(), strunaForces->getSize()*strunaForces->getElementSize(), stream);
    cuEventRecord(syncEvent, stream);
}

double CudaCalcStrunaForceKernel::addForces(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return 0;
    // Wait until executeOnWorkerThread() is finished.
    cu.getWorkThread().flush();
    cuStreamWaitEvent(cu.getCurrentStream(), syncEvent, 0);
    // Add in the forces.
    if (includeForces) {
        void* args[] = {&strunaForces->getDevicePointer(), &cu.getForce().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer()};
        cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
    }
    // Return plugin energy.
    sm_energy*=str2omm_e;
    return sm_energy;
}
