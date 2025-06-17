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
 * Contributors: V. Ovchinnikov (DYNAMO Plugin interface)                     *
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
#ifdef __DYNAMO_SUBSET
        if (atomlist_d) {cudaFree(atomlist_d); atomlist_d=NULL;}
#endif
    }
    cu.setAsCurrent();
    if (dynamoForces != NULL) delete dynamoForces;
#ifdef __DYNAMO_SUBSET
    if (dynamoSubForces != NULL) delete dynamoSubForces;
#endif
    cuStreamDestroy(stream);
    cuEventDestroy(syncEvent);
}

void CudaCalcDynamoForceKernel::initialize(const System& system, const DynamoForce& force) {
    natoms = system.getNumParticles();
    //
    cu.setAsCurrent();
    cuStreamCreate(&stream, CU_STREAM_NON_BLOCKING);
    cuEventCreate(&syncEvent, CU_EVENT_DISABLE_TIMING);
    elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
#ifdef __DYNAMO_DEBUG
    cout<<__STRNG(_WHOAMI)<<" OPENMM is using "<<(cu.getUseDoublePrecision()?"DOUBLE":"SINGLE")<<" PRECISION"<<endl;
#endif
    dynamoForces = new CudaArray(cu, 3*natoms, elementSize, "dynamoForces"); //allocate device force array
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaDynamoKernelSources::dynamoForce, defines);
#ifdef __DYNAMO_SUBSET
    dynamoSubForces = new CudaArray(cu, 3*natoms, elementSize, "dynamoSubForces"); //overdimensioned
    cudaMemset((void*)dynamoSubForces->getDevicePointer(), 0.0, 3*natoms*elementSize); // do not know whether needed
    expandSubForcesKernel = cu.getKernel(module, "expandSubForces");
#endif
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
    r=( _FLOAT *) calloc(3 * natoms, sizeof(_FLOAT));
    fr=(_FLOAT *) calloc(3 * natoms, sizeof(_FLOAT));
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
    _FLOAT *m=NULL; //mass
    _FLOAT *q=NULL; //charge
    m = (_FLOAT*) malloc(natoms * sizeof(_FLOAT)); // allocate memory
    for (int i = 0; i < natoms; i++)
        m[i] = system.getParticleMass(i);

    q = (_FLOAT*) calloc(natoms, sizeof(_FLOAT));
    // If there's a NonbondedForce, get charges from it (otherwise, they will remain zero)

    for (int j = 0; j < system.getNumForces(); j++) {
        const NonbondedForce* nonbonded = dynamic_cast<const NonbondedForce*>(&system.getForce(j));
        if (nonbonded != NULL) {
            double sigma, epsilon, qi;
            for (int i = 0; i < natoms; i++) {
                nonbonded->getParticleParameters(i, qi, sigma, epsilon);
                q[i]=qi;
            }
        }
    }
    // initialize dynamo
    // NOTE: here, we are finding out whether the plugin interface is compiled in double or single precision (NOT OPENMM!)
    // Note that the actual fortran plugin can _also_ be compiled in double/single prec
    int ierr = (sizeof(_FLOAT)==sizeof(double)) ? \
      master_init_plugin(natoms, 0, (double*)m, NULL, (double*)q, NULL, inputfile, ilen, logfile, llen, &atomlist, usesPeriodic, (double*)box, NULL) : \
      master_init_plugin(natoms, 1, NULL, (float*)m, NULL, (float*)q, inputfile, ilen, logfile, llen, &atomlist, usesPeriodic, NULL, (float*)box);

    if (atomlist!=NULL) { // atom indices provided; check their range
     int atomid;
     for ( int l = 0 ; l++ < atomlist[0]  ; ) { // first value in atomlist is the list length
      atomid = atomlist[l] - 1 ; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the plugin
      if ( atomid < 0 || atomid >= natoms ) {
       cerr<<__STRNG(_whoami)<<"Atom ID"<<atomid<<"is out of range"<<endl; ierr=1;
      } // atomid valid
     } // over list
#ifdef __DYNAMO_SUBSET
// allocate & populate device atomlist
     natoms_requested=atomlist[0];
     cudaError_t cuda_err=cudaMalloc((void**)&(atomlist_d),(natoms_requested+1)*sizeof(int));
     if (cuda_err!=cudaSuccess) {ierr=1;cerr<<__STRNG(_whoami)<<"Error allocating CUDA memory:"<<cudaGetErrorString(cuda_err)<<endl;}
     cudaMemcpy(atomlist_d,atomlist,(natoms_requested+1)*sizeof(int),cudaMemcpyHostToDevice);
#else
     natoms_requested=natoms;
#endif
    } else {
     natoms_requested=natoms;
    }
    free(m);
    free(q);
    pos.resize(natoms); // this means we are requesting all atoms to be communicated from the GPU at every time step, no ?
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
    _FLOAT * rptr; // pointer to coordinate array
    int* aptr; // pointer to atom index array
    int ii, i, j, ierr;
    // buffer for uploading forces to the device:
    bool qdble=cu.getUseDoublePrecision();
#ifdef __DYNAMO_DEBUG
    cout<<__STRNG(_WHOAMI)<<" OPENMM is using "<<(qdble?"DOUBLE":"SINGLE")<<" PRECISION"<<endl;
    CUresult res ; // for CUDA error status
#endif
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
    if (qdble) {
     double *frc = (double*) cu.getPinnedBuffer();
#include "dynamoRunPlugin.h"
    } else {
     float *frc = (float*) cu.getPinnedBuffer();
#include "dynamoRunPlugin.h"
    }
    cu.setAsCurrent();
//    cuMemcpyHtoDAsync(dynamoForces->getDevicePointer(), cu.getPinnedBuffer(), dynamoForces->getSize()*dynamoForces->getElementSize(), stream);
#ifdef __DYNAMO_SUBSET
    cuMemcpyHtoDAsync(dynamoSubForces->getDevicePointer(), cu.getPinnedBuffer(), 3*natoms_requested*dynamoSubForces->getElementSize(), stream);
#else
    cuMemcpyHtoDAsync(dynamoForces->getDevicePointer(), cu.getPinnedBuffer(), 3*natoms_requested*dynamoForces->getElementSize(), stream);
// DBG : print dynamoForces
//  double *frc = (double*) cu.getPinnedBuffer();
//  for ( int i=0 ; i<3*natoms ; i+=3) { cout<<i<<" frc ="<<"("<<frc[i]<<","<<frc[i+1]<<","<<frc[i+2]<<")"<<endl; }
#endif
//    if (res!=CUDA_SUCCESS) {cerr<<__STRNG(_whoami)<<"cuMemcpyHostToDAsync error:"<<res<<endl;}
    cuEventRecord(syncEvent, stream);
}

double CudaCalcDynamoForceKernel::addForces(bool includeForces, bool includeEnergy, int groups) {
// return 0;
    if ((groups&forceGroupFlag) == 0)
        return 0;
    // Wait until executeOnWorkerThread() is finished.
    cu.getWorkThread().flush();
    cuStreamWaitEvent(cu.getCurrentStream(), syncEvent, 0);
    // Add in the forces.

#ifdef __DYNAMO_DEBUG
//==== DBG check index array === NOTE: the inds are mangled after a few iterations !
    int * indexArray = (int*) malloc(natoms*sizeof(int));
    for (int i=0 ; i<natoms ; i++) { indexArray[i]=-i; }
//    cuMemcpyDtoH(&indexArray[0], cu.getAtomIndexArray().getDevicePointer(), natoms*sizeof(int));
    cudaError_t err=cudaMemcpy(indexArray, (void *)cu.getAtomIndexArray().getDevicePointer(), natoms*sizeof(int),cudaMemcpyDeviceToHost); // same as above
    if (err!=cudaSuccess) {cerr<<__STRNG(_whoami)<<"CUDA indexArr error:"<<cudaGetErrorString(err)<<endl;}
//    for (int i=0 ; i<natoms ; i++) { cout<<i<<"~="<<indexArray[i]<<endl; } // these are not equal, i.e. atoms migrate, which is why AtomIndexArray is needed !

// DBG : print dynamoForces
    double *fdebug = (double*) malloc(3*natoms*elementSize);
/*
    CUresult res=cuMemcpyDtoH(&fdebug[0], dynamoForces->getDevicePointer(), 3*natoms*elementSize); // does NOT work error #1 invalid argument
    if (res!=CUDA_SUCCESS) {cerr<<__STRNG(_whoami)<<"CUDA fdebug error:"<<CudaContext::getErrorString(res)<<endl;}
// check dynamoForces:
    cout<<"INIT:"<<dynamoForces->isInitialized()<<endl;
    cout<<"SIZE:"<<dynamoForces->getSize()<<endl;
    cout<<"ESIZE, DBLE, FLOAT:"<<dynamoForces->getElementSize()<<" "<<sizeof(double)<<" "<<sizeof(float)<<endl;
    cout<<"NAME:"<<dynamoForces->getName()<<endl;
    cout<<"PTR:"<<dynamoForces->getDevicePointer()<<endl;
    cout<<"USE DOUBLE?:"<<cu.getUseDoublePrecision()<<endl;
*/
    err=cudaMemcpy(fdebug, (void*)dynamoForces->getDevicePointer(), 3*natoms*elementSize,cudaMemcpyDeviceToHost);
    if (err!=cudaSuccess) {cerr<<__STRNG(_whoami)<<"CUDA fdebug error:"<<cudaGetErrorString(err)<<endl;}
//    err=cudaMemcpy(fdebug, (void*)cu.getForce().getDevicePointer(), 3*natoms*elementSize,cudaMemcpyDeviceToHost);
    for (int i=0 ; i<3*natoms ; i+=3) { cout<<i<<":="<<"("<<fdebug[i]<<","<<fdebug[i+1]<<","<<fdebug[i+2]<<")"<<endl; }  // get all 0s !
//    for (int i=0 ; i<3*natoms ; i+=3) { cout<<i<<":="<<"("<<fr[i]<<","<<fr[i+1]<<","<<fr[i+2]<<")"<<endl; }  // fine
    free(fdebug);
//==== DBG
#endif
    if (includeForces) {
#ifdef __DYNAMO_SUBSET
     if (atomlist!=NULL) {
        void* args[] = {&dynamoSubForces->getDevicePointer(), &dynamoForces->getDevicePointer(), &atomlist_d};
        cu.executeKernel(expandSubForcesKernel, args, atomlist[0]); // kernel, arguments, workunits, blocksize=-1 ;
//        cu.executeKernel(expandSubForcesKernel, args, natoms); // kernel, arguments, workunits, blocksize=-1 ;
     }
#endif
     void* args[] = {&dynamoForces->getDevicePointer(), &cu.getForce().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer()};
     cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
    }
    // Return plugin energy.
    master_energy*=str2omm_e;
    return (double)master_energy;
}
