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

#include "OpenCLStrunaKernels.h"
#include "OpenCLStrunaKernelSources.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/opencl/OpenCLBondedUtilities.h"
#include "openmm/opencl/OpenCLForceInfo.h"
#include <cstring>
#include <map>

using namespace StrunaPlugin;
using namespace OpenMM;
using namespace std;

class OpenCLCalcStrunaForceKernel::StartCalculationPreComputation : public OpenCLContext::ForcePreComputation {
public:
    StartCalculationPreComputation(OpenCLCalcStrunaForceKernel& owner) : owner(owner) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        owner.beginComputation(includeForces, includeEnergy, groups);
    }
    OpenCLCalcStrunaForceKernel& owner;
};

class OpenCLCalcStrunaForceKernel::ExecuteTask : public OpenCLContext::WorkTask {
public:
    ExecuteTask(OpenCLCalcStrunaForceKernel& owner) : owner(owner) {
    }
    void execute() {
        owner.executeOnWorkerThread();
    }
    OpenCLCalcStrunaForceKernel& owner;
};

class OpenCLCalcStrunaForceKernel::AddForcesPostComputation : public OpenCLContext::ForcePostComputation {
public:
    AddForcesPostComputation(OpenCLCalcStrunaForceKernel& owner) : owner(owner) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        return owner.addForces(includeForces, includeEnergy, groups);
    }
    OpenCLCalcStrunaForceKernel& owner;
};

OpenCLCalcStrunaForceKernel::~OpenCLCalcStrunaForceKernel() {
    if (strunaForces != NULL)
        delete strunaForces;
}

void OpenCLCalcStrunaForceKernel::initialize(const System& system, const StrunaForce& force) {
    queue = cl::CommandQueue(cl.getContext(), cl.getDevice());
    int elementSize = (cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    strunaForces = new OpenCLArray(cl, 3*system.getNumParticles(), elementSize, "strunaForces");
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLStrunaKernelSources::strunaForce, defines);
    addForcesKernel = cl::Kernel(program, "addForces");
    forceGroupFlag = (1<<force.getForceGroup());
    cl.addPreComputation(new StartCalculationPreComputation(*this));
    cl.addPostComputation(new AddForcesPostComputation(*this));

    natoms = system.getNumParticles();
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
    //
    usesPeriodic = system.usesPeriodicBoundaryConditions();
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
    int ierr=sm_init_plugin(natoms, m, q, inputfile, ilen, logfile, llen, &atomlist);
    free(m);
    free(q);
    hasInitialized = true;
}

double OpenCLCalcStrunaForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    // This method does nothing.  The actual calculation is started by the pre-computation, continued on
    // the worker thread, and finished by the post-computation.
    
    return 0;
}

void OpenCLCalcStrunaForceKernel::beginComputation(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return;
    contextImpl.getPositions(pos);
    // The actual force computation will be done on a different thread.
    cl.getWorkThread().addTask(new ExecuteTask(*this));
}

void OpenCLCalcStrunaForceKernel::executeOnWorkerThread() {
    int iteration = cl.getStepCount();
    double* rptr; // pointer to coordinate array
    int* aptr; // pointer to atom index array
    int i, j, ierr;
    // buffer for uploading forces to the device:
    bool qdble=cl.getUseDoublePrecision();
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
      *(rptr++) = pos[i][0]*nm2A; //units
      *(rptr++) = pos[i][1]*nm2A;
      *(rptr++) = pos[i][2]*nm2A;
     }
    // compute plugin forces and energy
     ierr=sm_dyna_plugin(iteration, r, fr, &sm_energy, &atomlist, usesPeriodic, box); // might return valid atomlist
    // copy plugin forces
//=============
     if (qdble) { // double precision version
      double *frc = (double*) cl.getPinnedBuffer();
      if (atomlist!=NULL) { // atom indices provided; use them for adding forces
       for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
        i=*aptr - 1; // for zero offset (e.g. first coordinate lives in r[0]
        j=3*i;
        frc[j]= fr[j++]*str2omm_f; //units
        frc[j]= fr[j++]*str2omm_f;
        frc[j]= fr[j]*str2omm_f;
       }
      } else { // no atomlist provided; loop over all atoms
       for (j=0 ; j < 3*natoms ; j++) {
        frc[j]= fr[j++]*str2omm_f; //units
       }
      } // atomlist
//=============
     } else { // single precision, identical code
      float *frc = (float*) cl.getPinnedBuffer();
      if (atomlist!=NULL) { // atom indices provided; use them for adding forces
       for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
        i=*aptr - 1; // for zero offset (e.g. first coordinate lives in r[0]
        j=3*i;
        frc[j]= (float) fr[j++]*str2omm_f; //units
        frc[j]= (float) fr[j++]*str2omm_f;
        frc[j]= (float) fr[j]*str2omm_f;
       }
      } else { // no atomlist provided; loop over all atoms
       for (j=0 ; j < 3*natoms ; j++) {
        frc[j]= fr[j++]*str2omm_f; //units
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
     ierr=sm_dyna_plugin(iteration, r, fr, &sm_energy, &atomlist, usesPeriodic, box); // atomlist should not be modified in this call
//
     if (qdble) { // double
      double *frc = (double*) cl.getPinnedBuffer();
      for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
       i=*aptr - 1; // zero offset (see above)
       j=3*i ;
       frc[j]= fr[j++]*str2omm_f; //units
       frc[j]= fr[j++]*str2omm_f;
       frc[j]= fr[j]*str2omm_f;
      }
     } else { // single
      float *frc = (float*) cl.getPinnedBuffer();
      for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
       i=*aptr - 1; // zero offset (see above)
       j=3*i ;
       frc[j]= (float) fr[j++]*str2omm_f; //units
       frc[j]= (float) fr[j++]*str2omm_f;
       frc[j]= (float) fr[j]*str2omm_f;
      }
     } // qdble
    } // atomlist == NULL
    // upload forces to device
    queue.enqueueWriteBuffer(strunaForces->getDeviceBuffer(), CL_FALSE, 0, strunaForces->getSize()*strunaForces->getElementSize(), cl.getPinnedBuffer(), NULL, &syncEvent);
}

double OpenCLCalcStrunaForceKernel::addForces(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return 0;
    // Wait until executeOnWorkerThread() is finished.
    cl.getWorkThread().flush();
    vector<cl::Event> events(1);
    events[0] = syncEvent;
    syncEvent = cl::Event();
    queue.enqueueWaitForEvents(events);
    // Add in the forces.
    if (includeForces) {
        addForcesKernel.setArg<cl::Buffer>(0, strunaForces->getDeviceBuffer());
        addForcesKernel.setArg<cl::Buffer>(1, cl.getForceBuffers().getDeviceBuffer());
        addForcesKernel.setArg<cl::Buffer>(2, cl.getAtomIndexArray().getDeviceBuffer());
        cl.executeKernel(addForcesKernel, cl.getNumAtoms());
    }
    // Return the energy.
    return sm_energy;
}
