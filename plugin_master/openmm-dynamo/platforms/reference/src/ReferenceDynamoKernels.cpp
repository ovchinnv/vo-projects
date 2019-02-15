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

#include "ReferenceDynamoKernels.h"
#include "DynamoForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include <cstring>
#include <cstdlib>

using namespace DynamoPlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static RealVec* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealVec*) data->periodicBoxVectors;
}

ReferenceCalcDynamoForceKernel::ReferenceCalcDynamoForceKernel(std::string name, \
               const OpenMM::Platform& platform, \
               OpenMM::ContextImpl& contextImpl) : CalcDynamoForceKernel(name, platform), \
               contextImpl(contextImpl), hasInitialized(false), atomlist(NULL), natoms(0), r(NULL), fr(NULL) {
}

ReferenceCalcDynamoForceKernel::~ReferenceCalcDynamoForceKernel() {
/*    if (hasInitialized) {
        master_done_plugin();
        free(r);
        free(fr);
    }*/
}

void ReferenceCalcDynamoForceKernel::initialize(const System& system, const DynamoForce& force) {
    //
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
    r=(_FLOAT*) calloc(3 * natoms, sizeof(double));
    fr=(_FLOAT*) calloc(3 * natoms, sizeof(double));
    // PBC flag
    usesPeriodic = system.usesPeriodicBoundaryConditions();
    OpenMM::Vec3 boxVectors[3];
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
    int ierr=master_init_plugin(natoms, m, q, inputfile, ilen, logfile, llen, &atomlist, usesPeriodic, box);
    free(m);
    free(q);
    hasInitialized = (ierr==0);
}

double ReferenceCalcDynamoForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    //
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    long int iteration = data->stepCount;
    _FLOAT* rptr; // pointer to positions array
    _FLOAT* fptr; // pointer to force array
    int* aptr; // pointer to atom index array
    int i, j, ierr;
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& frc = extractForces(context);
    _FLOAT master_energy;
    RealVec * boxVectors;
    if (usesPeriodic) {
     boxVectors = extractBoxVectors(context);
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
    //
    // copy coordinates :
    if (atomlist==NULL) { // atomlist is not defined; therefore, provide all coords
     for (i=0, rptr=r ; i < natoms ; i++) {
      *(rptr++) = pos[i][0]*nm2A; // include unit conversion
      *(rptr++) = pos[i][1]*nm2A;
      *(rptr++) = pos[i][2]*nm2A;
     }
    // compute plugin forces and energy
     ierr = (sizeof(_FLOAT)==sizeof(double)) ? \
      master_dyna_plugin(iteration, r, (double*)fr, NULL, 0, &master_energy, &atomlist, usesPeriodic, box) : \
      master_dyna_plugin(iteration, r, NULL, (float*)fr, 1, &master_energy, &atomlist, usesPeriodic, box)   ; // might return valid atomlist
    // copy plugin forces
     if (atomlist!=NULL) { // atom indices provided; use them for adding forces
      for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
       j=*aptr - 1; // for zero offset (e.g. first coordinate lives in r[0])
       fptr=fr + 3*j ;
       frc[j][0]+= (*(fptr++))*str2omm_f; // include unit conversion
       frc[j][1]+= (*(fptr++))*str2omm_f;
       frc[j][2]+= (*(fptr))*str2omm_f;
      }
     } else { // no atomlist provided; loop over all atoms
      for (i=0, fptr=fr ; i < natoms ; i++) {
       frc[i][0]+= *(fptr++)*str2omm_f;
       frc[i][1]+= *(fptr++)*str2omm_f;
       frc[i][2]+= *(fptr++)*str2omm_f;
      }
     } // atomlist
    } else { // atomlist not null : loop over only the desired indices
     for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
      j=*aptr - 1;
      rptr=r + 3*j;
      *(rptr++) = pos[j][0]*nm2A;
      *(rptr++) = pos[j][1]*nm2A;
      *(rptr)   = pos[j][2]*nm2A;
     }
//
     ierr = (sizeof(_FLOAT)==sizeof(double)) ? \
      master_dyna_plugin(iteration, r, (double*)fr, NULL, 0, &master_energy, &atomlist, usesPeriodic, box) : \
      master_dyna_plugin(iteration, r, NULL, (float*)fr, 1, &master_energy, &atomlist, usesPeriodic, box)   ; // atomlist should not be modified in this call
//
     for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
      j=*aptr - 1;
      fptr=fr + 3*j ;
      frc[j][0]+= *((fptr++))*str2omm_f;
      frc[j][1]+= *((fptr++))*str2omm_f;
      frc[j][2]+= *((fptr))*str2omm_f;
     }
    } // atomlist == NULL
    master_energy*=str2omm_e ; // convert units
    return master_energy;
}
