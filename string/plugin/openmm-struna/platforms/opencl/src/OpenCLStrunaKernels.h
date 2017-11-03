#ifndef OPENCL_STRUNA_KERNELS_H_
#define OPENCL_STRUNA_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "StrunaKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/opencl/OpenCLContext.h"
#include "openmm/opencl/OpenCLArray.h"
#include <stdbool.h>
#include "struna.h"
#include <vector>

namespace StrunaPlugin {

/**
 * This kernel is invoked by StrunaForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcStrunaForceKernel : public CalcStrunaForceKernel {
public:
    OpenCLCalcStrunaForceKernel(std::string name, const OpenMM::Platform& platform, OpenMM::ContextImpl& contextImpl, OpenMM::OpenCLContext& cl) :
            CalcStrunaForceKernel(name, platform), contextImpl(contextImpl), cl(cl), hasInitialized(false), strunaForces(NULL), r(NULL), fr(NULL), atomlist(NULL),\
            natoms(0) {
    }
    ~OpenCLCalcStrunaForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the StrunaForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const StrunaForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * The is called by the pre-computation to start the calculation running.
     */
    void beginComputation(bool includeForces, bool includeEnergy, int groups);
    /**
     * This is called by the worker thread to do the computation.
     */
    void executeOnWorkerThread();
    /**
     * This is called by the post-computation to add the forces to the main array.
     */
    double addForces(bool includeForces, bool includeEnergy, int groups);
private:
    class ExecuteTask;
    class StartCalculationPreComputation;
    class AddForcesPostComputation;
    bool hasInitialized; 
    _Bool usesPeriodic;
    OpenMM::ContextImpl& contextImpl;
    OpenMM::OpenCLContext& cl;
    OpenMM::OpenCLArray* strunaForces;
    cl::Kernel addForcesKernel;
    cl::CommandQueue queue;
    cl::Event syncEvent;
    int forceGroupFlag;
    std::vector<OpenMM::Vec3> pos;
    OpenMM::Vec3 boxVectors[3];
// VO 2017
    double sm_energy; // plugin energy
    int natoms;    // number of particles
    int* atomlist; // list of atom indices involved in restraints
    double *r, *fr; // positions and forces
    double box[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ; // box vectors
};

} // namespace StrunaPlugin

#endif /*OPENCL_STRUNA_KERNELS_H_*/
