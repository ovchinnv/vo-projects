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

/**
 * This tests the CUDA implementation of StrunaForce.
 */

#include "StrunaForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace StrunaPlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerStrunaOpenCLKernelFactories();

void testForce() {
    // Create a System that applies a force based on the distance between two atoms.

    const int numParticles = 4;
    System system;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i] = Vec3(i, 0.1*i, -0.3*i);
    }
    string script ="test.sm";
    string log ="test.log";
    // create script file below
    ofstream testfile;
    testfile.open(script);
    testfile << "output {\n";
    testfile << " minwarnlev=-1\n";
    testfile << "}\n";
    testfile << "smcv_init {\n";
    testfile << " defer_init_dyna=yes\n";
    testfile << " : init maxcv 1\n";
    testfile << " : add dist_com sele atomid=1 end sele atomid=3 force "<<2.0/418.4<<" \n"; // 2kJ/mol/nm^2
    testfile << " : set ind 1 cvz 0 col main \n";
    testfile << " : dyna rstr\n";
    testfile << "}\n";
    testfile.close();
    //
    StrunaForce* struna = new StrunaForce(script, log);
    system.addForce(struna);
    LangevinIntegrator integ(300.0, 1.0, 1.0);
    Platform& platform = Platform::getPlatformByName("OpenCL");
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Compute the forces and energy.

    State state = context.getState(State::Energy | State::Forces);
    Vec3 delta = positions[0]-positions[2];
    double energy = delta.dot(delta)
    double dist = sqrt(energy);
    ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(-delta/dist, state.getForces()[0], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(), state.getForces()[1], 1e-5);
    ASSERT_EQUAL_VEC(delta/dist, state.getForces()[2], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(), state.getForces()[3], 1e-5);
}

int main(int argc, char* argv[]) {
    try {
        registerStrunaOpenCLKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("OpenCL").setPropertyDefaultValue("OpenCLPrecision", string(argv[1]));
        testForce();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
