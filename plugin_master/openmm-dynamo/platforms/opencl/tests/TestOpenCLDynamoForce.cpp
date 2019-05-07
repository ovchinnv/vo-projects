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
 * This tests the CUDA implementation of DynamoForce.
 */

#include "DynamoForce.h"
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

using namespace DynamoPlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerDynamoOpenCLKernelFactories();

void testForce() {
    // Create a System that applies a force based on the distance between two atoms.

    const int numParticles = 22;
    double kf = 2 ; // kJ/mol/nm2
    System system;
    vector<Vec3> positions(numParticles);
    positions[0]=Vec3(-0.186, -1.490, -0.181);
    positions[1]=Vec3(-0.926, -2.447, -0.497);
    positions[2]=Vec3(0.008, -1.112, 0.725);
    positions[3]=Vec3(0.533, -0.396, 1.184);
    positions[4]=Vec3(-0.502, -1.632, 1.410);
    positions[5]=Vec3(-0.216, -2.590, 1.371);
    positions[6]=Vec3(-0.309, -1.255, 2.315);
    positions[7]=Vec3(-1.480, -1.560, 1.212);
    positions[8]=Vec3(-0.096, 2.144, -0.669);
    positions[9]=Vec3(0.871, 2.385, -0.588);
    positions[10]=Vec3(-0.565, 2.318, 0.197);
    positions[11]=Vec3(-0.520, 2.679, -1.400);
    positions[12]=Vec3(-0.172, 1.172, -0.892);
    positions[13]=Vec3(-1.139, 0.931, -0.973);
    positions[14]=Vec3(0.746, 0.780, -0.955);
    positions[15]=Vec3(1.713, 1.021, -0.873);
    positions[16]=Vec3(0.634, -0.654, -1.283);
    positions[17]=Vec3(0.099, -0.774, -2.218);
    positions[18]=Vec3(2.063, -1.223, -1.276);
    positions[19]=Vec3(2.670, -0.716, -2.057);
    positions[20]=Vec3(2.556, -1.051, -0.295);
    positions[21]=Vec3(2.070, -2.314, -1.490);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i]*=0.1 ; // to nm
    }

    string script ="test.sm";
    string log ="test.log";
    string structure_file ="test.pdb";
    // create script file below
    ofstream testfile;
    testfile.open(script.c_str());
    testfile << "output {\n";
    testfile << " minwarnlev=-1\n";
    testfile << "}\n";
    testfile << "structure_file="<<structure_file<<endl;
    testfile << "structure_filetype=PDB\n";
    testfile << "smcv_init {\n";
    testfile << " ! defer_init_dyna=yes\n";
    testfile << " : init maxcv 1\n";
    testfile << " : add dist_com sele atomid=1 end sele atomid=3 end force "<<kf/418.4<<" \n"; // kcal/mol/A2
    testfile << " : set ind 1 cvz 0 col main \n";
    testfile << " : copy main comp\n";
    testfile << " : copy main curr\n";
    testfile << " : fill col inst\n";
    testfile << " : print col inst name cv.ini\n";
    testfile << " : stat hist hnam test.hist hcol inst rene renm test.ene\n";
    testfile << " : hist add\n";
    testfile << " : stat\n";
    testfile << " : dyna rstr hisf 1 stat staf 1\n";
    testfile << "}\n";
    testfile << "smcv_done {\n";
    testfile << " : print col inst name cv.fin\n";
    testfile << "}\n";
    testfile.close();
    // create pdb file :
    testfile.open(structure_file.c_str());
    testfile << "ATOM      1  C   ALA     2      -0.186  -1.490  -0.181  4.00  4.00      DIA  C\n"; // NOTE coords change energy in Dynamo, but not in OPENMM !
    testfile << "ATOM      2  O   ALA     2      -0.926  -2.447  -0.497  0.00  0.00      DIA  O\n";
    testfile << "ATOM      3  NT  ALA     2       0.008  -1.112   0.725  5.00  5.00      DIA  N\n";
    testfile << "ATOM      4  HNT ALA     2       0.533  -0.396   1.184  0.00  0.00      DIA  H\n";
    testfile << "ATOM      5  CAT ALA     2      -0.502  -1.632   1.410  0.00  0.00      DIA  C\n";
    testfile << "ATOM      6  HT1 ALA     2      -0.216  -2.590   1.371  0.00  0.00      DIA  H\n";
    testfile << "ATOM      7  HT2 ALA     2      -0.309  -1.255   2.315  0.00  0.00      DIA  H\n";
    testfile << "ATOM      8  HT3 ALA     2      -1.480  -1.560   1.212  0.00  0.00      DIA  H\n";
    testfile << "ATOM      9  CAY ALA     2      -0.096   2.144  -0.669  0.00  0.00      DIA  C\n";
    testfile << "ATOM     10  HY1 ALA     2       0.871   2.385  -0.588  0.00  0.00      DIA  H\n";
    testfile << "ATOM     11  HY2 ALA     2      -0.565   2.318   0.197  0.00  0.00      DIA  H\n";
    testfile << "ATOM     12  HY3 ALA     2      -0.520   2.679  -1.400  0.00  0.00      DIA  H\n";
    testfile << "ATOM     13  CY  ALA     2      -0.172   1.172  -0.892  1.00  1.00      DIA  C\n";
    testfile << "ATOM     14  OY  ALA     2      -1.139   0.931  -0.973  0.00  0.00      DIA  O\n";
    testfile << "ATOM     15  N   ALA     2       0.746   0.780  -0.955  2.00  2.00      DIA  N\n";
    testfile << "ATOM     16  HN  ALA     2       1.713   1.021  -0.873  0.00  0.00      DIA  H\n";
    testfile << "ATOM     17  CA  ALA     2       0.634  -0.654  -1.283  3.00  3.00      DIA  C\n";
    testfile << "ATOM     18  HA  ALA     2       0.099  -0.774  -2.218  0.00  0.00      DIA  H\n";
    testfile << "ATOM     19  CB  ALA     2       2.063  -1.223  -1.276  0.00  0.00      DIA  C\n";
    testfile << "ATOM     20  HB1 ALA     2       2.670  -0.716  -2.057  0.00  0.00      DIA  H\n";
    testfile << "ATOM     21  HB2 ALA     2       2.556  -1.051  -0.295  0.00  0.00      DIA  H\n";
    testfile << "ATOM     22  HB3 ALA     2       2.070  -2.314  -1.490  0.00  0.00      DIA  H\n";
    testfile << "END\n";
    testfile.close();
    //
    DynamoForce* dynamo = new DynamoForce(script, log);
    system.addForce(dynamo);
    LangevinIntegrator integ(300.0, 1.0, 1.0);
    Platform& platform = Platform::getPlatformByName("OpenCL");
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Compute the forces and energy.

    State state = context.getState(State::Energy | State::Forces);
    Vec3 delta = positions[0]-positions[2];
    cout << delta <<endl;
    double energy = delta.dot(delta);
    double dist = sqrt(energy);
    cout << dist << endl;
    energy*=(0.5*kf) ;
    cout << energy << endl;
    cout << state.getPotentialEnergy() << endl;
    cout << state.getForces()[0] << endl;
    cout << -delta/dist*kf << endl;

    ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 1e-6);
    ASSERT_EQUAL_VEC(-delta*kf, state.getForces()[0], 1e-6);
    ASSERT_EQUAL_VEC(Vec3(), state.getForces()[1], 1e-6); // i.e. expect zero force
    ASSERT_EQUAL_VEC(delta*kf, state.getForces()[2], 1e-6);
    ASSERT_EQUAL_VEC(Vec3(), state.getForces()[3], 1e-6); // expect zero force
}

int main(int argc, char* argv[]) {
    try {
        registerDynamoOpenCLKernelFactories();
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
