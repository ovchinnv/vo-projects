#ifndef OPENMM_STRUNAFORCE_H_
#define OPENMM_STRUNAFORCE_H_

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

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <string>
#include "internal/windowsExportStruna.h"

namespace StrunaPlugin {

/**
 * This class implements a connection between OpenMM and STRUNA (http://www.victorovchinnikov.com/research/string-method).  It is a Force object that you
 * add to the System with addForce().  Its behavior is defined by a STRUNA input script file name, which you pass to the constructor
 * as a string.
 *
 * Be aware the STRUNA numbers atoms starting from 1, whereas OpenMM numbers them starting from 0.
 */

class OPENMM_EXPORT_STRUNA StrunaForce : public OpenMM::Force {
public:
    /**
     * Create a StrunaForce.
     *
     * @param script    STRUNA input script file name
     * @param log       STRUNA input log file name
     */
    StrunaForce(const std::string& script, const std::string& log);
    /**
     * STRUNA input script file name and log file name
     */
    const std::string& getScript() const;
    const std::string& getLog() const;
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    std::string script, log;
};

} // namespace StrunaPlugin

#endif /*OPENMM_STRUNAFORCE_H_*/
