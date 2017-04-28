%module openmmstruna


%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"

%{
#include "StrunaForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
%}

namespace StrunaPlugin {

class StrunaForce : public OpenMM::Force {
public:
    StrunaForce(const std::string& script);
    const std::string& getScript() const;
};

}
