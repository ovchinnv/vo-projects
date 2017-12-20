%module openmmdynamo


%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"

%{
#include "DynamoForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
%}

namespace DynamoPlugin {

class DynamoForce : public OpenMM::Force {
public:
    DynamoForce(const std::string& script, const std::string& log);
    const std::string& getScript() const;
    const std::string& getLog() const;
};

}
