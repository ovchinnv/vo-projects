#include "../plugin.H"

#ifndef GLOBALMASTERDY_H
#define GLOBALMASTERDY_H
#endif

class Molecule;
class SubmitReduction;

class GlobalMasterDY : public GlobalMaster {
 public:
  GlobalMasterDY();
  ~GlobalMasterDY();
 private:
  Bool initialized;
  void calculate();
  void destroy();
  Molecule *molecule;
  SimParameters * params;
  long int iteration; // not sure we need to keep a local counter
//  SubmitReduction *reduction;
  __CFLOAT sm_energy;
  __CINT *atomlist=NULL; // to maintain a list of atoms that are needed by plugin; memory managed by FORTRAN side
  __CFLOAT *r=NULL, *fr=NULL ; //positions and forces
  __CFLOAT box[12] ; // three unit cell vectors followed by cell origin
  __CBOOL usesPeriodic ;// whether the simulation uses PBC
};
