/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#if !defined(GLOBALMASTERFREEENERGY_H)
#define GLOBALMASTERFREEENERGY_H

#include "NamdTypes.h"
#include <map>

class Molecule;
class SimParameters;
class SubmitReduction;

class GlobalMasterFreeEnergy : public GlobalMaster {
 public:
  GlobalMasterFreeEnergy();
  ~GlobalMasterFreeEnergy();
 private:
  virtual void calculate();
  void initialize();
  void user_initialize();
  Molecule *molecule;
  SimParameters *simParams;
  char *config;
  SubmitReduction *reduction;
//
  ARestraintManager  m_RestraintManager;
  ALambdaManager     m_LambdaManager;
// VO 2013 : store a map of positions indexed by atomid
// kept public for simplicity
  std::map<AtomID, Position> positions;
//
public:
  // These all return -1 on error.
  int getAtomID(const char *segid, int resid, const char *aname);
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid);
  int requestAtom(int atomid);
  int addForce(int atomid, Force force);
  int getMoleculeSize(); // total number of atoms ; VO 2013
  void fetchPositions() {
   AtomIDList::const_iterator a_i = getAtomIdBegin();
   AtomIDList::const_iterator a_e = getAtomIdEnd();
   PositionList::const_iterator p_i = getAtomPositionBegin();
   for ( ; a_i != a_e; ++a_i, ++p_i ){ positions[*a_i] = *p_i; }
  }
// access to positions & molecule
// (kludge to work around design flaws)
  friend class ARestraint;
#ifdef FE_RESTRAINT_RMSD_FORTRAN
  friend class AnRMSDRestraint;
#endif
};
#endif
