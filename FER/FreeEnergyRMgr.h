/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

//-------------------------------------------------------------------------
// ARestraintManager contains a (potentially long) list of restraint POINTERS
// written by David Hurwitz, March to May 1998.
//-------------------------------------------------------------------------
#if !defined(RMGR_HPP)
  #define RMGR_HPP

//typedef ARestraint* pRestr;

// to start, there's room for this number of restraint pointers
// each time array size is exceeded, its size is increased by this many times.
const int kNumToStart = 1024;
const int kMultiplier = 4;

class GlobalMasterFreeEnergy;

class ARestraintManager {
private:
  ARestraint** Restraints; // list of restraint pointers
  int  NumRestraints;        // number of pointers in the list
  int  MaxNum;               // max num pointers without allocating more mem
  AFixedPosRestraint  Dummy; // for setting ARestraint statics

public:
  ARestraintManager();
  ~ARestraintManager();
  ARestraint*  operator[] (int Index);
  void   Add(ARestraint* pRestraint);
  int    GetNumRestraints() {return(NumRestraints);}
  double   AddForces(GlobalMasterFreeEnergy& CFE); // returns energy
  void   PrintRestraintInfo();
  void   Print_dU_dLambda();
  double Sum_dU_dLambdas();
  Bool_t ThereIsAForcingRestraint();
  void   PrintPreInfo(int Index);
  void   SetLambdaKf(double LambdaKf)   {Dummy.SetLambdaKf(LambdaKf);}
  void   SetLambdaRef(double LambdaRef) {Dummy.SetLambdaRef(LambdaRef);}
  void   SetLambdas(double LambdaKf, double LambdaRef)
  {
     Dummy.SetLambdaKf(LambdaKf);
     Dummy.SetLambdaRef(LambdaRef);
  }
};

#endif
