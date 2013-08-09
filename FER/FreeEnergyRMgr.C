/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <memory.h>
#include <string.h>
// #include <iomanip.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "Vector.h"
#include "FreeEnergyGroup.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"

#include "NamdTypes.h"
#include "GlobalMaster.h"
#include "GlobalMasterFreeEnergy.h"

#define DEBUGM
#include "Debug.h"

ARestraintManager::ARestraintManager() {
//------------------------------------------------------------------------
// allocate space for restraint POINTERS
//------------------------------------------------------------------------
  Restraints = new ARestraint*[kNumToStart];
  NumRestraints = 0;
  MaxNum = kNumToStart;
}


ARestraintManager::~ARestraintManager() {
//------------------------------------------------------------------------
// when the RestraintManager goes out of scope,
// free the space that was allocated for each restraint
// (see: ARestraint* GetRestraint(char*, int&)
// then, free the space that was allocated for the pointers
//------------------------------------------------------------------------
  for (int i=0; i<NumRestraints; i++) {
    delete Restraints[i];
  }
  delete []Restraints;
}


ARestraint* ARestraintManager::operator[] (int Index) {
//------------------------------------------------------------------------
// get a pointer
//------------------------------------------------------------------------
  ASSERT( (Index>=0) && (Index<NumRestraints) );
  return(Restraints[Index]);
}


void ARestraintManager::Add(ARestraint* pRestraint) {
//------------------------------------------------------------------------
// add a pointer to the list.  if there's not enough room, make room.
//------------------------------------------------------------------------
  ARestraint**  newRestraints;

  // if there's no room for a new pointer
  if (NumRestraints == MaxNum) {
    // create an array with more space
    MaxNum *= kMultiplier;
    newRestraints = new ARestraint*[MaxNum];
    // fast copy from the full array to the new one (memcpy(dest, src, bytes))
    memcpy(newRestraints, Restraints, sizeof(ARestraint*)*NumRestraints);
    // return the space used for the full array
    delete []Restraints;
    // point to the bigger array
    Restraints = newRestraints;
  }

  // add the int to the int array
  Restraints[NumRestraints] = pRestraint;
  NumRestraints++;
}

double ARestraintManager::AddForces(GlobalMasterFreeEnergy& CFE) {
//---------------------------------------------------------------------------
// each restraint is responsible for all steps of force computation (VO: 2013)
//---------------------------------------------------------------------------
  double total_energy=0.;
  for (int i=0; i<NumRestraints; i++) { 
#ifdef DEBUGM
   Restraints[i]->ApplyForce(CFE, 1, 0.00001); // run an finite difference test at each calculation
#else
   Restraints[i]->ApplyForce(CFE);
#endif
   total_energy += Restraints[i]->GetEnergy();
   return total_energy;
  }
}
//
double ARestraintManager::Sum_dU_dLambdas() {
//---------------------------------------------------------------------------
// sum up dU/dLambda from each forcing restraint
//---------------------------------------------------------------------------
  double Sum=0;

  for (int i=0; i<NumRestraints; i++) {
    Sum += Restraints[i]->Get_dU_dLambda();
  }
  return(Sum);
}


Bool_t ARestraintManager::ThereIsAForcingRestraint() {
//---------------------------------------------------------------------------
// return kTrue if there's at least one forcing restraint
//---------------------------------------------------------------------------
  for (int i=0; i<NumRestraints; i++) {
    if (Restraints[i]->IsMoving()) {
      return(kTrue);
    }
  }
  return(kFalse);
}

void ARestraintManager::PrintRestraintInfo() {
  char  Str[100];
  for (int i=0; i<NumRestraints; i++) {
    Restraints[i]->GetStr(Str);
    iout << "Free Energy Restraint #" << i+1 << " ("<< Str << "): ";
    Restraints[i]->PrintInfo();
    iout<< " Energy= "<< Restraints[i]->GetEnergy();
    iout << std::endl << endi;
  }
}

void ARestraintManager::Print_dU_dLambda() {
//---------------------------------------------------------------------------
// if restraint is a forcing restraint, print dU/dLambda.
//---------------------------------------------------------------------------
#if defined(_VERBOSE_PMF)
  for (int i=0; i<NumRestraints; i++) {
    if (Restraints[i]->IsMoving()) {
      PrintPreInfo(i);
      iout << "dU/dLambda = ";
      iout << Restraints[i]->Get_dU_dLambda() << std::endl << endi;
    }
  }
#endif
}
