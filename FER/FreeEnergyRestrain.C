/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.
// modifications by Victor Ovchinnikov 2013
#ifdef FE_RESTRAINT_RMSD_FORTRAN
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "bestfit.h" //external Fortran wrappers
#ifdef __cplusplus
}
#endif
#endif
//
#include <string.h>
#include <cstdlib>
#include <algorithm>
//#include <iomanip.h>
#include "common.h"
#include "InfoStream.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyEnums.h"
#include "Vector.h"
#include "FreeEnergyGroup.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"

#include "NamdTypes.h"
#include "GlobalMaster.h"
#include "GlobalMasterFreeEnergy.h"
#include "Molecule.h"

#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

// initialize static member variables
// (lambda is the same for all Moving restraints)
double  ARestraint::m_LambdaKf  = 1.0;
double  ARestraint::m_LambdaRef = 0.0;

ARestraint::ARestraint() {
//-----------------------------------------------------------------
// constructor for base class
//-----------------------------------------------------------------
  m_pGroups = NULL;
  m_pCOMs = NULL;
  m_NumGroups = 0;
}


ARestraint::~ARestraint() {
//-----------------------------------------------------------------
// free space that may have been allocated for Groups and COM's
//-----------------------------------------------------------------
  if (m_pGroups != NULL) {
    ASSERT(m_pCOMs != NULL);
    delete []m_pGroups;
    delete []m_pCOMs;
  }
}


void ARestraint::EarlyExit(const char* Str, int AtomID) {
//-----------------------------------------------------------------
// unrecoverable error
//-----------------------------------------------------------------
  char  NumStr[40];
  char  Msg[80];
  strcpy(Msg, Str);

  iout << "FreeEnergy: " << std::endl << endi;
  if (AtomID >= 0) {
   sprintf(NumStr, "%d", AtomID);
   strcat(Msg, " for AtomID: ");
   strcat(Msg, NumStr);
  }
  iout << "FreeEnergy: " << Msg;
  iout << std::endl << endi;
  NAMD_die("FreeEnergy: Fatal Error with Free Energy Restraints");
}

void ARestraint::SetGroup(AGroup& Group, int GroupIndex) {
//-----------------------------------------------------------------
// set one group of atoms
//-----------------------------------------------------------------
  ASSERT( (GroupIndex>=0) && (GroupIndex<m_NumGroups) );
  m_pGroups[GroupIndex] = Group;
}


void ARestraint::SetGroups(AGroup& Group1) {
//-----------------------------------------------------------------
// set one group of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 1);
  m_pGroups[0] = Group1;
}


void ARestraint::SetGroups(AGroup& Group1, AGroup& Group2) {
//-----------------------------------------------------------------
// set two groups of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 2);
  m_pGroups[0] = Group1;
  m_pGroups[1] = Group2;
}


void ARestraint::SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3) {
//-----------------------------------------------------------------
// set three groups of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 3);
  m_pGroups[0] = Group1;
  m_pGroups[1] = Group2;
  m_pGroups[2] = Group3;
}


void ARestraint::SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3, AGroup& Group4) {
//-----------------------------------------------------------------
// set four groups of atoms
//-----------------------------------------------------------------
  ASSERT(m_NumGroups >= 4);
  m_pGroups[0] = Group1;
  m_pGroups[1] = Group2;
  m_pGroups[2] = Group3;
  m_pGroups[3] = Group4;
}

void ARestraint::DistributeForce(int WhichGroup, Force force,
                                 GlobalMasterFreeEnergy& CFE) {
//----------------------------------------------------------------------
// Distribute Force among the group of atoms specified by WhichGroup
//
// note:  m_pGroups points to an array of Groups
//        m_pGroups[WhichGroup] references one of the Groups
//        m_pGroups[WhichGroup][i] returns an AtomID from the Group
//        (operator[] is defined to return an item from the Group)
//----------------------------------------------------------------------
  ASSERT( (WhichGroup>=0) && (WhichGroup<m_NumGroups) );

  int NumAtoms = m_pGroups[WhichGroup].GetSize();
  double const* Weights=m_pGroups[WhichGroup].GetWeights();
  double weight;

  // distribute Force according to mass of each atom in the group
  for (int i=0; i<NumAtoms; i++) {
    int aid = m_pGroups[WhichGroup][i];

    Vector scaledforce = Weights[i] * force ;
//
    DebugM(2,"Applying force "<<scaledforce<<" to atom "<<aid<<"\n");
    if ( (aid >= 0) && (aid < CFE.molecule->numAtoms ) ) {
     CFE.modifyForcedAtoms().add(aid);
     CFE.modifyAppliedForces().add( scaledforce );
    } else {
      EarlyExit( "Atom ID out of bounds:", aid);
    }//if
  }
}

void ARestraint::UpdateCOMs(GlobalMasterFreeEnergy& CFE) {
//-----------------------------------------------------------------
// calculate the center-of-mass of each group of atoms
//
// note:  m_pGroups points to an array of Groups
//        m_pGroups[i] references one of the Groups
//        m_pGroups[i][j] returns an AtomID from the Group
//        (operator[] is defined to return an item from the Group)
//-----------------------------------------------------------------
  int      i, j, aid;
  Position pos ;
  double const* Weights;
  double Weight;
//
  ASSERT(m_NumGroups > 0);
// fetch positions
  // for each group of atoms
  for (i=0; i<m_NumGroups; i++) {
    Weights=m_pGroups[i].GetWeights();
    Vector COM(0.,0.,0) ;
    // for each atom in the group
    for (j=0; j<m_pGroups[i].GetSize(); j++) {
      aid = m_pGroups[i][j]; // get atom ID
      Vector pos=CFE.positions[aid]; // fetch from map
      Weight = Weights[j];
      COM += pos * Weight;
    }
    m_pCOMs[i] = COM;
  }
}

void APosRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE) {
// for each center-of-mass
 Force force; // force vector
//
 for (int j=0; j<m_NumGroups; j++) {
  // update centers of mass
  UpdateCOMs(CFE);
  // apply restraining force in direction opposite to the gradient vector
  force = GetGradient(j) * (-1.0);
  DistributeForce(j, force, CFE);
 }
}

APosRestraint::APosRestraint() {
//-----------------------------------------------------------------
// each APosRestraint restrains 1 group of atoms to a location
//-----------------------------------------------------------------
  m_NumGroups = 1;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new Vector[m_NumGroups];
}


void APosRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the position for this position restraint
//--------------------------------------------------------------------
//
#if defined(_VERBOSE_PMF)
  iout << "Position = " << m_pCOMs[0];
  iout << "  Target = " << GetPosTarget();
  iout << "  Distance = " << GetDistance();
  iout << std::endl << endi;
#else
  iout<<"("<<m_pCOMs[0]<<")"<< GetPosTarget()<< "  "<<GetDistance()<< " | ";
#endif
}

double APosRestraint::GetE(Vector RefPos, double LambdaKf) {
//--------------------------------------------------------------------
// calculate and return the Energy for this position restraint.
//
//     E = (Kf/2) * (|ri - rref|)**2
//
// where |ri - rref| is the distance between a) the center-of-mass
// of the restrained atoms and b) the reference position.
//
// Note:  COM is calculated before this routine is called.
//--------------------------------------------------------------------
  return( ((m_Kf*LambdaKf)/2.0) * (m_pCOMs[0] - RefPos).length2() );
}


Vector APosRestraint::GetGrad(int /* WhichGroup */, Vector RefPos, double LambdaKf) {
//-------------------------------------------------------------------------
// calculate and return the gradient for this position restraint.
//     E = (Kf/2) * (|ri - rref|)**2
// return:  grad(E)
// Notes: COM is calculated before this routine is called.
//-------------------------------------------------------------------------
  return ( m_pCOMs[0] - RefPos ) * m_Kf * LambdaKf ;
}

double AFixedPosRestraint::GetEnergy() {
//--------------------------------------------------------------------
// return the Energy for this fixed position restraint.
//--------------------------------------------------------------------
  return(GetE(m_RefPos));
}


Vector AFixedPosRestraint::GetGradient(int WhichGroup) {
//-------------------------------------------------------------------------
// return the Gradient for this fixed position restraint.
//-------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefPos));
}


double ABoundPosRestraint::GetEnergy() {
//--------------------------------------------------------------------
// calculate and return the Energy for this bound position restraint.
//
// Note:  This is an exception because the form for the E term is
//        different from the other postion restraints.
//--------------------------------------------------------------------
  double  E, Dist, Diff;

  E = 0.0;
  Dist = (m_pCOMs[0] - (m_RefPos)).length();
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    Diff = Dist - m_RefDist;
    E = (m_Kf/2.0) * (Diff*Diff);
  }
  return(E);
}

Vector ABoundPosRestraint::GetGradient(int /* WhichGroup */) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this bound position restraint.
//
// Note:  This is an exception because the form for the E term is
//        different from the other postion restraints.
//---------------------------------------------------------------------------
  double  Dist;
  Vector Vec(0.,0.,0.);   // Vec is initialized to (0,0,0)

  // WhichGroup = 0;  // don't care -- there's only 1 atom restrained
  Dist = (m_pCOMs[0]-m_RefPos).length();
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    Vec = (m_pCOMs[0]- m_RefPos) * m_Kf * (Dist - m_RefDist) / __MAX(Dist, kALittle); //  prevent divide overflow
  }
  return(Vec);
}


double AMovingPosRestraint::GetEnergy() {
//--------------------------------------------------------------------
// return the Energy for this Moving position restraint.
//
// rref = lambda*r1 + (1-lambda)*r0.
// where r0 is the starting position and r1 is the final position
//--------------------------------------------------------------------
  Vector  RefPos;

  RefPos = m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef);
  return(GetE(RefPos, m_LambdaKf));
}


Vector AMovingPosRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving position restraint.
//
// rref = lambda*r1 + (1-lambda)*r0.
// where r0 is the starting position and r1 is the final position
//---------------------------------------------------------------------------
  Vector  RefPos;

  RefPos = m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefPos, m_LambdaKf));
}


double AMovingPosRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving position restraint
//---------------------------------------------------------------------------
  Vector  RefPos;
//  double   T1, T2, T3;

  RefPos = m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef);
//  T1 = (m_pCOMs[0][0] - RefPos[0]) * (m_StartPos[0] - m_StopPos[0]);
//  T2 = (m_pCOMs[0][1] - RefPos[1]) * (m_StartPos[1] - m_StopPos[1]);
//  T3 = (m_pCOMs[0][2] - RefPos[2]) * (m_StartPos[2] - m_StopPos[2]);
//  return( m_Kf * m_LambdaKf * (T1+T2+T3) );
    return m_Kf * m_LambdaKf * ( (m_pCOMs[0] - RefPos) * (m_StartPos - m_StopPos) );// replaced above by a scalar (dot) product
}


ADistRestraint::ADistRestraint() {
//-----------------------------------------------------------------------
// each ADistRestraint restrains the distance between 2 groups of atoms
//-----------------------------------------------------------------------
  m_NumGroups = 2;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new Vector[m_NumGroups];
}


void ADistRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the distance for this distance restraint
//--------------------------------------------------------------------
#if defined(_VERBOSE_PMF)
  iout << "Distance = "<< (m_pCOMs[0]-m_pCOMs[1]).length();
  iout << "  Target = ";<< GetDistTarget();
  iout << std::endl << endi;
#else
  iout <<(m_pCOMs[0]-m_pCOMs[1]).length()<< " "<< GetDistTarget()<< " | ";
#endif
}


double ADistRestraint::GetE(double RefDist, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the Energy for this distance restraint.
//
//     E = (Kf/2) * (di-dref)**2
//
// where di is the distance between 2 centers-of-mass of restrained atoms,
// and dref is the reference distance.
// Note:  COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  double Dist, Diff;

  Dist = (m_pCOMs[0]-m_pCOMs[1]).length();
  Diff = Dist - RefDist;
  return( ((m_Kf*LambdaKf)/2.0) * (Diff*Diff) );
}


Vector ADistRestraint::GetGrad(int WhichGroup,
                                double RefDist, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this distance restraint.
//
//     E = (Kf/2) * (di-dref)**2
//
// return:  grad(E)
//
// Notes: COM is calculated before this routine is called.
//        m_pCOMS[0 & 1] reference the COM's of each group of atoms
//---------------------------------------------------------------------------
  ASSERT( (WhichGroup==0) || (WhichGroup==1) );
//
  Vector Vec = m_pCOMs[1] - m_pCOMs[0];
  double Dist = Vec.length();
  return  (WhichGroup * 2. - 1.) * (  Vec * m_Kf * LambdaKf * (Dist - RefDist)/__MAX(Dist, kALittle) );// WG=0 => prefix=-1 ; WG=1=> prefix=+1
}

double AFixedDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this fixed distance restraint.
//---------------------------------------------------------------------------
  return(GetE(m_RefDist));
}


Vector AFixedDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed distance restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefDist));
}

double ABoundDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound distance restraint.
//---------------------------------------------------------------------------
  double Dist, E;
//
  E = 0.0;
  Dist = (m_pCOMs[0] - m_pCOMs[1]).length();
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    E = GetE(m_RefDist);
  }
  return(E);
}

Vector ABoundDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound distance restraint.
//---------------------------------------------------------------------------
  double  Dist;
  Vector Vec;

  Dist = (m_pCOMs[0] - m_pCOMs[1]).length();
  if (((m_Bound==kUpper) && (Dist>m_RefDist)) ||
      ((m_Bound==kLower) && (Dist<m_RefDist))) {
    Vec = GetGrad(WhichGroup, m_RefDist);
  }
  return(Vec);
}


double AMovingDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this Moving distance restraint.
//---------------------------------------------------------------------------
  double  RefDist;

  RefDist = m_StopDist*m_LambdaRef + m_StartDist*(1.0-m_LambdaRef);
  return(GetE(RefDist, m_LambdaKf));
}


Vector AMovingDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving distance restraint.
//---------------------------------------------------------------------------
  double  RefDist;
  
  RefDist = m_StopDist*m_LambdaRef + m_StartDist*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefDist, m_LambdaKf));
}


double AMovingDistRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving distance restraint
//---------------------------------------------------------------------------
  double  Dist;
  double  RefDist;
//
  Dist = (m_pCOMs[0] - m_pCOMs[1]).length();
  RefDist = m_StopDist*m_LambdaRef + m_StartDist*(1.0-m_LambdaRef);
  return( m_Kf * m_LambdaKf * (Dist-RefDist)*(m_StartDist-m_StopDist) );
}

void ADistRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE) {
// for each center-of-mass
 Vector Force;
 for (int j=0; j<m_NumGroups; j++) {
  // update centers of mass
  UpdateCOMs(CFE);
  // apply restraining force in opposite direction from gradient
  Force = GetGradient(j) * (-1.0);
  DistributeForce(j, Force, CFE);
 }
}

AnAngleRestraint::AnAngleRestraint() {
//-----------------------------------------------------------------------
// each AnAngleRestraint restrains the angle between 3 groups of atoms
//-----------------------------------------------------------------------
  m_NumGroups = 3;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new Vector[m_NumGroups];
}

void AnAngleRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the angle for this angle restraint
//--------------------------------------------------------------------
  double  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]) * (180/kPi);
#if defined(_VERBOSE_PMF)
  iout << "Angle = "<<Angle<< " degrees";
  iout << "  Target = "<<GetAngleTarget()*180/kPi<< " degrees"<< std::endl << endi;
#else
  iout <<Angle<<" "<<GetAngleTarget()*180/kPi<< " | ";
#endif
}


double AnAngleRestraint::GetE(double RefAngle, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the Energy for this angle restraint.
//
//     E = (Kf/2) * (Theta-ThetaRef)**2
//
// where Theta is the angle between 3 centers-of-mass of restrained atoms,
// m_pCOMs[0] -- m_pCOMs[1] -- m_pCOMs[2].
// ThetaRef is the reference angle.
//
// Note:  COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  double Angle, Diff;

  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  Diff = Angle - RefAngle;
  return( ((m_Kf*LambdaKf)/2.0) * (Diff*Diff) );
}

double AnAngleRestraint::GetAngle(Vector& A, Vector& B, Vector& C) {
//-----------------------------------------------------------------
// determine the angle formed by the points A-B-C
//-----------------------------------------------------------------
  double u;
  double a = (B-C).length();
  double b = (A-C).length();
  double c = (A-B).length();

  u = (a*a + c*c - b*b) / (2.0*a*c);
  // protect against acos(<-1.0) and acos(>1.0)
  if (u < -1.0) {u = -1.0;}
  if (u >  1.0) {u =  1.0;}
  return(acos(u));
}

Vector AnAngleRestraint::GetGrad(int WhichGroup,
                                  double RefAngle, double LambdaKf) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this angle restraint.
//
//     E = (Kf/2) * (Theta-ThetaRef)**2
//
// return:  grad(E)
//
// Notes: COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  Vector A, B, C;
  double  Angle;
  double  a, b, c;
  double  u;
  double  Const1, Const2, Const3, Const4, Const5;
  Vector Vec1, Vec2;

  ASSERT( (WhichGroup==0) || (WhichGroup==1) || (WhichGroup==2) );
  A = m_pCOMs[0];
  B = m_pCOMs[1];
  C = m_pCOMs[2];

  a = (B-C).length();
  b = (A-C).length();
  c = (A-B).length();

  u = (a*a + c*c - b*b) / (2.0*a*c);

  // protect against acos(<-1.0), acos(>1.0), sqrt(<0), and divide by 0
  if (u < -1.0)       {u = -1.0;}
  if (u > kAlmostOne) {u = kAlmostOne;}
  Angle = acos(u);

  Const1 = -1.0 / sqrt(1.0 - u*u);
  Const2 = m_Kf * LambdaKf * (Angle-RefAngle) * Const1;

  Const3 = -a/(2.0*c*c) + 1.0/(a+a) + (b*b)/(2.0*a*c*c);
  Const4 = -c/(2.0*a*a) + 1.0/(c+c) + (b*b)/(2.0*c*a*a);
  Const5 = -b/(a*c);

  if (WhichGroup == 0) {
    Vec1 = (A-C) * (Const5/b);
    Vec2 = (A-B) * (Const3/c);
  }
  else if (WhichGroup == 1) {
    Vec1 = (B-A) * (Const3/c);
    Vec2 = (B-C) * (Const4/a);
  }
  else {
    Vec1 = (C-A) * (Const5/b);
    Vec2 = (C-B) * (Const4/a);
  }
  return( (Vec1+Vec2)*Const2);
}

void AnAngleRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE) {
// for each center-of-mass
 Vector Force;
 for (int j=0; j<m_NumGroups; j++) {
  // update centers of mass
  UpdateCOMs(CFE);
  // apply restraining force in opposite direction from gradient
  Force = GetGradient(j) * (-1.0);
  DistributeForce(j, Force, CFE);
 }
}

double AFixedAngleRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this fixed angle restraint.
//---------------------------------------------------------------------------
  return(GetE(m_RefAngle));
}

Vector AFixedAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed angle restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefAngle));
}

double ABoundAngleRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound angle restraint.
//---------------------------------------------------------------------------
  double  E, Angle;

  E = 0.0;
  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  if (((m_Bound==kUpper) && (Angle>m_RefAngle)) ||
      ((m_Bound==kLower) && (Angle<m_RefAngle))) {
    E = GetE(m_RefAngle);
  }
  return(E);
}

Vector ABoundAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound angle restraint
//---------------------------------------------------------------------------
  double  Angle;
  Vector Vec;

  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  if (((m_Bound==kUpper) && (Angle>m_RefAngle)) ||
      ((m_Bound==kLower) && (Angle<m_RefAngle))) {
    Vec = GetGrad(WhichGroup, m_RefAngle);
  }
  return(Vec);
}


double AMovingAngleRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this Moving angle restraint.
//---------------------------------------------------------------------------
  double  RefAngle;

  RefAngle = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetE(RefAngle, m_LambdaKf));
}


Vector AMovingAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving angle restraint.
//---------------------------------------------------------------------------
  double  RefAngle;

  RefAngle = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefAngle, m_LambdaKf));
}


double AMovingAngleRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving angle restraint
//---------------------------------------------------------------------------
  double  Angle;
  double  RefAngle;

  Angle = GetAngle(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2]);
  RefAngle = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return( m_Kf * m_LambdaKf * (Angle-RefAngle)*(m_StartAngle-m_StopAngle) );
}


ADiheRestraint::ADiheRestraint() {
//----------------------------------------------------------------------------
// each ADiheRestraint restrains the dihedral angle between 4 groups of atoms
//----------------------------------------------------------------------------
  m_NumGroups = 4;
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new Vector[m_NumGroups];
}


void ADiheRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the dihedral angle for this dihedral restraint
//--------------------------------------------------------------------
  double  Dihedral;

  Dihedral = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]) * (180/kPi);
  Dihedral = GetDiheTarget1() * (180/kPi);
  Dihedral = GetDiheTarget2() * (180/kPi);

#if defined(_VERBOSE_PMF)
  iout << "Dihedral = "<<Dihedral<< " degrees"<<"  Target = "<<GetDiheTarget1() * (180/kPi)<< " degrees";
  if (TwoTargets()) {
    iout <<" to "<<GetDiheTarget2() * (180/kPi)<< " degrees";
  }
  iout<< std::endl << endi;
#else
  iout << Dihedral <<"  "<<GetDiheTarget1() * (180/kPi);
  if (TwoTargets()) {
    iout << ", "<<GetDiheTarget2() * (180/kPi);
  }
  iout << " | ";
#endif
}


double ADiheRestraint::GetE(double RefAngle, double Const) {
//---------------------------------------------------------------------------
// calculate and return the Energy for this angle restraint.
//
//    E = (E0/2) * (1 - cos(Chi - ChiRef))
//
// where Chi is the dihedral angle between 4 centers-of-mass of restrained atoms,
// m_pCOMs[0] -- m_pCOMs[1] -- m_pCOMs[2] -- m_pCOMs[3].
// ChiRef is the reference angle.
//
// Note:  COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  double  Angle;

  Angle = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  return( (Const/2.0) * (1.0 - cos(Angle-RefAngle)) );
}


Vector ADiheRestraint::GetGrad(int WhichGroup,
                                double RefAngle, double Const) {
//---------------------------------------------------------------------------
// calculate and return the gradient for this dihedral angle restraint.
//
//    E = (E0/2) * (1 - cos(Chi - ChiRef))
//
// return:  grad(E)
//
// Notes: COM's are calculated before this routine is called.
//---------------------------------------------------------------------------
  Vector A, B, C, D;

  ASSERT((WhichGroup==0)||(WhichGroup==1)||(WhichGroup==2)||(WhichGroup==3));

  if ((WhichGroup==0) || (WhichGroup==1)) {
    A = m_pCOMs[0];
    B = m_pCOMs[1];
    C = m_pCOMs[2];
    D = m_pCOMs[3];
  }
  // re-state the problem so the gradient is solved for either atoms 0 or 1
  else {
    A = m_pCOMs[3];
    B = m_pCOMs[2];
    C = m_pCOMs[1];
    D = m_pCOMs[0];
    if (WhichGroup==3) {WhichGroup=0;}
    if (WhichGroup==2) {WhichGroup=1;}
  }

  Vector CD(D - C);
  Vector CB(B - C);
  Vector BC(C - B);
  Vector BA(A - B);
  Vector AC(C - A);
  Vector CDxCB, BCxBA;
  Vector Vec;
  double  phi;
  double  top, bot, u;

  CDxCB = cross(CD,CB);
  BCxBA = cross(BC,BA);

  top = CDxCB*BCxBA; // dotpr
  bot = CDxCB.length() * BCxBA.length();

  u = top/bot;
  // protect against acos(<-1.0), acos(>1.0), sqrt(<0), and divide by 0
  if (u < kAlmostMinusOne) {u = kAlmostMinusOne;}
  if (u > kAlmostOne)      {u = kAlmostOne;}

  // get dihedral using atan
  phi = GetDihe(A,B,C,D);

  ASSERT((WhichGroup==0) || (WhichGroup==1));
  if (WhichGroup==0) {
    Vector dP1( 0,      0,      0    );
    Vector dP2( 0,      0,      0    );
    Vector dP3( 0,      0,      0    );
    Vector dP4( 0,     -BC[2],  BC[1]);
    Vector dP5( BC[2],  0,     -BC[0]);
    Vector dP6(-BC[1],  BC[0],  0    );
    Vec = gradU(CDxCB, BCxBA, dP1, dP2, dP3, dP4, dP5, dP6);
  }
  else {
    Vector dP1( 0,     -CD[2],  CD[1]);
    Vector dP2( CD[2],  0,     -CD[0]);
    Vector dP3(-CD[1],  CD[0],  0    );
    Vector dP4( 0,      AC[2], -AC[1]);
    Vector dP5(-AC[2],  0,      AC[0]);
    Vector dP6( AC[1], -AC[0],  0    );
    Vec = gradU(CDxCB, BCxBA, dP1, dP2, dP3, dP4, dP5, dP6);
  }
//
  Vec *= (Const/2.0) * sin(phi-RefAngle) * (-1.0/sqrt(1.0 - u*u));

  // flip gradient for negative angles
  if (phi < 0) {
    Vec *= -1.0;
  }
//
  return(Vec);
}


Vector ADiheRestraint::gradU(Vector& P1P2P3, Vector& P4P5P6,
                              Vector& dP1,    Vector& dP2,    Vector& dP3,
                              Vector& dP4,    Vector& dP5,    Vector& dP6) {
//----------------------------------------------------------------
// calculate the gradient for ((P1P2P3.dot.P4P5P6)/(mag(P1P2P3)*mag(P4P5P6)))
// P1P2P3 = (P1)i + (P2)j + (P3)k
// P4P5P6 = (P4)i + (P5)j + (P6)k
// dP1 = (d(P1)/dx)i + (d(P1)/dy)j +(d(P1)/dz)k
// dP2 = (d(P2)/dx)i + (d(P2)/dy)j +(d(P2)/dz)k
// dP3 = (d(P3)/dx)i + (d(P3)/dy)j +(d(P3)/dz)k
// dP4 = (d(P4)/dx)i + (d(P4)/dy)j +(d(P4)/dz)k
// dP5 = (d(P5)/dx)i + (d(P5)/dy)j +(d(P5)/dz)k
// dP6 = (d(P6)/dx)i + (d(P6)/dy)j +(d(P6)/dz)k
//----------------------------------------------------------------
  double  Mag123, Mag456, Dot;
  double  Const1, Const2, Const3;
  double  P1, P2, P3, P4, P5, P6;
  Vector RetVec;

  P1 = P1P2P3[0];  P2 = P1P2P3[1];  P3 = P1P2P3[2];
  P4 = P4P5P6[0];  P5 = P4P5P6[1];  P6 = P4P5P6[2];

  Mag123 = P1P2P3.length();
  Mag456 = P4P5P6.length();
  Dot    = P1P2P3.dot(P4P5P6);

  Const1 =         1.0 / (Mag123*Mag456);
  Const2 = -Dot * (1.0 / (Mag123*Mag456*Mag456*Mag456));
  Const3 = -Dot * (1.0 / (Mag456*Mag123*Mag123*Mag123));

  RetVec = (dP4*P1 + dP1*P4 + dP5*P2 + dP2*P5 + dP6*P3 + dP3*P6) * Const1 +
           (dP4*P4 + dP5*P5 + dP6*P6)                            * Const2 +
           (dP1*P1 + dP2*P2 + dP3*P3)                            * Const3;

  return(RetVec);
}

double ADiheRestraint::GetDihe(Vector& A, Vector& B, Vector& C, Vector& D) {
//-----------------------------------------------------------------
// determine the dihedral angle formed by the points A-B-C-D
//-----------------------------------------------------------------
  Vector CD(D - C);
  Vector CB(B - C);
  Vector BC(C - B);
  Vector BA(A - B);
  Vector CDxCB, BCxBA;
  double  top, bot, cos_u, sin_u, Angle;
  Vector topVec;

  CDxCB = cross(CD,CB);
  BCxBA = cross(BC,BA);

  top = CDxCB.dot(BCxBA);
  bot = CDxCB.length() * BCxBA.length();
  cos_u = top/bot;

  // protect against acos(<-1.0) and acos(>1.0)
  if (cos_u < -1.0) {cos_u = -1.0;}
  if (cos_u >  1.0) {cos_u =  1.0;}

  topVec = cross(CDxCB,BCxBA);
  sin_u = (topVec/bot).dot(CB/CB.length());

  // protect against asin(<-1.0) and asin(>1.0)
  if (sin_u < -1.0) {sin_u = -1.0;}
  if (sin_u >  1.0) {sin_u =  1.0;}

  Angle = atan2(sin_u, cos_u);
  return(Angle);
}

void ADiheRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE) {
// for each center-of-mass
 Vector Force;
 for (int j=0; j<m_NumGroups; j++) {
  // update centers of mass
  UpdateCOMs(CFE);
  // apply restraining force in opposite direction from gradient
  Force = GetGradient(j) * (-1.0);
  DistributeForce(j, Force, CFE);
 }
}

double AFixedDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this fixed dihedral angle restraint.
//---------------------------------------------------------------------------
  return(GetE(m_RefAngle, m_Kf));
}


Vector AFixedDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed dihedral angle restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, m_RefAngle, m_Kf));
}


double ABoundDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound dihedral angle restraint.
//---------------------------------------------------------------------------
  double  E, Dihe, Const;

  Const = m_Kf / (1.0 - cos(m_IntervalAngle));
  Dihe = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  // dihedral angle is between LowerAngle and UpperAngle
  if ( (Dihe>m_LowerAngle) && (Dihe<m_UpperAngle) ) {
    E = 0.0;
  }
  // dihedral angle is between LowerAngle and LowerAngle-IntervalAngle
  else if ( (Dihe<m_LowerAngle) && (Dihe>(m_LowerAngle-m_IntervalAngle)) ) {
    E = GetE(m_LowerAngle, Const);
  }
  // dihedral angle is between UpperAngle and UpperAngle+IntervalAngle
  else if ( (Dihe>m_UpperAngle) && (Dihe<(m_UpperAngle+m_IntervalAngle)) ) {
    E = GetE(m_UpperAngle, Const);
  }
  // dihedral angle is more than UpperAngle or less than LowerAngle
  else {
    E = Const;
  }
  return(E);
}


Vector ABoundDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound dihedral angle restraint.
//---------------------------------------------------------------------------
  Vector Vec;
  double  Dihe, Const;

  Const = m_Kf / (1.0 - cos(m_IntervalAngle));
  Dihe = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  // dihedral angle is between LowerAngle and LowerAngle-IntervalAngle
  if ( (Dihe<m_LowerAngle) && (Dihe>(m_LowerAngle-m_IntervalAngle)) ) {
    Vec = GetGrad(WhichGroup, m_LowerAngle, Const);
  }
  // dihedral angle is between UpperAngle and UpperAngle+IntervalAngle
  else if ( (Dihe>m_UpperAngle) && (Dihe<(m_UpperAngle+m_IntervalAngle)) ) {
    Vec = GetGrad(WhichGroup, m_UpperAngle, Const);
  }
  return(Vec);
}


double AMovingDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this Moving dihedral angle restraint.
//---------------------------------------------------------------------------
  double  RefDihe;

  RefDihe = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetE(RefDihe, m_Kf*m_LambdaKf));
}


Vector AMovingDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving dihedral angle restraint.
//---------------------------------------------------------------------------
  double  RefDihe;
  
  RefDihe = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return(GetGrad(WhichGroup, RefDihe, m_Kf*m_LambdaKf));
}


double AMovingDiheRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving dihedral angle restraint
//---------------------------------------------------------------------------
  double  Dihe;
  double  RefDihe;
//
  Dihe = GetDihe(m_pCOMs[0], m_pCOMs[1], m_pCOMs[2], m_pCOMs[3]);
  RefDihe = m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef);
  return((m_Kf/2)*m_LambdaKf * sin(Dihe-RefDihe) * (m_StartAngle-m_StopAngle));
}
//-------------------------------------------------------------------------
#ifdef FE_RESTRAINT_RMSD_FORTRAN
//------------------------RMSD Classes-------------------------------------
AnRMSDRestraint::AnRMSDRestraint() {
//
  m_NumGroups=2; // orientation group and Moving group
  m_pGroups = new AGroup[m_NumGroups];
  m_pCOMs = new Vector[1] ; // there is only one center-of-mass :  that of the orientation group (1)
  ugrad=NULL;
  qdiffrot=1;
  u[0]=1.;
  u[1]=0.;
  u[2]=0.;
  u[3]=0.;
  u[4]=1.;
  u[5]=0.;
  u[6]=0.;
  u[7]=0.;
  u[8]=1.;
}
//
AnRMSDRestraint::~AnRMSDRestraint() {
//derived class destructor
  if (ugrad !=NULL) delete []ugrad;
}
//
void AnRMSDRestraint::PrintInfo() {
#if defined(_VERBOSE_PMF)
  iout << "RMS Distance = "<<instRMSD<<"    Target = "<<refRMSD<< std::endl << endi;
#else
  iout << instRMSD << "  " << refRMSD << " | ";
#endif
}
//
void AnRMSDRestraint::qDiffRotCheck() {
//
 qdiffrot = 1;
//
 if (m_NumGroups==2) {
   int size1 = m_pGroups[0].GetSize();
   int size2 = m_pGroups[1].GetSize();
   if ( size1==size2 ) { //same size
    int const* inds1 = m_pGroups[0].GetInds();
    int const* inds2 = m_pGroups[1].GetInds();
    qdiffrot=0 ; // assume that the indices are the same, check below
    for (int i=0 ; i < size1 ; i++) {
     if ( inds1[i]!=inds2[i] ) { // compare indices
      qdiffrot = 1; // they are different, exit
      break;
     } //if
    } //for
   } //size
 }//numgroups
 else {
  EarlyExit("AnRMSDRestraint::qDiffRotCheck() : incorrect number of groups");
 }
//
}//fn

void AnRMSDRestraint::ApplyForce(GlobalMasterFreeEnergy &CFE) {
//
 int i, j, k, p, q, aid ;
//
 const _Bool qswapdim=1; // fortran-compatible boolean
 const int norient = m_pGroups[0].GetSize();                            //number of orientation atoms
 const int nforced = m_pGroups[1].GetSize();                            //number of forced atoms
 const int* iatom_f = m_pGroups[1].GetInds();                           //orientation index array
 const int* iatom_o = NULL;//forced index array
//
 double* rcurrent_f= (double*) malloc(3*norient*sizeof(double))    ;    // orientation coordinates
 double* forces_f  = (double*) malloc(3*norient*sizeof(double))    ;    // forces on orientation coordinates
 double* rcurrent_o=NULL;
 double* forces_o=NULL;
//
 const double* forcedWeights = m_pGroups[1].GetWeights();
 const double* orientWeights = NULL;
// reference coordinates
 const double* rtarget_f= m_pGroups[0].GetCoords();                            //orientation coordinates of the target structure
 const double* rtarget_o= NULL;
 double* rtarget_rot_f=NULL;
//
//
 if (qdiffrot) {
  iatom_o = m_pGroups[0].GetInds() ; //forced index array
  rcurrent_o = (double*) malloc(3*norient*sizeof(double)) ; // forced coordinates
  forces_o = (double*) malloc(3*norient*sizeof(double)) ; // forced coordinates
  orientWeights = m_pGroups[0].GetWeights() ;
  rtarget_o= m_pGroups[0].GetCoords();  //forced coordinates of the target structure
 } else { // point all orientation data to the forced data
  iatom_o = iatom_f;
  orientWeights=forcedWeights;
  rcurrent_o=rcurrent_f;
  rtarget_o=rtarget_f;
  forces_o=forces_f;
 }
//
//      load coordinates
 for (i=0, j=0 ; i<nforced ; i++) {// forced
  aid=iatom_f[i];
  Vector pos = CFE.positions[aid];
  rcurrent_f[j++]=pos.x;          rcurrent_f[j++]=pos.y;         rcurrent_f[j++]=pos.z;
 };
//
 if (qorient && qdiffrot) { // orientation
  for (i=0, j=0 ; i<norient ; i++) {
   aid=iatom_o[i];
   Vector pos = CFE.positions[aid];
   rcurrent_o[j++]=pos.x;          rcurrent_o[j++]=pos.y;         rcurrent_o[j++]=pos.z;
  };
 }
//
 if (qorient) {
  double* COM = (double*) com(rcurrent_o,orientWeights,norient,qswapdim);   // pointer returned is of type void, cast to double
  for (i=0;i<3*norient;){rcurrent_o[i++]-=COM[0];rcurrent_o[i++]-=COM[1];rcurrent_o[i++]-=COM[2]; }
  if (qdiffrot) { // also subtract COM from forced coordinates
   for (i=0;i<3*nforced;){rcurrent_f[i++]-=COM[0];rcurrent_f[i++]-=COM[1];rcurrent_f[i++]-=COM[2]; }
  }
  free(COM);
 }
//
// compute orientation and derivatives
//
 if (qorient) {
  if (qdiffrot) {
    double* ugrad=(double *) malloc(27*norient*sizeof(double)); //allocate gradient array
    RMSBestFitGrad(rtarget_o,rcurrent_o,orientWeights,norient,u,ugrad,1,norient,qswapdim); // compute matrix and gradients u: targ -> curr
    rtarget_rot_f=(double*) matmul(u, rtarget_f, 3, 3, nforced);
//
//    compute forces on orientation atoms
//    r1 = matmul(transpose(v%rcurrent_f-v%rtarget_rot_f), v%forcedWeights // shorthand way
    double r1[3] ; r1[0]=0.; r1[1]=0.; r1[2]=0.;
    double r2[3] ;
//
    for (i=0, k=0; i<nforced ; i++) {
     double w = forcedWeights[i];
     r1[0]=r1[0] + w *( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++;
     r1[1]=r1[1] + w *( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++;
     r1[2]=r1[2] + w *( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++;
    }
//     compute COM forces on orientation atoms:
    for (i=0, j=0; i<norient; i++) {
     double w = orientWeights[i];
     forces_o[j++] = -r1[0] * w;
     forces_o[j++] = -r1[1] * w;
     forces_o[j++] = -r1[2] * w;
    }
//     compute quadratic using ugrad:
//
    const double* r3=rtarget_f;
    for (j=0, k=0 ; j<nforced ; j++, r3+=3) {
//
      double w = forcedWeights[j];
//
//    also add forces on the forcing atoms here:
      forces_f[k] =  w * ( rcurrent_f[k]-rtarget_rot_f[k] ); //x
      r2[0] = w * rcurrent_f[k];
      k++;
//
      forces_f[k] =  w * ( rcurrent_f[k]-rtarget_rot_f[k] ); //y
      r2[1] = w * rcurrent_f[k]; 
      k++;
//
      forces_f[k] =  w * ( rcurrent_f[k]-rtarget_rot_f[k] ); //z
      r2[2] = w * rcurrent_f[k]; 
      k++;
//
      double* f=forces_o;
      double* ugradi=ugrad;
      for ( i=0 ; i<norient ; i++ , f+=3, ugradi+=27) {
//
       double * ugradx = ugradi ; 
       double * ugrady = ugradx + 9; 
       double * ugradz = ugrady + 9;

       int ind;
       for (q=0, ind=0 ; q<3; q++) {
        for (p=0; p<3; p++, ind++) { // p varies the fastest ;  ind goes from 0 to 9
         double d=r2[p] * r3[q];
//
         f[0]-= d * ugradx[ind];
         f[1]-= d * ugrady[ind];
         f[2]-= d * ugradz[ind];
        }; // q
       }; // p
      }; // i (orientation atoms)
    }; // j (forced atoms)
    free(ugrad);
//=======================================================================================
  } else { // qorient but not qdiffrot -- regular TMD as we know it
     RMSBestFit(rtarget_o,rcurrent_o,orientWeights,norient,u,qswapdim); // compute matrix only, no gradients u: targ -> curr
//   transform target structure to overlap with current
     rtarget_rot_f=(double*) matmul(u, rtarget_f, 3, 3, nforced); // will deallocate below
  } //qdiffrot
 } else { // not qorient
     rtarget_rot_f = (double*) rtarget_f;
 } // qorient
//
//   add forces on the forcing atoms:
//
 for (j=0, k=0; j<nforced; j++) {
  double w = forcedWeights[j];
  forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++ ;
 } // for
//  compute RMSD
 instRMSD = rmsd(rcurrent_f, rtarget_rot_f, forcedWeights, nforced, qswapdim);
//  if ( refRMSD == -1.0 ) refRMSD = instRMSD; // compute target RMS from first coordinate
//  force prefactor
 double pref = ComputeForcePrefactor(); // overloaded in derived classes
//
// send forces to the CFE object
// forcing atoms
 double* f=forces_f;
 for (i=0; i<nforced; i++, f+=3) {
  int aid = iatom_f[i];
  Vector force(f[0], f[1], f[2]);
  force*=(-pref);
//
  DebugM(2,"Applying force "<<force<<" to forced atom "<<aid<<"\n");
  CFE.modifyForcedAtoms().add(aid);
  CFE.modifyAppliedForces().add( force );
 }
// orientation atoms
 if (qdiffrot) {
  f=forces_o;
  for (i=0; i<norient; i++, f+=3) {
   int aid = iatom_o[i];
   Vector force(f[0], f[1], f[2]);
   force*=(-pref);
//
   DebugM(2,"Applying force "<<force<<" to orientation atom "<<aid<<"\n");
   CFE.modifyForcedAtoms().add(aid);
   CFE.modifyAppliedForces().add( force );
  }
 }
//
// ----------------------v clean up memory
//
 if ( qdiffrot || qorient ) free(rtarget_rot_f);
 free(forces_f);
 free(rcurrent_f);
 if(qdiffrot) {
  free(forces_o);
  free(rcurrent_o);
 }
// ----------------------^ clean up memory
} //ApplyForce
#endif //RMSD restraints
