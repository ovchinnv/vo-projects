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
//
#include "NamdTypes.h"
#include "GlobalMaster.h"
#include "GlobalMasterFreeEnergy.h"
#include "Molecule.h"
//
#include "Debug.h"

// initialize static member variables
// (lambda is the same for all Moving restraints)
double  ARestraint::LambdaKf  = 1.0;
double  ARestraint::LambdaRef = 0.0;

ARestraint::ARestraint() {
//-----------------------------------------------------------------
// constructor for base class
//-----------------------------------------------------------------
  Groups = NULL;
  COMs = NULL;
  NumGroups = 0;
}


ARestraint::~ARestraint() {
//-----------------------------------------------------------------
// free space that may have been allocated for Groups and COM's
//-----------------------------------------------------------------
  if (Groups != NULL) {
    ASSERT(COMs != NULL);
    delete []Groups;
    delete []COMs;
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
  ASSERT( (GroupIndex>=0) && (GroupIndex<NumGroups) );
  Groups[GroupIndex] = Group;
}


void ARestraint::SetGroups(AGroup* groups, int n) {
//-----------------------------------------------------------------
// set groups
//-----------------------------------------------------------------
  ASSERT(NumGroups >= n);
  for (int i=0; i < n ; i++) { Groups[i] = groups[i] ; }
}

void ARestraint::DistributeForce(int WhichGroup, Force force,
                                 GlobalMasterFreeEnergy& CFE) {
//----------------------------------------------------------------------
// Distribute Force among the group of atoms specified by WhichGroup
//
// note:  Groups points to an array of Groups
//        Groups[WhichGroup] references one of the Groups
//        Groups[WhichGroup][i] returns an AtomID from the Group
//        (operator[] is defined to return an item from the Group)
//----------------------------------------------------------------------
  ASSERT( (WhichGroup>=0) && (WhichGroup<NumGroups) );

  int NumAtoms = Groups[WhichGroup].GetSize();
  double const* Weights=Groups[WhichGroup].GetWeights();
  double weight;

  // distribute Force according to mass of each atom in the group
  for (int i=0; i<NumAtoms; i++) {
    int aid = Groups[WhichGroup][i];

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
// note:  Groups points to an array of Groups
//        Groups[i] references one of the Groups
//        Groups[i][j] returns an AtomID from the Group
//        (operator[] is defined to return an item from the Group)
//-----------------------------------------------------------------
  int      i, j, aid;
  Position pos ;
  double const* Weights;
  double Weight;
//
  ASSERT(NumGroups > 0);
// fetch positions
  // for each group of atoms
  for (i=0; i<NumGroups; i++) {
    Weights=Groups[i].GetWeights();
    Vector COM(0.,0.,0) ;
    // for each atom in the group
    for (j=0; j<Groups[i].GetSize(); j++) {
      aid = Groups[i][j]; // get atom ID
      Vector pos=CFE.positions[aid]; // fetch from map
      Weight = Weights[j];
      COM += pos * Weight;
    }
    COMs[i] = COM;
  }
}

void APosRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE, bool qfdtest, double dh) {
// for each center-of-mass
 Force force; // force vector
//
 for (int j=0; j<NumGroups; j++) {
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
  NumGroups = 1;
  Groups = new AGroup[NumGroups];
  COMs = new Vector[NumGroups];
}

void APosRestraint::PrintInfo() {
  iout << " Position= " << COMs[0];
  iout << " Target= "   << GetPosTarget();
  iout << " Distance= " << GetDistance();
}

double APosRestraint::GetE(Vector RefPos, double LambdaKf) {
//--------------------------------------------------------------------
//     E = (Kf/2) * (|ri - rref|)**2
// where |ri - rref| is the distance between a) the center-of-mass
// of the restrained atoms and b) the reference position.
// Note:  COM is calculated before this routine is called.
//--------------------------------------------------------------------
  return( ((Kf*LambdaKf)/2.0) * (COMs[0] - RefPos).length2() );
}


Vector APosRestraint::GetGrad(int /* WhichGroup */, Vector RefPos, double LambdaKf) {
//-------------------------------------------------------------------------
//     E = (Kf/2) * (|ri - rref|)**2
// return:  grad(E)
//-------------------------------------------------------------------------
  return ( COMs[0] - RefPos ) * Kf * LambdaKf ;
}

double AFixedPosRestraint::GetEnergy() {return(GetE(RefPos));}
Vector AFixedPosRestraint::GetGradient(int WhichGroup) { return(GetGrad(WhichGroup, RefPos));}


double ABoundPosRestraint::GetEnergy() {
  double  E, Dist, Diff;

  E = 0.0;
  Dist = (COMs[0] - (RefPos)).length();
  if (((Bound==kUpper) && (Dist>RefDist)) ||
      ((Bound==kLower) && (Dist<RefDist))) {
    Diff = Dist - RefDist;
    E = (Kf/2.0) * (Diff*Diff);
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
  Dist = (COMs[0]-RefPos).length();
  if (((Bound==kUpper) && (Dist>RefDist)) ||
      ((Bound==kLower) && (Dist<RefDist))) {
    Vec = (COMs[0]- RefPos) * Kf * (Dist - RefDist) / __MAX(Dist, kALittle); //  prevent divide overflow
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

  RefPos = StopPos*LambdaRef + StartPos*(1.0-LambdaRef);
  return(GetE(RefPos, LambdaKf));
}


Vector AMovingPosRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving position restraint.
//
// rref = lambda*r1 + (1-lambda)*r0.
// where r0 is the starting position and r1 is the final position
//---------------------------------------------------------------------------
  Vector  RefPos;

  RefPos = StopPos*LambdaRef + StartPos*(1.0-LambdaRef);
  return(GetGrad(WhichGroup, RefPos, LambdaKf));
}


double AMovingPosRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving position restraint
//---------------------------------------------------------------------------
  Vector  RefPos;
//  double   T1, T2, T3;

  RefPos = StopPos*LambdaRef + StartPos*(1.0-LambdaRef);
//  T1 = (COMs[0][0] - RefPos[0]) * (StartPos[0] - StopPos[0]);
//  T2 = (COMs[0][1] - RefPos[1]) * (StartPos[1] - StopPos[1]);
//  T3 = (COMs[0][2] - RefPos[2]) * (StartPos[2] - StopPos[2]);
//  return( Kf * LambdaKf * (T1+T2+T3) );
    return Kf * LambdaKf * ( (COMs[0] - RefPos) * (StartPos - StopPos) );// replaced above by a scalar (dot) product
}


ADistRestraint::ADistRestraint() {
//-----------------------------------------------------------------------
// each ADistRestraint restrains the distance between 2 groups of atoms
//-----------------------------------------------------------------------
  NumGroups = 2;
  Groups = new AGroup[NumGroups];
  COMs = new Vector[NumGroups];
}


void ADistRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the distance for this distance restraint
//--------------------------------------------------------------------
  iout << " Distance= "<< (COMs[0]-COMs[1]).length();
  iout << " Target= "<< GetDistTarget();
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

  Dist = (COMs[0]-COMs[1]).length();
  Diff = Dist - RefDist;
  return( ((Kf*LambdaKf)/2.0) * (Diff*Diff) );
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
//        COMS[0 & 1] reference the COM's of each group of atoms
//---------------------------------------------------------------------------
  ASSERT( (WhichGroup==0) || (WhichGroup==1) );
//
  Vector Vec = COMs[1] - COMs[0];
  double Dist = Vec.length();
  return  (WhichGroup * 2. - 1.) * (  Vec * Kf * LambdaKf * (Dist - RefDist)/__MAX(Dist, kALittle) );// WG=0 => prefix=-1 ; WG=1=> prefix=+1
}

double AFixedDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this fixed distance restraint.
//---------------------------------------------------------------------------
  return(GetE(RefDist));
}


Vector AFixedDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed distance restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, RefDist));
}

double ABoundDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound distance restraint.
//---------------------------------------------------------------------------
  double Dist, E;
//
  E = 0.0;
  Dist = (COMs[0] - COMs[1]).length();
  if (((Bound==kUpper) && (Dist>RefDist)) ||
      ((Bound==kLower) && (Dist<RefDist))) {
    E = GetE(RefDist);
  }
  return(E);
}

Vector ABoundDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound distance restraint.
//---------------------------------------------------------------------------
  double  Dist;
  Vector Vec;

  Dist = (COMs[0] - COMs[1]).length();
  if (((Bound==kUpper) && (Dist>RefDist)) ||
      ((Bound==kLower) && (Dist<RefDist))) {
    Vec = GetGrad(WhichGroup, RefDist);
  }
  return(Vec);
}


double AMovingDistRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this Moving distance restraint.
//---------------------------------------------------------------------------
  double  RefDist;

  RefDist = StopDist*LambdaRef + StartDist*(1.0-LambdaRef);
  return(GetE(RefDist, LambdaKf));
}


Vector AMovingDistRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving distance restraint.
//---------------------------------------------------------------------------
  double  RefDist;
  
  RefDist = StopDist*LambdaRef + StartDist*(1.0-LambdaRef);
  return(GetGrad(WhichGroup, RefDist, LambdaKf));
}


double AMovingDistRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving distance restraint
//---------------------------------------------------------------------------
  double  Dist;
  double  RefDist;
//
  Dist = (COMs[0] - COMs[1]).length();
  RefDist = StopDist*LambdaRef + StartDist*(1.0-LambdaRef);
  return( Kf * LambdaKf * (Dist-RefDist)*(StartDist-StopDist) );
}

void ADistRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE, bool qfdtest, double dh) {
// for each center-of-mass
 Vector Force;
 for (int j=0; j<NumGroups; j++) {
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
  NumGroups = 3;
  Groups = new AGroup[NumGroups];
  COMs = new Vector[NumGroups];
}

void AnAngleRestraint::PrintInfo() {
  double  Angle = GetAngle(COMs[0], COMs[1], COMs[2]) * (180/kPi);
  iout << " Angle(deg)= "<<Angle;
  iout << " Target= "<<GetAngleTarget()*180/kPi;
}


double AnAngleRestraint::GetE(double RefAngle, double LambdaKf) {
//     E = (Kf/2) * (Theta-ThetaRef)**2
//
// where Theta is the angle between 3 centers-of-mass of restrained atoms,
// COMs[0] -- COMs[1] -- COMs[2].
// ThetaRef is the reference angle.
//---------------------------------------------------------------------------
  double Angle, Diff;

  Angle = GetAngle(COMs[0], COMs[1], COMs[2]);
  Diff = Angle - RefAngle;
  return( ((Kf*LambdaKf)/2.0) * (Diff*Diff) );
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
//     E = (Kf/2) * (Theta-ThetaRef)**2
// return:  grad(E)
//---------------------------------------------------------------------------
  Vector A, B, C;
  double  Angle;
  double  a, b, c;
  double  u;
  double  Const1, Const2, Const3, Const4, Const5;
  Vector Vec1, Vec2;

  ASSERT( (WhichGroup==0) || (WhichGroup==1) || (WhichGroup==2) );
  A = COMs[0];
  B = COMs[1];
  C = COMs[2];

  a = (B-C).length();
  b = (A-C).length();
  c = (A-B).length();

  u = (a*a + c*c - b*b) / (2.0*a*c);

  // protect against acos(<-1.0), acos(>1.0), sqrt(<0), and divide by 0
  if (u < -1.0)       {u = -1.0;}
  if (u > kAlmostOne) {u = kAlmostOne;}
  Angle = acos(u);

  Const1 = -1.0 / sqrt(1.0 - u*u);
  Const2 = Kf * LambdaKf * (Angle-RefAngle) * Const1;

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

void AnAngleRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE, bool qfdtest, double dh) {
// for each center-of-mass
 Vector Force;
 for (int j=0; j<NumGroups; j++) {
  // update centers of mass
  UpdateCOMs(CFE);
  // apply restraining force in opposite direction from gradient
  Force = GetGradient(j) * (-1.0);
  DistributeForce(j, Force, CFE);
 }
}

double AFixedAngleRestraint::GetEnergy() {return(GetE(RefAngle));}
Vector AFixedAngleRestraint::GetGradient(int WhichGroup) {return(GetGrad(WhichGroup, RefAngle));}

double ABoundAngleRestraint::GetEnergy() {
  double  E, Angle;

  E = 0.0;
  Angle = GetAngle(COMs[0], COMs[1], COMs[2]);
  if (((Bound==kUpper) && (Angle>RefAngle)) ||
      ((Bound==kLower) && (Angle<RefAngle))) {
    E = GetE(RefAngle);
  }
  return(E);
}

Vector ABoundAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this bound angle restraint
//---------------------------------------------------------------------------
  double  Angle;
  Vector Vec;

  Angle = GetAngle(COMs[0], COMs[1], COMs[2]);
  if (((Bound==kUpper) && (Angle>RefAngle)) ||
      ((Bound==kLower) && (Angle<RefAngle))) {
    Vec = GetGrad(WhichGroup, RefAngle);
  }
  return(Vec);
}


double AMovingAngleRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this Moving angle restraint.
//---------------------------------------------------------------------------
  double  RefAngle;

  RefAngle = StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef);
  return(GetE(RefAngle, LambdaKf));
}


Vector AMovingAngleRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving angle restraint.
//---------------------------------------------------------------------------
  double  RefAngle;

  RefAngle = StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef);
  return(GetGrad(WhichGroup, RefAngle, LambdaKf));
}


double AMovingAngleRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving angle restraint
//---------------------------------------------------------------------------
  double  Angle;
  double  RefAngle;

  Angle = GetAngle(COMs[0], COMs[1], COMs[2]);
  RefAngle = StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef);
  return( Kf * LambdaKf * (Angle-RefAngle)*(StartAngle-StopAngle) );
}


ADiheRestraint::ADiheRestraint() {
//----------------------------------------------------------------------------
// each ADiheRestraint restrains the dihedral angle between 4 groups of atoms
//----------------------------------------------------------------------------
  NumGroups = 4;
  Groups = new AGroup[NumGroups];
  COMs = new Vector[NumGroups];
}


void ADiheRestraint::PrintInfo() {
//--------------------------------------------------------------------
// print the dihedral angle for this dihedral restraint
//--------------------------------------------------------------------
  double  Dihedral;
  Dihedral = GetDihe(COMs[0], COMs[1], COMs[2], COMs[3]) * (180/kPi);
//
  iout << " Dihedral(deg)= "<<Dihedral;
  iout << " Target= "<<GetDiheTarget1() * (180/kPi);
  if (TwoTargets()) iout <<" to "<<GetDiheTarget2() * (180/kPi);
}

double ADiheRestraint::GetE(double RefAngle, double Const) {
//    E = (E0/2) * (1 - cos(Chi - ChiRef))
//
// where Chi is the dihedral angle between 4 centers-of-mass of restrained atoms,
// COMs[0] -- COMs[1] -- COMs[2] -- COMs[3].
// ChiRef is the reference angle.
//---------------------------------------------------------------------------
  double  Angle;

  Angle = GetDihe(COMs[0], COMs[1], COMs[2], COMs[3]);
  return( (Const/2.0) * (1.0 - cos(Angle-RefAngle)) );
}

Vector ADiheRestraint::GetGrad(int WhichGroup,
                                double RefAngle, double Const) {
//    E = (E0/2) * (1 - cos(Chi - ChiRef))
//
// return:  grad(E)
//---------------------------------------------------------------------------
  Vector A, B, C, D;

  ASSERT((WhichGroup==0)||(WhichGroup==1)||(WhichGroup==2)||(WhichGroup==3));

  if ((WhichGroup==0) || (WhichGroup==1)) {
    A = COMs[0];
    B = COMs[1];
    C = COMs[2];
    D = COMs[3];
  }
  // re-state the problem so the gradient is solved for either atoms 0 or 1
  else {
    A = COMs[3];
    B = COMs[2];
    C = COMs[1];
    D = COMs[0];
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

void ADiheRestraint::ApplyForce(GlobalMasterFreeEnergy& CFE, bool qfdtest, double dh) {
// for each center-of-mass
 Vector Force;
 for (int j=0; j<NumGroups; j++) {
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
  return(GetE(RefAngle, Kf));
}


Vector AFixedDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this fixed dihedral angle restraint.
//---------------------------------------------------------------------------
  return(GetGrad(WhichGroup, RefAngle, Kf));
}


double ABoundDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this bound dihedral angle restraint.
//---------------------------------------------------------------------------
  double  E, Dihe, Const;

  Const = Kf / (1.0 - cos(IntervalAngle));
  Dihe = GetDihe(COMs[0], COMs[1], COMs[2], COMs[3]);
  // dihedral angle is between LowerAngle and UpperAngle
  if ( (Dihe>LowerAngle) && (Dihe<UpperAngle) ) {
    E = 0.0;
  }
  // dihedral angle is between LowerAngle and LowerAngle-IntervalAngle
  else if ( (Dihe<LowerAngle) && (Dihe>(LowerAngle-IntervalAngle)) ) {
    E = GetE(LowerAngle, Const);
  }
  // dihedral angle is between UpperAngle and UpperAngle+IntervalAngle
  else if ( (Dihe>UpperAngle) && (Dihe<(UpperAngle+IntervalAngle)) ) {
    E = GetE(UpperAngle, Const);
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

  Const = Kf / (1.0 - cos(IntervalAngle));
  Dihe = GetDihe(COMs[0], COMs[1], COMs[2], COMs[3]);
  // dihedral angle is between LowerAngle and LowerAngle-IntervalAngle
  if ( (Dihe<LowerAngle) && (Dihe>(LowerAngle-IntervalAngle)) ) {
    Vec = GetGrad(WhichGroup, LowerAngle, Const);
  }
  // dihedral angle is between UpperAngle and UpperAngle+IntervalAngle
  else if ( (Dihe>UpperAngle) && (Dihe<(UpperAngle+IntervalAngle)) ) {
    Vec = GetGrad(WhichGroup, UpperAngle, Const);
  }
  return(Vec);
}


double AMovingDiheRestraint::GetEnergy() {
//---------------------------------------------------------------------------
// return the Energy for this Moving dihedral angle restraint.
//---------------------------------------------------------------------------
  double  RefDihe;

  RefDihe = StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef);
  return(GetE(RefDihe, Kf*LambdaKf));
}


Vector AMovingDiheRestraint::GetGradient(int WhichGroup) {
//---------------------------------------------------------------------------
// return the gradient for this Moving dihedral angle restraint.
//---------------------------------------------------------------------------
  double  RefDihe;
  
  RefDihe = StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef);
  return(GetGrad(WhichGroup, RefDihe, Kf*LambdaKf));
}


double AMovingDiheRestraint::Get_dU_dLambda() {
//---------------------------------------------------------------------------
// return dU/dLambda for this Moving dihedral angle restraint
//---------------------------------------------------------------------------
  double  Dihe;
  double  RefDihe;
//
  Dihe = GetDihe(COMs[0], COMs[1], COMs[2], COMs[3]);
  RefDihe = StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef);
  return((Kf/2)*LambdaKf * sin(Dihe-RefDihe) * (StartAngle-StopAngle));
}
//-------------------------------------------------------------------------
#ifdef FE_RESTRAINT_RMSD_FORTRAN
//------------------------RMSD Classes-------------------------------------
AnRMSDRestraint::AnRMSDRestraint() {
//
  NumGroups=2; // orientation group and Moving group
  Groups = new AGroup[NumGroups];
  COMs = new Vector[1] ; // there is only one center-of-mass :  that of the orientation group (1)
  ugrad=NULL;
  qdiffrot=1;
// qorient, qtrans set on input
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
  iout << " RMS Distance= "<<instRMSD;
  iout << " Target= "<<refRMSD;
}
//
void AnRMSDRestraint::qDiffRotCheck() {
//
 qdiffrot = 1;
//
 if (NumGroups==2) {
  int norient = Groups[0].GetSize();
  int nforced = Groups[1].GetSize();
  int const *iatom_o=Groups[0].GetInds();
  int const *iatom_f=Groups[1].GetInds();
  const double* orientWeights = Groups[0].GetWeights() ;
  const double* forcedWeights = Groups[1].GetWeights() ;
//
  if (qorient && !qtrans) { // sanity check
    EarlyExit("AnRMSDRestraint::qDiffRotCheck() : Orientation requires translation to be enabled. Aborting.");
  }
//
  if (nforced==0) {
   EarlyExit("AnRMSDRestraint::qDiffRotCheck() : No forcing atoms found (nothing to to). Aborting.");
  } else if (qorient) {
   if (norient == 0) {
    iout << " Warning (AnRMSDRestraint::qDiffRotCheck) : No orientation atoms found. Orientation is off." << std::endl<<endi;
    qorient=0 ;
    qtrans=0;
    Groups[0].Clear() ; // orientation group
   } else if (norient < 2) {
    iout << " Warning (AnRMSDRestraint::qDiffRotCheck) : Fewer than two (2) orientation atoms found. Rotation is off." << std::endl<<endi;
    qorient=0 ;
   }
  } else if (qtrans && (norient==0)) {
    iout << " Warning (AnRMSDRestraint::qDiffRotCheck) : No orientation atoms found. Translation is off." << std::endl<<endi;
    qorient=0;
    qtrans=0;
    Groups[0].Clear() ; // orientation group
  }
//
  if ( (qtrans) && (norient==nforced)) {
   qdiffrot=0 ; // assume that the indices are the same, check below
   std::map<int, double > imap; // maps aid to o/f atom weights
// first store forced atom indices
   for (int i=0; i<nforced; i++) {
    int aid = iatom_f[i];
    imap[aid]=forcedWeights[i]; //store forced atom weight
   }
// now loop over orientation atom indices
   for (int i=0; i<norient; i++) {
    int aid = iatom_o[i];
    std::map<int, double > ::iterator ind = imap.find(aid) ; // check whether id already in the map
    if (ind == imap.end()) {// not found
     qdiffrot=1; break;
    } else if (imap[aid] != orientWeights[i]) {
     qdiffrot=1; break;
    }
   } //i
  }// qorient && norient = nforced
// check for the special case with two orientation atoms, make sure forced atoms are contained within the orientation group
  if (norient==2 && qdiffrot && qorient) { // OK for translation with no rotation
   for (int i=0 ; i < nforced ; i++) {
    if  ( !( (iatom_f[i] == iatom_o[0] ) || (iatom_f[i] == iatom_o[1] ) ) ) { // compare indices
     EarlyExit("AnRMSDRestraint::qDiffRotCheck() : With two orientation atoms, \
the orientation group must contain the forcing group. Aborting.");
    break;
    }
   }
  }
 }//numgroups
 else {
  EarlyExit("AnRMSDRestraint::qDiffRotCheck() : RMSD restraints require two groups");
 }
DebugM(9,"AnRMSDREstraint::qDiffRotCheck() : q_orient="<<qorient<<", q_translate="<<qtrans<<"\n");
}

void AnRMSDRestraint::SetGroups(AGroup* groups, int n) {
//-----------------------------------------------------------------
// set two groups of atoms
//-----------------------------------------------------------------
  ASSERT(n >= 2);
  Groups[0] = groups[0];
  Groups[1] = groups[1];
  qDiffRotCheck();
//
// remove COM from target atoms if needed
//
 if (qtrans) {
  bool qswapdim=1;
  const int norient = Groups[0].GetSize();                              //number of orientation atoms
  const int nforced = Groups[1].GetSize();                             //number of forced atoms
  double* rtarget_o= Groups[0].GetCoords();
  const double* orientWeights = Groups[0].GetWeights() ;
  double* COM = (double*) com(rtarget_o,orientWeights,norient,qswapdim);
// remove COM from forcing group
  double* rtarget_f= Groups[1].GetCoords();
  for (int i=0;i<3*nforced;){rtarget_f[i++]-=COM[0];rtarget_f[i++]-=COM[1];rtarget_f[i++]-=COM[2]; }
// from orientation group (optional)
  for (int i=0;i<3*norient;){rtarget_o[i++]-=COM[0];rtarget_o[i++]-=COM[1];rtarget_o[i++]-=COM[2]; }
DebugM(9, "Add RMSD groups : Removed COM from target atoms: "<<COM[0]<<" "<<COM[1]<<" "<<COM[2]<<"\n");
  free(COM);
 }
}
//

void AnRMSDRestraint::ApplyForce(GlobalMasterFreeEnergy &CFE, bool fdtest, double dh) {
//
 int i, j, k, p, q, aid ;
//
 const _Bool qswapdim=1; // fortran-compatible boolean
 const int norient = Groups[0].GetSize();                            //number of orientation atoms
 const int nforced = Groups[1].GetSize();                            //number of forced atoms
 const int* iatom_f = Groups[1].GetInds();                           //forced index array
 const int* iatom_o = NULL;//orientation index array
//
 double* rcurrent_f= (double*) malloc(3*nforced*sizeof(double)); __MEMOK(rcurrent_f)   ;    // forced coordinates
 double* forces_f  = (double*) malloc(3*nforced*sizeof(double)); __MEMOK(forces_f)     ;    // forces on forced coordinates
 double* rcurrent_o=NULL;
 double* forces_o=NULL;
//
 const double* forcedWeights = Groups[1].GetWeights();
 const double* orientWeights = NULL;
// reference coordinates
 const double* rtarget_f= Groups[1].GetCoords();                            //forced coordinates of the target structure
 const double* rtarget_o= NULL;
 double* rtarget_rot_f=NULL;
 double* COM=NULL;
//
DebugM(9, "RMSD Force calculation : qdiffrot, qorient, qtrans: "<<qdiffrot<<" "<<qorient<<" "<<qtrans<<"\n");
DebugM(9, "RMSD Force calculation : nforced, norient: "<<nforced<<" "<<norient<<"\n");
//
 if (qdiffrot && qtrans) {
  iatom_o = Groups[0].GetInds() ; //orientation index array
  rcurrent_o = (double*) malloc(3*norient*sizeof(double)) ;__MEMOK(rcurrent_o) ;
  forces_o = (double*) malloc(3*norient*sizeof(double)) ; __MEMOK(forces_o)  ;
  orientWeights = Groups[0].GetWeights() ;
  rtarget_o= Groups[0].GetCoords();  //orientation coordinates of the target structure
 } else { // point all orientation data to the forced data
  iatom_o = iatom_f;
  orientWeights=forcedWeights;
  rcurrent_o=rcurrent_f;
  rtarget_o=rtarget_f;
  forces_o=forces_f;
 }
//
//      load coordinates
DebugM(9,"RMSD force : loading forced atom coordinatess\n");
 for (i=0, j=0 ; i<nforced ; i++) {// forced
  aid=iatom_f[i];
  Vector pos = CFE.positions[aid];
  rcurrent_f[j++]=pos.x;          rcurrent_f[j++]=pos.y;         rcurrent_f[j++]=pos.z;
 };
//
 if (qtrans && qdiffrot) { // orientation
DebugM(9,"RMSD force :  loading orientation atom coordinates\n");
  for (i=0, j=0 ; i<norient ; i++) {
   aid=iatom_o[i];
   Vector pos = CFE.positions[aid];
   rcurrent_o[j++]=pos.x;          rcurrent_o[j++]=pos.y;         rcurrent_o[j++]=pos.z;
  };
 }
//
// aa -- debug v
// DebugM(1,rcurrent_f<< " "<< rcurrent_o <<" "<<norient<<" "<<nforced<<" "<<qdiffrot<<" "<<qorient<<"\n"); 
// for (i=0;i<3*nforced;) { DebugM(1,"RF : "<<rcurrent_f[i] << " : "<<rcurrent_f[i+1] << " : "<<rcurrent_f[i+2] <<"\n"); i+=3; }
// for (i=0;i<3*norient;) { DebugM(1,"RO : "<<rcurrent_o[i] << " : "<<rcurrent_o[i+1] << " : "<<rcurrent_o[i+2] <<"\n"); i+=3; }
// double* COM = (double*) com(rtarget_o,orientWeights,norient,qswapdim);
// DebugM(1, "RMSD Force calculation : TCOM: "<<COM[0]<<" "<<COM[1]<<" "<<COM[2]<<"\n");free(COM);
// aa Debug ^
//
 if (qtrans) { // subtract COM from forced coordinates
  COM = (double*) com(rcurrent_o,orientWeights,norient,qswapdim);   // pointer returned is of type void, cast to double
  for (i=0;i<3*nforced;){rcurrent_f[i++]-=COM[0];rcurrent_f[i++]-=COM[1];rcurrent_f[i++]-=COM[2]; }
//  if (qdiffrot) {
//   for (i=0;i<3*norient;){rcurrent_o[i++]-=COM[0];rcurrent_o[i++]-=COM[1];rcurrent_o[i++]-=COM[2]; } // not needed
//  }
DebugM(9, "RMSD Force calculation : COM: "<<COM[0]<<" "<<COM[1]<<" "<<COM[2]<<"\n");
 }
//
// compute orientation and derivatives
//
 if (qorient) {
  if (qdiffrot) {
    double* ugrad=(double *) malloc(27*norient*sizeof(double)); __MEMOK(ugrad)  ;//allocate gradient array
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
     r1[0]+= w *( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++;
     r1[1]+= w *( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++;
     r1[2]+= w *( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++;
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
      for (j=0, k=0; j<nforced; j++) {//   add forces on the forcing atoms:
       double w = forcedWeights[j];
       forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++ ;
       forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++ ;
       forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++ ;
      } // for
  } //qdiffrot
 } else { // not qorient (qtrans=1 covered here : in that case, COM has been subtracted)
    rtarget_rot_f=(double*)rtarget_f ; // point to static target structure w/o rotation ; need dirty cast b/c rtarget_f const
    if (qdiffrot && qtrans) { // note additional force on COM atoms
     double r1[3] ; r1[0]=0.; r1[1]=0.; r1[2]=0.;
     for (j=0, k=0; j<nforced; j++) {//   add forces on the forcing atoms:
      double w = forcedWeights[j];
      forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; r1[0]+=forces_f[k] ; k++ ;
      forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; r1[1]+=forces_f[k] ; k++ ;
      forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; r1[2]+=forces_f[k] ; k++ ;
     } // j-for
//     compute COM forces on orientation atoms:
     for (i=0, j=0; i<norient; i++) {
      double w = orientWeights[i];
      forces_o[j++] = -r1[0] * w;
      forces_o[j++] = -r1[1] * w;
      forces_o[j++] = -r1[2] * w;
     }
    } else { // not qdiffrot & qtrans
     for (j=0, k=0; j<nforced; j++) {//   add forces on the forcing atoms:
      double w = forcedWeights[j];
      forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++ ;
      forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++ ;
      forces_f[k] = w * ( rcurrent_f[k] - rtarget_rot_f[k] ) ; k++ ;
     } // for
    } //qdiffrot
 } // qorient
//
//======== done with gradient computation ^ ; note that gradients may be in two parts corresponding to orientation & forcing terms
//  compute RMSD
 instRMSD = rmsd(rcurrent_f, rtarget_rot_f, forcedWeights, nforced, qswapdim);
 ComputeRefRMSD(); // does (almost) nothing for fixed restraint; needed for moving restraint
//
 DebugM(1,"Current RMSD : "<<instRMSD<<" : Reference RMSD : "<<refRMSD<<"\n");
 DebugM(1,"LambdaRef : "<<LambdaRef<<" : LambdaKf : "<<LambdaKf<<"\n");
//
//  force prefactor
 double pref = ComputeForcePrefactor(); // overloaded in derived classes
 DebugM(9,"Force prefactor is "<<pref<<" (force constant was "<<Kf<<")\n");
//
//==================================================================================================================================
//========== finite difference test v (if requested); putting this test here is the simplest (though inelegant) solution
 if (fdtest) {
  double oriRMSD = instRMSD; // save correct rmsd
//
  std::map<int, Vector> fmap;      // maps aid to forces
  std::map<int, std::vector<int> > imap; // maps aid to o/f atom indices
  double* f=forces_f;
  for (i=0; i<nforced; i++, f+=3) {
   int aid = iatom_f[i];
   Vector force(f[0], f[1], f[2]); // will multiply below by pref to get gradient
   fmap[aid]=force; // store force for comparison
   imap[aid]=std::vector<int> (1,i); //store forced atom index
  }
  if (qdiffrot) {
   f=forces_o;
   for (i=0; i<norient; i++, f+=3) {
    int aid = iatom_o[i];
    Vector force(f[0], f[1], f[2]);
    std::map<int, Vector>::iterator ind = fmap.find(aid) ; // check whether id already in the map
    if (ind==fmap.end()) {// not found
     fmap[aid]=force; // force
     imap[aid]=std::vector<int> (1,-1); // indicate that forced (first) index undefined (-1)
    } else {
     fmap[aid]+=force; // add force contribution from o atoms
    }
    imap[aid].push_back(i); // add orientation index
   }
  } // diffrot
// now iterate over all positions, perturb coordinates, and compute finite differences
// 
  double err, maxerr=0;
// remove COM from orientation atoms (otherwise, COM would be subtracted twice from r_f)
  if (qtrans) {
   for (i=0;i<3*norient;){rcurrent_o[i++]-=COM[0];rcurrent_o[i++]-=COM[1];rcurrent_o[i++]-=COM[2]; }
  }
//
  std::map<int, Vector>::iterator fc=fmap.begin(); // does not want to go inside loop
  for (std::map<int, std::vector<int> >::iterator i=imap.begin(); fc!=fmap.end(); ++fc, ++i ) {
   int aid = i->first;
   double* rf =                 (i->second[0] > -1) ? (rcurrent_f + 3*(i->second[0]) ) : NULL ; // pointer to forcing coordinate
   double* ro = (qdiffrot && (i->second.size()> 1)) ? (rcurrent_o + 3*(i->second[1]) ) : NULL ; // pointer to orientation coordinate
//
//DebugM(1, "RF : " << rf << "RO : " << ro<<"\n" );
   for (int k=0; k<3 ; k++ ) { // over the components
    double saved ; if (rf) saved = rf[k] ; else saved = ro[k] ;
//
    double de = 0. ; // finite difference
    for (int l=1; l > -2 ; l-=2) { // central finite-diff 
     if (rf) rf[k] += dh*l;
     if (ro) ro[k] += dh*l;
 // recompute rmsd:

     if (qtrans) {
      free(COM); COM = (double*) com(rcurrent_o,orientWeights,norient,qswapdim); // note : make sure COM has been subtracted above 
      u[0]=1.; u[1]=0.;  u[2]=0.;  u[3]=0.; u[4]=1.;  u[5]=0.;  u[6]=0.;  u[7]=0.; u[8]=1.;
      if (qorient) {
       RMSBestFit(rtarget_o,rcurrent_o,orientWeights,norient,u,qswapdim);
       free(rtarget_rot_f); 
      } //qorient
      rtarget_rot_f=(double*) matmul(u, rtarget_f, 3, 3, nforced); // rotated target structure
 // instead of subtracting COM from rcurrent, add it to rtarget (to avoid modding rcurrent_f, which is reused in the loop)
      for (int i=0;i<3*nforced;){rtarget_rot_f[i++]+=COM[0];rtarget_rot_f[i++]+=COM[1];rtarget_rot_f[i++]+=COM[2]; }
 //compute RMSD
      instRMSD = rmsd(rcurrent_f, rtarget_rot_f, forcedWeights, nforced, qswapdim);
     } else { //not qtrans
      instRMSD = rmsd(rcurrent_f, rtarget_f, forcedWeights, nforced, qswapdim);
     } // qtrans
     de += GetEnergy() * l;
// restore correct coordinates
     if (rf) rf[k] = saved ;
     if (ro) ro[k] = saved ;
    } //l

//DebugM(1, "AtomID : "<< (i->first)<<" analytical : " << pref * fc-> second[k]<< "  FD: "<< de/(2.*dh)<<"\n") ;
DebugM(9, "AtomID : "<< (i->second[0])<<" analytical : " << pref * fc-> second[k]<< "  FD: "<< de/(2.*dh)<<"\n") ;
//DebugM(1, "AtomID : "<< (i->second[0])<<" analytical : " << pref * fc-> second[k]<< "  FD: "<< de/(2.)<<"\n") ;
//    fc->second[k] = pref * (fc->second[k]) ;
//
    err = fabs ( pref * fc->second[k] - de/(2.*dh) ) ; // analytical derivative minus the finite difference
    if ( err  > maxerr ) maxerr=err;
   } // over components (k)
  } // over map indices
  instRMSD=oriRMSD; // restore correct rmsd
//
  DebugM(1, "RMSD ApplyForce : The Maximum difference between analytical and FD force is " <<maxerr<<"\n" );
 }
// ========================== end of Finite Difference test ^ ================================
//
// ====================== compute forces from gradients and communicate them to Master =======
//
// send forces to the CFE object
// forcing atoms
 double* f=forces_f;
 for (i=0; i<nforced; i++, f+=3) {
  int aid = iatom_f[i];
  Vector force(f[0], f[1], f[2]);
  force*=(-pref);
//
  DebugM(9,"Applying force "<<force<<" to forced atom "<<aid<<"\n");
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
   DebugM(9,"Applying force "<<force<<" to orientation atom "<<aid<<"\n");
   CFE.modifyForcedAtoms().add(aid);
   CFE.modifyAppliedForces().add( force );
  }
 }
//
// ----------------------v clean up memory
//
 if (qtrans) {
  free(COM);
  if (qorient || fdtest) {
   free(rtarget_rot_f);
  }
  if(qdiffrot) {
   free(forces_o);
   free(rcurrent_o);
  }
 } //qtrans
 free(forces_f);
 free(rcurrent_f);
// ----------------------^ clean up memory
} //ApplyForce
#endif //RMSD restraints
