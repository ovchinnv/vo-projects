/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
//
// written by David Hurwitz, March to May 1998.
//
// modifications by Victor Ovchinnikov 2013
//
//#define  DEBUGM // debugging output
#define MIN_DEBUG_LEVEL 1
#define MAX_DEBUG_LEVEL 10

#if !defined(RESTRAINT_HPP)
  #define RESTRAINT_HPP

class GlobalMasterFreeEnergy;

//****************************************************************************
//  ARestraint:
//****************************************************************************
class ARestraint {
//-----------------------------------------------------------------------------
// ARestraint is the base class for all restraints :
//
//             -----------------ARestraint------------------------------------\
//           /                /           \                 \                  \
//         /                /               \                 \                 \
//        |                 |                |                 |                 \
//   APosRestraint   ADistRestraint   AnAngleRestraint   ADiheRestraint       AnRMSDRestraint
//        |                 |                |                 |                  |
//        |                 |                |                 |                  |
//   AFixedPosRestraint     |         AFixedAngleRestraint     |              AFixedRMSDRestraint
//   ABoundPosRestraint     |         ABoundAngleRestraint     |              ABoundRMSDRestraint
//   AMovingPosRestraint   |         AMovingAngleRestraint   |              AMovingRMSDRestraint
//                          |                                  |
//                   AFixedDistRestraint                 AFixedDiheRestraint
//                   ABoundDistRestraint                 ABoundDiheRestraint
//                   AMovingDistRestraint               AMovingDiheRestraint
//-----------------------------------------------------------------------------
protected:
  double    Kf;
  int       NumGroups;
  AGroup*   Groups;
  Vector*   COMs;

  // lambda is the same for all Moving restraints
  static double  LambdaKf;
  static double  LambdaRef;

public:
  ARestraint();
  virtual ~ARestraint();
  int     GetNumGroups()    {return(NumGroups);}
  void    SetKf(double k)  {Kf=k;}
  double  GetKf()           {return(Kf);}
  void    SetLambdaKf(double l)   {LambdaKf=l;}
  void    SetLambdaRef(double l) {LambdaRef=l;}
  double  GetLambdaKf()                  {return(LambdaKf);}
  double  GetLambdaRef()                 {return(LambdaRef);}
  virtual void    SetGroup(AGroup& Group, int GroupIndex);
  virtual void    SetGroups(AGroup* groups, int ngroup);
//---------------------------------------------//
  // pure virtual functions
  virtual void     ApplyForce(GlobalMasterFreeEnergy& CFE, bool test_fd=0, double dh=0.00001)=0; // apply forces; test analytical forces against FD
  virtual double   GetEnergy() = 0;
  virtual void     GetStr(char* Str) = 0;
  virtual void     PrintInfo() = 0;
  virtual double   Get_dU_dLambda() {return(0.0);}
  virtual Bool_t   IsMoving() {return(kFalse);}
protected:
  void    UpdateCOMs(GlobalMasterFreeEnergy& CFE);
  void    DistributeForce(int WhichGroup, Vector Force, GlobalMasterFreeEnergy& CFE);
  void    EarlyExit(const char* Str, int AtomID = -999); //VO 2013 : invalid default value provided for a generic error message
};
//
// Positional restraints:
//****************************************************************************
class APosRestraint : public ARestraint {
//---------------------------------------------------------------------------
// APosRestraint is a derived class of ARestraint
// this is a virtual class because it does not implement GetEnergy amd GetStr
//---------------------------------------------------------------------------
protected:
  Vector RefPos;
public:
  APosRestraint();
  void    PrintInfo();
  void    SetRefPos(Vector Pos) {RefPos=Pos;}
  double  GetDistance()  { return  (RefPos - COMs[0]).length() ;}
  Vector  GetPosTarget()  { return(RefPos);}
protected:
  void    ApplyForce(GlobalMasterFreeEnergy& CFE, bool fd, double dh);
  double  GetE(Vector RefPos, double LambdaKf=1.0);
  Vector GetGrad(int WhichGroup, Vector RefPos, double LambdaKf=1.0);
  Vector virtual GetGradient(int)=0;
};
//****************************************************************************
//  AFixedPosRestraint, ABoundPosRestraint, AMovingPosRestraint:
//****************************************************************************
class AFixedPosRestraint : public APosRestraint {
//-------------------------------------------------------------------
// AFixedPosRestraint is derived from APosRestraint
//-------------------------------------------------------------------
public:
  double  GetEnergy(); // inherited from top
  void    GetStr(char* Str) {strcpy(Str, "Fixed Position Restraint");} // return my name
protected:
  Vector GetGradient(int WhichGroup); // local
};

class ABoundPosRestraint : public APosRestraint {
//-------------------------------------------------------------------
// ABoundPosRestraint is a derived class of APosRestraint
//-------------------------------------------------------------------
private:
  double  RefDist;
  Bound_t Bound;
public:
  void    SetRefDist(double Dist) {RefDist=Dist;}
  void    SetBound(Bound_t b) {Bound=b;}
  double  GetRefDist()            {return(RefDist);}
  Bound_t GetBound()              {return(Bound);}
  double  GetEnergy();
  void    GetStr(char* Str) {
    strcpy(Str, "Bound Position Restraint");
  }
protected:
  Vector GetGradient(int WhichGroup);
};

class AMovingPosRestraint : public APosRestraint {
//-------------------------------------------------------------------
// AMovingPosRestraint is a derived class of APosRestraint
//-------------------------------------------------------------------
private:
  Vector StartPos;
  Vector StopPos;
public:
  void    SetStartPos(Vector Pos)       {StartPos=Pos;}
  void    SetStopPos(Vector Pos)        {StopPos=Pos;}
  Vector GetStartPos()                  {return(StartPos);}
  Vector GetStopPos()                   {return(StopPos);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Position Restraint");
  }
  Vector GetPosTarget() {
    return(StopPos*LambdaRef + StartPos*(1.0-LambdaRef));
  }
  double  GetDistance() {
    Vector RefPos = StopPos*LambdaRef + StartPos*(1.0-LambdaRef);
    return (RefPos-COMs[0]).length();
  }
protected:
  Vector GetGradient(int WhichGroup);
};

//
// Distance restraints:
//****************************************************************************
//
class ADistRestraint : public ARestraint {
//-------------------------------------------------------------------
// ADistRestraint is a derived class of ARestraint 
// this is an abstract class
//-------------------------------------------------------------------
protected:
  double  RefDist;
public:
  ADistRestraint();
  void    PrintInfo();
  void    SetRefDist(double Dist)  {RefDist=Dist;}
  double  GetDistTarget()          {return(RefDist);}
protected:
  void    ApplyForce(GlobalMasterFreeEnergy& CFE, bool fd, double dh);
  double  GetE(double RefDist, double LambdaKf=1.0);
  Vector GetGrad(int WhichGroup, double RefDist, double LambdaKf=1.0);
  Vector virtual GetGradient(int)=0;
};
//****************************************************************************
//  AFixedDistRestraint, ABoundDistRestraint, AMovingDistRestraint:
//****************************************************************************
class AFixedDistRestraint : public ADistRestraint {
//-------------------------------------------------------------------
// AFixedDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
public:
  double  GetEnergy();
  void    GetStr(char* Str) {
    strcpy(Str, "Fixed Distance Restraint");
  }
protected:
  Vector GetGradient(int WhichGroup); // to be called only from class methods
};

class ABoundDistRestraint : public ADistRestraint {
//-------------------------------------------------------------------
// ABoundDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  Bound_t Bound;
public:
  void    SetBound(Bound_t b)  {Bound=b;}
  Bound_t GetBound()           {return(Bound);}
  double  GetEnergy();
  void    GetStr(char* Str) {
    strcpy(Str, "Bound Distance Restraint");
  }
protected:
  Vector GetGradient(int WhichGroup);
};

class AMovingDistRestraint : public ADistRestraint {
//-------------------------------------------------------------------
// AMovingDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  double  StartDist;
  double  StopDist;
public:
  void    SetStartDist(double Dist)      {StartDist=Dist;}
  void    SetStopDist(double Dist)       {StopDist=Dist;}
  double  GetStartDist()                 {return(StartDist);}
  double  GetStopDist()                  {return(StopDist);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Distance Restraint");
  }
  double  GetDistTarget() {
    return(StopDist*LambdaRef + StartDist*(1.0-LambdaRef));
  }
protected:
  Vector GetGradient(int WhichGroup);
};
//
// Angle restraints:
//****************************************************************************
//
class AnAngleRestraint : public ARestraint {
//-------------------------------------------------------------------
// AnAngleRestraint is a derived class of ARestraint
//-------------------------------------------------------------------
protected:
  double  RefAngle;     // in radians
public:
  AnAngleRestraint();
  void    PrintInfo();
  void    SetRefAngle(double Angle)  {RefAngle=Angle;}
  double  GetAngleTarget()           {return(RefAngle);}
protected:
  void    ApplyForce(GlobalMasterFreeEnergy& CFE, bool fd, double dh);
  double  GetE(double RefAngle, double LambdaKf=1.0);
  Vector GetGrad(int WhichGroup, double RefAngle, double LambdaKf=1.0);
  double  GetAngle(Vector& A, Vector& B, Vector& C);
  Vector virtual GetGradient(int)=0;
};
//****************************************************************************
//  AFixedAngleRestraint, ABoundAngleRestraint, AMovingAngleRestraint:
//****************************************************************************
class AFixedAngleRestraint : public AnAngleRestraint {
//-------------------------------------------------------------------
// AFixedAngleRestraint is a derived class of AnAngleRestraint
//-------------------------------------------------------------------
public:
  double  GetEnergy();
  Vector GetGradient(int WhichGroup);
  void    GetStr(char* Str) {
    strcpy(Str, "Fixed AngleRestraint");
  }
};

class ABoundAngleRestraint : public AnAngleRestraint {
//-------------------------------------------------------------------
// ABoundAngleRestraint is a derived class of AnAngleRestraint
//-------------------------------------------------------------------
private:
  double  RefAngle;     // in radians
  Bound_t Bound;
public:
  void    SetRefAngle(double Angle)  {RefAngle=Angle;}
  void    SetBound(Bound_t b)    {Bound=b;}
  Bound_t GetBound()                 {return(Bound);}
  double  GetEnergy();
  void    GetStr(char* Str) {
    strcpy(Str, "Bound Angle Restraint");
  }
protected:
  Vector GetGradient(int WhichGroup);
};

class AMovingAngleRestraint : public AnAngleRestraint {
//-------------------------------------------------------------------
// AMovingAngleRestraint is a derived class of AnAngleRestraint
//-------------------------------------------------------------------
private:
  double  StartAngle;     // in radians
  double  StopAngle;      // in radians
public:
  void    SetStartAngle(double Angle)    {StartAngle=Angle;}
  void    SetStopAngle(double Angle)     {StopAngle=Angle;}
  double  GetStartAngle()                {return(StartAngle);}
  double  GetStopAngle()                 {return(StopAngle);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Angle Restraint");
  }
  double  GetAngleTarget() {
    return(StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef));
  }
protected:
  Vector GetGradient(int WhichGroup);
};
//
// Dihedral restraints:
//****************************************************************************
//
class ADiheRestraint : public ARestraint {
//-------------------------------------------------------------------
// ADiheRestraint is a derived class of ARestraint
//-------------------------------------------------------------------
protected:
  double  RefAngle;     // in radians
public:
  ADiheRestraint();
  void    PrintInfo();
  Bool_t  TwoTargets()      {return(kFalse);}
  double  GetDiheTarget1()  {return(RefAngle);}
  double  GetDiheTarget2()  {return(0);}
protected:
  void    ApplyForce(GlobalMasterFreeEnergy& CFE, bool fd, double dh);
  double  GetE(double RefDihe, double Const);
  Vector GetGrad(int WhichGroup, double RefDihe, double Const);
  Vector gradU(Vector& P1P2P3, Vector& P4P5P6,
                Vector& dP1,    Vector& dP2,    Vector& dP3,
                Vector& dP4,    Vector& dP5,    Vector& dP6);
  double  GetDihe(Vector& A, Vector& B, Vector& C, Vector& D);
  Vector virtual GetGradient(int)=0;
};

//****************************************************************************
//  AFixedDiheRestraint, ABoundDiheRestraint, AMovingDiheRestraint:
//****************************************************************************
class AFixedDiheRestraint : public ADiheRestraint {
//-------------------------------------------------------------------
// AFixedDiheRestraint is a derived class of ADiheRestraint
//-------------------------------------------------------------------
public:
  void    SetRefAngle(double Angle)  {RefAngle=Angle;}
  double  GetEnergy();
  void    GetStr(char* Str) {
    strcpy(Str, "Fixed   Dihedral Restraint");
  }
protected:
  Vector GetGradient(int WhichGroup);
};

class ABoundDiheRestraint : public ADiheRestraint {
//-------------------------------------------------------------------
// ABoundDiheRestraint is a derived class of ADiheRestraint
//-------------------------------------------------------------------
private:
  double  LowerAngle;     // radians (between 0 and 2pi)
  double  UpperAngle;     // radians (between 0 and 2pi)
  double  IntervalAngle;  // radians
public:
  void    SetLowerAngle(double Angle)     {LowerAngle=Angle;}
  void    SetUpperAngle(double Angle)     {UpperAngle=Angle;}
  void    SetIntervalAngle(double Angle)  {IntervalAngle=Angle;}
  double  GetLowerAngle()                 {return(GetDiheTarget1());} // not used -- remove ?
  double  GetUpperAngle()                 {return(GetDiheTarget2());} // not used -- remove ?
  double  GetIntervalAngle()              {return(IntervalAngle);}
  double  GetEnergy();
  void    GetStr(char* Str) {
    strcpy(Str, "Bound Dihedral Restraint");
  }
  Bool_t  TwoTargets()      {return(kTrue);}
  double  GetDiheTarget1()  {return(LowerAngle);}
  double  GetDiheTarget2()  {return(UpperAngle);}
//
protected:
  Vector GetGradient(int WhichGroup);
};

class AMovingDiheRestraint : public ADiheRestraint {
//-------------------------------------------------------------------
// AMovingDiheRestraint is a derived class of ADiheRestraint
//-------------------------------------------------------------------
private:
  double  StartAngle;
  double  StopAngle;
public:
  void    SetStartAngle(double Angle)    {StartAngle=Angle;}
  void    SetStopAngle(double Angle)     {StopAngle=Angle;}
  double  GetStartAngle()                {return(StartAngle);}
  double  GetStopAngle()                 {return(StopAngle);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Dihedral Restraint");
  }
  Bool_t  TwoTargets()     {return(kFalse);}
  double  GetDiheTarget1() {
    return(StopAngle*LambdaRef + StartAngle*(1.0-LambdaRef));
  }
protected:
  Vector GetGradient(int WhichGroup);
};
//
#ifdef FE_RESTRAINT_RMSD_FORTRAN
// RMSD restraints:
//****************************************************************************
class AnRMSDRestraint : public ARestraint {
//---------------------------------------------------------------------------
// AnRMSDRestraint is a derived class of ARestraint
// this is a virtual class because it does not implement GetEnergy amd GetStr
//---------------------------------------------------------------------------
protected:
  double refRMSD;   // target RMS distance
  double instRMSD;  //instantaneous RMS distance
  double u[9];
  double* ugrad;    // derivatives of rotation matrix
  bool qorient  ;   // whether RMSD is to be computed after best-fit orientation
  bool qdiffrot ;   // whether orientation atoms are the same as Moving (RMSD) atoms
                    // if yes, gradient computation does not require rotation matrix derivatives (simpler)
  bool qtrans ;     // if not qorient, translate to COM before applying forces (thus computing orientational best-fit only)
public:
  AnRMSDRestraint();
  ~AnRMSDRestraint(); // modified constructor to free ugrad
  void    PrintInfo();
  void    SetRefRMSD(double v)  {refRMSD=v;}
  void    SetQorient(bool q)  {qorient=q;}
  void    SetQtrans(bool q)  {qtrans=q;}
  virtual void    ComputeRefRMSD() {
   if ( refRMSD < 0.0 ) { // this may be useful if one does not care to specify the initial RMSD
//   DebugM(1,"Setting Reference RMSD to : "<<instRMSD<<"\n");
    refRMSD = instRMSD; // assuming instRMSD has been computed
   }
  }
//  double  GetRMSD() {return instRMSD;} // RMSD to target structure
//  double  GetRMSDTarget() {return refRMSD;} // desired RMSD to target structure
  void    qDiffRotCheck(); // Compute qdiffrot
  void    SetGroups(AGroup *, int );
protected:
  void    ApplyForce(GlobalMasterFreeEnergy& CFE, bool fd, double dh);
  virtual double ComputeForcePrefactor()=0;
};
//****************************************************************************
//  Derived classes : AFixedRMSDRestraint, ABoundRMSDRestraint, AMovingRMSDRestraint:
//****************************************************************************
class AFixedRMSDRestraint : public AnRMSDRestraint {
//-------------------------------------------------------------------
public:
  AFixedRMSDRestraint() : AnRMSDRestraint() {;}
//
  double  GetEnergy() { double e = (instRMSD - refRMSD); return ( 0.5 * Kf * e * e ); }
  void    GetStr(char* Str) { strcpy(Str, "Fixed RMSD Restraint"); }
protected:
  double ComputeForcePrefactor() {
   double RMS = __MAX(instRMSD, kALittle); //  prevent divide overflow
   double pref = Kf * (1.0 - refRMSD/RMS);
   return pref;
  }
};

class ABoundRMSDRestraint : public AnRMSDRestraint {
//-------------------------------------------------------------------
// ABoundDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  Bound_t Bound;
public:
  ABoundRMSDRestraint() : AnRMSDRestraint() { Bound = kUpper ;}
  void    SetBound(Bound_t b)  {Bound=b;}
  Bound_t GetBound()               {return(Bound);}
  double  GetEnergy() { double e = (instRMSD - refRMSD); if (Bound==kLower) e = -e; return ( e>0. ? 0.5 * Kf * e * e :  0.0 ); }
  void    GetStr(char* Str) {  strcpy(Str, "Bound RMSD Restraint"); }
protected:
  double ComputeForcePrefactor() {
   double dRMSD = (instRMSD - refRMSD); 
   double RMS = __MAX(instRMSD, kALittle); //  prevent divide overflow
   if (Bound==kUpper) return ( dRMSD > 0. ? Kf * dRMSD/RMS : 0.0); //Upper bonud
                      return ( dRMSD < 0. ? Kf * dRMSD/RMS : 0.0); //Lower bound
  }
};

class AMovingRMSDRestraint : public AnRMSDRestraint {
//-------------------------------------------------------------------
// AMovingDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  double  startRMSD;
  double  stopRMSD;
public:
  AMovingRMSDRestraint() : AnRMSDRestraint() {;}
  void    SetStartRMSD(double v)      {startRMSD=v;}
  void    SetStopRMSD(double v)       {stopRMSD =v;}
  double  GetStartRMSD()              {return(startRMSD);}
  double  GetStopRMSD()               {return(stopRMSD);}
  double  GetEnergy() { double e = (instRMSD - refRMSD); return ( 0.5 * Kf * e * e ); }
  virtual void ComputeRefRMSD() {
   if ( startRMSD < 0.0 ) startRMSD = instRMSD;
   refRMSD = stopRMSD*LambdaRef + startRMSD*(1.0-LambdaRef);
  }
  double  Get_dU_dLambda() { return ( Kf * (instRMSD - refRMSD) * (stopRMSD - startRMSD) ) ; }
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) { strcpy(Str, "Moving RMSD Restraint"); }
protected:
  double ComputeForcePrefactor() {
   double RMS = __MAX(instRMSD, kALittle); //  prevent divide overflow
   double pref = Kf * (1.0 - refRMSD/RMS);
   return pref;
  }
};
#endif //RMSD
#endif
