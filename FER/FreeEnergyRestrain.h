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
  double    m_Kf;
  int       m_NumGroups;
  AGroup*   m_pGroups;
  Vector*   m_pCOMs;

  // lambda is the same for all Moving restraints
  static double  m_LambdaKf;
  static double  m_LambdaRef;

public:
  ARestraint();
  virtual ~ARestraint();
  int     GetNumGroups()    {return(m_NumGroups);}
  void    SetKf(double Kf)  {m_Kf=Kf;}
  double  GetKf()           {return(m_Kf);}
  void    SetLambdaKf(double LambdaKf)   {m_LambdaKf=LambdaKf;}
  void    SetLambdaRef(double LambdaRef) {m_LambdaRef=LambdaRef;}
  double  GetLambdaKf()                  {return(m_LambdaKf);}
  double  GetLambdaRef()                 {return(m_LambdaRef);}
  void    SetGroup(AGroup& Group, int GroupIndex);
//the four functions below should really be in the derived classes (which know how many groups they need)
  void    SetGroups(AGroup& Group1);
  void    SetGroups(AGroup& Group1, AGroup& Group2);
  void    SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3);
  void    SetGroups(AGroup& Group1, AGroup& Group2, AGroup& Group3, AGroup& Group4);
//---------------------------------------------//

  // (pure) virtual functions
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
  Vector m_RefPos;
public:
  APosRestraint();
  void    PrintInfo();
  void    SetRefPos(Vector Pos) {m_RefPos=Pos;}
  double  GetDistance()          { return  (m_RefPos - m_pCOMs[0]).length() ;}
  Vector  GetPosTarget()          { return(m_RefPos);}
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
  double  m_RefDist;
  Bound_t m_Bound;
public:
  void    SetRefDist(double Dist) {m_RefDist=Dist;}
  void    SetBound(Bound_t Bound) {m_Bound=Bound;}
  double  GetRefDist()            {return(m_RefDist);}
  Bound_t GetBound()              {return(m_Bound);}
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
  Vector m_StartPos;
  Vector m_StopPos;
public:
  void    SetStartPos(Vector Pos)       {m_StartPos=Pos;}
  void    SetStopPos(Vector Pos)        {m_StopPos=Pos;}
  Vector GetStartPos()                  {return(m_StartPos);}
  Vector GetStopPos()                   {return(m_StopPos);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Position Restraint");
  }
  Vector GetPosTarget() {
    return(m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef));
  }
  double  GetDistance() {
    Vector RefPos = m_StopPos*m_LambdaRef + m_StartPos*(1.0-m_LambdaRef);
    return (RefPos-m_pCOMs[0]).length();
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
  double  m_RefDist;
public:
  ADistRestraint();
  void    PrintInfo();
  void    SetRefDist(double Dist)  {m_RefDist=Dist;}
  double  GetDistTarget()          {return(m_RefDist);}
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
  Bound_t m_Bound;
public:
  void    SetBound(Bound_t Bound)  {m_Bound=Bound;}
  Bound_t GetBound()               {return(m_Bound);}
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
  double  m_StartDist;
  double  m_StopDist;
public:
  void    SetStartDist(double Dist)      {m_StartDist=Dist;}
  void    SetStopDist(double Dist)       {m_StopDist=Dist;}
  double  GetStartDist()                 {return(m_StartDist);}
  double  GetStopDist()                  {return(m_StopDist);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Distance Restraint");
  }
  double  GetDistTarget() {
    return(m_StopDist*m_LambdaRef + m_StartDist*(1.0-m_LambdaRef));
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
  double  m_RefAngle;     // in radians
public:
  AnAngleRestraint();
  void    PrintInfo();
  void    SetRefAngle(double Angle)  {m_RefAngle=Angle;}
  double  GetAngleTarget()           {return(m_RefAngle);}
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
  double  m_RefAngle;     // in radians
  Bound_t m_Bound;
public:
  void    SetRefAngle(double Angle)  {m_RefAngle=Angle;}
  void    SetBound(Bound_t Bound)    {m_Bound=Bound;}
  Bound_t GetBound()                 {return(m_Bound);}
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
  double  m_StartAngle;     // in radians
  double  m_StopAngle;      // in radians
public:
  void    SetStartAngle(double Angle)    {m_StartAngle=Angle;}
  void    SetStopAngle(double Angle)     {m_StopAngle=Angle;}
  double  GetStartAngle()                {return(m_StartAngle);}
  double  GetStopAngle()                 {return(m_StopAngle);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Angle Restraint");
  }
  double  GetAngleTarget() {
    return(m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef));
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
  double  m_RefAngle;     // in radians
public:
  ADiheRestraint();
  void    PrintInfo();
  Bool_t  TwoTargets()      {return(kFalse);}
  double  GetDiheTarget1()  {return(m_RefAngle);}
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
  void    SetRefAngle(double Angle)  {m_RefAngle=Angle;}
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
  double  m_LowerAngle;     // radians (between 0 and 2pi)
  double  m_UpperAngle;     // radians (between 0 and 2pi)
  double  m_IntervalAngle;  // radians
public:
  void    SetLowerAngle(double Angle)     {m_LowerAngle=Angle;}
  void    SetUpperAngle(double Angle)     {m_UpperAngle=Angle;}
  void    SetIntervalAngle(double Angle)  {m_IntervalAngle=Angle;}
  double  GetLowerAngle()                 {return(GetDiheTarget1());} // not used -- remove ?
  double  GetUpperAngle()                 {return(GetDiheTarget2());} // not used -- remove ?
  double  GetIntervalAngle()              {return(m_IntervalAngle);}
  double  GetEnergy();
  void    GetStr(char* Str) {
    strcpy(Str, "Bound Dihedral Restraint");
  }
  Bool_t  TwoTargets()      {return(kTrue);}
  double  GetDiheTarget1()  {return(m_LowerAngle);}
  double  GetDiheTarget2()  {return(m_UpperAngle);}
//
protected:
  Vector GetGradient(int WhichGroup);
};

class AMovingDiheRestraint : public ADiheRestraint {
//-------------------------------------------------------------------
// AMovingDiheRestraint is a derived class of ADiheRestraint
//-------------------------------------------------------------------
private:
  double  m_StartAngle;
  double  m_StopAngle;
public:
  void    SetStartAngle(double Angle)    {m_StartAngle=Angle;}
  void    SetStopAngle(double Angle)     {m_StopAngle=Angle;}
  double  GetStartAngle()                {return(m_StartAngle);}
  double  GetStopAngle()                 {return(m_StopAngle);}
  double  GetEnergy();
  double  Get_dU_dLambda();
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) {
    strcpy(Str, "Moving Dihedral Restraint");
  }
  Bool_t  TwoTargets()     {return(kFalse);}
  double  GetDiheTarget1() {
    return(m_StopAngle*m_LambdaRef + m_StartAngle*(1.0-m_LambdaRef));
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
public:
  AnRMSDRestraint();
  AnRMSDRestraint(AGroup&, AGroup&);
  ~AnRMSDRestraint(); // modified constructor to free ugrad
  void    PrintInfo();
  void    SetRefRMSD(double v)  {refRMSD=v;}
  void    ComputeRefRMSD() {;} // nothing: restraint is fixed by default
//  double  GetRMSD() {return instRMSD;} // RMSD to target structure
//  double  GetRMSDTarget() {return refRMSD;} // desired RMSD to target structure
  void    qDiffRotCheck(); // Compute qdiffrot
  void    SetGroups(AGroup&, AGroup&);
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
  AFixedRMSDRestraint(AGroup& g1, AGroup& g2) : AnRMSDRestraint(g1, g2) {;}
//
  double  GetEnergy() { double e = (instRMSD - refRMSD); return ( 0.5 * m_Kf * e * e ); }
  void    GetStr(char* Str) { strcpy(Str, "Fixed RMSD Restraint"); }
protected:
  double ComputeForcePrefactor() {
   double RMS = __MAX(instRMSD, kALittle); //  prevent divide overflow
   double pref = m_Kf * (1.0 - refRMSD/RMS);
   return pref;
  }
};

class ABoundRMSDRestraint : public AnRMSDRestraint {
//-------------------------------------------------------------------
// ABoundDistRestraint is a derived class of ADistRestraint
//-------------------------------------------------------------------
private:
  Bound_t m_Bound;
public:
  ABoundRMSDRestraint() : AnRMSDRestraint() { m_Bound = kUpper ;}
  ABoundRMSDRestraint(AGroup& g1, AGroup& g2) : AnRMSDRestraint(g1, g2) { m_Bound = kUpper ;}
  void    SetBound(Bound_t Bound)  {m_Bound=Bound;}
  Bound_t GetBound()               {return(m_Bound);}
  double  GetEnergy() { double e = (instRMSD - refRMSD); if (m_Bound==kLower) e = -e; return ( e>0. ? 0.5 * m_Kf * e * e :  0.0 ); }
  void    GetStr(char* Str) {  strcpy(Str, "Bound RMSD Restraint"); }
protected:
  double ComputeForcePrefactor() {
   double dRMSD = (instRMSD - refRMSD); 
   double RMS = __MAX(instRMSD, kALittle); //  prevent divide overflow
   if (m_Bound==kUpper) return ( dRMSD > 0. ? m_Kf * dRMSD/RMS : 0.0); //Upper bonud
                        return ( dRMSD < 0. ? m_Kf * dRMSD/RMS : 0.0); //Lower bound
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
  AMovingRMSDRestraint(AGroup& g1, AGroup& g2) : AnRMSDRestraint(g1, g2) {;}
  void    SetStartRMSD(double v)      {startRMSD=v;}
  void    SetStopRMSD(double v)       {stopRMSD =v;}
  double  GetStartRMSD()              {return(startRMSD);}
  double  GetStopRMSD()               {return(stopRMSD);}
  double  GetEnergy() { double e = (instRMSD - refRMSD); return ( 0.5 * m_Kf * e * e ); }
  void    ComputeRefRMSD() { refRMSD = stopRMSD*m_LambdaRef + startRMSD*(1.0-m_LambdaRef); }
  double  Get_dU_dLambda() { return ( m_Kf * (instRMSD - refRMSD) * (stopRMSD - startRMSD) ) ; }
  Bool_t  IsMoving() {return(kTrue);}
  void    GetStr(char* Str) { strcpy(Str, "Moving RMSD Restraint"); }
protected:
  double ComputeForcePrefactor() {
   double RMS = __MAX(instRMSD, kALittle); //  prevent divide overflow
   double pref = m_Kf * (1.0 - refRMSD/RMS);
   return pref;
  }
};
#endif //RMSD
#endif
