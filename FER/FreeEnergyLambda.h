/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#if !defined(LAMBDA_HPP)
  #define LAMBDA_HPP

class ALambdaControl {
//private:
public:

  // don't forget to change operator= if member variables change
  int     NumSteps;         // for pmf block
  int     NumEquilSteps;    // for mcti block
  int     NumAccumSteps;    // "
  int     NumRepeats;       // "
  int     NumPrintSteps;    // for pmf & mcti blocks
  int     StartStep;        // "
  int     StopStep;         // "
  double  LambdaKf;         // "
  double  LambdaRef;        // "
  feptask_t  Task;             // "
  double  Sum_dU_dLambda;   // for accumulating dU/dLambda
  int     Num_dU_dLambda;   // number averaged
  double  MCTI_Integration; // for accumulating <dU/dLambda> * dLambda

  static int  CurrStep;     // for all pmf & mcti blocks

public:
  ALambdaControl();
  void    Init(ALambdaControl& PriorBlock);
  double  GetLambdaKf();
  double  GetLambdaRef();
  bool    IsActive() { return ( (CurrStep>=StartStep) && (CurrStep<=GetLastStep()) ) ; } 
  Bool_t  IsTimeToPrint();
  bool    IsFirstStep()               { return ( CurrStep == StartStep ); }
  Bool_t  IsTimeToPrint_dU_dLambda();
  Bool_t  IsTimeToClearAccumulator();
  Bool_t  IsEndOf_MCTI_Step();
  Bool_t  IsEndOf_MCTI();
  void    PrintHeader(double dT);
  void    PrintLambdaHeader(double dT);
  void    IncCurrStep() {CurrStep++;}
  ALambdaControl&  operator= (ALambdaControl& PmfBlock);
  void    GetTaskStr(char* Str);
  void    GetPaddedTaskStr(char* Str);
  void    Integrate_MCTI();
  void    Accumulate(double dU_dLambda);
  double  GetIntegration();
  double  GetAccumulation();
  void    ZeroAccumulator() {
    Sum_dU_dLambda = 0.0;
    Num_dU_dLambda = 0;
  }

  int    GetNumSteps()                 {return( (GetLastStep() - StartStep + 1) );}
  int    GetNumStepsSoFar()            {return( CurrStep - StartStep + 1);}
  int    GetNumAccumStepsSoFar();
  int    GetNum_dU_dLambda()           {return(Num_dU_dLambda);}
  void   SetNumSteps(int Steps)        {NumSteps=Steps;}
  void   SetNumEquilSteps(int Steps)   {NumEquilSteps=Steps;}
  void   SetNumAccumSteps(int Steps)   {NumAccumSteps=Steps;}
  void   SetNumPrintSteps(int Steps)   {NumPrintSteps=Steps;}
  void   SetNumRepeats(int Repeats)    {NumRepeats=Repeats;}
  void   SetStartStep(int Step)        {StartStep=Step;}
  void   SetStopStep(int Step)         {StopStep=Step;}
  void   SetLambdaKf(double l)         {LambdaKf=l;}
  void   SetLambdaRef(double l)        {LambdaRef=l;}
  void   SetTask(feptask_t t)          {Task=t;}
  feptask_t GetTask()                     {return(Task);}

private:
  bool   IsLastStep() { return (CurrStep==GetLastStep()); }
  int    GetLastStep();
};

#endif

