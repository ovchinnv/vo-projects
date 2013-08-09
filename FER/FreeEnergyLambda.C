/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <string.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyLambda.h"
#include "Vector.h"

//#define DEBUGM
#include "Debug.h"

// initialize static member variables
int  ALambdaControl::CurrStep = 0;

ALambdaControl::ALambdaControl() {
//------------------------------------------------------------------------
// initialize member variables
//------------------------------------------------------------------------
  // initialize these for the first block.
  StartStep = 0;
  LambdaKf  = 1.0;         // applied to Kf
  LambdaRef = 0.0;         // applied to reference (pos, dist, angle, dihe)
  Sum_dU_dLambda = 0.0;    // for free energy measurement
  Num_dU_dLambda = 0;      // "
  MCTI_Integration = 0.0;  // "

  // set the rest to illegal values.
  Task =           kUnknownTask;
  NumSteps =      -1;
  NumEquilSteps = -1;
  NumAccumSteps = -1;
  NumPrintSteps = -1;
  NumRepeats =    -1;
  StopStep =      -1;
}

void ALambdaControl::Init(ALambdaControl& PriorBlock) {
//------------------------------------------------------------------------
// initialize this object using the settings for the prior object.
//------------------------------------------------------------------------
  Task =          PriorBlock.Task;
  NumSteps =      PriorBlock.NumSteps;
  NumEquilSteps = PriorBlock.NumEquilSteps;
  NumAccumSteps = PriorBlock.NumAccumSteps;
  NumPrintSteps = PriorBlock.NumPrintSteps;
  NumRepeats =    PriorBlock.NumRepeats;
  switch (PriorBlock.Task) {

    case kUp:       LambdaKf=1.0;  LambdaRef=1.0; break;
    case kStepUp:   LambdaKf=1.0;  LambdaRef=1.0; break;
    case kDown:     LambdaKf=1.0;  LambdaRef=0.0; break;
    case kStepDown: LambdaKf=1.0;  LambdaRef=0.0; break;
    case kStop:     LambdaKf=1.0;  LambdaRef=PriorBlock.LambdaRef; break;

    case kGrow:     LambdaKf=1.0;  LambdaRef=PriorBlock.LambdaRef; break;
    case kStepGrow: LambdaKf=1.0;  LambdaRef=PriorBlock.LambdaRef; break;
    case kFade:     LambdaKf=0.0;  LambdaRef=PriorBlock.LambdaRef; break;
    case kStepFade: LambdaKf=0.0;  LambdaRef=PriorBlock.LambdaRef; break;
    case kNoGrow:   LambdaKf =  PriorBlock.LambdaKf;
                    LambdaRef = PriorBlock.LambdaRef;  break;
    default:        ASSERT(kFalse); break;  //should never get here

  }
  StartStep = PriorBlock.GetLastStep();
}


ALambdaControl& ALambdaControl::operator= (ALambdaControl& PmfBlock) {
//------------------------------------------------------------------------
// copy everything from PmfBlock to this block.
//------------------------------------------------------------------------
  NumSteps =      PmfBlock.NumSteps;
  NumEquilSteps = PmfBlock.NumEquilSteps;
  NumAccumSteps = PmfBlock.NumAccumSteps;
  NumRepeats =    PmfBlock.NumRepeats;
  NumPrintSteps = PmfBlock.NumPrintSteps;
  StartStep =     PmfBlock.StartStep;
  StopStep =      PmfBlock.StopStep;
  LambdaKf =      PmfBlock.LambdaKf;
  LambdaRef =     PmfBlock.LambdaRef;
  Task =          PmfBlock.Task;
  return(*this);
}


void ALambdaControl::Accumulate(double dU_dLambda) {
//------------------------------------------------------------------------
// if lambda is constant, sum dU/dLambda
// if lambda is changing, sum dU/dLambda * dLambda
//------------------------------------------------------------------------

  Num_dU_dLambda++;
  switch (Task) {
    // lambda is constant
    case kStop:      case kNoGrow:    case kStepUp:
    case kStepDown:  case kStepGrow:  case kStepFade:
      Sum_dU_dLambda += dU_dLambda;
      break;
    // lambda is increasing
    case kUp:
    case kGrow:
      Sum_dU_dLambda += dU_dLambda / (double)NumSteps;
      break;
    // lambda is decreasing
    case kDown:
    case kFade:
      Sum_dU_dLambda -= dU_dLambda / (double)NumSteps;
      break;
    // should never get here
    default:
      ASSERT(kFalse);
      break;
  }
}


void ALambdaControl::Integrate_MCTI() {
//------------------------------------------------------------------------
// integrate MCTI:  <dU/dLambda> * dLambda
//------------------------------------------------------------------------
  ASSERT(Task==kStepUp   || Task==kStepDown || 
         Task==kStepGrow || Task==kStepFade);

  switch (Task) {
    // lambda is increasing
    case kStepUp:
    case kStepGrow:
      MCTI_Integration += GetAccumulation() / (double)NumRepeats;
      break;
    // lambda is decreasing
    case kStepDown:
    case kStepFade:
      MCTI_Integration -= GetAccumulation() / (double)NumRepeats;
      break;
    // should never get here
    default:
      ASSERT(kFalse);
      break;
  }
}


double ALambdaControl::GetIntegration() {
//------------------------------------------------------------------------
// get MCTI integral.  integral(dU/dLambda> * dLambda
//------------------------------------------------------------------------
  ASSERT(Task==kStepUp   || Task==kStepDown || 
         Task==kStepGrow || Task==kStepFade);

  return(MCTI_Integration);
}


double ALambdaControl::GetAccumulation() {
//------------------------------------------------------------------------
// if lambda is constant, return (sum/N)
// if lambda is changing, return (sum)
//------------------------------------------------------------------------
  switch (Task) {
    // lambda is constant
    case kStop:    case kNoGrow:
    case kStepUp:  case kStepDown:  case kStepGrow:  case kStepFade:
      ASSERT(Num_dU_dLambda != 0);
      return(Sum_dU_dLambda / (double)Num_dU_dLambda);
    // lambda is changing
    case kUp:  case kDown:  case kGrow:  case kFade:
      return(Sum_dU_dLambda);
    // should never get here
    default:
      ASSERT(kFalse);
      return(0);
  }
}


int ALambdaControl::GetNumAccumStepsSoFar() {
//------------------------------------------------------------------------
// return the total number of steps dU/dLambda has been accumulated
// this is only called during MCTI
//------------------------------------------------------------------------
  ASSERT(Task==kStepUp   || Task==kStepDown || 
         Task==kStepGrow || Task==kStepFade);

  int Count, Maybe;

  Count = (CurrStep-StartStep) / (NumAccumSteps+NumEquilSteps);
  Count *= NumAccumSteps;
  Maybe = (CurrStep-StartStep) % (NumAccumSteps+NumEquilSteps);
  Maybe -= NumEquilSteps;
  if (Maybe > 0) {
    Count += Maybe;
  }
  return(Count);
}


int ALambdaControl::GetNumSteps() {
//------------------------------------------------------------------------
// get the number of steps needed for this pmf or mcti block
//------------------------------------------------------------------------
  // make sure StopStep is calculated
  GetLastStep();
  
  return( (StopStep - StartStep) + 1 );
}


Bool_t ALambdaControl::IsLastStep() {
//------------------------------------------------------------------------
// return true if we're on the last step of this pmf or mcti block
//------------------------------------------------------------------------
	if (CurrStep == GetLastStep()) {
    return(kTrue);
  }
  else {
    return(kFalse);
  }
}


int ALambdaControl::GetLastStep() {
//------------------------------------------------------------------------
// get the last step of this task
//------------------------------------------------------------------------
  // if it's already calculated, just return it
  if (StopStep > 0) {
    return(StopStep);
  }
  // otherwise calculate it
  switch (Task) {
    case kStepUp:
    case kStepDown:
    case kStepGrow:
    case kStepFade:
      StopStep = StartStep +
                  (NumAccumSteps+NumEquilSteps) * NumRepeats;
      break;
    default:
      StopStep = StartStep + NumSteps;
      break;
  }
  // and return it
  return(StopStep);
}


void ALambdaControl::GetPaddedTaskStr(char* Str) {
//------------------------------------------------------------------------
// get a string that describes this task
//------------------------------------------------------------------------
  switch (Task) {
    case kUp:        strcpy(Str, "      Up");      break;
    case kDown:      strcpy(Str, "    Down");      break;
    case kStop:      strcpy(Str, "    Stop");      break;
    case kGrow:      strcpy(Str, "    Grow");      break;
    case kFade:      strcpy(Str, "    Fade");      break;
    case kNoGrow:    strcpy(Str, "  NoGrow");      break;
    case kStepUp:    strcpy(Str, "  StepUp");      break;
    case kStepDown:  strcpy(Str, "StepDown");      break;
    case kStepGrow:  strcpy(Str, "StepGrow");      break;
    case kStepFade:  strcpy(Str, "StepFade");      break;
    default:         strcpy(Str, "Bug Alert!!!");  break;
  }
}


void ALambdaControl::GetTaskStr(char* Str) {
//------------------------------------------------------------------------
// get a string that describes this task
//------------------------------------------------------------------------
  switch (Task) {
    case kUp:        strcpy(Str, "Up");            break;
    case kDown:      strcpy(Str, "Down");          break;
    case kStop:      strcpy(Str, "Stop");          break;
    case kGrow:      strcpy(Str, "Grow");          break;
    case kFade:      strcpy(Str, "Fade");          break;
    case kNoGrow:    strcpy(Str, "NoGrow");        break;
    case kStepUp:    strcpy(Str, "StepUp");        break;
    case kStepDown:  strcpy(Str, "StepDown");      break;
    case kStepGrow:  strcpy(Str, "StepGrow");      break;
    case kStepFade:  strcpy(Str, "StepFade");      break;
    default:         strcpy(Str, "Bug Alert!!!");  break;
  }
}


Bool_t ALambdaControl::IsTimeToClearAccumulator() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's time
// to start accumulating dU/dLambda from zero.
// (clear the accumulator on the 1st step)
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  switch (Task) {
    // for pmf blocks, clear accumulator on first step
    case kUp:      case kDown:    case kStop:
    case kGrow:    case kFade:    case kNoGrow:
      if (CurrStep == StartStep) {
        RetVal = kTrue;
      }
      break;
    // for mcti blocks, clear accumulator after each equilibration
    case kStepUp:  case kStepDown:  case kStepGrow:  case kStepFade:
      if ( (CurrStep-(StartStep+NumEquilSteps)) % 
             (NumEquilSteps+NumAccumSteps) == 1) {
        RetVal = kTrue;
      }
      break;
    // should never get here
    default:
      ASSERT(kFalse);
      break;
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsFirstStep() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's the
// first step of the control object.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  if (CurrStep == (StartStep+1)) {
    RetVal = kTrue;
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsTimeToPrint() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's time
// to print out restraint information.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if printing is required
  if (NumPrintSteps > 0) {
    // if number-of-steps from StartStep is an even multiple of NumPrintSteps
    // or if it's the last step of this pmf or mcti block,
    // then it's time to print
    if ( IsLastStep() || (((CurrStep-StartStep)%NumPrintSteps)==0) ) {
      RetVal = kTrue;
    }
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsEndOf_MCTI() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if this is
// the last time step of an mcti block.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if an MCTI block is currently active
  if ( (Task==kStepUp)   || (Task==kStepDown) ||
       (Task==kStepGrow) || (Task==kStepFade) ) {
    // if this is the last step
    if (CurrStep == GetLastStep()) {
      RetVal = kTrue;
    }
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsEndOf_MCTI_Step() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if this is
// the last time step of an mcti step.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if an MCTI block is currently active
  if ( (Task==kStepUp)   || (Task==kStepDown) ||
       (Task==kStepGrow) || (Task==kStepFade) ) {
    // if this is the last step of accumulation
    if (((CurrStep-StartStep)%(NumEquilSteps+NumAccumSteps))==0) {
      // then this is the last time step of this mcti step
      RetVal = kTrue;
    }
  }
  return(RetVal);
}


Bool_t ALambdaControl::IsTimeToPrint_dU_dLambda() {
//------------------------------------------------------------------------
// ASSUMING that this object is currently active, decide if it's time
// to print out dU/dLambda information.
//------------------------------------------------------------------------
  Bool_t  RetVal=kFalse;

  // if printing is required
  if (NumPrintSteps > 0) {
    // if number-of-steps from StartStep is an even multiple of NumPrintSteps
    // or if it's the last step of this pmf or mcti block,
    // then it might be time to print du/dLambda
    if ( IsLastStep() || (((CurrStep-StartStep)%NumPrintSteps)==0) ) {
      // for mcti blocks
      if ((Task==kStepUp)   || (Task==kStepDown) ||
          (Task==kStepGrow) || (Task==kStepFade)) {
        // it's only time to print if we're no longer equilibrating
        if ( ((CurrStep-StartStep-1) % (NumEquilSteps+NumAccumSteps))
             > NumEquilSteps) {
          RetVal = kTrue;
        }
      }
      // for other blocks (up, down, grow, fade) it's time
      else {
        RetVal = kTrue;
      }
    }
  }
  return(RetVal);
}


void ALambdaControl::PrintLambdaHeader(double dT) {
//----------------------------------------------------------------------------
// ASSUMING that this object is currently active, write out the header for
// a new lambda control object
//----------------------------------------------------------------------------
  // double Time;
  
  // calculate current time in femto-seconds
  // Time = (double)CurrStep * dT;

iout << "FreeEnergy: ";
#if !defined(_VERBOSE_PMF)
  iout << "nstep  time(ps)  ";
  iout << "    task  lambdaKf  lambdaRef     delta-G  #steps  n*{value  target |}" << std::endl;
  iout << "FreeEnergy: -----  --------  ";
  iout << "--------  --------  ---------  ----------  ------  -------------------" << std::endl;
  iout << endi;
#endif
}


void ALambdaControl::PrintHeader(double dT) {
//----------------------------------------------------------------------------
// ASSUMING that this object is currently active, write out the current time
//----------------------------------------------------------------------------
  double  Time;
  char    Str[100], Str2[100];

  // calculate current time in femto-seconds
  Time = (double)CurrStep * dT;

#if defined(_VERBOSE_PMF)
  iout << "FreeEnergy: " << std::endl << endi;
  iout << "FreeEnergy: ";
  iout << "Time Step = "  << CurrStep            <<    ",  ";
  iout << "Time = ";
  // write out time in ps
  iout << Time/1000.0     << " ps,  ";
  iout << "Lambda_Kf = "  << LambdaKf            <<    ",  ";
  iout << "Lambda_Ref = " << LambdaRef           <<     "  ";
  GetTaskStr(Str);
  iout << "(" << Str << ")";
  iout << std::endl << endi;
  iout << "FreeEnergy: ";
  iout << "------------------------------------------------";
  iout << "-------------------";
  iout << std::endl << endi;
#else
  sprintf(Str, "%5d", CurrStep);
  // write out time in ps
  sprintf(Str2, "%8.3f", Time/1000.0);
  iout << "FreeEnergy: ";
  iout << Str << "  " << Str2 << "  ";
  GetPaddedTaskStr(Str);
  iout << Str << "  ";
  sprintf(Str, "%8.5f", LambdaKf);
  iout << Str << "  ";
  sprintf(Str, "%9.5f", LambdaRef);
  iout << Str << "  ";
#endif
}


Bool_t ALambdaControl::IsActive() {
//------------------------------------------------------------------------
// determine if this object is currently active
//------------------------------------------------------------------------
  if ( (CurrStep>=StartStep) && (CurrStep<=GetLastStep()) ) {
    return(kTrue);
  }
  else {
    return(kFalse);
  }
}


double ALambdaControl::GetLambdaKf() {
//------------------------------------------------------------------------
// calculate LambdaKf for grow, fade, stepgrow, and stepfade.
// LambdaKf=1.0, for up, down, stepup, stepdown, and stop.
// for nogrow, LambdaKf is Lambda from the config file.
//------------------------------------------------------------------------
  if (IsActive()) {
    int N;
    switch (Task) {
      case kGrow:
        LambdaKf = (1.0*(CurrStep-StartStep)) / NumSteps;
        break;
      case kFade:
        LambdaKf = 1.0 - (1.0*(CurrStep-StartStep)) / NumSteps;
        break;
      case kStepGrow:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        LambdaKf = (1.0 * (N+1)) / NumRepeats;
        break;
      case kStepFade:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        LambdaKf = 1.0 - (1.0 * (N+1)) / NumRepeats;
        break;
      case kNoGrow:
        break;              // return prior setting of LambdaKf
      default:
        LambdaKf=1.0;
    }
  }
  else {
    LambdaKf=1.0;
  }
  return(LambdaKf);
}


double ALambdaControl::GetLambdaRef() {
//------------------------------------------------------------------------
// calculate LambdaRef for up, down, stepup, and stepdown.
// for stop, LambdaRef is Lambda from the config file.
// for grow, fade, stepgrow, stepfade, and nogrow,
//   LambdaRef is LambdaT from the config file.
//------------------------------------------------------------------------
  int     N;
  if (IsActive()) {
    switch (Task) {
      case kUp:
        LambdaRef = (1.0*(CurrStep-StartStep))/NumSteps;
        break;
      case kDown:
        LambdaRef = 1.0 - (1.0*(CurrStep-StartStep))/NumSteps;
        break;
      case kStepUp:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        LambdaRef = (1.0 * (N+1))/NumRepeats;
        break;
      case kStepDown:
        N = (int) ( (CurrStep-StartStep-1) / (NumEquilSteps+NumAccumSteps) );
        LambdaRef = 1.0 - (1.0 * (N+1))/NumRepeats;
      default: 
        break;             // return prior setting of LambdaRef
    }
  }
  else {
    LambdaRef=0.0;
  }
//
DebugM(1,"GetLambdaRef: Active "<<IsActive()<<" Task "<<Task<<kUp<<"\n");
DebugM(1,"GetLambdaRef: CurrStep "<<CurrStep<<" StartStep "<<StartStep<<" NumSteps " <<NumSteps<<" LambdaRef "<<LambdaRef<<"\n");
//
  return(LambdaRef);
}
