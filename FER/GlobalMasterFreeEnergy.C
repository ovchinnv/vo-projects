/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include <string.h>
#include "InfoStream.h"
#include "NamdTypes.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"
#include "Vector.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"
#include "GlobalMaster.h"
#include "GlobalMasterFreeEnergy.h"
#include "FreeEnergyParse.h"

#include "Node.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "SimParameters.h"

#include <stdio.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

void GlobalMasterFreeEnergy::calculate() {
  /* zero out the forces */
  modifyForcedAtoms().resize(0);
  modifyAppliedForces().resize(0);

  /* XXX is this line needed at all? */
  modifyGroupForces().resize(getGroupMassEnd() - getGroupMassBegin());
  modifyGroupForces().setall(Vector(0,0,0));

//  iout << iDEBUG << "Free energy perturbation - calculate()\n" << endi; 
  m_LambdaManager.IncCurrStep();
//-----------------------------------------------------------------
// get lambdas from LambdaManager, inform RestraintManager.
// calculate gradients for each center-of-mass of each restraint,
// and apply the forces to the atoms involved
//-----------------------------------------------------------------
  double  LambdaKf, LambdaRef;
  double  Sum_dU_dLambdas;

  if (m_LambdaManager.GetLambdas(LambdaKf, LambdaRef)) {

    // stuff that's done every time step
    m_RestraintManager.SetLambdas(LambdaKf, LambdaRef);
    // fetch all positions (I actually think this should be done at the GlobalMaster level)
    fetchPositions();
    // restraint manager calls individual restraints to add their forces
    m_RestraintManager.AddForces(*this);

    if (m_LambdaManager.IsTimeToClearAccumulator()) {
      m_LambdaManager.ZeroAccumulator();
    }
    Sum_dU_dLambdas = m_RestraintManager.Sum_dU_dLambdas();
    m_LambdaManager.Accumulate(Sum_dU_dLambdas);

    // for integrating all the MCTI averages
    if (m_LambdaManager.IsEndOf_MCTI_Step()) {
      m_LambdaManager.Integrate_MCTI();
    }

    // stuff that's done when it's time to print
    if (m_LambdaManager.IsFirstStep()) {
      m_LambdaManager.PrintLambdaHeader(simParams->dt);
    }
    if (m_LambdaManager.IsTimeToPrint()) {
      m_LambdaManager.PrintHeader(simParams->dt);
      if (m_LambdaManager.IsTimeToPrint_dU_dLambda()) {
        m_RestraintManager.Print_dU_dLambda_Info();
        if (m_RestraintManager.ThereIsAForcingRestraint()) {
          m_LambdaManager.Print_dU_dLambda_Summary(Sum_dU_dLambdas);
        }
      }
      else {
        m_LambdaManager.PrintSomeSpaces();
      }
      m_RestraintManager.PrintEnergyInfo();
      m_RestraintManager.PrintRestraintInfo();
      if (m_LambdaManager.IsEndOf_MCTI()) {
        m_LambdaManager.Print_MCTI_Integration();
      }
    }
  }
}


void GlobalMasterFreeEnergy::user_initialize() {
//-----------------------------------------------------------------
// read all the input from config
//-----------------------------------------------------------------

  iout << iINFO << "  FREE ENERGY PERTURBATION CONFIG\n"; 
  iout << iINFO << "***********************************\n"; 
  int config_len = strlen(config);
  if ( config_len < 10000 ) {
    iout << config;
  } else {
    char *new_config = new char[10000 + 10];
    strncpy(new_config,config,10000);
    new_config[10000] = 0;
    strcat(new_config,"\n...\n");
    iout << new_config;
    delete [] new_config;
  }
  iout << iINFO << "***********************************\n" << endi; 

  ReadInput(config, m_RestraintManager, m_LambdaManager, *this, simParams->dt);

  // exit if there aren't enough steps to complete all pmf & mcti blocks
  int Total = m_LambdaManager.GetTotalNumSteps();
  if (Total > simParams->N) {
    iout << "FreeEnergy: Not enough steps to complete pfm & mcti blocks" << std::endl;
    iout << "FreeEnergy:   Num Steps Needed =    " << Total << std::endl;
    iout << "FreeEnergy:   Num Steps Requested = " << simParams->N << std::endl << endi;
    NAMD_die("FreeEnergy: Fatal Run-Time Error");
  }
}
//---------------------------------------------------------------------
int GlobalMasterFreeEnergy::getAtomID  (const char *segid, int resid, const char *aname){return molecule->get_atom_from_name(segid,resid,aname);}
int GlobalMasterFreeEnergy::getNumAtoms(const char* segid, int resid){return molecule->get_residue_size(segid,resid);} // 0 on error
int GlobalMasterFreeEnergy::getAtomID  (const char *segid, int resid, int index){return molecule->get_atom_from_index_in_residue(segid,resid,index);}
double GlobalMasterFreeEnergy::getMass (int atomid){return ((atomid < 0)||(atomid >= molecule->numAtoms) ? -1.: molecule->atommass(atomid));}//-1 on failure
int GlobalMasterFreeEnergy::getMoleculeSize() { return molecule->numAtoms; }; // total number of atoms ; VO 2013

int GlobalMasterFreeEnergy::requestAtom(int atomid)
{
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  modifyRequestedAtoms().add(atomid);
  return 0;  // success
}

// no need for this v (VO 2013)
/*int GlobalMasterFreeEnergy::addForce(int atomid, Force force)
{
  DebugM(2,"Forcing "<<atomid<<" with "<<force<<"\n");
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  modifyForcedAtoms().add(atomid);
  modifyAppliedForces().add(force);
  return 0;  // success
}
*/

GlobalMasterFreeEnergy::GlobalMasterFreeEnergy() : GlobalMaster() {
  DebugM(3,"Constructing GlobalMasterFreeEnergy\n");
  molecule = Node::Object()->molecule;
  simParams = Node::Object()->simParameters;
  // now set up the free energy stuff
  initialize();
}

GlobalMasterFreeEnergy::~GlobalMasterFreeEnergy() {
  DebugM(3,"Destructing GlobalMasterFreeEnergy\n");
  delete config;
}


void GlobalMasterFreeEnergy::initialize() {
  DebugM(4,"Initializing master\n");

  // Get our script
  StringList *script = Node::Object()->configList->find("freeEnergyConfig");

  config = new char[1];
  config[0] = '\0';

  for ( ; script; script = script->next) {
    if ( strstr(script->data,"\n") ) {
      size_t add_len = strlen(script->data);
      size_t config_len = 0;
      config_len = strlen(config);
      char *new_config = new char[config_len + add_len + 2];
      strcpy(new_config,config);
      strcat(new_config,script->data);
      strcat(new_config,"\n");  // just to be safe
      delete [] config;
      config = new_config;
    } else {
      FILE *infile = fopen(script->data,"r");
      if ( ! infile ) {
	char errmsg[256];
	sprintf(errmsg,"Error trying to read file %s!\n",script->data);
	NAMD_die(errmsg);
      }
      fseek(infile,0,SEEK_END);
      size_t add_len = ftell(infile);
      size_t config_len = 0;
      config_len = strlen(config);
      char *new_config = new char[config_len + add_len + 3];
      strcpy(new_config,config);
      delete [] config;
      config = new_config;
      new_config += config_len;
      rewind(infile);
      fread(new_config,sizeof(char),add_len,infile);
      new_config += add_len;
      new_config[0] = '\n';
      new_config[1] = '\0';
      fclose(infile);
    }
  }

  // iout << iDEBUG << "Free energy perturbation - initialize()\n" << endi; 
  user_initialize();
}

