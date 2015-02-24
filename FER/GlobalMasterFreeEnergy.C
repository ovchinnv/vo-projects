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
// get lambdas from LambdaManager
// calculate gradients for each center-of-mass of each restraint,
// and apply the forces to the atoms involved
//-----------------------------------------------------------------

  double  LambdaKf, LambdaRef;
  double  Sum_dU_dLambdas;

  if (m_RestraintManager.GetNumRestraints() > 0) { // no need to worry if nothing is defined

   bool LambdasActive=(m_LambdaManager.GetLambdas(LambdaKf, LambdaRef));
   if (LambdasActive) m_RestraintManager.SetLambdas(LambdaKf, LambdaRef); // set new lamdas on RestraintManager
// fetch all positions (I actually think this should be done at the GlobalMaster level)
   fetchPositions();
// restraint manager calls individual restraints to add their forces
   double total_restraint_energy = m_RestraintManager.AddForces(*this);

   if (LambdasActive) {
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
    if (m_LambdaManager.IsTimeToPrint()) {
      m_LambdaManager.PrintHeader(simParams->dt);
      if (m_LambdaManager.IsTimeToPrint_dU_dLambda()) {
        m_RestraintManager.Print_dU_dLambda();
        if (m_RestraintManager.ThereIsAForcingRestraint()) {
          m_LambdaManager.Print_dU_dLambda_Summary(Sum_dU_dLambdas);
        }
      }
      m_RestraintManager.PrintRestraintInfo();
      if (m_LambdaManager.IsEndOf_MCTI()) {
        m_LambdaManager.Print_MCTI_Integration();
      }
    }
   }
// add restraint energy to MISC column
   reduction->item(REDUCTION_MISC_ENERGY) += total_restraint_energy;
   reduction->submit();
  } // numRestraints > 0
}


void GlobalMasterFreeEnergy::user_initialize() {
//-----------------------------------------------------------------
// read all the input from config
//-----------------------------------------------------------------

  iout << iINFO << "=================================================\n"; 
  iout << iINFO << "  FREE ENERGY CONFIG SCRIPT (COMMENTS REMOVED)\n"; 
  iout << iINFO << "=================================================\n"; 
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
  iout << iINFO << "=================================================\n"<<endi; 
//
  ReadInput(config, m_RestraintManager, m_LambdaManager, *this, simParams->dt);

  // exit if there aren't enough steps to complete all pmf & mcti blocks
  int Total = m_LambdaManager.GetTotalNumSteps();
  if (Total > simParams->N) {
    iout << "===================================================================" << std::endl;
    iout << "FreeEnergy: WARNING : Not enough steps to complete PMF/MCTI blocks" << std::endl;
    iout << "FreeEnergy:   Steps Needed   : "<< Total << std::endl;
    iout << "FreeEnergy:   Steps Available: "<< simParams->N << std::endl << endi;
    iout << "===================================================================" << std::endl;
// continue anyway : maybe user wants to accumulate stats for as long as possible 
//    NAMD_die("FreeEnergy: Fatal Run-Time Error");
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


GlobalMasterFreeEnergy::GlobalMasterFreeEnergy() : GlobalMaster() {
  DebugM(3,"Constructing GlobalMasterFreeEnergy\n");
  molecule = Node::Object()->molecule;
  simParams = Node::Object()->simParameters;
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  // now set up the free energy stuff
  initialize();
}

GlobalMasterFreeEnergy::~GlobalMasterFreeEnergy() {
  DebugM(3,"Destructing GlobalMasterFreeEnergy\n");
  delete config;
  delete reduction;
}


void GlobalMasterFreeEnergy::initialize() {
  DebugM(4,"Initializing master\n");

  // Get our script
  StringList *script = Node::Object()->configList->find("freeEnergyConfig");

  config = new char[1];
  config[0] = '\0';
 //
 // functionality to remove comments: (VO 2013)
  int comlen = 1;
  const char* (comment[comlen]);
  comment[0]="#";
//  comment[1]="//"; // keep this because may refer to file paths with extra slashes
//
  for ( ; script; script = script->next) {
    char* filescript=NULL;
    char* istr = script->data; // beginning
    char* eol = strstr(istr,"\n"); // check for newlines in script
//============================================= external file script ================================
    if (!eol) { // no newlines found above -- will assume this is a filename spec
      FILE *infile = fopen(script->data,"r");
      if ( ! infile ) {
	char errmsg[256];
	sprintf(errmsg,"Error trying to read file %s!\n",script->data);
	NAMD_die(errmsg);
      }
      fseek(infile,0,SEEK_END); // go to end of file
      size_t filesize = ftell(infile); // get position in the file (bytes)
      filescript = new char[filesize + 2];
      rewind(infile);
      fread(filescript,sizeof(char),filesize,infile);
      filescript[filesize++] = '\n';
      filescript[filesize]   = '\0';
      fclose(infile);
      istr=filescript;
    }
//============================================= embedded script =====================================
    while (char* eol = strstr(istr,"\n")) { // proceed until there are no more endlines
// check for comments
     char* estr = eol; // start assuming there are no comments
     for (int j = 0 ; j < comlen ; j++) {
      char* comm = strstr( istr, comment[j] );
      if (comm) { //found comment
       if ( ( estr - comm ) > 0 ) estr=comm; // comment occurs before current end of command
      }
     } // all comment marks
     int datalen = estr - istr;
     if (datalen>0) { // do not bother with blank or all-comment lines
      size_t add_len = datalen;
      size_t config_len = 0;
      config_len = strlen(config);
      char *new_config = new char[config_len + add_len + 2]; // include null
      strcpy(new_config,config);
      strncat(new_config,istr, datalen);
      strcat(new_config,"\n");
      delete [] config;
      config = new_config;
     }
     istr=eol+1;  //skip to end of line and skip newline char
    } // while
    if (filescript) delete[] filescript;
  } //for

  // iout << iDEBUG << "Free energy perturbation - initialize()\n" << endi; 
  user_initialize();
}

