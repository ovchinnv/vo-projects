#ifndef __BINDC
#define __BINDC
#endif

#include "InfoStream.h"
#include "NamdTypes.h"
#include "Node.h"
#include "Molecule.h"
#include "GlobalMaster.h"
#include "SimParameters.h"
#include "ReductionMgr.h"
#include "GlobalMasterSMCV.h"


GlobalMasterSMCV::GlobalMasterSMCV() {
   __CCHAR *inputfile = NULL, *logfile = NULL;
   __CINT ilen=0, llen=0, ierr=0;
   __CCHAR* deflogfile = "smcv.log" ;
   SimParameters *params = Node::Object()->simParameters;
   molecule=Node::Object()->molecule ;
   __CINT natoms=molecule->numAtoms;
//
   initialized=0;
   smcv_energy=0.;
   CkPrintf("# SMCV PLUGIN: Initializing ...\n");
//
// find input file
//
   inputfile = params->SMCVConfigFileName ;
   if (!inputfile) {
    CkPrintf("# SMCV PLUGIN: input file not specified\n");
    CkExit();
   } else {
    CkPrintf("# SMCV PLUGIN: input file is '"); 
    CkPrintf(inputfile);
    CkPrintf("'\n");
   }
//
// find log file
//
   logfile = params->SMCVLogFileName ;
   if (!logfile) {
    CkPrintf("# SMCV PLUGIN: log file not specified, will write to '");
    logfile=deflogfile;
    CkPrintf(logfile);
    CkPrintf("'\n");
   } else {
     CkPrintf("# SMCV PLUGIN: log file is '"); 
     CkPrintf(logfile);
     CkPrintf("'\n");
   }
//
   __CFLOAT *mass=NULL, *charge=NULL;
//
   mass =   (__CFLOAT *) calloc(molecule->numAtoms, sizeof(__CFLOAT));
   charge = (__CFLOAT *) calloc(molecule->numAtoms, sizeof(__CFLOAT));
//
   for (int i=0; i<molecule->numAtoms; i++){
    mass[i]=molecule->atommass(i);
    charge[i]=molecule->atomcharge(i);
   }
//
   ilen=strlen(inputfile);
   llen=strlen(logfile);
   ierr=smcv_init_from_namd(natoms, mass, charge, inputfile, ilen, logfile, llen, &atomlist) ;
// mass and charge no longer needed
   free(mass);
   free(charge);
//
   int atomid;
   if (atomlist!=NULL) { // atom indices provided; add them
    // first element gives list size:
    for ( int l = 0 ; l < atomlist[0]  ; ++l) { // first value in atomlist is the list length
     atomid = atomlist[l] - 1 ; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the SMCV plugin
     if ( atomid >= 0 && atomid < molecule->numAtoms ) {
      modifyRequestedAtoms().add(atomid); // note: atomlist will be deallocated in Fortran
     } else {
      CkPrintf("# SMCV PLUGIN: Atom ID ");
      CkPrintf("%d", atomid);
      CkPrintf("is out of range\n");
      CkExit();
     } // atomid valid
    } // over list
    r  =     (__CFLOAT *) calloc(3*molecule->numAtoms, sizeof(__CFLOAT));
    fr =     (__CFLOAT *) calloc(3*molecule->numAtoms, sizeof(__CFLOAT));
   } else {
    CkPrintf("# SMCV PLUGIN: No restraint atoms found (nothing to do)\n");
    CkExit();
   } // atomlist
//
   iteration = params->firstTimestep;
   reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
} // GlobalMasterSMCV


GlobalMasterSMCV::~GlobalMasterSMCV() {
    smcv_done_from_namd();
    delete reduction;
    if (r) free(r);
    if (fr) free(fr);
} // ~GlobalMasterSMCV


void GlobalMasterSMCV::calculate() {
 // load coordinates
  int atomid;
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();
  Vector p ;
  for ( x=r, y=x+molecule->numAtoms, z=y+molecule->numAtoms ; a_i != a_e; ++a_i, ++p_i ) {
    p = *p_i;
    x[*a_i] = p.x; // a_i dereferences to atomid ; p_i dereferences to position (instance of Vector)
    y[*a_i] = p.y; // however, in principle, access can be random (at the very least, nonconsecutive)
    z[*a_i] = p.z;
  }
  //call smcv
  __CINT ierr=smcv_dyna_from_namd(iteration++, r, fr, &smcv_energy);
  // apply forces
  if (atomlist!=NULL) { // atom indices provided; add them
    // first element gives list size:
   for ( int l = 0 ; l < atomlist[0]  ; ++l ) { // first value in atomlist is the list length
     atomid = atomlist[l] - 1 ; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the SMCV plugin
     modifyForcedAtoms().add(atomid);
     p.x=fx[atomid];
     p.y=fy[atomid];
     p.z=fz[atomid];
     modifyAppliedForces().add(p);
   } // for
  } // atomlist
} // calculate

