#ifndef __BINDC
#define __BINDC
#endif

#include "InfoStream.h"
#include "NamdTypes.h"
#include "Node.h"
#include "Molecule.h"
#include "GlobalMaster.h"
#include "SimParameters.h"
//#include "ReductionMgr.h"
#include "GlobalMasterSMCV.h"
//#define DEBUGM
#include "Debug.h"

GlobalMasterSMCV::GlobalMasterSMCV() {
   __CCHAR *inputfile = NULL, *logfile = NULL;
   __CINT ilen=0, llen=0, ierr=0;
   __CCHAR* deflogfile = "smcv.log" ;
   params = Node::Object()->simParameters;
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
    for ( int l = 0 ; l++ < atomlist[0]  ; ) { // first value in atomlist is the list length
     atomid = atomlist[l] - 1 ; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the SMCV plugin
     if ( atomid >= 0 && atomid < molecule->numAtoms ) {
      DebugM(1,"Adding atom "<<atomid<<"\n");
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
//   reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
} // GlobalMasterSMCV


GlobalMasterSMCV::~GlobalMasterSMCV() {
 params=NULL;
 destroy();
} // ~GlobalMasterSMCV


void GlobalMasterSMCV::calculate() {
 // load coordinates
  int atomid;
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();
  Vector p ;
  for ( ; a_i != a_e; ++a_i, ++p_i ) { // this loop is nly over atoms that have been requested
    atomid = (*a_i);
#ifdef DEBUGM
    CkPrintf("%d\n",atomid);
#endif
    atomid*=3; // index into coordinate array
    p = *p_i;
    r[atomid++] = p.x; // a_i dereferences to atomid ; p_i dereferences to position (instance of Vector)
    r[atomid++] = p.y; // however, in principle, access can be random (at the very least, nonconsecutive)
    r[atomid]   = p.z;
  }
  //call smcv
  __CINT ierr=smcv_dyna_from_namd(iteration++, r, fr, &smcv_energy);

#ifdef DEBUGM
// check that forces are received correctly
  for (int i=0;i<3*molecule->numAtoms;i++){
   CkPrintf("%15.10f %15.10f\n",r[i],fr[i]);
  }
#endif
  // apply forces
  modifyAppliedForces().resize(0);
  modifyForcedAtoms().resize(0);
  //
  if (atomlist!=NULL) { // atom indices provided; add them
    // first element gives list size:
   for ( int l = 0 ; l++ < atomlist[0]  ; ) { // first value in atomlist is the list length
     atomid = atomlist[l] - 1; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the SMCV plugin
     modifyForcedAtoms().add(atomid);
     atomid*=3 ; // index into force array
     p.x=fr[atomid++];
     p.y=fr[atomid++];
     p.z=fr[atomid];
     DebugM(1,"Adding force "<< p << " to atom "<<atomid<<"\n");
     modifyAppliedForces().add(p);
   } // for
  } // atomlist
  // destructor is never called, so finalize after run length exceeded
  if (iteration>params->N) destroy();
} // calculate

void GlobalMasterSMCV::destroy(){
  CkPrintf("# SMCV PLUGIN: Finalizing...\n");
  smcv_done_from_namd();
//    delete reduction;
  if (r) free(r);
  if (fr) free(fr);
}
