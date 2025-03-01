#ifndef __BINDC
#define __BINDC
#endif

#include "InfoStream.h"
#include "NamdTypes.h"
#include "Node.h"
#include "Molecule.h"
#include "Lattice.h"
#include "GlobalMaster.h"
#include "SimParameters.h"
//#include "ReductionMgr.h"
#include <stdbool.h>
#include "GlobalMasterSM.h"
//#define DEBUGM
#include "Debug.h"

GlobalMasterSM::GlobalMasterSM() {
   __CCHAR *inputfile = NULL, *logfile = NULL;
   __CINT ilen=0, llen=0, ierr=0;
   const __CCHAR* deflogfile = "struna.log" ;
   params = Node::Object()->simParameters;
   molecule=Node::Object()->molecule ;
   __CINT natoms=molecule->numAtoms;
//
   initialized=0;
   sm_energy=0.;
   CkPrintf("# STRUNA PLUGIN: Initializing ...\n");
//
// find input file
//
   inputfile = params->SMConfigFileName ;
   if (!inputfile) {
    CkPrintf("# STRUNA PLUGIN: input file not specified (syntax : StrunaConfigFile <inputfile>)\n");
    CkExit();
   } else {
    CkPrintf("# STRUNA PLUGIN: input file is '"); 
    CkPrintf(inputfile);
    CkPrintf("'\n");
   }
//
// find log file
//
   logfile = params->SMLogFileName ;
   if (!logfile) {
    CkPrintf("# STRUNA PLUGIN: log file not specified (syntax : StrunaLogFile <logfile>), will write to '");
    logfile=(__CCHAR *) deflogfile;
    CkPrintf(logfile);
    CkPrintf("'\n");
   } else {
     CkPrintf("# STRUNA PLUGIN: log file is '"); 
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
// PBC
//   CkPrintf("# STRUNA PLUGIN: checking if simulation uses periodic boundary conditions\n"); 
   Lattice &lattice = params->lattice;
   usesPeriodic = ( lattice.a_p() && lattice.b_p() && lattice.c_p() ) ;
   if (usesPeriodic) {
    CkPrintf("# STRUNA PLUGIN: Simulation uses periodic boundary conditions\n"); 
    Vector const a = lattice.a() ;
    Vector const b = lattice.b() ;
    Vector const c = lattice.c() ;
    Position const o = lattice.origin() ;
    box[0]=a.x;
    box[1]=a.y;
    box[2]=a.z;
    box[3]=b.x;
    box[4]=b.y;
    box[5]=b.z;
    box[6]=c.x;
    box[7]=c.y;
    box[8]=c.z;
    box[9]=o.x;
    box[10]=o.y;
    box[11]=o.z;
   } else {
    for ( int i=0 ; i < 12 ; i++ ) { box[i]=0.0 ; }
   }
// initialize struna
   ierr=sm_init_plugin(natoms, mass, charge, inputfile, ilen, logfile, llen, &atomlist, usesPeriodic, box);
// mass and charge no longer needed
   free(mass);
   free(charge);
//
   int atomid;
   if (atomlist!=NULL) { // atom indices provided; add them
    // first element gives list size:
    for ( int l = 0 ; l++ < atomlist[0]  ; ) { // first value in atomlist is the list length
     atomid = atomlist[l] - 1 ; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the SM plugin
     if ( atomid >= 0 && atomid < molecule->numAtoms ) {
      DebugM(1,"Adding atom "<<atomid<<"\n");
      modifyRequestedAtoms().add(atomid); // note: atomlist will be deallocated in Fortran
     } else {
      CkPrintf("# STRUNA PLUGIN: Atom ID ");
      CkPrintf("%d", atomid);
      CkPrintf("is out of range\n");
      CkExit();
     } // atomid valid
    } // over list
   } else { // no indices provided, assume that initialization is deferred until dynamics; request all coordinates
    for ( atomid = 0 ; atomid < molecule->numAtoms ; atomid++) {
      DebugM(1,"Adding atom "<<atomid<<"\n");
      modifyRequestedAtoms().add(atomid);
    }
   } // atomlist
   r  = (__CFLOAT *) calloc(3*molecule->numAtoms, sizeof(__CFLOAT));
   fr = (__CFLOAT *) calloc(3*molecule->numAtoms, sizeof(__CFLOAT));
//
//   iteration = params->firstTimestep; // not sure we need to keep a local counter
//   reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
} // GlobalMasterSM


GlobalMasterSM::~GlobalMasterSM() {
 params=NULL;
 destroy();
} // ~GlobalMasterSM


void GlobalMasterSM::calculate() {
 // load coordinates
  int atomid;
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();
  Vector p ;
  for ( ; a_i != a_e; ++a_i, ++p_i ) { // this loop is only over atoms that have been requested
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
//
#ifdef DEBUGM
// check that forces are received correctly
  for (int i=0;i<3*molecule->numAtoms;i++){
   CkPrintf("%15.10f %15.10f\n",r[i],fr[i]);
  }
#endif
  if (usesPeriodic) {
   Vector const a = lattice->a() ;
   Vector const b = lattice->b() ;
   Vector const c = lattice->c() ;
   Position const o = lattice->origin() ;
   box[0]=a.x;
   box[1]=a.y;
   box[2]=a.z;
   box[3]=b.x;
   box[4]=b.y;
   box[5]=b.z;
   box[6]=c.x;
   box[7]=c.y;
   box[8]=c.z;
   box[9]=o.x;
   box[10]=o.y;
   box[11]=o.z;
  }
  // apply forces
  modifyAppliedForces().resize(0);
  modifyForcedAtoms().resize(0);
  //
  if (atomlist!=NULL) { // atom indices already provided; compute and add forces
    //call sm
    __CINT ierr=sm_dyna_plugin(step, r, fr, NULL, 0, &sm_energy, &atomlist, usesPeriodic, box);
    // first element gives list size:
   for ( int l = 0 ; l++ < atomlist[0]  ; ) { // first value in atomlist is the list length
     atomid = atomlist[l] - 1; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the SM plugin
     modifyForcedAtoms().add(atomid);
     atomid*=3 ; // index into force array
     p.x=fr[atomid++];
     p.y=fr[atomid++];
     p.z=fr[atomid];
     DebugM(1,"Adding force "<< p << " to atom "<<atomid<<"\n");
     modifyAppliedForces().add(p);
   } // for
  } else { // atomlist is NULL, but we expect it to be populated by the plugin at the first dynamics step
    //call sm
    __CINT ierr=sm_dyna_plugin(step, r, fr, NULL, 0, &sm_energy, &atomlist, usesPeriodic, box);
    //
    if (atomlist!=NULL) { // atom indices provided; add them
     modifyRequestedAtoms().resize(0);
    // first element gives list size:
     for ( int l = 0 ; l++ < atomlist[0]  ; ) { // first value in atomlist is the list length
      atomid = atomlist[l] - 1 ; // subtract one because atom indices are offset from 0 in NAMD, but from 1 in the SM plugin
      if ( atomid >= 0 && atomid < molecule->numAtoms ) {
       DebugM(1,"Adding atom "<<atomid<<"\n");
       modifyRequestedAtoms().add(atomid); // note: atomlist will be deallocated in Fortran
       // now add forces :
       modifyForcedAtoms().add(atomid);
       atomid*=3 ; // index into force array
       p.x=fr[atomid++];
       p.y=fr[atomid++];
       p.z=fr[atomid];
       DebugM(1,"Adding force "<< p << " to atom "<<atomid<<"\n");
       modifyAppliedForces().add(p);
      } else {
       CkPrintf("# STRUNA PLUGIN: Atom ID ");
       CkPrintf("%d", atomid);
       CkPrintf("is out of range\n");
       CkExit();
      } // atomid valid
     } // over list
   } else { // atom list still undefined
    CkPrintf("# STRUNA PLUGIN: No restraint atoms found (nothing to do)\n");
    CkExit(); 
   } // atomlist
  }
  // destructor is not always called, so finalize after run length exceeded
  if (step>params->N) destroy(); // this is problematic for multiple runs per script, e.g. minimize, then run MD
} // calculate

void GlobalMasterSM::destroy(){
  CkPrintf("# STRUNA PLUGIN: Finalizing...\n");
  sm_done_plugin();
// Keep in case there are multiple runs
//  if (r) free(r);
//  if (fr) free(fr);
}
