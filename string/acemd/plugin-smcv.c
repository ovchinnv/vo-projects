// SMCV plugin for ACEMD
// C, not Fortran!
#ifdef int
#undef int
#endif

#include "tcl.h"
#include "aceplug.h"
#include "string.h"
#include "smcv.h"
#include <stdlib.h>

aceplug_err_t aceplug_init(struct aceplug_sim_t *s, int argc, char **argkey, char **argval) {
 __CCHAR *inputfile = NULL, *logfile = NULL;
 __CINT ilen=0, llen=0, ierr=0;
 __CCHAR* deflogfile = "smcv.log" ;
 __CINT * atomlist = NULL ; // to maintain a list of atoms that are needed by plugin; other coordinates to be omitted
//
 printf("# SMCV PLUGIN: Initializing ...\n");
//
 for (int i=0; i<argc; i++){
  if ( !(strcmp(argkey[i], "input")) || !(strcmp(argkey[i], "inputfile"))) {
   inputfile=argval[i];
  } else if ( !(strcmp(argkey[i], "output")) || !(strcmp(argkey[i], "outputfile")) || !(strcmp(argkey[i], "log")) || !(strcmp(argkey[i], "logfile"))) {
   logfile=argval[i];
  }
 } // for
 if (inputfile == NULL) {
  printf("# SMCV PLUGIN: input file not specified (syntax : input <inputfile> log <logfile>)\n");
  exit(1);
 } else {
  printf("# SMCV PLUGIN: input file is '"); 
  printf(inputfile);
  printf("'\n");
 }
//
 if (logfile == NULL) {
  printf("# SMCV PLUGIN: log file not specified (syntax : input <inputfile> log <logfile>); using '");
  logfile=deflogfile;
  printf(logfile);
  printf("'\n");
 } else {
  printf("# SMCV PLUGIN: log file is '"); 
  printf(logfile);
  printf("'\n");
 }

 __CINT natoms=0;
 __CFLOAT *m=NULL, *q=NULL;

 natoms = s->natoms ;
//
 s->privdata = NULL; // initialize private pointer
//
 m = (__CFLOAT *) calloc(natoms, sizeof(__CFLOAT));
 q = (__CFLOAT *) calloc(natoms, sizeof(__CFLOAT));
 for (int i=0; i<natoms; i++){
  m[i]=s->mass[i];
  q[i]=s->charge[i];
 }
//
 ilen=strlen(inputfile);
 llen=strlen(logfile);
 ierr=smcv_init_from_acemd(natoms, m, q, inputfile, ilen, logfile, llen, &atomlist) ;
 if (atomlist!=NULL) { // atom indices provided; store them as private data
  s->privdata = atomlist ;
 }
 free(m);
 free(q);
//
 return (ierr>0) ? ACEPLUG_ERR : ACEPLUG_OK ;
}
//==========================================================
aceplug_err_t aceplug_calcforces(struct aceplug_sim_t *s) {
 int n = s->natoms ;
 int i, j;
 __CINT ierr;
 __CINT * atomlist = NULL, *m ; // list of atom indices required by plugin
 __CFLOAT *r = malloc( n * 3 * sizeof(__CFLOAT) ) ; //positions
 float *fr   = calloc( n * 3 , sizeof(float) ) ; //forces
 __CFLOAT *rptr;
 float *fptr;
 __CFLOAT smcv_energy;
 long int iteration =s->step ;
//
 if ( s-> plugin_load_positions() ) {return ACEPLUG_ERR;}
 if ( s-> plugin_load_forces() ) {return ACEPLUG_ERR;}
//
 if (s->privdata==NULL) { // provide all coords
  for (i=0, rptr=r ; i< s->natoms ; i++) {
   *(rptr++) = s -> pos[i].x;
   *(rptr++) = s -> pos[i].y;
   *(rptr++) = s -> pos[i].z;
  }
// compute plugin forces
  ierr=smcv_dyna_from_acemd(iteration, r, fr, &smcv_energy, &atomlist); // might return valid atomlist
// apply plugin forces
  if (atomlist!=NULL) { // atom indices provided; store them as private data and use for adding forces
   s->privdata = atomlist++ ; // point to atomlist and go to first atom index
   for (m=atomlist+atomlist[-1] ; atomlist<m ; atomlist++) { // iterate until atomlist points to the last index
    j=*atomlist ;
    fptr=fr + 3*j ;
    s -> frc[j].x+= *(fptr++);
    s -> frc[j].y+= *(fptr++);
    s -> frc[j].z+= *(fptr);
   }
  } else { // no atomlist provided; loop over all atoms
   for (i=0, fptr=fr ; i< s->natoms ; i++) {
    s -> frc[i].x+= *(fptr++);
    s -> frc[i].y+= *(fptr++);
    s -> frc[i].z+= *(fptr++);
   }
  } // atomlist
 } else { // privdata not null : loop over only the desired indices
  atomlist = (__CINT*)s->privdata+1 ;// first entry skipped; it is the number of atoms (see line below)
  for (m=atomlist + atomlist[-1]  ; atomlist<m ; atomlist++) { // iterate until atomlist points to the last index
   j=*atomlist ;
   rptr=r + 3*j ;
   *(rptr++) = s -> pos[j].x;
   *(rptr++) = s -> pos[j].y;
   *(rptr++) = s -> pos[j].z;
  }
//
  ierr=smcv_dyna_from_acemd(iteration, r, fr, &smcv_energy, &atomlist); // atomlist should not be modified in this call
//
  atomlist = (__CINT*)s->privdata+1 ;// first entry skipped; it is the number of atoms (see line below)
  for (m=atomlist+atomlist[-1] ; atomlist<m ; atomlist++) { // iterate until atomlist points to the last index
   j=*atomlist ;
   fptr=fr + 3*j ;
   s -> frc[j].x+= *(fptr++);
   s -> frc[j].y+= *(fptr++);
   s -> frc[j].z+= *(fptr);
  }

 } // s-> privdata
//
 free(r);
 free(fr);
 if ( s-> plugin_update_forces() ) { return ACEPLUG_ERR;}
 if ( s-> plugin_add_energy(ENERGY_EXTERNAL,smcv_energy) ) { return ACEPLUG_ERR;}
//
 return (ierr>0) ? ACEPLUG_ERR : ACEPLUG_OK ;
}
//==========================================================
aceplug_err_t aceplug_endstep(struct aceplug_sim_t *s) {
 return ACEPLUG_OK;
}
//==========================================================
aceplug_err_t aceplug_terminate(struct aceplug_sim_t *s) {
 printf("# SMCV PLUGIN: Finalizing...\n");
 smcv_done_from_acemd();
 return ACEPLUG_OK;
}

