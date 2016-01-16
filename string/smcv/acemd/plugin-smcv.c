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
 m = (__CFLOAT *) calloc(natoms, sizeof(__CFLOAT));
 q = (__CFLOAT *) calloc(natoms, sizeof(__CFLOAT));
 for (int i=0; i<natoms; i++){
  m[i]=s->mass[i];
  q[i]=s->charge[i];
 }
//
 ilen=strlen(inputfile);
 llen=strlen(logfile);
 ierr=smcv_init_from_acemd(natoms, m, q, inputfile, ilen, logfile, llen) ;
 free(m);
 free(q);
//
 return (ierr>0) ? ACEPLUG_ERR : ACEPLUG_OK ;
}
//==========================================================
aceplug_err_t aceplug_calcforces(struct aceplug_sim_t *s) {
 int n = s->natoms ;
 __CFLOAT *r = malloc( n * 3 * sizeof(__CFLOAT) ) ; //positions
 float *fr   = calloc( n * 3 , sizeof(float) ) ; //forces
 __CFLOAT smcv_energy;
 long int iteration =s->step ;
//
 if ( s-> plugin_load_positions() ) {return ACEPLUG_ERR;}
 if ( s-> plugin_load_forces() ) {return ACEPLUG_ERR;}
 for (int i=0, j=0 ; i< s->natoms ; i++) {
  r[j++] = s -> pos[i].x;
  r[j++] = s -> pos[i].y;
  r[j++] = s -> pos[i].z;
 }
//
 smcv_master_from_acemd(iteration,r,fr, &smcv_energy);
//
 for (int i=0, j=0 ; i< s->natoms ; i++) {
  s -> frc[i].x+=fr[j++];
  s -> frc[i].y+=fr[j++];
  s -> frc[i].z+=fr[j++];
 }
 free(r);
 free(fr);
 if ( s-> plugin_update_forces() ) { return ACEPLUG_ERR;}
 if ( s-> plugin_add_energy(ENERGY_EXTERNAL,smcv_energy) ) { return ACEPLUG_ERR;}
//
 return ACEPLUG_OK;
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

