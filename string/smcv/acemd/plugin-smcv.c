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
 const __CCHAR* deflogfile = "smcv.log" ;
 for (int i=0; i<argc; i++){
  if ( !(strcmp(argkey[i], "input")) || !(strcmp(argkey[i], "inputfile"))) {
   inputfile=argval[i];
  } else if ( !(strcmp(argkey[i], "output")) || !(strcmp(argkey[i], "outputfile")) || !(strcmp(argkey[i], "log")) || !(strcmp(argkey[i], "logfile"))) {
   logfile=argval[i]; 
  }
 } // for
 if (inputfile == NULL) {
  printf("#SMCV PLUGIN: input file not specified (syntax : input <inputfile> log <logfile>)\n");
  exit(1);
 }
//
 if (logfile == NULL) {
  printf("#SMCV PLUGIN: log file not specified (syntax : input <inputfile> log <logfile>); uisng ");
  printf(deflogfile);
  printf("\n");
  logfile=deflogfile;
 }

 __CINT natoms=0;
 __CFLOAT *m=NULL, *q=NULL;

 natoms = s->natoms ;
 m = (__CFLOAT *) calloc(natoms, sizeof(__CFLOAT));
 q = (__CFLOAT *) calloc(natoms, sizeof(__CFLOAT));
 for (int i=0; i<natoms; i++){
  m[i]=s->mass[i];
  q[i]=s->charge[i];
 }

 smcv_init_from_acemd(natoms, m, q, inputfile, logfile) ;
//
 return 0;
}
