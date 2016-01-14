// SMCV plugin for ACEMD
#include "tcl.h"
#include "aceplug.h"
#include "string.h"
#include "smcv.h"
#include <stdlib.h>

aceplug_err_t aceplug_init(struct aceplug_sim_t *s, int argc, char **argkey, char **argval) {
 char* inputfile = NULL ;
 for (int i=0; i<argc; i++){
  if ( !(strcmp(argkey[i], "input")) || !(strcmp(argkey[i], "inputfile"))) {
   inputfile=argval[i];
   break; // take first argument
  } //if
 } // for
 if (inputfile == NULL) {
  printf("#SMCV PLUGIN: input file not specified (syntax : [input|inputfile] <inputfile>)\n");
  exit(1);
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

 smcv_init_from_acemd(natoms, m, q, inputfile) ;
//
 return 0;
}
