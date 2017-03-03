// plugin for ACEMD
// C, not Fortran!
#ifdef int
#undef int
#endif

#include "tcl.h"
#include "aceplug.h"
#include "plugin.h"
#include "string.h"
#include <stdlib.h>

// structure for plugin private data :
typedef struct {
 __CINT *atomlist ; // to maintain a list of atoms that are needed by plugin; other coordinates to be omitted
 __CFLOAT *r ; // coordinate array
 float *fr ; // forces array
} pdata ;

aceplug_err_t aceplug_init(struct aceplug_sim_t *s, int argc, char **argkey, char **argval) {
 __CCHAR *inputfile = NULL, *logfile = NULL;
 __CINT ilen=0, llen=0, ierr=0;
 __CCHAR *deflogfile = "plugin.log" ;
// for plugin private data :
 pdata *private_data = malloc(sizeof(pdata)) ;
//
 printf("# MASTER PLUGIN: Initializing ...\n");
//
 int i;
 for (i=0; i<argc; i++){
  if ( !(strcmp(argkey[i], "input")) || !(strcmp(argkey[i], "inputfile")) || !(strcmp(argkey[i], "INPUT")) || !(strcmp(argkey[i], "INPUTFILE")) ) {
   inputfile=argval[i];
  } else if ( !(strcmp(argkey[i], "output")) || !(strcmp(argkey[i], "outputfile")) || !(strcmp(argkey[i], "log")) || !(strcmp(argkey[i], "logfile")) ||\
              !(strcmp(argkey[i], "OUTPUT")) || !(strcmp(argkey[i], "OUTPUTFILE")) || !(strcmp(argkey[i], "LOG")) || !(strcmp(argkey[i], "LOGFILE")) ) {
   logfile=argval[i];
  } //if
 } // for
 if (inputfile == NULL) {
  printf("# MASTER PLUGIN: input file not specified (syntax : input <inputfile>)\n");
  exit(1);
 } else {
  printf("# MASTER PLUGIN: input file is '"); 
  printf(inputfile);
  printf("'\n");
 }
//
 if (logfile == NULL) {
  printf("# MASTER PLUGIN: log file not specified (syntax : log <logfile>); using '");
  logfile=deflogfile;
  printf(logfile);
  printf("'\n");
 } else {
  printf("# MASTER PLUGIN: log file is '"); 
  printf(logfile);
  printf("'\n");
 }
//
 __CINT natoms=0;
 natoms = s->natoms ;
//
 private_data->r  = (__CFLOAT *) calloc(3 * natoms, sizeof(__CFLOAT));
 private_data->fr = (float *) calloc(3 * natoms, sizeof(float));
 private_data->atomlist = NULL ; // NOTE that atomlist will be allocated in the FORTRAN code (which knows how to do it)
 s->privdata = private_data; // initialize private data pointer
//
 ilen=strlen(inputfile);
 llen=strlen(logfile);
 ierr=__NM(init_from_plugin)(inputfile, ilen, logfile, llen, &(private_data->atomlist));
//
 return (ierr>0) ? ACEPLUG_ERR : ACEPLUG_OK ;
}
//==========================================================
aceplug_err_t aceplug_calcforces(struct aceplug_sim_t *s) {
 int n = s->natoms ;
 int i, j;
 __CINT ierr;
// for accessing the plugin private data :
 pdata *private_data = s->privdata ;
 __CINT *atomlist = private_data->atomlist ;
 __CFLOAT *r = private_data->r ;
 float *fr = private_data->fr ;
//
 __CINT *aptr ;
 __CFLOAT *rptr;
 float *fptr;
 __CFLOAT __NM(energy);
 long int iteration =s->step ;
//
//
 if ( s-> plugin_load_positions() ) {return ACEPLUG_ERR;}
 if ( s-> plugin_load_forces() ) {return ACEPLUG_ERR;}
//
 if (atomlist==NULL) { // atomlist is not defined; therefore, provide all coords
  for (i=0, rptr=r ; i< s->natoms ; i++) {
   *(rptr++) = s -> pos[i].x;
   *(rptr++) = s -> pos[i].y;
   *(rptr++) = s -> pos[i].z;
  }
// compute plugin forces
//  ierr=__NM(dyna_from_acemd) (iteration, r, fr, &__NM(energy)); // might return valid atomlist
  ierr=__NM(dyna_from_plugin) (iteration, r, NULL, fr, &__NM(energy), 1); // single-prec;
  private_data->atomlist=atomlist ;
// apply plugin forces
  if (atomlist!=NULL) { // atom indices provided; store them as private data and use for adding forces
   for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
    j=*aptr - 1; // for zero offset (e.g. first coordinate lives in r[0]
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
 } else { // atomlist not null : loop over only the desired indices
  for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
   j=*aptr - 1;
   rptr=r + 3*j ;
   *(rptr++) = s -> pos[j].x;
   *(rptr++) = s -> pos[j].y;
   *(rptr)   = s -> pos[j].z;
  }
//
//  ierr=__NM(dyna_from_acemd)(iteration, r, fr, &__NM(energy));
  ierr=__NM(dyna_from_plugin)(iteration, r, NULL, fr, &__NM(energy), 1); // single-prec flag
//
  for (aptr=atomlist+1 ; aptr<atomlist + 1 + (*atomlist) ; aptr++) { // iterate until atomlist points to the last index
   j=*aptr - 1;
   fptr=fr + 3*j ;
   s -> frc[j].x+= *(fptr++);
   s -> frc[j].y+= *(fptr++);
   s -> frc[j].z+= *(fptr);
  }
//
 } // atomlist == NULL
//
 if ( s-> plugin_update_forces() ) { return ACEPLUG_ERR;}
 if ( s-> plugin_add_energy(ENERGY_EXTERNAL, __NM(energy)) ) { return ACEPLUG_ERR;}
// we know that energy is computed, but it does not seem to make into the ACEMD energy output
//printf("%12.5f\n", __NM(energy));
//
 return (ierr>0) ? ACEPLUG_ERR : ACEPLUG_OK ;
}
//==========================================================
aceplug_err_t aceplug_endstep(struct aceplug_sim_t *s) {
 return ACEPLUG_OK;
}
//==========================================================
aceplug_err_t aceplug_terminate(struct aceplug_sim_t *s) {
 pdata *private_data = s->privdata ;
 printf("# MASTER PLUGIN: Finalizing...\n");
 __NM(done_from_plugin)();
 free(private_data->r);
 free(private_data->fr);
 private_data->atomlist=NULL ; // this memory is managed by FORTRAN
 free(private_data);
 return ACEPLUG_OK;
}
