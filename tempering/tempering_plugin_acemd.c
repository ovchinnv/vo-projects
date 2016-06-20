// Tempering Plugin for ACEMD
#include "tcl.h"
#include "aceplug.h"
#include "tempering.h"
#include "string.h"
#include <stdlib.h>
#include <math.h>

aceplug_err_t aceplug_init(struct aceplug_sim_t *s, int argc, char **argkey, char **argval) {
 __CCHAR *inputfile = NULL, *logfile = NULL;
 __CINT ilen=0, llen=0, smlen=0, ierr=0;
 __CCHAR *deflogfile = "tempering.log" ;
//
 printf("# TEMPERING PLUGIN: Initializing ...\n");
//
 int i;
 for (i=0; i<argc; i++){
  if ( !(strcmp(argkey[i], "input")) || !(strcmp(argkey[i], "inputfile")) || !(strcmp(argkey[i], "INPUT")) || !(strcmp(argkey[i], "INPUTFILE")) ) {
   inputfile=argval[i];
  } else if ( !(strcmp(argkey[i], "output")) || !(strcmp(argkey[i], "outputfile")) || !(strcmp(argkey[i], "log")) || !(strcmp(argkey[i], "logfile")) ||\
              !(strcmp(argkey[i], "OUTPUT")) || !(strcmp(argkey[i], "OUTPUTFILE")) || !(strcmp(argkey[i], "LOG")) || !(strcmp(argkey[i], "LOGFILE")) ) {
   logfile=argval[i];
  }
 } // for
 if (inputfile == NULL) {
  printf("# TEMPERING PLUGIN: input file not specified (syntax : input <inputfile>)\n");
  exit(1);
 } else {
  printf("# TEMPERING PLUGIN: input file is '"); 
  printf(inputfile);
  printf("'\n");
 }
//
 if (logfile == NULL) {
  printf("# TEMPERING PLUGIN: log file not specified (syntax : log <logfile>); using '");
  logfile=deflogfile;
  printf(logfile);
  printf("'\n");
 } else {
  printf("# TEMPERING PLUGIN: log file is '"); 
  printf(logfile);
  printf("'\n");
 }
//
 ilen=strlen(inputfile);
 llen=strlen(logfile);
 ierr=tempering_init_from_acemd(inputfile, ilen, logfile, llen);
//
 return (ierr>0) ? ACEPLUG_ERR : ACEPLUG_OK ;
}
//==========================================================
aceplug_err_t aceplug_endstep(struct aceplug_sim_t *s) {

 double edata[ENERGY_MAX] ;      // energy data array
// double *ke = edata + ENERGY_KE; // pointer to kinetic energy value
//double *inst_temperature = edata + TEMPERATURE; // pointer to the value of kinetic temperature
 double *pe = edata + ENERGY_PE; // pointer to potential energy value
 double current_temperature ;    // current temperature of thermostat
 double new_temperature ;        // new temperature computed by tempering code
 double velocity_scale ;         // velocity scaling factor
 long int iteration ;
 __CINT ierr ; // error code returned from FORTRAN side
// get thermostat temp
 if ( s-> plugin_get_temperature(&current_temperature) ) {return ACEPLUG_ERR;}
#ifdef __DEBUG
 printf("# TEMPERING PLUGIN: got thermostat temperature :"); 
 printf("%12.5f\n", current_temperature);
#endif

// get energy
 if ( s-> plugin_get_energy_temp(edata) ) {return ACEPLUG_ERR;}
#ifdef __DEBUG
 printf("# TEMPERING PLUGIN: got energy :"); 
 printf("%12.5f\n", *pe);
#endif

// obtain new temperature from tempering routine
 ierr=tempering_dyna_from_acemd(iteration, pe, current_temperature, &new_temperature);
#ifdef __DEBUG
 printf("# TEMPERING PLUGIN: got new temperature :"); 
 printf("%12.5f\n", new_temperature);
#endif
// rescale velocities
 velocity_scale = sqrt(new_temperature/current_temperature) ;
 s->plugin_scale_velocities(velocity_scale);
// set new thermostat temperature
 if ( s-> plugin_set_temperature(new_temperature) ) { return ACEPLUG_ERR;}

 return (ierr>0) ? ACEPLUG_ERR : ACEPLUG_OK ;
}
//==========================================================
aceplug_err_t aceplug_terminate(struct aceplug_sim_t *s) {
 printf("# TEMPERING PLUGIN: Finalizing...\n");
 tempering_done_from_acemd();
 return ACEPLUG_OK;
}
