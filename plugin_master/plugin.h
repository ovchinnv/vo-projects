/* C prototypes */
// C, not Fortran!
extern int master_init_from_plugin (const char * , const int, const char *, const int, int **);
extern int master_dyna_from_plugin ( const long int, const double *, double *, float *, double *, const _Bool );
extern void master_done_from_plugin ();
