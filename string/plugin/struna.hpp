/* C prototypes for plugin*/
// C, not Fortran!
extern "C" int sm_init_plugin(const int, const double *, const double *, const char * , const int, const char *, const int, int **, const _Bool, const double *);
extern "C" int sm_dyna_plugin( const long int, const double *, double *, float *, const _Bool, double *, int **, const _Bool, const double * );
extern "C" void sm_done_plugin();
