/* C prototypes */
// C, not Fortran!
extern int smcv_init_from_acemd (const int, const double *,const double *, const char * , const int, const char *, const int, int **);
extern void smcv_done_from_acemd ();
extern void smcv_dyna_from_acemd ( const long int, const double *, float *, double *, int ** );
