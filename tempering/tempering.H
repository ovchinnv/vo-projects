/* C prototypes */
// C, not Fortran!
#ifdef int
#undef int
#endif

#ifdef __CPLUSPLUS
#define __EXTERN extern "C"
#else
#define __EXTERN extern
#endif
__EXTERN __CINT tempering_init_from_acemd ( const __CCHAR * , const __CINT, const __CCHAR *, const __CINT );
__EXTERN __CINT tempering_dyna_from_acemd ( const long int, const __CFLOAT*, const __CFLOAT, __CFLOAT * );
__EXTERN   void tempering_done_from_acemd ();
