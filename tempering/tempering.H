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
__EXTERN __CINT tempering_init_from_plugin ( const __CCHAR * , const __CINT, const __CCHAR *, const __CINT );
__EXTERN __CINT tempering_dyna_from_plugin ( const long int, const __CFLOAT *INPUT, const __CFLOAT, __CFLOAT *OUTPUT);
__EXTERN   void tempering_done_from_plugin ();
