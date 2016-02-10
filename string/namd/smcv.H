/* C prototypes */
// C, not Fortran!
#ifdef int
#undef int
#endif
#ifdef float
#undef float
#endif

#ifdef __CPLUSPLUS
#define __EXTERN extern "C"
#else
#define __EXTERN extern
#endif
__EXTERN __CINT smcv_init_from_namd (const __CINT, const __CFLOAT *,const  __CFLOAT *, const __CCHAR * , const __CINT, const __CCHAR *, const __CINT,  __CINT **);
__EXTERN   void smcv_done_from_namd ();
__EXTERN __CINT smcv_dyna_from_namd ( const long int, const __CFLOAT *, float *, __CFLOAT *, __CINT ** );
