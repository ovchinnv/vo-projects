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
__EXTERN __CINT sm_init_from_acemd (const __CINT, const __CFLOAT *,const  __CFLOAT *, \
                                    const __CCHAR * , const __CINT, \
                                    const __CCHAR *, const __CINT, \
                                    __CINT **, \
                                    const __CCHAR *, const __CINT);
__EXTERN __CINT sm_dyna_from_acemd ( const long int, const __CFLOAT *, float *, __CFLOAT *, __CINT ** );
__EXTERN   void sm_done_from_acemd ();
