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

__EXTERN __CINT __NM(init_plugin)(const __CINT, const __CFLOAT *, const  __CFLOAT *, \
                                  const __CCHAR *, const __CINT, \
                                  const __CCHAR *, const __CINT, \
                                        __CINT **, const __CBOOL, const __CFLOAT *);
__EXTERN __CINT __NM(dyna_plugin)(const long int, const __CFLOAT *, double *, float *, const __CBOOL, \
                                        __CFLOAT *, __CINT **, const __CBOOL, const __CFLOAT * );
__EXTERN   void __NM(done_plugin)();

