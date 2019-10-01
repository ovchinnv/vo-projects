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

__EXTERN __CINT __NM(init_plugin)(const __CINT, const __CBOOL, const double *, const float *, const double *, const float *, \
                                  const __CCHAR *, const __CINT, \
                                  const __CCHAR *, const __CINT, \
                                        __CINT **, const __CBOOL, const double *, const float *);
__EXTERN __CINT __NM(dyna_plugin)(const long int, const __CBOOL, const double *, const float *, double *, float *, \
                                        double *, float *, __CINT **, const __CBOOL, const double *, const float * );
__EXTERN   void __NM(done_plugin)();

