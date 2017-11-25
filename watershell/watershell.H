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
__EXTERN __CINT watershell_init_from_plugin (const __CCHAR * , const __CINT, \
                                             const __CCHAR *, const __CINT, \
                                             __CINT **);
__EXTERN __CINT watershell_dyna_from_plugin ( const long int, const __CFLOAT *, float *, __CFLOAT * );
__EXTERN   void watershell_done_from_plugin ();
