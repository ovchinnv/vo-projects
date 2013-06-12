/* C prototypes */
#ifdef __CPLUSPLUS
#define __EXTERN "extern C"
//#elif defined(__SWIG)
//#define __EXTERN extern
#else
#define __EXTERN extern
#endif
//best-fit routines
__EXTERN void RMSBestFit        (__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CBOOL);
__EXTERN void RMSBestFitEval    (__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CFLOAT *, __CBOOL);
__EXTERN void RMSBestFitGrad    (__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CFLOAT *, __CINT, __CINT, __CBOOL);
__EXTERN void RMSBestFitGradEval(__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CFLOAT *, __CINT, __CINT, __CFLOAT *,__CBOOL);
//compute center of mass
__EXTERN void* com(__CFLOAT*, __CFLOAT*,__CINT,__CBOOL);
//multiply two matrices
__EXTERN void* matmul(__CFLOAT*, __CFLOAT*, __CINT, __CINT, __CINT);
//compute rmsd between two vectors
__EXTERN __CFLOAT rmsd(__CFLOAT*, __CFLOAT*, __CFLOAT*, __CINT, __CBOOL);
