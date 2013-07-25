/* C prototypes */
#ifdef __CPLUSPLUS
#define __EXTERN "extern C"
//#elif defined(__SWIG)
//#define __EXTERN extern
#else
#define __EXTERN extern
#endif
//best-fit routines
__EXTERN void RMSBestFit        (const __CFLOAT *,const __CFLOAT *,const  __CFLOAT *,const __CINT, __CFLOAT *,const __CBOOL);
__EXTERN void RMSBestFitEval    (const __CFLOAT *,const __CFLOAT *,const  __CFLOAT *,const __CINT, __CFLOAT *, __CFLOAT *,const __CBOOL);
__EXTERN void RMSBestFitGrad    (const __CFLOAT *,const __CFLOAT *,const  __CFLOAT *,const __CINT, __CFLOAT *, __CFLOAT *,\
                                 const  __CINT, const __CINT, const __CBOOL);
__EXTERN void RMSBestFitGradEval(const __CFLOAT *,const __CFLOAT *,const  __CFLOAT *,const __CINT, __CFLOAT *, __CFLOAT *,\
                                 const __CINT,const __CINT, __CFLOAT *,const __CBOOL);
//compute center of mass
__EXTERN void* com(const __CFLOAT*,const __CFLOAT*,const __CINT,const __CBOOL);
//multiply two matrices
__EXTERN void* matmul(const __CFLOAT*,const __CFLOAT*,const __CINT,const __CINT,const __CINT);
//compute rmsd between two vectors
__EXTERN __CFLOAT rmsd(const __CFLOAT*,const __CFLOAT*,const __CFLOAT*,const __CINT,const __CBOOL);
