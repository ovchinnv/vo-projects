/* C prototypes */
#ifdef __CPLUSPLUS
#define __EXTERNC "extern C"
#else 
#define __EXTERNC
#endif
//best-fit
__EXTERNC void RMSBestFit        (__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CBOOL);
__EXTERNC void RMSBestFitEval    (__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CFLOAT *, __CBOOL);
__EXTERNC void RMSBestFitGrad    (__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CFLOAT *, __CINT, __CINT, __CBOOL);
__EXTERNC void RMSBestFitGradEval(__CFLOAT *,__CFLOAT *, __CFLOAT *, __CINT, __CFLOAT *, __CFLOAT *, __CINT, __CINT, __CFLOAT *,__CBOOL);
//com
__EXTERNC void* com(__CFLOAT*, __CFLOAT*,__CINT,__CBOOL);
//matmul
__EXTERNC void* matmul(__CFLOAT*, __CFLOAT*, __CINT, __CINT, __CINT);
//rmsd
__EXTERNC __CFLOAT rmsd(__CFLOAT*, __CFLOAT*, __CFLOAT*, __CINT, __CBOOL);
