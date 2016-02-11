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
__EXTERN __CINT smcv_dyna_from_namd ( const long int, const __CFLOAT *, __CFLOAT *, __CFLOAT * );

#ifndef GLOBALMASTERSMCV_H
#define GLOBALMASTERSMCV_H
#endif

class Molecule;
class SubmitReduction;

class GlobalMasterSMCV : public GlobalMaster {
 public:
  GlobalMasterSMCV();
  ~GlobalMasterSMCV();
 private:
  Bool initialized;
  void calculate();
  void destroy();
  Molecule *molecule;
  SimParameters * params;
//  SubmitReduction *reduction;
  long int iteration;
  __CFLOAT smcv_energy;
  __CINT *atomlist=NULL; // to maintain a list of atoms that are needed by plugin; memory managed by FORTRAN side
  __CFLOAT *r=NULL, *fr=NULL ; //positions and forces
};
