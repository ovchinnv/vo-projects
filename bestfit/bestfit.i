#include "source.defs"
#undef __CPLUSPLUS
#define __SWIG
// prototype
%module bestfit
%{
#include "bestfit.h"
%}
// parse header file to generate wrappers
%include "bestfit.h"
