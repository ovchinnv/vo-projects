#!/bin/csh -f

set this_prog = $0

if ($#argv < 1) then
  echo "Usage : $this_prog extracts the memory usage of a cuda program."
  echo "Syntax: $this_prog <.cu file>"
  echo ""
  echo "Example: $this_prog matrixmul.cu"
  echo ""
  exit 1
endif

set cu_file = $1
#set cuda_sdk_inc = cuda_sdk/inc
set cuda_sdk_inc = ../../common/inc

echo "nvcc --ptxas-options -v -I. -I$cuda_sdk_inc -I/usr/local/cuda/include -DUNIX  -o $cu_file.cubin -cubin $cu_file"
nvcc --ptxas-options -v -I. -I$cuda_sdk_inc -I/usr/local/cuda/include -DUNIX  -o $cu_file.cubin -cubin $cu_file

