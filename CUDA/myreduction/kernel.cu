#ifdef __TEX
#define in( _I ) tex1Dfetch(tex_reduce, _I )
#else
#define in( _I ) in[ _I ]
#endif

// this is the best kernel for the chromnebook
 __global__ void reduction_cuda(__CTYPE *in, __CTYPE *out, int n) {

  unsigned int tx = threadIdx.x;
  unsigned int ix = blockIdx.x * (2*_BDIM) + tx;
  unsigned int stride = 2*_BDIM*gridDim.x;
//
  __shared__ volatile float local[_BDIM]; // necessary to avoid shmem=>register optimization, to avoid syncthreads within warp
  local[tx]=(__CTYPE)_REDINI; // initializer
//
// serial part
  while (ix<n-_BDIM){
   _REDOP(local[tx],local[tx],in(ix));
   _REDOP(local[tx],local[tx],in(ix+_BDIM));
   ix+=stride;
  }
  while (ix<n){
   _REDOP(local[tx],local[tx],in(ix));
   ix+=stride;
  }
//  local[tx]=v;
  __syncthreads() ;
//parallel part:

// loop only if the block size is beyond a certain size
#define _REDUCE(_STEP) if (tx<_STEP) _REDOP(local[tx],local[tx],local[tx + _STEP]) ; __syncthreads();
#if _BDIM>=1024
  for (unsigned int step=_BDIM>>1 ; step > 256 ; step>>=1) {
   _REDUCE(step)
  }
#endif
#if _BDIM >= 512
  _REDUCE(256)
#endif
#if _BDIM >= 256
  _REDUCE(128)
#endif
#if _BDIM >= 128
  _REDUCE(64)
#endif
#undef _REDUCE
#define _REDUCE(_STEP) _REDOP(local[tx],local[tx],local[tx + _STEP]) ;
  if (tx<32) {
#if _BDIM >= 64
   _REDUCE(32)
#endif
#if _BDIM >= 32
   _REDUCE(16)
#endif
#if _BDIM >= 16
   _REDUCE(8)
#endif
#if _BDIM >= 8
   _REDUCE(4)
#endif
#if _BDIM >= 4
   _REDUCE(2)
#endif
#if _BDIM >= 2
   _REDUCE(1)
#endif
  }
  if(tx==0) out[blockIdx.x]=local[0];
 }

/****************************************************************************/
 __global__ void reduction_cuda_simple(__CTYPE *in, __CTYPE *out, int n) {

  unsigned int tx=threadIdx.x;
  unsigned int ix=2*_BDIM*blockIdx.x + tx ;
  __shared__ __CTYPE local[2*_BDIM];

//note that I could do with less storage space here
  if (ix+_BDIM < n) local[tx+_BDIM]=in[ix+_BDIM];
  if (ix<n)         local[tx]=in[ix];
//  local[tx+_BDIM]=in[ix+_BDIM];
//  local[tx]=in[ix];

// NOTE : block size must be a power of two
// otherwise the loop below will not work
  for (unsigned int step=_BDIM ; step > 0 ; step/=2) {
   if (tx < step) local[tx]+=local[tx+step];
   __syncthreads();
  }
  if (tx==0) out[blockIdx.x]=local[0];
 }

/****************************************************************************/
// slighy faster
 __global__ void reduction_cuda_simple2(__CTYPE *in, __CTYPE *out, int n) {

  unsigned int tx=threadIdx.x;
  unsigned int ix=2*_BDIM*blockIdx.x + tx ;
  __shared__ __CTYPE local[_BDIM];

//note that I could do with less storage space here
  local[tx]=(__CTYPE)0;
  if (ix+_BDIM < n) local[tx]+=in[ix+_BDIM];
  if (ix<n)         local[tx]+=in[ix];
  __syncthreads();
//  local[tx+_BDIM]=in[ix+_BDIM];
//  local[tx]=in[ix];

// NOTE : block size must be a power of two
// otherwise the loop below will not work
  for (unsigned int step=_BDIM>>1 ; step > 0 ; step>>=1) {
   if (tx < step) local[tx]+=local[tx+step];
   __syncthreads();
  }
  if (tx==0) out[blockIdx.x]=local[0];
 }

/****************************************************************************/
// significantly faster
// further optimization would require doing more work per thread, as in the first kernel
 __global__ void reduction_cuda_simple3(__CTYPE *in, __CTYPE *out, int n) {

  unsigned int tx=threadIdx.x;
  unsigned int ix=2*_BDIM*blockIdx.x + tx ;
  __shared__  volatile __CTYPE local[_BDIM]; // volatile needed to avoid syncthreads

//note that I could do with less storage space here
  local[tx]=(__CTYPE)0;
  if (ix+_BDIM < n) local[tx]+=in[ix+_BDIM];
  if (ix<n)         local[tx]+=in[ix];
  __syncthreads();
//  local[tx+_BDIM]=in[ix+_BDIM];
//  local[tx]=in[ix];

// NOTE : block size must be a power of two
// otherwise the loop below will not work
// unroll
//
#ifdef _REDUCE
#undef _REDUCE
#endif

#define _REDUCE(_STEP) if (tx<_STEP) local[tx] += local[tx + _STEP] ; __syncthreads();

// loop only if the block size is beyond a certain size
#if _BDIM>=1024
  for (unsigned int step=_BDIM>>1 ; step > 256 ; step>>=1) {
   __REDUCE(step);
  }
#endif

#if _BDIM >= 512
  _REDUCE(256)

#endif

#if _BDIM >= 256
  _REDUCE(128)

#endif

#if _BDIM >= 128
  _REDUCE(64)

#endif

#undef _REDUCE
#define _REDUCE(_STEP) local[tx] += local[tx + _STEP] ; // __syncthreads();

  if (tx<32) {
#if _BDIM >= 64
    _REDUCE(32)
#endif
#if _BDIM >= 32
    _REDUCE(16)
#endif
#if _BDIM >= 16
    _REDUCE(8)
#endif
#if _BDIM >= 8
    _REDUCE(4)
#endif
#if _BDIM >= 4
    _REDUCE(2)
#endif
#if _BDIM >= 2
    _REDUCE(1)
#endif
  }
  if (tx==0) out[blockIdx.x]=local[0];
 }
