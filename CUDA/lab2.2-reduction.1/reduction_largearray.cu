/*
 * Copyright 1993-2006 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.   
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.  This source code is a "commercial item" as 
 * that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer software" and "commercial computer software 
 * documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein.
 */

#ifdef _WIN32
#  define NOMINMAX 
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil.h>

// includes, kernels
#include <reduction_largearray_kernel.cu>  

#define DEFAULT_NUM_ELEMENTS 512 
#define MAX_RAND 3


////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest( int argc, char** argv);

void WriteFile(float*, char* file_name, int size);

extern "C" 
unsigned int compare( const float* reference, const float* data, 
                     const unsigned int len);
extern "C" 
void computeGold( float* reference, float* idata, const unsigned int len);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int 
main( int argc, char** argv) 
{
    runTest( argc, argv);
    return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
//! Run a scan test for CUDA
////////////////////////////////////////////////////////////////////////////////
void
runTest( int argc, char** argv) 
{
    CUT_DEVICE_INIT(argc, argv);

    float device_time;
    float host_time;
    int num_elements = 0; // Must support large, non-power-of-2 arrays

    // allocate host memory to store the input data
    unsigned int mem_size = sizeof( float) * num_elements;
    float* h_data = (float*) malloc( mem_size);

    // initialize the input data on the host to be integer values
    // between 0 and 1000
    // Use DEFAULT_NUM_ELEMENTS num_elements
    num_elements = DEFAULT_NUM_ELEMENTS;

    // allocate host memory to store the input data
    mem_size = sizeof( float) * num_elements;
    h_data = (float*) malloc( mem_size);

    // initialize the input data on the host
    for( unsigned int i = 0; i < num_elements; ++i) 
    {
	h_data[i] = (int)(rand() % MAX_RAND);
    }
   
    
    unsigned int timer;
    CUT_SAFE_CALL(cutCreateTimer(&timer));

      
    // compute reference solution
    float* reference = (float*) malloc( mem_size);  
	cutStartTimer(timer);
    computeGold( reference, h_data, num_elements);
	cutStopTimer(timer);
    printf("\n\n**===-------------------------------------------------===**\n");
    printf("Processing %d elements...\n", num_elements);
    printf("Host CPU Processing time: %f (ms)\n", cutGetTimerValue(timer));
    host_time = cutGetTimerValue(timer);
    CUT_SAFE_CALL(cutDeleteTimer(timer));


    // allocate device memory input and output arrays
    float* d_idata = NULL;
    float* d_odata = NULL;

    CUDA_SAFE_CALL( cudaMalloc( (void**) &d_idata, mem_size));
    CUDA_SAFE_CALL( cudaMalloc( (void**) &d_odata, mem_size));
    
    // copy host memory to device input array
    CUDA_SAFE_CALL( cudaMemcpy( d_idata, h_data, mem_size, cudaMemcpyHostToDevice) );
    // initialize all the other device arrays to be safe
    CUDA_SAFE_CALL( cudaMemcpy( d_odata, h_data, mem_size, cudaMemcpyHostToDevice) );

    // **===-----------------------------------------------------------===**

    // Run just once to remove startup overhead for more accurate performance 
    // measurement
    reductionArray(d_odata, d_idata, 16);

    // Run the reduction 
    // **===-------- Lab 2.2 - Code segment 1.1:
    //       Modify here to allocate timer and start timer
    //       for the kernel this function -----------===**
    
    // **===-------- Lab 2.2 - Modify the body of this function -----------===**
    reductionArray(d_odata, d_idata, num_elements);
    // **===---------Lab 2.2 - Code segment 1.2: -----===**
    //               call CUDA thread synchronization function
    //               to wait until kernel finishes

    // **===-------- Lab 2.2 - Code segment 1.2: -----===**
    //		     stop timer and delete timer here

    // There's no need to change anything below
    printf("CUDA Processing time: %f (ms)\n", cutGetTimerValue(timer));
    device_time = cutGetTimerValue(timer);
    printf("Speedup: %fX\n", host_time/device_time);

    // **===-----------------------------------------------------------===**


    // copy result from device to host
    CUDA_SAFE_CALL(cudaMemcpy( h_data, d_odata, sizeof(float) * num_elements, 
                               cudaMemcpyDeviceToHost));

    // We can use an epsilon of 0 since values are integral and in a range 
    // that can be exactly represented
    float epsilon = 0.0f, result = fabs(*h_data - *reference); 
    bool result_regtest = result <= epsilon ;
    printf( "Test %s\n", result_regtest ? "PASSED" : "FAILED");
    printf( "device: %f  host: %f\n", h_data[0], *reference);

    // cleanup memory
    cutDeleteTimer(timer);
    free( h_data);
    free( reference);
    cudaFree( d_odata);
    cudaFree( d_idata);
}

