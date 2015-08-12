#ifndef BLELLOCH_TEMPLATES
#define BLELLOCH_TEMPLATES

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Includes CUDA
#include <cuda_runtime.h>

// Utilities and timing functions
#include <helper_functions.h>    // includes cuda.h and cuda_runtime_api.h

// CUDA helper functions
#include <helper_cuda.h>         // helper functions for CUDA error check

template <class T>
__global__ void blellochReduceKernel(T *input, int n, int size)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    const int dId = size - (tId*2*n + 1);
    if (dId < 0) return;
    if (dId - n < 0)
    {
        return;
    }
    input[dId] += input[dId-n];
}

template <class T>
__global__ void blellochDownsweepKernel(T *input, int n, int size)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    const int dId = size - (tId*2*n + 1);
    T left;
    T right;
    if (dId < 0) return;
    if (dId - n < 0)
    {
        return;
    }
    right = input[dId];
    left = input[dId-n];
    input[dId-n] = right;
    input[dId] = left+right;
}

template <class T>
__global__ void setLastValueToZero(T *input, int size)
{
    input[size-1] = 0;
}

template <class T>
__global__ void printTable(T *input, int size)
{
    printf("[%d",(int)input[0]);
    for (int i=1; i<size; i++)
    {
        printf(",%d",(int)input[i]);
    }
    printf("]\n");
}

template <class T>
void invokeBlellochKernel(T *input, int size)
{
    dim3 blockSize(512,1,1);
    dim3 gridSize(1,1,1);
    int n;

    for (n=1; n < size; n = n*2)
    {
        gridSize.x = ((size+2*n-1)/(2*n) + blockSize.x -1) / blockSize.x;
        //printf("Will run %d blocks\n",gridSize.x);
        blellochReduceKernel<T><<<gridSize, blockSize>>>(input, n, size);
        cudaDeviceSynchronize();
    }

    setLastValueToZero<T><<<1,1>>>(input, size);
    cudaDeviceSynchronize();

    for (n = n/2; n >= 1; n = n/2)
    {
        gridSize.x = ((size+2*n-1)/(2*n) + blockSize.x -1) / blockSize.x;
        //printf("Will run %d blocks\n",gridSize.x);
        blellochDownsweepKernel<T><<<gridSize, blockSize>>>(input, n, size);
        cudaDeviceSynchronize();
    }
}


#endif // BLELLOCH_TEMPLATES