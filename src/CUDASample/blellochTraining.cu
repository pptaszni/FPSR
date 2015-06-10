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


__global__ void blellochReduceKernel(float *input, int n, int size)
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

__global__ void blellochDownsweepKernel(float *input, int n, int size)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    const int dId = size - (tId*2*n + 1);
    float left;
    float right;
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

__global__ void setLastValueToZero(float *input, int size)
{
    input[size-1] = 0;
}

__global__ void printTable(float *input, int size)
{
    printf("[%d",input[0]);
    for (int i=1; i<size; i++)
    {
        printf(",%d",(int)input[i]);
    }
    printf("]\n");
}

void invokeBlellochKernel(float *input, int size)
{
    dim3 blockSize(512,1,1);
    dim3 gridSize(1,1,1);
    int n;

    for (n=1; n < size; n = n*2)
    {
        gridSize.x = ((size+2*n-1)/(2*n) + blockSize.x -1) / blockSize.x;
        printf("Will run %d blocks\n",gridSize.x);
        blellochReduceKernel<<<gridSize, blockSize>>>(input, n, size);
        cudaDeviceSynchronize();
        printTable<<<1,1>>>(input,size);
        cudaDeviceSynchronize();
    }

    setLastValueToZero<<<1,1>>>(input, size);
    cudaDeviceSynchronize();
    printTable<<<1,1>>>(input,size);
    cudaDeviceSynchronize();

    for (n = n/2; n >= 1; n = n/2)
    {
        gridSize.x = ((size+2*n-1)/(2*n) + blockSize.x -1) / blockSize.x;
        printf("Will run %d blocks\n",gridSize.x);
        blellochDownsweepKernel<<<gridSize, blockSize>>>(input, n, size);
        cudaDeviceSynchronize();
        printTable<<<1,1>>>(input,size);
        cudaDeviceSynchronize();
    }
}



void blellochSample()
{

    printf("Blelloch sample starts ...\n");
    dim3 blockSize(1,1,1);
    dim3 gridSize(1,1,1);
    unsigned int tableSize = 50;
    float *d_tableToScan;
    float *h_tableToScan;


    h_tableToScan = (float*) malloc(sizeof(float)*tableSize);
    checkCudaErrors(cudaMalloc(&d_tableToScan, sizeof(float)*tableSize));
    for (int i=0;i<tableSize;i++) h_tableToScan[i]=i;
    checkCudaErrors(cudaMemcpy(d_tableToScan, h_tableToScan, sizeof(float)*tableSize, cudaMemcpyHostToDevice));

    invokeBlellochKernel(d_tableToScan, tableSize);


}
