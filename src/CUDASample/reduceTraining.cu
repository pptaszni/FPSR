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

typedef float(*reduceFPtr)(float,float);

__device__ float d_Max(float a, float b)
{
    return max(a, b);
}

__device__ float d_Add(float a, float b)
{
    return a+b;
}

__global__ void reduceFun(float *input, unsigned int s, float *output, reduceFPtr func, bool last)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= s) return;
    if (!last && tId < s)
    {
        input[tId] = (*func)(input[tId], input[tId+s]);
    }
    if (last) input[s] = input[s+s];
    if (tId == 0) *output = input[0];
}

void invokeReduceFun(float *input, int inputSize, float *output, reduceFPtr func)
{
    dim3 blockSize(1,1,1);
    dim3 gridSize(1,1,1);
    bool odd = false;
    unsigned int s = inputSize;
    while (s > 1)
    {
        odd = (s%2);
        s >>= 1;
        blockSize.x = 512;
        gridSize.x = (s + blockSize.x - 1) / blockSize.x;
        reduceFun<<<gridSize, blockSize>>>(input, s, output, func, false);
        cudaDeviceSynchronize();
        if (odd)
        {
            reduceFun<<<1, 1>>>(input, s, output, func, true);
            cudaDeviceSynchronize();
            s++;
        }
    }

}

__device__ reduceFPtr d_p1 = d_Add;

void reduceSample()
{
    printf("Starting reduce sample ...\n");
    dim3 blockSize(1,1,1);
    dim3 gridSize(1,1,1);
    unsigned int tableSize = 1000;
    float *d_tableToReduce;
    float *h_tableToReduce;
    float *d_result;
    float *h_result;
    reduceFPtr h_p1;
    checkCudaErrors(cudaMemcpyFromSymbol(&h_p1, d_p1, sizeof(reduceFPtr)));

    h_tableToReduce = (float*) malloc(sizeof(float)*tableSize);
    h_result = (float*) malloc(sizeof(float));
    checkCudaErrors(cudaMalloc(&d_tableToReduce, sizeof(float)*tableSize));
    checkCudaErrors(cudaMalloc(&d_result, sizeof(float)));

    for (unsigned int i=0; i<tableSize; i++) h_tableToReduce[i]=(float)i;

    checkCudaErrors(cudaMemcpy(d_tableToReduce, h_tableToReduce, sizeof(float)*tableSize, cudaMemcpyHostToDevice));
    cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
    invokeReduceFun(d_tableToReduce, tableSize, d_result, h_p1);
    cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());

    checkCudaErrors(cudaMemcpy(h_result, d_result, sizeof(float), cudaMemcpyDeviceToHost));

    printf("Obtained result: %f\n",*h_result);
    checkCudaErrors(cudaFree(d_tableToReduce));
    checkCudaErrors(cudaFree(d_result));
    free(h_tableToReduce);
    free(h_result);
}


