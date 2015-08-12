#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>

// Includes CUDA
#include <cuda_runtime.h>

// Utilities and timing functions
#include <helper_functions.h>    // includes cuda.h and cuda_runtime_api.h

// CUDA helper functions
#include <helper_cuda.h>         // helper functions for CUDA error check

#include "blellochTemplates.hpp"


struct keyValuePair
{
    unsigned int key;
    unsigned int value;
};

template <>
__global__ void printTable<keyValuePair>(keyValuePair *input, int size)
{
    printf("[%x",input[0].key);
    for (int i=1; i<size; i++)
    {
        printf(",%x",(int)input[i].key);
    }
    printf("]\n");
}

__global__
void computeHistogramAndRelativeOffsets(keyValuePair* input, unsigned int* output,
    unsigned int* offsets, int inputSize, unsigned int mask, unsigned int iter)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= inputSize) return;
    unsigned int binNum = (input[tId].key & mask) >> iter;
    offsets[binNum*inputSize + tId] = 1;
    atomicAdd(&(output[binNum]),1);
}

__global__
void scatterToAddresses(keyValuePair* input, keyValuePair* output, unsigned int* globAddr,
    unsigned int* locAddr, int inputSize, int numBins, unsigned int mask, unsigned int iter)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= inputSize) return;
    unsigned int binNum = (input[tId].key & mask) >> iter;
    unsigned int relOff = locAddr[binNum*inputSize+tId];
    output[globAddr[binNum]+relOff] = input[tId];
}

void runRadixSort(keyValuePair* h_input, keyValuePair* h_output, int size)
{
    int numBits = 8;
    int numBins = 1 << numBits;
    int blockSize = 1;
    int gridSize = 1;
    keyValuePair* d_input;
    keyValuePair* d_output;
    unsigned int* d_absAddr;
    unsigned int* d_offsets;
    unsigned int mask;

    checkCudaErrors(cudaMalloc(&d_input, sizeof(keyValuePair)*size));
    checkCudaErrors(cudaMalloc(&d_output, sizeof(keyValuePair)*size));
    checkCudaErrors(cudaMalloc(&d_absAddr, sizeof(unsigned int)*numBins));
    checkCudaErrors(cudaMalloc(&d_offsets, sizeof(unsigned int)*numBins*size));
    checkCudaErrors(cudaMemcpy(d_input, h_input,
        sizeof(keyValuePair)*size, cudaMemcpyHostToDevice));

    for (unsigned int i=0; i < 8*sizeof(unsigned int); i+=numBits)
    {
        mask = (numBins - 1) << i;
        //printf("maska: %x, iter: %d\n", mask, i);
        blockSize = (size < 512) ? size : 512;
        gridSize = (size+blockSize-1) / blockSize;
        checkCudaErrors(cudaMemset(d_absAddr, 0, sizeof(unsigned int)*numBins));
        checkCudaErrors(cudaMemset(d_offsets, 0, sizeof(unsigned int)*numBins*size));
        cudaDeviceSynchronize();
        computeHistogramAndRelativeOffsets<<<gridSize, blockSize>>>(
            d_input, d_absAddr, d_offsets, size, mask, i);
        cudaDeviceSynchronize();
        invokeBlellochKernel<unsigned int>(d_absAddr, numBins);
        cudaDeviceSynchronize();
        //printTable<keyValuePair><<<1,1>>>(d_input, size);
        //cudaDeviceSynchronize();
        for (int j=0; j<numBins; j++)
        {
            invokeBlellochKernel<unsigned int>(d_offsets+size*j, size);
            cudaDeviceSynchronize();
        }
        scatterToAddresses<<<gridSize, blockSize>>>(d_input, d_output, d_absAddr,
            d_offsets, size, numBins, mask, i);
        cudaDeviceSynchronize();
        //printTable<keyValuePair><<<1,1>>>(d_output, size);
        //cudaDeviceSynchronize();
        std::swap(d_input,d_output);
    }

    checkCudaErrors(cudaMemcpy(h_output, d_input,
        sizeof(keyValuePair)*size, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(d_input));
    checkCudaErrors(cudaFree(d_output));
    checkCudaErrors(cudaFree(d_absAddr));
    checkCudaErrors(cudaFree(d_offsets));
}

void initializeTable(keyValuePair* tab, int size)
{
    int seed = time(NULL);
    srand(seed);
    for (int i=0; i<size; i++)
    {
        tab[i].key = rand()%1000;
        tab[i].value = i;
    }
}

void radixSortSample()
{

    printf("RadixSort sample starts ...\n");
    dim3 blockSize(1,1,1);
    dim3 gridSize(1,1,1);
    unsigned int tableSize = 150;
    keyValuePair *h_tableToSort;
    keyValuePair *h_result;

    h_tableToSort = (keyValuePair*) malloc(sizeof(keyValuePair)*tableSize);
    h_result = (keyValuePair*) malloc(sizeof(keyValuePair)*tableSize);
    initializeTable(h_tableToSort, tableSize);
    runRadixSort(h_tableToSort, h_result, tableSize);

    printf("Results:\n");
    printf("[%d",h_result[0].key);
    for (int i=1; i<tableSize; i++)
    {
        printf(",%d",h_result[i].key);
    }
    printf("]\n");

    free(h_tableToSort);
    free(h_result);
}
