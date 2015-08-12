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

#include "blellochTemplates.hpp"



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
    for (int i=0;i<tableSize;i++)
    {
        h_tableToScan[i]=i;
    }
    h_tableToScan[0] = 3;
    h_tableToScan[1] = 7;
    checkCudaErrors(cudaMemcpy(d_tableToScan, h_tableToScan,
        sizeof(float)*tableSize, cudaMemcpyHostToDevice));

    invokeBlellochKernel<float>(d_tableToScan, tableSize);

    checkCudaErrors(cudaMemcpy(h_tableToScan, d_tableToScan,
        sizeof(float)*tableSize, cudaMemcpyDeviceToHost));

    printf("[%d",(int)h_tableToScan[0]);
    for (int i=1; i<tableSize; i++)
    {
        printf(",%d",(int)h_tableToScan[i]);
    }
    printf("]\n");

    free(h_tableToScan);
    cudaFree(d_tableToScan);


}
