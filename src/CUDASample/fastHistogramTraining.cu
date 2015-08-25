#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <fstream>

// CUDA helper functions
#include <helper_cuda.h>         // helper functions for CUDA error check

// Normal distribution generators
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>

#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/device_ptr.h>


void initializeValues(unsigned int *data,
        const unsigned int numElems,
        const unsigned int numBins)
{
    srand(time(NULL));
    unsigned int mean = rand() % 100 + 462;
    const float stddev = 100.f;
    std::cout << "Mean value: " << mean << std::endl;
    thrust::minstd_rand rng;
    thrust::random::normal_distribution<float> normalDist((float)mean, stddev);
    for (size_t i = 0; i < numElems; ++i){
        data[i] = min(max((int)normalDist(rng), 0), numBins - 1);
    }
}

void storeResults(unsigned int* histo,
        const unsigned int numBins)
{
    std::ofstream outFile;
    outFile.open("out.hist", std::ios::trunc);
    for(int i=0; i<numBins; i++)
    {
        outFile << histo[i] << std::endl;
    }
}



__global__
void divideValuesByKeys(unsigned long* const d_vals,
        unsigned long* const d_keys,
        const unsigned int numBins)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= numBins) return;
    if (d_keys[tId] != 0) d_vals[tId] = d_vals[tId]/d_keys[tId];
}

__global__
void convertUI2UL(const unsigned int* const in,
        unsigned long* const out,
        const unsigned int numElems)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= numElems) return;
    out[tId] = static_cast<unsigned long>(in[tId]);
}

__global__
void convertUL2UI(unsigned long* const in,
        unsigned int* const out,
        const unsigned int numElems)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= numElems) return;
    out[tId] = static_cast<unsigned int>(in[tId]);
}

__global__
void saveToHistogram(unsigned long* const vals,
        unsigned long* const keys,
        unsigned int* const hist,
        const unsigned int numBins)
{
    int j=0;
    for (int i=0; i<numBins; i++)
    {
        if( i==keys[j] )
        {
            hist[i] = static_cast<unsigned int>(vals[j]);
            j++;
        }
        else
        {
            hist[i] = 0;
        }
    }
}

void computeHistogram(const unsigned int* const d_vals,
        unsigned int* const d_histo,
        const unsigned int numBins,
        const unsigned int numElems)
{
    unsigned int blockSize;
    unsigned int gridSize;
    unsigned long *d_values;
    unsigned long *d_keys;
    unsigned long *d_outValues;
    unsigned long *d_outKeys;
    
    blockSize = 1024;
    gridSize = (numElems+blockSize-1) / blockSize;
    checkCudaErrors(cudaMalloc(&d_values, sizeof(unsigned long)*numElems));
    checkCudaErrors(cudaMalloc(&d_keys, sizeof(unsigned long)*numElems));
    checkCudaErrors(cudaMalloc(&d_outValues, sizeof(unsigned long)*numElems));
    checkCudaErrors(cudaMalloc(&d_outKeys, sizeof(unsigned long)*numElems));
    convertUI2UL<<<gridSize, blockSize>>>(d_vals, d_values, numElems);
    convertUI2UL<<<gridSize, blockSize>>>(d_vals, d_keys, numElems);
    cudaDeviceSynchronize();
    thrust::device_ptr<unsigned long> devptr_keys(d_keys);
    thrust::device_ptr<unsigned long> devptr_values(d_values);
    thrust::device_ptr<unsigned long> devptr_outKeys(d_outKeys);
    thrust::device_ptr<unsigned long> devptr_outValues(d_outValues);
    thrust::sort_by_key(devptr_keys, devptr_keys+numElems, devptr_values);
    cudaDeviceSynchronize();
    thrust::reduce_by_key(devptr_keys, devptr_keys+numElems, devptr_values, devptr_outKeys, devptr_outValues);
    cudaDeviceSynchronize();
    divideValuesByKeys<<<gridSize, blockSize>>>(d_outValues, d_outKeys, numBins);
    cudaDeviceSynchronize();
    saveToHistogram<<<1, 1>>>(d_outValues, d_outKeys, d_histo, numBins);
    cudaDeviceSynchronize();
    cudaFree(d_values);
    cudaFree(d_keys);
    cudaFree(d_outValues);
    cudaFree(d_outKeys);
}

int AddFunc(int a, int b)
{
    a=1;
    b=1;
    thrust::plus<int> oper;
    return oper(a,b);
}

void fastHistogramSample()
{

    printf("FastHistogram sample starts ...\n");
    dim3 blockSize(1,1,1);
    dim3 gridSize(1,1,1);
    const unsigned int numBins = 1024;
    const unsigned int numElems = 10000 * numBins;
    unsigned int *h_vals = new unsigned int[numElems];
    unsigned int *h_histo = new unsigned int[numBins];
    unsigned int *d_vals;
    unsigned int *d_histo;

    initializeValues(h_vals, numElems, numBins);
    checkCudaErrors(cudaMalloc(&d_vals, sizeof(unsigned int)*numElems));
    checkCudaErrors(cudaMalloc(&d_histo, sizeof(unsigned int)*numBins));
    checkCudaErrors(cudaMemset(d_histo, 0, sizeof(unsigned int)*numBins));
    checkCudaErrors(cudaMemcpy(d_vals, h_vals, sizeof(unsigned int)*numElems, cudaMemcpyHostToDevice));
    computeHistogram(d_vals, d_histo, numBins, numElems);
    checkCudaErrors(cudaMemcpy(h_histo, d_histo, sizeof(unsigned int)*numBins, cudaMemcpyDeviceToHost));
    storeResults(h_histo, numBins);

    delete h_vals;
    delete h_histo;
    cudaFree(d_vals);
    cudaFree(d_histo);
    printf("DONE!\n");

}
