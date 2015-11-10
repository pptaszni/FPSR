#include <stdlib.h>
#include <iostream>
#include <string>
#include <opencv2/opencv.hpp>

// CUDA helper functions
#include <helper_cuda.h>         // helper functions for CUDA error check

// remember to delete the allocated imagePtr in the end
void loadImageRGBA(const std::string &filename,
        uchar4 **imagePtr,
        size_t *numRows, size_t *numCols)
{
    cv::Mat image = cv::imread(filename.c_str(), CV_LOAD_IMAGE_COLOR);
    if (image.empty()) {
        std::cerr << "Couldn't open file: " << filename << std::endl;
        exit(1);
    }

    if (image.channels() != 3) {
        std::cerr << "Image must be color!" << std::endl;
        exit(1);
    }

    if (!image.isContinuous()) {
        std::cerr << "Image isn't continuous!" << std::endl;
        exit(1);
    }
    
    cv::Mat imageRGBA;
    cv::cvtColor(image, imageRGBA, CV_BGR2RGBA);
    
    *imagePtr = new uchar4[image.rows * image.cols];
    
    unsigned char *cvPtr = imageRGBA.ptr<unsigned char>(0);
    for (size_t i = 0; i < image.rows * image.cols; ++i) {
        (*imagePtr)[i].x = cvPtr[4 * i + 0];
        (*imagePtr)[i].y = cvPtr[4 * i + 1];
        (*imagePtr)[i].z = cvPtr[4 * i + 2];
        (*imagePtr)[i].w = cvPtr[4 * i + 3];
    }
    
    *numRows = image.rows;
    *numCols = image.cols;
}

void saveImageRGBA(const uchar4* const image,
        const size_t numRows, const size_t numCols,
        const std::string &output_file)
{
    int sizes[2];
    sizes[0] = numRows;
    sizes[1] = numCols;
    cv::Mat imageRGBA(2, sizes, CV_8UC4, (void *)image);
    cv::Mat imageOutputBGR;
    cv::cvtColor(imageRGBA, imageOutputBGR, CV_RGBA2BGR);
    //output the image
    cv::imwrite(output_file.c_str(), imageOutputBGR);
}

__global__
void computeMask(const uchar4* const sourceImg,
        bool* const mask,
        const size_t size)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= size) return;
    unsigned int pixSum = sourceImg[tId].x + sourceImg[tId].y + sourceImg[tId].z;
    if (pixSum < 3*255)
    {
        mask[tId] = true;
    }
    else
    {
        mask[tId] = false;
    }
}

__global__
void computeBorderPixAndInteriorPix(const bool* const mask,
        bool* const interiorPixels,
        bool* const borderPixels,
        const size_t numRows,
        const size_t numCols)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= numRows*numCols) return;
    if (!mask[tId])
    {
        interiorPixels[tId] = false;
        borderPixels[tId] = false;
        return;
    }
    if (mask[tId - numCols] && mask[tId + numCols]
            && mask[tId-1] && mask[tId+1])
    {
        interiorPixels[tId] = true;
        borderPixels[tId] = false;
    }
    else
    {
        interiorPixels[tId] = false;
        borderPixels[tId] = true;
    }
}

__global__
void computeG(const uchar4* const sourceImg,
        const bool* const interiorPixels,
        float3* const g,
        const size_t numRows,
        const size_t numCols)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= numRows*numCols) return;
    if (!interiorPixels[tId]) return;
    float sumX = 4.f*sourceImg[tId].x;
    float sumY = 4.f*sourceImg[tId].y;
    float sumZ = 4.f*sourceImg[tId].z;
    sumX -= (float)sourceImg[tId-numCols].x + (float)sourceImg[tId+numCols].x
        + (float)sourceImg[tId-1].x + (float)sourceImg[tId+1].x;
    sumY -= (float)sourceImg[tId-numCols].y + (float)sourceImg[tId+numCols].y
        + (float)sourceImg[tId-1].y + (float)sourceImg[tId+1].y;
    sumZ -= (float)sourceImg[tId-numCols].z + (float)sourceImg[tId+numCols].z
        + (float)sourceImg[tId-1].z + (float)sourceImg[tId+1].z;
    g[tId].x = sumX;
    g[tId].y = sumY;
    g[tId].z = sumZ;
}

__global__
void copySourceImgToBlendedVals(const uchar4* const sourceImg,
        float3* const blendedVals1,
        float3* const blendedVals2,
        const size_t size)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= size) return;
    blendedVals1[tId].x = (float)sourceImg[tId].x;
    blendedVals1[tId].y = (float)sourceImg[tId].y;
    blendedVals1[tId].z = (float)sourceImg[tId].z;
    blendedVals2[tId].x = (float)sourceImg[tId].x;
    blendedVals2[tId].y = (float)sourceImg[tId].y;
    blendedVals2[tId].z = (float)sourceImg[tId].z;
}

__global__
void computeIteration(const uchar4* const destImg,
        const bool* const interiorPixels,
        const bool* const borderPixels,
        const float3* const blendedVals1,
        const float3* const g,
        float3* const blendedVals2,
        const size_t numRows,
        const size_t numCols)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= numRows*numCols) return;
    if (!interiorPixels[tId]) return;
    float blendedSumX = 0.f;
    float blendedSumY = 0.f;
    float blendedSumZ = 0.f;
    float borderSumX = 0.f;
    float borderSumY = 0.f;
    float borderSumZ = 0.f;
    if (interiorPixels[tId-1])
    {
        blendedSumX += blendedVals1[tId-1].x;
        blendedSumY += blendedVals1[tId-1].y;
        blendedSumZ += blendedVals1[tId-1].z;
    }
    else
    {
        borderSumX += destImg[tId-1].x;
        borderSumY += destImg[tId-1].y;
        borderSumZ += destImg[tId-1].z;
    }

    if (interiorPixels[tId+1])
    {
        blendedSumX += blendedVals1[tId+1].x;
        blendedSumY += blendedVals1[tId+1].y;
        blendedSumZ += blendedVals1[tId+1].z;
    }
    else
    {
        borderSumX += destImg[tId+1].x;
        borderSumY += destImg[tId+1].y;
        borderSumZ += destImg[tId+1].z;
    }
    if (interiorPixels[tId-numCols])
    {
        blendedSumX += blendedVals1[tId-numCols].x;
        blendedSumY += blendedVals1[tId-numCols].y;
        blendedSumZ += blendedVals1[tId-numCols].z;
    }
    else
    {
        borderSumX += destImg[tId-numCols].x;
        borderSumY += destImg[tId-numCols].y;
        borderSumZ += destImg[tId-numCols].z;
    }
    if (interiorPixels[tId+numCols])
    {
        blendedSumX += blendedVals1[tId+numCols].x;
        blendedSumY += blendedVals1[tId+numCols].y;
        blendedSumZ += blendedVals1[tId+numCols].z;
    }
    else
    {
        borderSumX += destImg[tId+numCols].x;
        borderSumY += destImg[tId+numCols].y;
        borderSumZ += destImg[tId+numCols].z;
    }

    float next_valX = (blendedSumX+borderSumX+g[tId].x)/4.f;
    float next_valY = (blendedSumY+borderSumY+g[tId].y)/4.f;
    float next_valZ = (blendedSumZ+borderSumZ+g[tId].z)/4.f;
    blendedVals2[tId].x = min(255.f, max(0.f, next_valX));
    blendedVals2[tId].y = min(255.f, max(0.f, next_valY));
    blendedVals2[tId].z = min(255.f, max(0.f, next_valZ));
}

__global__
void copyBlendedValsToOutput(const float3* const blendedVals,
        const bool* const interiorPixels,
        uchar4* const blendedImg,
        const size_t size)
{
    const int tId = blockIdx.x*blockDim.x + threadIdx.x;
    if (tId >= size) return;
    if (!interiorPixels[tId]) return;
    blendedImg[tId].x = (unsigned char)blendedVals[tId].x;
    blendedImg[tId].y = (unsigned char)blendedVals[tId].y;
    blendedImg[tId].z = (unsigned char)blendedVals[tId].z;
}

void debugPrint(uchar4* data, size_t rows, size_t cols)
{
    std::ofstream outFile;
    outFile.open("log", std::ios::trunc);
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            outFile << "(" << (int)data[i*cols+j].x << ", "
                << (int)data[i*cols+j].y << ", "
                << (int)data[i*cols+j].z << ", "
                << (int)data[i*cols+j].w << ") ";
        }
        outFile << std::endl;
    }
    outFile.close();
}

void your_blend(const uchar4* const h_sourceImg,  //IN
        const size_t numRowsSource, const size_t numColsSource,
        const uchar4* const h_destImg, //IN
        uchar4* const h_blendedImg) //OUT
{
    const size_t imgSize = numRowsSource*numColsSource;
    unsigned int blockSize = 1024;
    unsigned int gridSize = (imgSize+blockSize-1) / blockSize;
    bool *d_mask;
    uchar4 *d_sourceImg;
    uchar4 *d_destImg;
    uchar4 *d_blendedImg;
    bool *d_interiorPixels;
    bool *d_borderPixels;
    float3 *d_g;
    float3 *d_blendedVals1;
    float3 *d_blendedVals2;

    // device mem allocs
    checkCudaErrors(cudaMalloc(&d_mask, sizeof(bool)*imgSize));
    checkCudaErrors(cudaMalloc(&d_sourceImg, sizeof(uchar4)*imgSize));
    checkCudaErrors(cudaMalloc(&d_destImg, sizeof(uchar4)*imgSize));
    checkCudaErrors(cudaMalloc(&d_blendedImg, sizeof(uchar4)*imgSize));
    checkCudaErrors(cudaMalloc(&d_interiorPixels, sizeof(bool)*imgSize));
    checkCudaErrors(cudaMalloc(&d_borderPixels, sizeof(bool)*imgSize));
    checkCudaErrors(cudaMalloc(&d_g, sizeof(float3)*imgSize));
    checkCudaErrors(cudaMalloc(&d_blendedVals1, sizeof(float3)*imgSize));
    checkCudaErrors(cudaMalloc(&d_blendedVals2, sizeof(float3)*imgSize));

    // memcpy to device
    checkCudaErrors(cudaMemcpy(d_sourceImg, h_sourceImg, sizeof(uchar4)*imgSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_destImg, h_destImg, sizeof(uchar4)*imgSize, cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();

    // data preprocess
    computeMask<<<gridSize, blockSize>>>(d_sourceImg, d_mask, imgSize);
    cudaDeviceSynchronize();
    computeBorderPixAndInteriorPix<<<gridSize, blockSize>>>(d_mask, d_interiorPixels,
            d_borderPixels, numRowsSource, numColsSource);
    cudaDeviceSynchronize();
    computeG<<<gridSize, blockSize>>>(d_sourceImg, d_interiorPixels, d_g,
            numRowsSource, numColsSource);
    copySourceImgToBlendedVals<<<gridSize, blockSize>>>(d_sourceImg, d_blendedVals1, d_blendedVals2, imgSize);
    cudaDeviceSynchronize();
    // start iterations
    for (int i=0; i<8000; i++)
    {
        computeIteration<<<gridSize, blockSize>>>(d_destImg, d_interiorPixels, d_borderPixels,
                d_blendedVals1, d_g, d_blendedVals2, numRowsSource, numColsSource);
        cudaDeviceSynchronize();
        std::swap(d_blendedVals1, d_blendedVals2); // output goes to d_blendedVals1
    }

    // copy dtsImg to outputImg
    checkCudaErrors(cudaMemcpy(d_blendedImg, d_destImg, sizeof(uchar4)*imgSize, cudaMemcpyDeviceToDevice));
    cudaDeviceSynchronize();
    copyBlendedValsToOutput<<<gridSize, blockSize>>>(d_blendedVals1, d_interiorPixels, d_blendedImg, imgSize);
    cudaDeviceSynchronize();
    checkCudaErrors(cudaMemcpy(h_blendedImg, d_blendedImg, sizeof(uchar4)*imgSize, cudaMemcpyDeviceToHost));
    cudaDeviceSynchronize();

    // free device memory
    checkCudaErrors(cudaFree(d_mask));
    checkCudaErrors(cudaFree(d_sourceImg));
    checkCudaErrors(cudaFree(d_destImg));
    checkCudaErrors(cudaFree(d_blendedImg));
    checkCudaErrors(cudaFree(d_interiorPixels));
    checkCudaErrors(cudaFree(d_borderPixels));
    checkCudaErrors(cudaFree(d_g));
    checkCudaErrors(cudaFree(d_blendedVals1));
    checkCudaErrors(cudaFree(d_blendedVals2));
}

void poissonBlendingSample()
{
    std::string input_source_file("datadropbox/input_source.png");
    std::string input_dest_file("datadropbox/input_dest.png");
    std::string output_file("datadropbox/output.png");
    uchar4 *h_sourceImg, *h_destImg, *h_blendedImg;
    size_t numRowsSource, numColsSource, numRowsDest, numColsDest;

    loadImageRGBA(input_source_file, &h_sourceImg, &numRowsSource, &numColsSource);
    loadImageRGBA(input_dest_file, &h_destImg, &numRowsDest, &numColsDest);

    //debugPrint(h_sourceImg, numRowsSource, numColsSource);

    assert(numRowsSource == numRowsDest);
    assert(numColsSource == numColsDest);
    h_blendedImg = new uchar4[numRowsSource*numColsSource];

    printf("Poisson Blending sample starts ... \n");

    your_blend(h_sourceImg, numRowsSource, numColsSource, h_destImg, h_blendedImg);

    printf("DONE!\n");

    saveImageRGBA(h_blendedImg, numRowsDest, numColsDest, output_file);

    delete h_sourceImg;
    delete h_destImg;
    delete h_blendedImg;
}
