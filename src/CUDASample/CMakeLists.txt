find_package(CUDA QUIET REQUIRED)

SET(CUDASample_SRC simpleSurfaceWrite.cu reduceTraining.cu blellochTraining.cu
    radixSortTraining.cu fastHistogramTraining.cu poissonBlending.cu)

set(GENCODE -gencode=arch=compute_32,code=sm_32)
set(GENCODE ${GENCODE} -gencode=arch=compute_52,code=sm_52)

cuda_add_library(CUDASample ${CUDASample_SRC} STATIC OPTIONS ${GENCODE})
cuda_add_cublas_to_target(CUDASample)
