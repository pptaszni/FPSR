find_package(CUDA QUIET REQUIRED)

SET(CUDASample_SRC simpleSurfaceWrite.cu reduceTraining.cu blellochTraining.cu
	radixSortTraining.cu)

set(GENCODE -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=compute_35)
set(GENCODE ${GENCODE} -gencode=arch=compute_20,code=sm_20)

cuda_add_library(CUDASample ${CUDASample_SRC} STATIC OPTIONS ${GENCODE})
cuda_add_cublas_to_target(CUDASample)
