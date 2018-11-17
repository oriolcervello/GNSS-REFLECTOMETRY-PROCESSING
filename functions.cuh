#ifndef LIBRARIES
#define LIBRARIES
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cufft.h>
#include <ctime>
#include <chrono>
#include <npps.h>
#define PI 3.14159265
using namespace std;
#endif

#ifndef FUNC
#define FUNC



//ERROR CHECK
#define CUDA_ERROR_CHECK

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )
#define cufftSafeCall(err)      __cufftSafeCall(err, __FILE__, __LINE__)

static const char *_cudaGetErrorEnum(cufftResult error)
{
	switch (error)
	{
	case CUFFT_SUCCESS:
		return "CUFFT_SUCCESS";

	case CUFFT_INVALID_PLAN:
		return "CUFFT_INVALID_PLAN";

	case CUFFT_ALLOC_FAILED:
		return "CUFFT_ALLOC_FAILED";

	case CUFFT_INVALID_TYPE:
		return "CUFFT_INVALID_TYPE";

	case CUFFT_INVALID_VALUE:
		return "CUFFT_INVALID_VALUE";

	case CUFFT_INTERNAL_ERROR:
		return "CUFFT_INTERNAL_ERROR";

	case CUFFT_EXEC_FAILED:
		return "CUFFT_EXEC_FAILED";

	case CUFFT_SETUP_FAILED:
		return "CUFFT_SETUP_FAILED";

	case CUFFT_INVALID_SIZE:
		return "CUFFT_INVALID_SIZE";

	case CUFFT_UNALIGNED_DATA:
		return "CUFFT_UNALIGNED_DATA";
	}

	return "<unknown>";
}

inline void __cufftSafeCall(cufftResult err, const char *file, const int line)
{
#ifdef CUDA_ERROR_CHECK
	if (CUFFT_SUCCESS != err) {
		fprintf(stderr, "CUFFT error in file '%s', line %d\n : %s\n", __FILE__, __LINE__, _cudaGetErrorEnum(err));
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
#endif

	return;
}

inline void __cudaSafeCall(cudaError err, const char *file, const int line)
{
#ifdef CUDA_ERROR_CHECK
	if (cudaSuccess != err)
	{
		fprintf(stderr, "CudaSafeCall error in file '%s' in line %i : %s.\n",
			__FILE__, __LINE__, cudaGetErrorString(err));
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
#endif

	return;
}

inline void __cudaCheckError(const char *file, const int line)
{
#ifdef CUDA_ERROR_CHECK
	cudaError err = cudaGetLastError();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",
			__FILE__, __LINE__, cudaGetErrorString(err));
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
#endif

	return;
}
//END ERROR CHECK

void readConfig(const char *, int, int *, int *, int *, int *, int*, int *, int *, int*
	,  int *, int *, int *, string *, string *,int *, int *);
void checkInputConfig(int , const char **, int , int , int , int , int , int, int, int, int
	,  int *, int *, int *, string *, string , int , int );
void planifftFunction(int , int , int , cufftHandle *);
void planfftFunction(int , int , int , cufftHandle *);
void readdata(int,int , cufftComplex *, string);
void readRealData(int, int,int, char *, string);
void writedata(int, cufftComplex *, string);
void writeMaxs(int, Npp32f *, int *,Npp32f *,int, string);
void writeIncoh(int, cuComplex *, string);
void writetime(int, string, long long *, long long *, long long *, long long *, long long *, long long *,
	long long *, long long *, long long *, long long *, long long *, long long *, long long *);
void maxCompute(int, Npp32f *, int,  Npp32f *, int *, Npp8u *);
void stdCompute(int, Npp32f *, int, Npp32f *, int *, Npp8u *, int);

__global__ void multip(int , cufftComplex *, cufftComplex *,int);
__global__ void extendRefSignal(int, cufftComplex *, int);
__global__ void applyDoppler(int, cufftComplex *, float, float, unsigned long long);
__global__ void inchoerentSum(int, cufftComplex *, Npp32f *, int, int );
__global__ void scale(int, cufftComplex *, int);
__global__ void maskAndShift(char *, cuComplex *, int);
__global__ void savePeak(int, cufftComplex *, cufftComplex *, int, int, int, int *);




#endif