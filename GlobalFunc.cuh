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


#ifndef GLOBFUNC
#define GLOBFUNC
__global__ void multip(int, cufftComplex *, cufftComplex *, int, bool);
__global__ void extendRefSignal(int, cufftComplex *, int);
__global__ void applyDoppler(int, cufftComplex *, float, float, unsigned long long, int, int, int, int);
__global__ void inchoerentSum(int, cufftComplex *, Npp32f *, int, int);
__global__ void scale(int, cufftComplex *, int);
__global__ void maskAndShift(char *, cuComplex *, int);
__global__ void savePeak(int, cufftComplex *, cufftComplex *, int, int, int, int *, int);
__global__ void selectMaxs(int, int, int, int *, Npp32f *);

#endif
