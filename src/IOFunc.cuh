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







#ifndef IOFUNC
#define IOFUNC
void readdata(int, int, cufftComplex *, string);
void readRealData(int, int, int, char *, string);
void writedata(int, cufftComplex *, string);
void writeMaxs(int, Npp32f *, int *, Npp32f *, Npp32f *, int, string, int, int, int, int);
void writeIncoh(int, cuComplex *, string);
void writetime(int, string, long long *, long long *, long long *, long long *, long long *, long long *,
	long long *, long long *, long long *, long long *, long long *, long long *, long long *);


#endif
