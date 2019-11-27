//==========================================================================
// Author: Oriol Cervelló (oriol.cn [at] protonmail.com) 
//==========================================================================
// License: GNU GPLv3.0
// Copyright (C) 2019  Oriol Cervelló
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//==========================================================================
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
__global__ void applyDoppler(int, cufftComplex *, float, float, unsigned long long int, int, int, int, int);
__global__ void inchoerentSum(int, cufftComplex *, Npp32f *, int, int);
__global__ void scale(int, cufftComplex *, int);
__global__ void maskAndShift(char *, cuComplex *, int);
__global__ void savePeak(int, cufftComplex *, cufftComplex *, int, int, int, int *, int);
__global__ void selectMaxs(int, int, int, int *, Npp32f *);
__global__ void copyInt2Float(__int16 *, cuComplex *, int);

#endif
