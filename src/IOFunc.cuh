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







#ifndef IOFUNC
#define IOFUNC
void readdata(int, int, cufftComplex *, string);
void readRealData(int, int, int, char *, string);
void readRealData2files(int, int, int, int, int, char *, string, string);
void readdataInt(int , int , __int16 *, string );
void writedata(int, cufftComplex *, string);
void writeMaxs(int, Npp32f *, int *, Npp32f *, Npp32f *, float, string, int, int, int, int);
void writeIncoh(int, cuComplex *, string);
void writetime(int, string, long long *, long long *, long long *, long long *, long long *, long long *,
	long long *, long long *, long long *, long long *, long long *, long long *, long long *);


#endif
