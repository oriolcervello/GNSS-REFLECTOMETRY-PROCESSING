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
#include <fstream>
#include <iostream>
#include <cuda.h>
#include <stdio.h>
using namespace std;



int main() {
	//-------
	//DEVICE MANAGEMENT

	cudaDeviceReset();
	int count;
	cudaGetDeviceCount(&count);
	fprintf(stderr, "count devices: %i \n", count);
	for (int i = 0; i < count; i++) {
		cudaDeviceProp aa;
		cudaSetDevice(i);
		cudaGetDeviceProperties(&aa, i);

		
		
		fprintf(stderr, "Device %i ", i);
		fprintf(stderr, ":\n  Name: %s \n", aa.name);
		fprintf(stderr, "  maxThreadsPerBlock: %i \n", aa.maxThreadsPerBlock);
		fprintf(stderr, "  max dim of block of x: %i \n", aa.maxThreadsDim[0]);
		fprintf(stderr, "  max dim of block of y: %i \n", aa.maxThreadsDim[1]);
		
		
		size_t freeMem, totalMem;

		cudaMemGetInfo(&freeMem, &totalMem);

		fprintf(stderr, "  Memory: \n");
		fprintf(stderr, "   Free = %zu, Total = %zu\n", freeMem, totalMem);
	}
	
}
