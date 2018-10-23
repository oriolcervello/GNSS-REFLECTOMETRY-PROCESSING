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
using namespace std;
#endif

#include"functions.cuh"

int main() {
	cudaDeviceReset();//reset device
	
	
	int fftsize = 32768;
	int numOfBytes = fftsize* 31;//bytes to read
	int samplesOfSignal = 4* numOfBytes;//complex samples
	int overlap = 32; //samples of overlaping 
	bool readbinary = 0, writebinary = 1;
	string* fileNames;
	fileNames = new string[2]{ "result.bin", "sin124xN.txt" };//names of files
	
	
	int numofFFTs = samplesOfSignal / (fftsize - overlap);
	int samples= numofFFTs * fftsize;//total samples needed
	if(samplesOfSignal > samples){samples = samplesOfSignal;}
	char *datainBytes, *deviceDatainBytes;
	cufftComplex *deviceDataFile1, *deviceDataFile2;
	cufftComplex *hostDataFile1, *hostDataFile2;

	int blockSize = 1024;
	int numBlocks;

	int const iterations = 5;
	long long read_elapsed_secs[iterations];
	long long fft_elapsed_secs[iterations];
	long long mul_elapsed_secs[iterations];
	long long ifft_elapsed_secs[iterations];
	long long write_elapsed_secs[iterations];
	long long elapsed_secs[iterations];
	long long shift_elapsed_secs[iterations];
	
	//ALLOCATE
	datainBytes = (char *)malloc(sizeof(char) * numOfBytes);
	CudaSafeCall(cudaMalloc((void **)&deviceDatainBytes, sizeof(char)*numOfBytes));
	hostDataFile1 = (cufftComplex *)malloc(sizeof(cufftComplex) * samples);
	hostDataFile2 = (cufftComplex *)malloc(sizeof(cufftComplex) * fftsize);
	CudaSafeCall(cudaMalloc(&deviceDataFile1, sizeof(cufftComplex)*samples));
	CudaSafeCall(cudaMalloc(&deviceDataFile2, sizeof(cufftComplex)*fftsize));
	cudaDeviceSynchronize();
	
	//FFT&IFFT PLANS
	cufftHandle plan;
	cufftHandle planref;
	cufftHandle inverseplan;
	planfftFunction(fftsize, numofFFTs, overlap, &plan);
	planfftFunction(fftsize, 1, 0, &planref);
	planifftFunction(fftsize, numofFFTs, overlap, &inverseplan);

	//LOOP
	int i;
	for (i = 0; i < iterations; i++) {
		auto Begin = std::chrono::high_resolution_clock::now();

		//READ DATA
		auto readdataBeg = chrono::high_resolution_clock::now();
		//readdata(samplesOfSignal, hostDataFile1, fileNames[0], readbinary);
		readrealdata(numOfBytes, datainBytes, fileNames[0]);
		
		auto elapsed_read = chrono::high_resolution_clock::now() - readdataBeg;
		readdata(fftsize - overlap, hostDataFile2, fileNames[1], readbinary);
		

		//CHECK: READED DATA 
		//cout << "read done\n";
		//writedata(fftsize, hostDataFile1, "rawsin.txt", writebinary);
		//writedata(fftsize, hostDataFile2, "rawsin2.txt", writebinary);

		//MEMORY FROM HOST TO DEVICE
		//CudaSafeCall(cudaMemcpy(deviceDataFile1, hostDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyHostToDevice));
		CudaSafeCall(cudaMemcpy(deviceDatainBytes, datainBytes, sizeof(char)*numOfBytes, cudaMemcpyHostToDevice));
		CudaSafeCall(cudaMemcpy(deviceDataFile2, hostDataFile2, sizeof(cufftComplex)*(fftsize - overlap), cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();

		//MASK AND SHIFT & EXTEND REFERENCE SIGNAL
		auto shiftBeg = chrono::high_resolution_clock::now();
		numBlocks = (numOfBytes + blockSize - 1) / blockSize;
		maskandshift << <numBlocks, blockSize >> > (deviceDatainBytes, deviceDataFile1, numOfBytes);
		CudaCheckError();
		

		numBlocks = (fftsize + blockSize - 1) / blockSize;
		
		extendRefSignal << <numBlocks, blockSize >> > (fftsize, deviceDataFile2, fftsize - overlap);
		CudaCheckError();
		cudaDeviceSynchronize();
		auto elapsed_shift = chrono::high_resolution_clock::now() - shiftBeg;

		//FFT
		auto fftBeg = chrono::high_resolution_clock::now();
		cufftSafeCall(cufftExecC2C(plan, deviceDataFile1, deviceDataFile1, CUFFT_FORWARD));
		cufftSafeCall(cufftExecC2C(planref, deviceDataFile2, deviceDataFile2, CUFFT_FORWARD));
		cudaDeviceSynchronize();
		auto elapsed_fft = chrono::high_resolution_clock::now() - fftBeg;

		//CHECK: MEMORY FROM DEVICE TO HOST (only for printing fft)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samples, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();

		//CHECK: FFT (only for printing fft)
		//fprintf(stderr, "%d FFt done of elements %d each\n", numofFFTs,fftsize);
		//writedata(fftsize, hostDataFile1, "fft.txt", writebinary);

		//COMPLEX CONJUGATE AND MULTIPLICATION
		numBlocks = (samples + blockSize - 1) / blockSize;
		auto mulBeg = chrono::high_resolution_clock::now();
		multip << <numBlocks, blockSize >> > (samples, deviceDataFile1, deviceDataFile2, fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();
		auto elapsed_mul = chrono::high_resolution_clock::now() - mulBeg;

		//CHECK: MEMORY FROM DEVICE TO HOST (only for printing multiplication result)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samples, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();

		//CHECK: MULTIPLICATION (only for printing multiplication result)
		//cout << "multiplication done\n";
		//writedata(fftsize, hostDataFile1, "mult.txt", writebinary);


		//IFFT (To obtain original again it has to be devided for the # of elements)
		auto ifftBeg = chrono::high_resolution_clock::now();
		cufftSafeCall(cufftExecC2C(inverseplan, deviceDataFile1, deviceDataFile1, CUFFT_INVERSE));
		cudaDeviceSynchronize();
		auto elapsed_ifft = chrono::high_resolution_clock::now() - ifftBeg;

		//MEMORY FROM HOST TO DEVICE FOR OUTPUT
		CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		
		//OUTPUT
		//cout << "IFFT done\n";
		auto writeBeg = chrono::high_resolution_clock::now();
		writedata(fftsize-overlap, hostDataFile1, "output.txt", writebinary);

		//ELAPSED TIME
		auto elapsed_write = chrono::high_resolution_clock::now() - writeBeg;
		auto elapsed_total = chrono::high_resolution_clock::now() - Begin;

		read_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_read).count();
		shift_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_shift).count();
		fft_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_fft).count();
		mul_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_mul).count();
		ifft_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_ifft).count();
		write_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_write).count();
		elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_total).count();
	}

	writetime(iterations, "Times_op3.txt", read_elapsed_secs, shift_elapsed_secs, fft_elapsed_secs,
		mul_elapsed_secs, ifft_elapsed_secs, write_elapsed_secs, elapsed_secs);

	cufftSafeCall(cufftDestroy(plan));
	cufftSafeCall(cufftDestroy(planref));
	cufftSafeCall(cufftDestroy(inverseplan));
	free(hostDataFile1);
	free(hostDataFile2);
	cudaFree(deviceDataFile1);
	cudaFree(deviceDataFile2);
	cudaDeviceReset();
	delete[] fileNames;
	
	return 0;
}