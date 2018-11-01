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

#include"functions.cuh"

int main(int argc, const char* argv[]) {
	cudaDeviceReset();//reset device
	
	//INPUTS
	int fftsize, fSampling, numofFFTs, overlap, quantofAverageIncoherent;
	bool readbinary = 1, writebinary = 0;
	int const numofDataLines = atoi(argv[2]);//substitut d'iterations
	string *fileNames;
	int *dataOffsetBeg, *dataOffsetEnd;
	int *doppler;

	fileNames = new string[numofDataLines];
	dataOffsetBeg = new int[numofDataLines];
	dataOffsetEnd = new int[numofDataLines];
	doppler = new int[numofDataLines];
	


	readConfig(argv[1], numofDataLines, &fftsize, &numofFFTs, &overlap, &fSampling, &quantofAverageIncoherent, &readbinary, &writebinary, dataOffsetBeg, dataOffsetEnd, doppler, fileNames);
	checkInputConfig(argc, argv, numofDataLines, fftsize, numofFFTs, overlap, fSampling, quantofAverageIncoherent, readbinary, writebinary, dataOffsetBeg, dataOffsetEnd, doppler, fileNames);




	//OTHER DECLARATIONS
	int samplesOfSignal = (numofFFTs * (fftsize-overlap))+overlap;//samples of data
	int samplesWithOverlap= numofFFTs * fftsize;//total samples needed
	if(samplesOfSignal > samplesWithOverlap){ samplesWithOverlap = samplesOfSignal;}
	int inchoerentNumofFFT = numofFFTs/ quantofAverageIncoherent;

	int blockSize = 1024;//threads per block
	int numBlocks, nBufferSize, samplePhaseMantain,i;
	
	
	int *devicearrayPos;
	cufftComplex *deviceDataFile1, *deviceDataFile2, *hostDataFile1, *hostDataFile2;
	Npp32f *deviceIncoherentSum, *devicearrayMaxs, *devicearrayStd;
	Npp8u * pDeviceBuffer;
	
	int const maxIter = 5;
	long long read_elapsed_secs[maxIter], fft_elapsed_secs[maxIter], mul_elapsed_secs[maxIter], ifft_elapsed_secs[maxIter],write_elapsed_secs[maxIter], elapsed_secs[maxIter], shift_elapsed_secs[maxIter];
	
	//ALLOCATE
	int *hostarrayPos = new int[inchoerentNumofFFT];
	Npp32f *hostarrayMaxs = new Npp32f[inchoerentNumofFFT];
	Npp32f *hostarrayStd = new Npp32f[inchoerentNumofFFT];
	hostDataFile1 = (cufftComplex *)malloc(sizeof(cufftComplex) * samplesWithOverlap);
	hostDataFile2 = (cufftComplex *)malloc(sizeof(cufftComplex) * fftsize);
	CudaSafeCall(cudaMalloc(&deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap));
	CudaSafeCall(cudaMalloc(&deviceDataFile2, sizeof(cufftComplex)*fftsize));
	CudaSafeCall(cudaMalloc(&deviceIncoherentSum, sizeof(Npp32f)*inchoerentNumofFFT*fftsize));
	CudaSafeCall(cudaMalloc(&devicearrayPos, sizeof(int)*inchoerentNumofFFT));
	CudaSafeCall(cudaMalloc(&devicearrayMaxs, sizeof(Npp32f)*inchoerentNumofFFT));
	CudaSafeCall(cudaMalloc(&devicearrayStd, sizeof(Npp32f)*inchoerentNumofFFT));
	nppsSumGetBufferSize_32f(fftsize, &nBufferSize);
	CudaSafeCall(cudaMalloc((void **)(&pDeviceBuffer), nBufferSize));
	cudaDeviceSynchronize();
	
	//FFT&IFFT PLANS
	cufftHandle plan;
	cufftHandle planref;
	cufftHandle inverseplan;
	planfftFunction(fftsize, numofFFTs, overlap, &plan);
	planfftFunction(fftsize, 1, 0, &planref);
	planifftFunction(fftsize, numofFFTs, 0, &inverseplan);

	//LOOP
	for (i = 0; i < numofDataLines; i++) {
		auto Begin = std::chrono::high_resolution_clock::now();

		//READ DATA
		auto readdataBeg = chrono::high_resolution_clock::now();
		//readdata(samplesOfSignal, hostDataFile1, fileNames[i], readbinary);
		readdatabinary(dataOffsetEnd[i]-dataOffsetBeg[i], dataOffsetBeg[i], hostDataFile1, fileNames[i]);
		readdata(fftsize - overlap, hostDataFile2, "prn_L1CA_32_100.bin", readbinary);
		auto elapsed_read = chrono::high_resolution_clock::now() - readdataBeg;

		//CHECK: READED DATA 
		//writedata(samplesOfSignal/2, hostDataFile1, "rawsin.txt", writebinary);
		//writedata(fftsize- overlap, hostDataFile2, "rawsin2.txt", writebinary);

		//MEMORY FROM HOST TO DEVICE
		CudaSafeCall(cudaMemcpy(deviceDataFile1, hostDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyHostToDevice));
		CudaSafeCall(cudaMemcpy(deviceDataFile2, hostDataFile2, sizeof(cufftComplex)*(fftsize - overlap), cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();

		//MULTIPLY BY DOPPLER
		samplePhaseMantain = (i * fftsize*numofFFTs)%fSampling;
		numBlocks = (samplesOfSignal + blockSize - 1) / blockSize;
		applyDoppler << <numBlocks, blockSize >> > (samplesOfSignal, deviceDataFile1, doppler[i], fSampling, samplePhaseMantain);
		CudaCheckError();
		cudaDeviceSynchronize();
	
		//CHECK: doppler (only for printing doppler)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesOfSignal/2, hostDataFile1, "dopplerout.txt", writebinary);
		
		//EXTEND REFERENCE SIGNAL
		numBlocks = (fftsize + blockSize - 1) / blockSize;
		auto shiftBeg = chrono::high_resolution_clock::now();
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

		//CHECK: FFT (only for printing fft)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "fft.txt", writebinary);

		//COMPLEX CONJUGATE AND MULTIPLICATION
		numBlocks = (samplesWithOverlap + blockSize - 1) / blockSize;
		auto mulBeg = chrono::high_resolution_clock::now();
		multip << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, deviceDataFile2, fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();
		auto elapsed_mul = chrono::high_resolution_clock::now() - mulBeg;

		//CHECK: MULTIPLICATION (only for printing multiplication result)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "mult.txt", writebinary);
		
		//IFFT (To obtain original again it has to be devided for the # of elements)
		auto ifftBeg = chrono::high_resolution_clock::now();
		cufftSafeCall(cufftExecC2C(inverseplan, deviceDataFile1, deviceDataFile1, CUFFT_INVERSE));
		cudaDeviceSynchronize();
		auto elapsed_ifft = chrono::high_resolution_clock::now() - ifftBeg;
		
		//SCALE
		numBlocks = (samplesWithOverlap + blockSize - 1) / blockSize;		
		scale << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, fftsize);
		cudaDeviceSynchronize();

		//INCOHERENT SUM
		numBlocks = (inchoerentNumofFFT*fftsize + blockSize - 1) / blockSize;
		inchoerentSum << <numBlocks, blockSize >> > (inchoerentNumofFFT*fftsize, deviceDataFile1, deviceIncoherentSum, quantofAverageIncoherent, fftsize);
		cudaDeviceSynchronize();

		//MAXIMUM AND STD
		maxAndStd(inchoerentNumofFFT, deviceIncoherentSum, fftsize, devicearrayMaxs, devicearrayStd,devicearrayPos, pDeviceBuffer);

		//CHECK: IFFT OR incho (not both at the same time)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceIncoherentSum, sizeof(Npp32f)*inchoerentNumofFFT*fftsize, cudaMemcpyDeviceToHost)); //TO PRINT INCHO SUM
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost)); //TO PRINT IFFT RESULT
		//cudaDeviceSynchronize();
		//writeIncohtxt(inchoerentNumofFFT*fftsize, hostDataFile1, "incoh.txt");//TO PRINT INCHO SUM
		//writedata(samplesWithOverlap, hostDataFile1,  "result.txt", writebinary);//TO PRINT IFFT RESULT 


		//MEMORY FROM HOST TO DEVICE FOR OUTPUT
		CudaSafeCall(cudaMemcpy(hostarrayMaxs, devicearrayMaxs, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(hostarrayStd, devicearrayStd, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(hostarrayPos, devicearrayPos, sizeof(int)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		
		//OUTPUT
		//cout<< hostDataFile1[0].x << " incho\n";
		auto writeBeg = chrono::high_resolution_clock::now();
		writeMaxstxt(inchoerentNumofFFT, hostarrayMaxs, hostarrayPos, hostarrayStd, "Maximums.txt");

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

	writetime(numofDataLines, "Times_op3.txt", read_elapsed_secs, shift_elapsed_secs, fft_elapsed_secs,
		mul_elapsed_secs, ifft_elapsed_secs, write_elapsed_secs, elapsed_secs);

	//FREE MEMORY
	cufftSafeCall(cufftDestroy(plan));
	cufftSafeCall(cufftDestroy(planref));
	cufftSafeCall(cufftDestroy(inverseplan));
	cudaFree(deviceDataFile1);
	cudaFree(deviceDataFile2);
	cudaFree(deviceIncoherentSum);
	cudaFree(devicearrayPos);
	cudaFree(devicearrayMaxs);
	cudaFree(pDeviceBuffer);
	cudaFree(devicearrayStd);
	cudaDeviceReset();
	delete[] fileNames;
	delete[] hostarrayPos;
	delete[] hostarrayMaxs;
	delete[] hostarrayStd;
	delete[] hostDataFile2;
	delete[] hostDataFile1;
	delete[] dataOffsetBeg;
	delete[] dataOffsetEnd;
	delete[] doppler;

	return 0;
}