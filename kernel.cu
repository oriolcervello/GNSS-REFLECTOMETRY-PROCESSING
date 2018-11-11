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
	
	//READ CONFIG FILE
	int fftsize, fSampling, numofFFTs, overlap, quantofAverageIncoherent, blockSize, peakRangeStd, peakSamplesToSave;
	int const numofDataLines = atoi(argv[2]);//substitut d'iterations
	string *fileDataNames, *fileRefNames;
	int *dataOffsetBeg, *dataOffsetEnd;
	int *doppler;

	fileDataNames = new string[numofDataLines];
	fileRefNames = new string[numofDataLines];
	dataOffsetBeg = new int[numofDataLines];
	dataOffsetEnd = new int[numofDataLines];
	doppler = new int[numofDataLines];

	readConfig(argv[1], numofDataLines, &fftsize, &numofFFTs, &overlap, &fSampling, &blockSize, &peakRangeStd, &peakSamplesToSave, &quantofAverageIncoherent, dataOffsetBeg, dataOffsetEnd, doppler, fileDataNames, fileRefNames);
	checkInputConfig(argc, argv, numofDataLines, fftsize, numofFFTs, overlap, fSampling, blockSize, peakRangeStd, peakSamplesToSave, quantofAverageIncoherent, dataOffsetBeg, dataOffsetEnd, doppler, fileDataNames, fileRefNames);

	//OTHER DECLARATIONS
	int samplesOfSignal = (numofFFTs * (fftsize-overlap))+overlap;//samples of complex data
	int bytesToRead = samplesOfSignal/4;
	if (samplesOfSignal % 4 != 0) { cout << "Warning bytesToRead rounded toward negative infinity: samplesOfSignal%4!=0 \n "; }
	int samplesWithOverlap= numofFFTs * fftsize;//total samples needed
	if(samplesOfSignal > samplesWithOverlap){ samplesWithOverlap = samplesOfSignal;}
	int inchoerentNumofFFT = numofFFTs/ quantofAverageIncoherent;
	if (numofFFTs % quantofAverageIncoherent != 0) {
		cout << "Error: numofFFTs / quantofAverageIncoherent != 0\n ";
		exit(-1);
	}

	string outputName;
	int numBlocks, nBufferSize,i;
	unsigned long long samplePhaseMantain;

	char *hostBytesOfData, *deviceBytesOfData;
	int *devicearrayPos,*hostarrayPos;
	cufftComplex *deviceDataFile1, *deviceDataFile2, *hostDataFile1, *hostDataFile2, *deviceDataToSave;
	Npp32f *deviceIncoherentSum, *devicearrayMaxs, *devicearrayStd, *hostarrayMaxs, *hostarrayStd;
	Npp8u * pDeviceBuffer;
	
	long long *read_elapsed_secs,*write_elapsed_secs, *elapsed_secs;
	
	//ALLOCATE
	read_elapsed_secs = new long long[numofDataLines];
	write_elapsed_secs = new long long[numofDataLines];
	elapsed_secs = new long long[numofDataLines];
	hostBytesOfData = (char *)malloc(sizeof(char) * bytesToRead);
	hostarrayPos = new int[inchoerentNumofFFT];
	hostarrayMaxs = new Npp32f[inchoerentNumofFFT];
	hostarrayStd = new Npp32f[inchoerentNumofFFT];
	hostDataFile1 = (cufftComplex *)malloc(sizeof(cufftComplex) * samplesWithOverlap);
	hostDataFile2 = (cufftComplex *)malloc(sizeof(cufftComplex) * fftsize);
	CudaSafeCall(cudaMalloc(&deviceBytesOfData, sizeof(char)*bytesToRead));
	CudaSafeCall(cudaMalloc(&deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap));
	CudaSafeCall(cudaMalloc(&deviceDataToSave, sizeof(cufftComplex)*peakSamplesToSave*numofFFTs));
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
	cudaDeviceSynchronize();

	//LOOP
	for (i = 0; i < numofDataLines; i++) {
		
		auto Begin = std::chrono::high_resolution_clock::now();
		//READ DATA
		readdata(dataOffsetEnd[i]-dataOffsetBeg[i], dataOffsetBeg[i], hostDataFile1, fileDataNames[i]);
		readdata(fftsize - overlap,0, hostDataFile2, fileRefNames[i]);

		/*if (dataOffsetEnd[i] - dataOffsetBeg[i] > bytesToRead) { cout << "Indices of reading in config file exceed bytesToRead decleared"; }
		if ((dataOffsetEnd[i] - dataOffsetBeg[i])%(fftsize*quantofAverageIncoherent)!=0) {
		cout << "Warning length of data won't complete last incho sum in DATALINE: "<<i<<"\n"; }
		readRealData(dataOffsetEnd[i] - dataOffsetBeg[i], dataOffsetBeg[i],bytesToRead, hostBytesOfData, fileNames[i]);*/
		auto elapsed_read = chrono::high_resolution_clock::now() - Begin;

		CudaSafeCall(cudaMemcpy(deviceDataFile1, hostDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyHostToDevice));
		CudaSafeCall(cudaMemcpy(deviceDataFile2, hostDataFile2, sizeof(cufftComplex)*(fftsize - overlap), cudaMemcpyHostToDevice));
		/*CudaSafeCall(cudaMemcpy(deviceBytesOfData, hostBytesOfData, sizeof(char)*bytesToRead, cudaMemcpyHostToDevice));*/
		cudaDeviceSynchronize();
		
		//CHECK: READED DATA 
		//writedata(samplesOfSignal/2, hostDataFile1, "rawsin.txt");
		//writedata(fftsize- overlap, hostDataFile2, "rawsin2.txt");

		//MASK AND SHIFT
		/*numBlocks = (bytesToRead + blockSize - 1) / blockSize;
		maskAndShift << <numBlocks, blockSize >> > (deviceBytesOfData, deviceDataFile1, bytesToRead);
		CudaCheckError();
		cudaDeviceSynchronize();
		*/
		//MULTIPLY BY DOPPLER
		samplePhaseMantain = (i * fftsize*numofFFTs);// %fSampling;----
		numBlocks = (samplesOfSignal + blockSize - 1) / blockSize;
		applyDoppler << <numBlocks, blockSize >> > (samplesOfSignal, deviceDataFile1, doppler[i], fSampling, samplePhaseMantain);
		CudaCheckError();
		cudaDeviceSynchronize();
	
		//CHECK: doppler (only for printing doppler)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesOfSignal/2, hostDataFile1, "dopplerout.txt");
		
		//EXTEND REFERENCE SIGNAL
		numBlocks = (fftsize + blockSize - 1) / blockSize;
		extendRefSignal << <numBlocks, blockSize >> > (fftsize, deviceDataFile2, fftsize - overlap);
		CudaCheckError();
		cudaDeviceSynchronize();

		//FFT
		cufftSafeCall(cufftExecC2C(plan, deviceDataFile1, deviceDataFile1, CUFFT_FORWARD));
		cufftSafeCall(cufftExecC2C(planref, deviceDataFile2, deviceDataFile2, CUFFT_FORWARD));
		cudaDeviceSynchronize();

		//CHECK: FFT (only for printing fft)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "fft.txt");

		//COMPLEX CONJUGATE AND MULTIPLICATION
		numBlocks = (samplesWithOverlap + blockSize - 1) / blockSize;
		multip << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, deviceDataFile2, fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();

		//CHECK: MULTIPLICATION (only for printing multiplication result)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "mult.txt");
		
		//IFFT
		cufftSafeCall(cufftExecC2C(inverseplan, deviceDataFile1, deviceDataFile1, CUFFT_INVERSE));
		cudaDeviceSynchronize();
		
		//SCALE (To take back original signal it has to be devided for the fftsize)
		numBlocks = (samplesWithOverlap + blockSize - 1) / blockSize;		
		scale << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();
		
		//CHECK: IFFT 
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost)); 
		//cudaDeviceSynchronize();	
		//writedata(samplesWithOverlap, hostDataFile1,  "IFFT-result.bin");

		//INCOHERENT SUM
		numBlocks = (inchoerentNumofFFT*fftsize + blockSize - 1) / blockSize;
		inchoerentSum << <numBlocks, blockSize >> > (inchoerentNumofFFT*fftsize, deviceDataFile1, deviceIncoherentSum, quantofAverageIncoherent, fftsize);
		CudaCheckError(); 
		cudaDeviceSynchronize();

		//CHECK: INCOHERENT
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceIncoherentSum, sizeof(Npp32f)*inchoerentNumofFFT*fftsize, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writeIncoh(inchoerentNumofFFT*fftsize, hostDataFile1, "incoh.bin");
	
		//MAXIMUM
		maxCompute(inchoerentNumofFFT, deviceIncoherentSum, fftsize, devicearrayMaxs, devicearrayPos, pDeviceBuffer);
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostarrayPos, devicearrayPos, sizeof(int)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(hostarrayMaxs, devicearrayMaxs, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		
		//SAVE PEAKS
		numBlocks = (numofFFTs*peakSamplesToSave + blockSize - 1) / blockSize;
		savePeak << <numBlocks, blockSize >> > (numofFFTs, deviceDataFile1, deviceDataToSave, peakSamplesToSave, quantofAverageIncoherent, fftsize, devicearrayPos);
		CudaCheckError();
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataToSave, sizeof(cuComplex)*numofFFTs*peakSamplesToSave, cudaMemcpyDeviceToHost));
		
		//STD
		stdCompute(inchoerentNumofFFT, deviceIncoherentSum, fftsize, devicearrayStd, hostarrayPos, pDeviceBuffer, peakRangeStd);
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostarrayStd, devicearrayStd, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		
		//OUTPUT
		auto writeBeg = chrono::high_resolution_clock::now();
		writeMaxs(inchoerentNumofFFT, hostarrayMaxs, hostarrayPos, hostarrayStd, "results/Maximums.txt");
		outputName = "results/PeaksIteration"+ to_string(i);
		outputName = outputName + ".bin";
		cout << outputName << "\n";
		writedata(numofFFTs*peakSamplesToSave, hostDataFile1, outputName);
	
		//ELAPSED TIME
		auto elapsed_write = chrono::high_resolution_clock::now() - writeBeg;
		auto elapsed_total = chrono::high_resolution_clock::now() - Begin;

		read_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_read).count();
		write_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_write).count();
		elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_total).count();
	}

	writetime(numofDataLines, "results/Times_op3.txt", read_elapsed_secs, write_elapsed_secs, elapsed_secs);

	//FREE MEMORY
	cufftSafeCall(cufftDestroy(plan));
	cufftSafeCall(cufftDestroy(planref));
	cufftSafeCall(cufftDestroy(inverseplan));
	cudaFree(deviceDataFile1);
	cudaFree(deviceDataFile2);
	cudaFree(deviceIncoherentSum);
	cudaFree(devicearrayPos);
	cudaFree(deviceBytesOfData);
	cudaFree(devicearrayMaxs);
	cudaFree(deviceDataToSave);
	cudaFree(pDeviceBuffer);
	cudaFree(devicearrayStd);
	cudaDeviceReset();
	delete[] fileDataNames;
	delete[] fileRefNames;
	delete[] hostBytesOfData;
	delete[] hostarrayPos;
	delete[] hostarrayMaxs;
	delete[] hostarrayStd;
	delete[] hostDataFile2;
	delete[] hostDataFile1;
	delete[] dataOffsetBeg;
	delete[] dataOffsetEnd;
	delete[] doppler;
	delete[] read_elapsed_secs;
	delete[] write_elapsed_secs;
	delete[] elapsed_secs;
	
	return 0;
}