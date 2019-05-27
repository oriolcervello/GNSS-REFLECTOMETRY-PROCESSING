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

#include"HostFunc.cuh"
#include "GlobalFunc.cuh"
#include "IOFunc.cuh"

int main(int argc, const char* argv[]) {
	cudaDeviceReset();//reset device

	//READ CONFIG FILE
	int fftsize, fSampling, numofFFTs, overlap, quantofAverageIncoherent, blockSize, peakRangeStd, peakSamplesToSave,
		samplesAvoidMax,ddmRes, ddmQuant;
	int const numofDataLines = atoi(argv[2]);//substitut d'iterations
	bool interferometic,writeoutputs=1;
	string *fileDataNames, *fileRefName, resultDirectory;
	int *dataOffsetBeg, *dataOffsetEnd,*dataOffsetBegInterferometric, *typeOfDataline, *dataOffsetEndInterferometric;
	float *doppler;

	fileRefName = new string[numofDataLines];
	fileDataNames = new string[numofDataLines];
	dataOffsetBeg = new int[numofDataLines];
	dataOffsetBegInterferometric = new int[numofDataLines];
	dataOffsetEndInterferometric = new int[numofDataLines];
	dataOffsetEnd = new int[numofDataLines];
	typeOfDataline = new int[numofDataLines];
	doppler = new float[numofDataLines];

	readConfig(argv[1], numofDataLines, &fftsize, &numofFFTs, &overlap, &fSampling, &blockSize, &peakRangeStd, &peakSamplesToSave, &quantofAverageIncoherent, dataOffsetBeg, dataOffsetEnd, doppler, fileDataNames, fileRefName, &ddmRes, &ddmQuant,&interferometic,dataOffsetBegInterferometric, &samplesAvoidMax, &resultDirectory,&writeoutputs, typeOfDataline, dataOffsetEndInterferometric);
	checkInputConfig(argc, argv, numofDataLines, fftsize, numofFFTs, overlap, fSampling, blockSize, peakRangeStd, peakSamplesToSave, quantofAverageIncoherent, dataOffsetBeg, dataOffsetEnd, doppler, fileDataNames, fileRefName, ddmRes, ddmQuant, interferometic, dataOffsetBegInterferometric, samplesAvoidMax, resultDirectory, writeoutputs, typeOfDataline, dataOffsetEndInterferometric);

	//OTHER DECLARATIONS
	int  originalSamplesOfSignal = (numofFFTs * (fftsize - overlap)) + overlap;//samples of complex data
	int samplesOfSignal = originalSamplesOfSignal *ddmQuant;//samples of complex data
	int bytesToRead = originalSamplesOfSignal/4;
	if (originalSamplesOfSignal % 4 != 0) { cout << "Warning bytesToRead rounded toward negative infinity: samplesOfSignal%4!=0 \n "; }
	int samplesWithOverlap= (numofFFTs * fftsize)*ddmQuant;//total samples needed
	if(samplesOfSignal > samplesWithOverlap){ samplesWithOverlap = samplesOfSignal;}
	int inchoerentNumofFFT = (numofFFTs/ quantofAverageIncoherent)*ddmQuant;
	if (numofFFTs % quantofAverageIncoherent != 0) {
		cout << "Error: numofFFTs / quantofAverageIncoherent != 0\n ";
		exit(-1);
	}

	string outputName;
	int numBlocks, nMaxBufferSize,nStdBufferSize,i,k,typ2=0,samplesDoppler= samplesOfSignal;
	int stdLength = (fftsize / 2) - ((peakSamplesToSave) / 2) - 1;
	unsigned long long int samplePhaseMantain;
	if (ddmQuant > 1) {
		samplesDoppler = samplesWithOverlap;
	}
	

	char *hostBytesOfData, *deviceBytesOfData;
	int *devicearrayPos,*hostarrayPos;
	cufftComplex *deviceDataFile1, *deviceDataFile2, *hostDataFile1, *hostDataFile2, *deviceDataToSave;
	Npp32f *deviceIncoherentSum, *devicearrayMaxs, *devicearrayStd,*devicearrayMean,*hostarrayMean, *hostarrayMaxs, *hostarrayStd;
	Npp8u *pStdDeviceBuffer,*pMaxDeviceBuffer;

	cufftHandle plan;
	
	long long *read_elapsed_secs,*write_elapsed_secs, *elapsed_secs, *mask_elapsed_secs, *doppler_elapsed_secs, 
		 *fft_elapsed_secs, *mult_elapsed_secs,*ifft_elapsed_secs, *extenddop_elapsed_secs, *incho_elapsed_secs
		, *max_elapsed_secs, *savep_elapsed_secs, *std_elapsed_secs;
	chrono::nanoseconds elapsed_read, mask_elapsed, extenddop_elapsed;

	int device2quant;
	if (interferometic == true) {
		device2quant = samplesWithOverlap;
	}
	else
	{
		device2quant = fftsize;
	}


	//ALLOCATE
	read_elapsed_secs = new long long[numofDataLines];
	mask_elapsed_secs = new long long[numofDataLines];
	doppler_elapsed_secs = new long long[numofDataLines];
	fft_elapsed_secs = new long long[numofDataLines];
	mult_elapsed_secs = new long long[numofDataLines];
	ifft_elapsed_secs = new long long[numofDataLines];
	extenddop_elapsed_secs = new long long[numofDataLines];
	incho_elapsed_secs = new long long[numofDataLines];
	max_elapsed_secs = new long long[numofDataLines];
	savep_elapsed_secs = new long long[numofDataLines];
	std_elapsed_secs = new long long[numofDataLines];
	write_elapsed_secs = new long long[numofDataLines];
	elapsed_secs = new long long[numofDataLines];

	hostBytesOfData = (char *)malloc(sizeof(char) * bytesToRead);
	hostarrayPos = new int[inchoerentNumofFFT];
	hostarrayMaxs = new Npp32f[inchoerentNumofFFT];
	hostarrayStd = new Npp32f[inchoerentNumofFFT];
	hostarrayMean = new Npp32f[inchoerentNumofFFT];
	hostDataFile1 = (cufftComplex *)malloc(sizeof(cufftComplex) * samplesWithOverlap);
	hostDataFile2 = (cufftComplex *)malloc(sizeof(cufftComplex) * device2quant);
	CudaSafeCall(cudaMalloc(&deviceBytesOfData, sizeof(char)*bytesToRead));
	CudaSafeCall(cudaMalloc(&deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap));
	CudaSafeCall(cudaMalloc(&deviceDataToSave, sizeof(cufftComplex)*peakSamplesToSave*numofFFTs*ddmQuant));
	CudaSafeCall(cudaMalloc(&deviceDataFile2, sizeof(cufftComplex)*device2quant));
	CudaSafeCall(cudaMalloc(&deviceIncoherentSum, sizeof(Npp32f)*inchoerentNumofFFT*fftsize));
	CudaSafeCall(cudaMalloc(&devicearrayPos, sizeof(int)*inchoerentNumofFFT));
	CudaSafeCall(cudaMalloc(&devicearrayMaxs, sizeof(Npp32f)*inchoerentNumofFFT));
	CudaSafeCall(cudaMalloc(&devicearrayStd, sizeof(Npp32f)*inchoerentNumofFFT));
	CudaSafeCall(cudaMalloc(&devicearrayMean, sizeof(Npp32f)*inchoerentNumofFFT));
	nppsMeanStdDevGetBufferSize_32f(stdLength, &nStdBufferSize);
	CudaSafeCall(cudaMalloc((void **)(&pStdDeviceBuffer), nStdBufferSize));
	nppsMaxGetBufferSize_32f(fftsize, &nMaxBufferSize);
	CudaSafeCall(cudaMalloc((void **)(&pMaxDeviceBuffer), nMaxBufferSize));
	cudaDeviceSynchronize();
	
	//MEMORY INFO
	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);
	cout<< "\n-MEMORY: \n";
	cout<< "Total GPU mem: "<< totalMem <<" bytes\n";
	size_t planBuffer = planMemEstimate(fftsize, numofFFTs, overlap);
	long long allocatedMem = sizeof(char)*bytesToRead + sizeof(cufftComplex)*samplesWithOverlap + nMaxBufferSize +
		sizeof(cufftComplex)*peakSamplesToSave*numofFFTs + sizeof(cufftComplex)*device2quant + sizeof(Npp32f)*inchoerentNumofFFT*fftsize+
		sizeof(Npp32f)*inchoerentNumofFFT*fftsize + sizeof(int)*inchoerentNumofFFT + sizeof(Npp32f)*inchoerentNumofFFT + sizeof(Npp32f)*inchoerentNumofFFT + nStdBufferSize;
	cout << "GPU mem allocated: " << allocatedMem <<" bytes\n";
	cout << "GPU total aprox mem used: " << allocatedMem+ planBuffer <<" bytes\n";
	long long allocatedRAM = sizeof(cufftComplex) * device2quant+ sizeof(cufftComplex) * samplesWithOverlap + sizeof(Npp32f) *inchoerentNumofFFT*3 + sizeof(int) *inchoerentNumofFFT+
		sizeof(char) * bytesToRead+ numofDataLines*13* sizeof(long long) +sizeof(int)*26+ numofDataLines*4*sizeof(int)+sizeof(fileRefName)+sizeof(fileDataNames)+sizeof(bool);
	cout << "RAM total aprox mem used: " << allocatedRAM << " bytes\n\n";





	//READ, EXTEND AND FFT OF REF SIGNAL
	if (interferometic == false) {
		prepareReference(fftsize, overlap, blockSize, hostDataFile2, deviceDataFile2, fileRefName[0]);
		delete[] fileRefName;
		delete[] dataOffsetBegInterferometric;
	}

	//LOOP
	for (i = 0; i < numofDataLines; i++) {
		
		auto begin = std::chrono::high_resolution_clock::now();
		//READ, MASK, AND EXTEND
		prepareData(dataOffsetEnd, dataOffsetBeg, bytesToRead, hostBytesOfData, fileDataNames,
			deviceBytesOfData, blockSize, ddmQuant, samplesOfSignal, samplesWithOverlap, deviceDataFile1
			, numofFFTs, fftsize, hostDataFile1,&elapsed_read, &mask_elapsed
			, &extenddop_elapsed,typeOfDataline,i);

		if (interferometic == true) {
			chrono::nanoseconds elapsed_read_inter, mask_elapsed_inter, extenddop_elapsed_inter;
			prepareData( dataOffsetEndInterferometric, dataOffsetBegInterferometric,
				bytesToRead, hostBytesOfData, fileRefName,deviceBytesOfData, blockSize, ddmQuant, samplesOfSignal,
				samplesWithOverlap, deviceDataFile2, numofFFTs, fftsize, hostDataFile2,&elapsed_read_inter, &mask_elapsed_inter
				, &extenddop_elapsed_inter, typeOfDataline, i);
			elapsed_read += elapsed_read_inter;
			mask_elapsed += mask_elapsed_inter;
			extenddop_elapsed += extenddop_elapsed_inter;
		}
		if (typeOfDataline[i] == 2) {
			i++;
			typ2++;
		}

		//CHECK: RAW DATA 
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*(dataOffsetEnd[i] - dataOffsetBeg[i])*4, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata((dataOffsetEnd[i] - dataOffsetBeg[i])*4, hostDataFile1, "results/rawdata.bin");

		//MULTIPLY BY DOPPLER
		auto dopplerbeg = std::chrono::high_resolution_clock::now();
		
		samplePhaseMantain = (unsigned long long int(i) * unsigned long long int(fftsize)*unsigned long long int(numofFFTs));
		//cout << "phase: " << samplePhaseMantain << "\n";
		numBlocks = (samplesDoppler + blockSize - 1) / blockSize;
		applyDoppler << <numBlocks, blockSize >> > (samplesDoppler, deviceDataFile1, doppler[i], fSampling, samplePhaseMantain, fftsize * numofFFTs, ddmQuant, ddmRes, fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();
		
		auto doppler_elapsed = chrono::high_resolution_clock::now() - dopplerbeg;
		
		//CHECK: doppler (only for printing doppler)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "results/dopplerout.bin");

		//FFT
		auto fftbeg = std::chrono::high_resolution_clock::now();
		planfftFunction(fftsize, numofFFTs, overlap, &plan);
		cudaDeviceSynchronize();
		for (k = 0; k < ddmQuant; k++) {
			cufftSafeCall(cufftExecC2C(plan, &deviceDataFile1[k*(numofFFTs * fftsize)], &deviceDataFile1[k*(numofFFTs * fftsize)], CUFFT_FORWARD));
			cudaDeviceSynchronize();
			if (interferometic == true) {
				cufftSafeCall(cufftExecC2C(plan, &deviceDataFile2[k*(numofFFTs * fftsize)], &deviceDataFile2[k*(numofFFTs * fftsize)], CUFFT_FORWARD));
				cudaDeviceSynchronize();
			}
		}
		
		cufftSafeCall(cufftDestroy(plan));
		auto fft_elapsed = chrono::high_resolution_clock::now() - fftbeg;

		//CHECK: FFT (only for printing fft)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "results/fft.bin");

		//COMPLEX CONJUGATE AND MULTIPLICATION
		auto multbeg = std::chrono::high_resolution_clock::now();
		numBlocks = (samplesWithOverlap + blockSize - 1) / blockSize;
		multip << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, deviceDataFile2, fftsize,interferometic);
		CudaCheckError();
		cudaDeviceSynchronize();
		auto mult_elapsed = chrono::high_resolution_clock::now() - multbeg;
		
		//CHECK: MULTIPLICATION (only for printing multiplication result)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "mult.txt");
		
		//IFFT
		auto ifftbeg = std::chrono::high_resolution_clock::now();
		planifftFunction(fftsize, numofFFTs, 0, &plan);
		cudaDeviceSynchronize();
		for (k = 0; k < ddmQuant; k++) {
			cufftSafeCall(cufftExecC2C(plan, &deviceDataFile1[k*(numofFFTs * fftsize)], &deviceDataFile1[k*(numofFFTs * fftsize)], CUFFT_INVERSE));
			cudaDeviceSynchronize();
		}
		cufftSafeCall(cufftDestroy(plan));
		auto ifft_elapsed = chrono::high_resolution_clock::now() - ifftbeg;

		//SCALE (To take back original signal it has to be devided for the fftsize)
		/*auto scalebeg = std::chrono::high_resolution_clock::now();
		numBlocks = (samplesWithOverlap + blockSize - 1) / blockSize;		
		scale << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();
		auto scale_elapsed = chrono::high_resolution_clock::now() - scalebeg;*/

		//CHECK: IFFT 
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost)); 
		//cudaDeviceSynchronize();	
		//writedata(samplesWithOverlap, hostDataFile1,  "results/IFFT.bin");

		//INCOHERENT SUM
		auto incohbeg = std::chrono::high_resolution_clock::now();
		numBlocks = (inchoerentNumofFFT*fftsize + blockSize - 1) / blockSize;
		inchoerentSum << <numBlocks, blockSize >> > (inchoerentNumofFFT*fftsize, deviceDataFile1, deviceIncoherentSum, quantofAverageIncoherent, fftsize);
		CudaCheckError(); 
		cudaDeviceSynchronize();
		auto incho_elapsed = chrono::high_resolution_clock::now() - incohbeg;
		
		//CHECK: INCOHERENT
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceIncoherentSum, sizeof(Npp32f)*inchoerentNumofFFT*fftsize, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writeIncoh(inchoerentNumofFFT*fftsize, hostDataFile1, "results/incoh.bin");
		
		//MAXIMUM
		auto maxbeg = std::chrono::high_resolution_clock::now();
		maxCompute(inchoerentNumofFFT, deviceIncoherentSum, fftsize, devicearrayMaxs, devicearrayPos, pMaxDeviceBuffer,samplesAvoidMax);
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostarrayPos, devicearrayPos, sizeof(int)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(hostarrayMaxs, devicearrayMaxs, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		auto max_elapsed = chrono::high_resolution_clock::now() - maxbeg;

		//SAVE PEAKS
		auto savepbeg = std::chrono::high_resolution_clock::now();
		if (ddmQuant > 1) {
			//numBlocks = ((numofFFTs / quantofAverageIncoherent) + blockSize - 1) / blockSize;
			selectMaxs << <1, blockSize >> > (numofFFTs, quantofAverageIncoherent, ddmQuant, devicearrayPos, devicearrayMaxs);
			CudaCheckError();
			cudaDeviceSynchronize();

			//CHECK MAX
			//CudaSafeCall(cudaMemcpy(&checkMax, devicearrayPos, sizeof(int)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
			//cout << checkMax[0] << "\n";
		}
		numBlocks = (numofFFTs*peakSamplesToSave*ddmQuant + blockSize - 1) / blockSize;
		savePeak << <numBlocks, blockSize >> > (numofFFTs, deviceDataFile1, deviceDataToSave, peakSamplesToSave, quantofAverageIncoherent, fftsize, devicearrayPos,ddmQuant);
		CudaCheckError();
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataToSave, sizeof(cuComplex)*numofFFTs*peakSamplesToSave*ddmQuant, cudaMemcpyDeviceToHost));
		auto savep_elapsed = chrono::high_resolution_clock::now() - savepbeg;
		
		//STD
		auto stdbeg = std::chrono::high_resolution_clock::now();
		stdCompute(inchoerentNumofFFT, deviceIncoherentSum, fftsize, devicearrayStd, hostarrayPos, pStdDeviceBuffer, peakRangeStd,stdLength, devicearrayMean);
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostarrayStd, devicearrayStd, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(hostarrayMean, devicearrayMean, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		auto std_elapsed = chrono::high_resolution_clock::now() - stdbeg;

		//OUTPUT
		auto writeBeg = chrono::high_resolution_clock::now();
		outputName = resultDirectory+"Maximums" + to_string(i-typ2) + ".bin";
		writeMaxs(inchoerentNumofFFT, hostarrayMaxs, hostarrayPos, hostarrayStd, hostarrayMean,doppler[i], outputName,i, ddmRes,ddmQuant,numofFFTs / quantofAverageIncoherent);
		
		
			outputName = resultDirectory + "PeaksIteration" + to_string(i-typ2) + ".bin";
			cout << outputName << "\n";
		if (writeoutputs == 1) {
			writedata(numofFFTs*peakSamplesToSave*ddmQuant, hostDataFile1, outputName);
		}
		//ELAPSED TIME
		auto elapsed_write = chrono::high_resolution_clock::now() - writeBeg;
		auto elapsed_total = chrono::high_resolution_clock::now() - begin;

		read_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_read).count();
		write_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_write).count();
		elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(elapsed_total).count();
		mask_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(mask_elapsed).count();
		doppler_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(doppler_elapsed).count();
		fft_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(fft_elapsed).count();
		mult_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(mult_elapsed).count();
		ifft_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(ifft_elapsed).count();
		extenddop_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(extenddop_elapsed).count();
		incho_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(incho_elapsed).count();
		max_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(max_elapsed).count();
		savep_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(savep_elapsed).count();
		std_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(std_elapsed).count();
	}
	outputName = resultDirectory + "Times.txt";
	writetime(numofDataLines, outputName, read_elapsed_secs, write_elapsed_secs, elapsed_secs,
		mask_elapsed_secs, extenddop_elapsed_secs, doppler_elapsed_secs,
		fft_elapsed_secs, mult_elapsed_secs, ifft_elapsed_secs, incho_elapsed_secs
		, max_elapsed_secs, savep_elapsed_secs, std_elapsed_secs);
	
	//FREE MEMORY
	//cufftSafeCall(cufftDestroy(plan));
	//cufftSafeCall(cufftDestroy(inverseplan));
	cudaFree(deviceDataFile1);
	cudaFree(deviceDataFile2);
	cudaFree(deviceIncoherentSum);
	cudaFree(devicearrayPos);
	cudaFree(deviceBytesOfData);
	cudaFree(devicearrayMaxs);
	cudaFree(devicearrayMean);
	cudaFree(deviceDataToSave);
	cudaFree(pStdDeviceBuffer);
	cudaFree(pMaxDeviceBuffer);
	cudaFree(devicearrayStd);
	cudaDeviceReset();
	delete[] fileDataNames;
	if (interferometic == true) {
		delete[] fileRefName;
		delete[] dataOffsetBegInterferometric;
	}
	delete[] hostBytesOfData;
	delete[] hostarrayPos;
	delete[] hostarrayMaxs;
	delete[] hostarrayMean;
	delete[] hostarrayStd;
	delete[] hostDataFile2;
	delete[] hostDataFile1;
	delete[] dataOffsetBeg;
	delete[] dataOffsetEnd;
	delete[] typeOfDataline;
	delete[] dataOffsetEndInterferometric;
	delete[] doppler;
	delete[] read_elapsed_secs;
	delete[] write_elapsed_secs;
	delete[] elapsed_secs;
	delete[] mask_elapsed_secs;
	delete[] doppler_elapsed_secs ;
	delete[] mult_elapsed_secs ;
	delete[] fft_elapsed_secs;
	delete[] ifft_elapsed_secs;
	delete[] extenddop_elapsed_secs;
	delete[] incho_elapsed_secs ;
	delete[] max_elapsed_secs ;
	delete[] savep_elapsed_secs ;
	delete[] std_elapsed_secs;
	return 0;
}