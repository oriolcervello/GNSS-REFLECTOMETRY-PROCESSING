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
	int fftsize, fSampling, numofFFTs, overlap, quantofAverageIncoherent, blockSize, peakRangeStd, peakSamplesToSave, ddmRes, ddmQuant;
	int const numofDataLines = atoi(argv[2]);//substitut d'iterations
	string *fileDataNames, fileRefName;
	int *dataOffsetBeg, *dataOffsetEnd;
	int *doppler;

	fileDataNames = new string[numofDataLines];
	dataOffsetBeg = new int[numofDataLines];
	dataOffsetEnd = new int[numofDataLines];
	doppler = new int[numofDataLines];

	readConfig(argv[1], numofDataLines, &fftsize, &numofFFTs, &overlap, &fSampling, &blockSize, &peakRangeStd, &peakSamplesToSave, &quantofAverageIncoherent, dataOffsetBeg, dataOffsetEnd, doppler, fileDataNames, &fileRefName, &ddmRes, &ddmQuant);
	checkInputConfig(argc, argv, numofDataLines, fftsize, numofFFTs, overlap, fSampling, blockSize, peakRangeStd, peakSamplesToSave, quantofAverageIncoherent, dataOffsetBeg, dataOffsetEnd, doppler, fileDataNames, fileRefName, ddmRes, ddmQuant);

	//OTHER DECLARATIONS
	int  originalSamplesOfSignal = (numofFFTs * (fftsize - overlap)) + overlap;//samples of complex data
	int samplesOfSignal = originalSamplesOfSignal *ddmQuant;//samples of complex data
	int bytesToRead = originalSamplesOfSignal /4;
	if (originalSamplesOfSignal % 4 != 0) { cout << "Warning bytesToRead rounded toward negative infinity: samplesOfSignal%4!=0 \n "; }
	int samplesWithOverlap= (numofFFTs * fftsize)*ddmQuant;//total samples needed
	if(samplesOfSignal > samplesWithOverlap){ samplesWithOverlap = samplesOfSignal;}
	int inchoerentNumofFFT = (numofFFTs/ quantofAverageIncoherent)*ddmQuant;
	if (numofFFTs % quantofAverageIncoherent != 0) {
		cout << "Error: numofFFTs / quantofAverageIncoherent != 0\n ";
		exit(-1);
	}

	string outputName;
	int numBlocks, nBufferSize,i,k;
	unsigned long long samplePhaseMantain;

	char *hostBytesOfData, *deviceBytesOfData;
	int *devicearrayPos,*hostarrayPos;
	cufftComplex *deviceDataFile1, *deviceDataFile2, *hostDataFile1, *hostDataFile2, *deviceDataToSave;
	Npp32f *deviceIncoherentSum, *devicearrayMaxs, *devicearrayStd, *hostarrayMaxs, *hostarrayStd;
	Npp8u * pDeviceBuffer;

	cufftHandle plan, planref;
	
	long long *read_elapsed_secs,*write_elapsed_secs, *elapsed_secs, *mask_elapsed_secs, *doppler_elapsed_secs, 
		 *fft_elapsed_secs, *mult_elapsed_secs,*ifft_elapsed_secs, *scale_elapsed_secs, *incho_elapsed_secs
		, *max_elapsed_secs, *savep_elapsed_secs, *std_elapsed_secs;
	
	//ALLOCATE
	read_elapsed_secs = new long long[numofDataLines];
	mask_elapsed_secs = new long long[numofDataLines];
	doppler_elapsed_secs = new long long[numofDataLines];
	fft_elapsed_secs = new long long[numofDataLines];
	mult_elapsed_secs = new long long[numofDataLines];
	ifft_elapsed_secs = new long long[numofDataLines];
	scale_elapsed_secs = new long long[numofDataLines];
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
	hostDataFile1 = (cufftComplex *)malloc(sizeof(cufftComplex) * samplesWithOverlap);
	hostDataFile2 = (cufftComplex *)malloc(sizeof(cufftComplex) * fftsize);
	CudaSafeCall(cudaMalloc(&deviceBytesOfData, sizeof(char)*bytesToRead));
	CudaSafeCall(cudaMalloc(&deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap));
	CudaSafeCall(cudaMalloc(&deviceDataToSave, sizeof(cufftComplex)*peakSamplesToSave*numofFFTs*ddmQuant));
	CudaSafeCall(cudaMalloc(&deviceDataFile2, sizeof(cufftComplex)*fftsize));
	CudaSafeCall(cudaMalloc(&deviceIncoherentSum, sizeof(Npp32f)*inchoerentNumofFFT*fftsize));
	CudaSafeCall(cudaMalloc(&devicearrayPos, sizeof(int)*inchoerentNumofFFT));
	CudaSafeCall(cudaMalloc(&devicearrayMaxs, sizeof(Npp32f)*inchoerentNumofFFT));
	CudaSafeCall(cudaMalloc(&devicearrayStd, sizeof(Npp32f)*inchoerentNumofFFT));
	nppsSumGetBufferSize_32f(fftsize, &nBufferSize);
	CudaSafeCall(cudaMalloc((void **)(&pDeviceBuffer), nBufferSize));
	cudaDeviceSynchronize();
	
	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);
	cout<< "\n-MEMORY: \n";
	//fprintf(stderr, "   Free = %zu, Total = %zu\n", freeMem, totalMem);
	cout<< "Total GPU mem: "<< totalMem <<" bytes\n";
	size_t planBuffer = planMemEstimate(fftsize, numofFFTs, overlap);
	long long allocatedMem = sizeof(char)*bytesToRead + sizeof(cufftComplex)*samplesWithOverlap +
		sizeof(cufftComplex)*peakSamplesToSave*numofFFTs + sizeof(cufftComplex)*fftsize + sizeof(Npp32f)*inchoerentNumofFFT*fftsize
		+ sizeof(int)*inchoerentNumofFFT + sizeof(Npp32f)*inchoerentNumofFFT + sizeof(Npp32f)*inchoerentNumofFFT + nBufferSize;
	cout << "GPU mem allocated: " << allocatedMem <<" bytes\n";
	cout << "GPU total aprox mem used: " << allocatedMem+ planBuffer <<" bytes\n\n";
	
	



	//READ, EXTEND AND FFT OF REF SIGNAL
	readdata(fftsize - overlap, 0, hostDataFile2, fileRefName);

	//writedata(fftsize - overlap, hostDataFile2, "rawref2.bin");

	CudaSafeCall(cudaMemcpy(deviceDataFile2, hostDataFile2, sizeof(cufftComplex)*(fftsize - overlap), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	
	numBlocks = (fftsize + blockSize - 1) / blockSize;
	extendRefSignal << <numBlocks, blockSize >> > (fftsize, deviceDataFile2, fftsize - overlap);
	CudaCheckError();


	planfftFunction(fftsize, 1, 0, &planref);
	cudaDeviceSynchronize();
	cufftSafeCall(cufftExecC2C(planref, deviceDataFile2, deviceDataFile2, CUFFT_FORWARD));
	cudaDeviceSynchronize();
	cufftSafeCall(cufftDestroy(planref));



	//LOOP
	for (i = 0; i < numofDataLines; i++) {
		
		auto begin = std::chrono::high_resolution_clock::now();
		//READ DATA
		/*readdata(dataOffsetEnd[i]-dataOffsetBeg[i], dataOffsetBeg[i], hostDataFile1, fileDataNames[i]);*/
		readRealData(dataOffsetEnd[i] - dataOffsetBeg[i], dataOffsetBeg[i],bytesToRead, hostBytesOfData, fileDataNames[i]);
		
		/*CudaSafeCall(cudaMemcpy(deviceDataFile1, hostDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyHostToDevice));*/
		CudaSafeCall(cudaMemcpy(deviceBytesOfData, hostBytesOfData, sizeof(char)*bytesToRead, cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();
		auto elapsed_read = chrono::high_resolution_clock::now() - begin;

		//MASK AND SHIFT
		auto maskbeg = std::chrono::high_resolution_clock::now();
		numBlocks = (bytesToRead + blockSize - 1) / blockSize;
		maskAndShift << <numBlocks, blockSize >> > (deviceBytesOfData, deviceDataFile1, bytesToRead);
		CudaCheckError();
		cudaDeviceSynchronize();
		auto mask_elapsed = chrono::high_resolution_clock::now() - maskbeg;
		
		//EXTEND FOR DOPPLER
		auto scalebeg = std::chrono::high_resolution_clock::now();
		if (ddmQuant > 1) {
			numBlocks = (samplesOfSignal + blockSize - 1) / blockSize;
			extendRefSignal << <numBlocks, blockSize >> > (samplesOfSignal, deviceDataFile1, originalSamplesOfSignal);
			CudaCheckError();
			cudaDeviceSynchronize();
		}
		auto scale_elapsed = chrono::high_resolution_clock::now() - scalebeg;



		//CHECK: RAW DATA 
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*(dataOffsetEnd[i] - dataOffsetBeg[i])*4, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata((dataOffsetEnd[i] - dataOffsetBeg[i])*4, hostDataFile1, "rawdata.bin");

		//MULTIPLY BY DOPPLER
		auto dopplerbeg = std::chrono::high_resolution_clock::now();
		samplePhaseMantain = (i * fftsize*numofFFTs);// %fSampling;----
		numBlocks = (samplesOfSignal + blockSize - 1) / blockSize;
		applyDoppler << <numBlocks, blockSize >> > (samplesOfSignal, deviceDataFile1, doppler[i], fSampling, samplePhaseMantain,originalSamplesOfSignal,ddmQuant,ddmRes,fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();
		auto doppler_elapsed = chrono::high_resolution_clock::now() - dopplerbeg;
		//CHECK: doppler (only for printing doppler)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesOfSignal/2, hostDataFile1, "dopplerout.txt");

		//FFT
		auto fftbeg = std::chrono::high_resolution_clock::now();
		planfftFunction(fftsize, numofFFTs, overlap, &plan);
		cudaDeviceSynchronize();
		for (k = 0; k < ddmQuant; k++) {
			cufftSafeCall(cufftExecC2C(plan, &deviceDataFile1[k*originalSamplesOfSignal], &deviceDataFile1[k*originalSamplesOfSignal], CUFFT_FORWARD));
		}
		cudaDeviceSynchronize();
		cufftSafeCall(cufftDestroy(plan));
		auto fft_elapsed = chrono::high_resolution_clock::now() - fftbeg;

		//CHECK: FFT (only for printing fft)
		//CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataFile1, sizeof(cufftComplex)*samplesWithOverlap, cudaMemcpyDeviceToHost));
		//cudaDeviceSynchronize();
		//writedata(samplesWithOverlap, hostDataFile1, "fft.txt");

		//COMPLEX CONJUGATE AND MULTIPLICATION
		auto multbeg = std::chrono::high_resolution_clock::now();
		numBlocks = (samplesWithOverlap + blockSize - 1) / blockSize;
		multip << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, deviceDataFile2, fftsize);
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
			cufftSafeCall(cufftExecC2C(plan, &deviceDataFile1[k*originalSamplesOfSignal], &deviceDataFile1[k*originalSamplesOfSignal], CUFFT_INVERSE));
		}
		cudaDeviceSynchronize();
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
		//writedata(samplesWithOverlap, hostDataFile1,  "IFFT-result.bin");

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
		//writeIncoh(inchoerentNumofFFT*fftsize, hostDataFile1, "incoh.bin");
		
		//MAXIMUM
		auto maxbeg = std::chrono::high_resolution_clock::now();
		maxCompute(inchoerentNumofFFT, deviceIncoherentSum, fftsize, devicearrayMaxs, devicearrayPos, pDeviceBuffer);
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostarrayPos, devicearrayPos, sizeof(int)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		CudaSafeCall(cudaMemcpy(hostarrayMaxs, devicearrayMaxs, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		auto max_elapsed = chrono::high_resolution_clock::now() - maxbeg;
		cudaDeviceSynchronize();
		
	


		//SAVE PEAKS
		auto savepbeg = std::chrono::high_resolution_clock::now();
		if (ddmQuant > 1) {
			numBlocks = ((numofFFTs / quantofAverageIncoherent) + blockSize - 1) / blockSize;
			selectMaxs << <numBlocks, blockSize >> > (numofFFTs, quantofAverageIncoherent, ddmQuant, devicearrayPos, devicearrayMaxs);
			CudaCheckError();
			cudaDeviceSynchronize();
		}
		




		numBlocks = (numofFFTs*peakSamplesToSave*ddmQuant + blockSize - 1) / blockSize;
		savePeak << <numBlocks, blockSize >> > (numofFFTs, deviceDataFile1, deviceDataToSave, peakSamplesToSave, quantofAverageIncoherent, fftsize, devicearrayPos,ddmQuant);
		CudaCheckError();
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostDataFile1, deviceDataToSave, sizeof(cuComplex)*numofFFTs*peakSamplesToSave*ddmQuant, cudaMemcpyDeviceToHost));
		auto savep_elapsed = chrono::high_resolution_clock::now() - savepbeg;
		
		//STD
		auto stdbeg = std::chrono::high_resolution_clock::now();
		stdCompute(inchoerentNumofFFT, deviceIncoherentSum, fftsize, devicearrayStd, hostarrayPos, pDeviceBuffer, peakRangeStd);
		cudaDeviceSynchronize();
		CudaSafeCall(cudaMemcpy(hostarrayStd, devicearrayStd, sizeof(Npp32f)*inchoerentNumofFFT, cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		auto std_elapsed = chrono::high_resolution_clock::now() - stdbeg;

		//OUTPUT
		auto writeBeg = chrono::high_resolution_clock::now();
		writeMaxs(inchoerentNumofFFT, hostarrayMaxs, hostarrayPos, hostarrayStd, doppler[i],"results/Maximums.txt");
		outputName = "results/PeaksIteration"+ to_string(i);
		outputName = outputName + ".bin";
		cout << outputName << "\n";
		writedata(numofFFTs*peakSamplesToSave*ddmQuant, hostDataFile1, outputName);
	
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
		scale_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(scale_elapsed).count();
		incho_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(incho_elapsed).count();
		max_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(max_elapsed).count();
		savep_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(savep_elapsed).count();
		std_elapsed_secs[i] = chrono::duration_cast<chrono::microseconds>(std_elapsed).count();
	}

	writetime(numofDataLines, "results/Times.txt", read_elapsed_secs, write_elapsed_secs, elapsed_secs,
		mask_elapsed_secs, doppler_elapsed_secs,
		fft_elapsed_secs, mult_elapsed_secs, ifft_elapsed_secs, scale_elapsed_secs, incho_elapsed_secs
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
	cudaFree(deviceDataToSave);
	cudaFree(pDeviceBuffer);
	cudaFree(devicearrayStd);
	cudaDeviceReset();
	delete[] fileDataNames;
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
	delete[] mask_elapsed_secs;
	delete[] doppler_elapsed_secs ;
	delete[] mult_elapsed_secs ;
	delete[] fft_elapsed_secs;
	delete[] ifft_elapsed_secs;
	delete[] scale_elapsed_secs;
	delete[] incho_elapsed_secs ;
	delete[] max_elapsed_secs ;
	delete[] savep_elapsed_secs ;
	delete[] std_elapsed_secs;
	return 0;
}