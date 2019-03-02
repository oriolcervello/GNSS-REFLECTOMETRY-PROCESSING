#include"Hostfunc.cuh"
#include "TextParser.cuh"


//INPUT CONFIG PARSER FUNCTIONS

void readConfig(const char *configFileName, int numofDataLines, int *fftsize, int *numofFFts, int *overlap, int *fSampling, int *blockSize, int *peakRangeStd, int *peakSamplesToSave,
	int* quantOfAverIncoh, int *dataOffsetBeg, int *dataOffsetEnd, int *doppler, string *fileNames,string *fileRefNames, int *ddmRes, int *ddmQuant,bool *interfer
,int *dataOffsetBegInterferometric,int *samplesAvoidMaxs,string *resultDirectory) {

	TextParser t(configFileName);

	TextParserSafeCall(t.seek("*FFTSIZE"));
	*fftsize = t.getint();
	TextParserSafeCall(t.seek("*NUMOFFFTS"));
	*numofFFts = t.getint();
	TextParserSafeCall(t.seek("*QUANTINCOHAVER"));
	*quantOfAverIncoh = t.getint();
	TextParserSafeCall(t.seek("*OVERLAP"));
	*overlap = t.getint();
	TextParserSafeCall(t.seek("*FSAMPLING"));
	*fSampling = t.getint();
	TextParserSafeCall(t.seek("*BLOCKSIZE"));
	*blockSize = t.getint();
	TextParserSafeCall(t.seek("*INTERFEROMETIC"));
	*interfer = t.getint();
	TextParserSafeCall(t.seek("*SAMPLESAVOIDMAX"));
	*samplesAvoidMaxs = t.getint();
	TextParserSafeCall(t.seek("*PEAKRANGESTD"));
	*peakRangeStd = t.getint();
	TextParserSafeCall(t.seek("*PEAKSAMPLESTOSAVE"));
	*peakSamplesToSave = t.getint();
	if (*interfer == false) {
		TextParserSafeCall(t.seek("*REFFILENAME"));
		fileRefNames[0] = t.getword(); 
	}
	TextParserSafeCall(t.seek("*RESULTSDIRECTORY"));
	*resultDirectory = t.getword();
	TextParserSafeCall(t.seek("*DDMFREQRES"));
	*ddmRes = t.getint();
	TextParserSafeCall(t.seek("*DDMNUMQUANT"));
	*ddmQuant = t.getint();

	if (*ddmQuant % 2 != 1) {
		cout << "ERROR: DDM QUANT has to be odd: 1(original)+2n(symethric)\n";
		exit(1);

	}


	TextParserSafeCall(t.seek("*QUANTDATALINES"));
	if (t.getint() != numofDataLines) {
		cout << "diferent num of Data lines in the file than declared on arguments \n  To execute enter arguments: NameconfigFile.ASE NumofDataLines\n";
		exit(1);
	}

	for (int i = 0; i < numofDataLines; i++) {
		TextParserSafeCall(t.seek("*DATALINE"));

		fileNames[i] = t.getword();
		dataOffsetBeg[i] = t.getint();
		dataOffsetEnd[i] = t.getint();
		doppler[i] = t.getint();
		if (*interfer == true) {
			fileRefNames[i] = t.getword();
			dataOffsetBegInterferometric[i]= t.getint();
		}
	}
}

void checkInputConfig(int argc, const char **argv, int numofDataLines, int fftsize, int numofFFts, int overlap, int fSampling,  int blockSize, int peakRangeStd, int peakSamplesToSave,
	int quantOfAverIncoh,  int *dataOffsetBeg, int *dataOffsetEnd, int *doppler, string *fileNames, string *fileRefNames, int ddmRes, int ddmQuant, bool interfer
, int *dataOffsetBegInterferometric, int samplesAvoidMaxs,string resultDirectory) {

	if (argc != 3) {
		cout << "Error: Wrong number of arguments\n"; 
		exit(0);
	}

	cout << "\n" << "-ARGS: " << argc << "\n";
	cout << "First: " << argv[0] << "\n";
	cout << "Second: " << argv[1] << "\n";
	cout << "Third: " << argv[2] << "\n\n";

	cout << "-INPUTS:\n";
	cout << "FFT Size: " << fftsize << "\n";
	cout << "Num. of FFT: " << numofFFts << "\n";
	cout << "Overlap: " << overlap << "\n";
	cout << "FSampling: " << fSampling << "\n";
	cout << "Quant of averg Inch.: " << quantOfAverIncoh << "\n";
	cout << "Blok Size: " << blockSize << "\n";
	cout << "Interferometric: " << interfer << "\n";
	cout << "Samples avoid MAxs: " << samplesAvoidMaxs << "\n";
	cout << "Peak samples for the std: " << peakRangeStd << "\n";
	cout << "Peak samples to save: " << peakSamplesToSave << "\n";
	if (interfer == false) {
		cout << "Ref File Name: " << fileRefNames[0] << "\n";
	}
	cout << "result directory: " << resultDirectory << "\n";
	cout << "DDM Res: " << ddmRes << "\n";
	cout << "DDM Quant: " << ddmQuant << "\n";


	cout << "Num of data lines: " << numofDataLines << "\n";
	cout << "Data lines: \n";
	for (int i = 0; i < numofDataLines; i++) {
		cout << fileNames[i] << "  ";
		cout << dataOffsetBeg[i] << "  ";
		cout << dataOffsetEnd[i] << "  ";
		cout << doppler[i] << "  ";
		if (interfer == true) {
			cout << fileRefNames[i] << "  ";
			cout << dataOffsetBegInterferometric[i] << "\n";
		}
		else {
			cout << "\n";
		}
	}

}

//PREPARE DATA FUNCTIONS

void prepareReference(int fftsize, int overlap,int blockSize ,cufftComplex *hostDataFile2, cufftComplex *deviceDataFile2,string fileRefName) {
	
	readdata(fftsize - overlap, 0, hostDataFile2, fileRefName);
	cufftHandle planref;
	CudaSafeCall(cudaMemcpy(deviceDataFile2, hostDataFile2, sizeof(cufftComplex)*(fftsize - overlap), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	if (overlap > 0) {
		int numBlocks = (fftsize + blockSize - 1) / blockSize;
		extendRefSignal << <numBlocks, blockSize >> > (fftsize, deviceDataFile2, fftsize - overlap);
		CudaCheckError();
	}
	planfftFunction(fftsize, 1, 0, &planref);
	cudaDeviceSynchronize();
	cufftSafeCall(cufftExecC2C(planref, deviceDataFile2, deviceDataFile2, CUFFT_FORWARD));
	cudaDeviceSynchronize();
	cufftSafeCall(cufftDestroy(planref));

}


void prepareData( int dataOffsetEnd,int dataOffsetBeg, int bytesToRead, char *hostBytesOfData, string fileDataNames,
	char *deviceBytesOfData, int blockSize, int ddmQuant, int samplesOfSignal, int samplesWithOverlap, cufftComplex *deviceDataFile1
     ,int numofFFTs, int fftsize, cufftComplex *hostDataFile1, chrono::nanoseconds *elapsed_read, chrono::nanoseconds *mask_elapsed
	,chrono::nanoseconds *extenddop_elapsed) {
	
	auto begin = std::chrono::high_resolution_clock::now();
	//READ DATA
	//readdata(dataOffsetEnd[i]-dataOffsetBeg[i], dataOffsetBeg[i], hostDataFile1, fileDataNames[i]);
	readRealData(dataOffsetEnd - dataOffsetBeg, dataOffsetBeg, bytesToRead, hostBytesOfData, fileDataNames);

	//CudaSafeCall(cudaMemcpy(deviceDataFile1, hostDataFile1, sizeof(cufftComplex)*samplesOfSignal, cudaMemcpyHostToDevice));
	CudaSafeCall(cudaMemcpy(deviceBytesOfData, hostBytesOfData, sizeof(char)*bytesToRead, cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	*elapsed_read = chrono::high_resolution_clock::now() - (begin);

	//MASK AND SHIFT
	auto maskbeg = std::chrono::high_resolution_clock::now();
	int numBlocks = (bytesToRead + blockSize - 1) / blockSize;
	maskAndShift << <numBlocks, blockSize >> > (deviceBytesOfData, deviceDataFile1, bytesToRead);
	CudaCheckError();
	cudaDeviceSynchronize();
	*mask_elapsed = chrono::high_resolution_clock::now() - maskbeg;

	//EXTEND FOR DOPPLER
	auto extenddopbeg = std::chrono::high_resolution_clock::now();
	if (ddmQuant > 1) {
		numBlocks = (samplesOfSignal + blockSize - 1) / blockSize;
		extendRefSignal << <numBlocks, blockSize >> > (samplesWithOverlap, deviceDataFile1, numofFFTs * fftsize);
		CudaCheckError();
		cudaDeviceSynchronize();
	}
	*extenddop_elapsed = chrono::high_resolution_clock::now() - extenddopbeg;


}

//FFT PLANS FUNCTIONS

void planfftFunction(int fftsize, int numofFFTs, int overlap, cufftHandle *plan) {

	int rank = 1;                           // --- 1D FFTs
	int n[] = { fftsize };                 // --- Size of the Fourier transform
	int istride = 1, ostride = 1;           // --- Distance between two successive input/output elements
	int idist = fftsize - overlap, odist = fftsize;// (DATASIZE / 2 + 1); // --- Distance between batches
	int inembed[] = { 0 };                  // --- Input size with pitch (ignored for 1D transforms)
	int onembed[] = { 0 };                  // --- Output size with pitch (ignored for 1D transforms)
	int batch = numofFFTs;// numofFFTs;                      // --- Number of batched executions
	cufftSafeCall(cufftPlanMany(plan, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch));

}

void planifftFunction(int fftsize, int numofFFTs, int overlap, cufftHandle *plan) {
	
	int rank = 1;                           // --- 1D FFTs
	int n[] = { fftsize };                 // --- Size of the Fourier transform
	int istride = 1, ostride = 1;           // --- Distance between two successive input/output elements
	int idist = fftsize, odist = fftsize - overlap;// (DATASIZE / 2 + 1); // --- Distance between batches
	int inembed[] = { 0 };                  // --- Input size with pitch (ignored for 1D transforms)
	int onembed[] = { 0 };                  // --- Output size with pitch (ignored for 1D transforms)
	int batch = numofFFTs;// numofFFTs;                      // --- Number of batched executions
	cufftSafeCall(cufftPlanMany(plan, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch));

}

size_t planMemEstimate(int fftsize, int numofFFTs, int overlap) {

	int rank = 1;                           // --- 1D FFTs
	int n[] = { fftsize };                 // --- Size of the Fourier transform
	int istride = 1, ostride = 1;           // --- Distance between two successive input/output elements
	int idist = fftsize, odist = fftsize - overlap;// (DATASIZE / 2 + 1); // --- Distance between batches
	int inembed[] = { 0 };                  // --- Input size with pitch (ignored for 1D transforms)
	int onembed[] = { 0 };                  // --- Output size with pitch (ignored for 1D transforms)
	int batch = numofFFTs;// numofFFTs;                      // --- Number of batched executions
	size_t workSize;
	cufftSafeCall(cufftEstimateMany( rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch,&workSize));

	cout << "cufft plan aprox buffer: " << workSize<< " bytes\n";
	return workSize;
}

//STATISTICS FUNCTIONS

void maxCompute(int numofIncoherentSums, Npp32f *deviceDataIncoherentSum, int fftsize, Npp32f *deviceArrayMaxs,
	 int *deviceArrayPos, Npp8u * pDeviceBuffer, int samplesAvoidMax) {

	for (int i = 0; i < numofIncoherentSums; i++) {

		nppsMaxIndx_32f(&deviceDataIncoherentSum[i*fftsize+(samplesAvoidMax)], fftsize+ samplesAvoidMax, &deviceArrayMaxs[i], &deviceArrayPos[i], pDeviceBuffer);
	}
}

void stdCompute(int numofIncoherentSums, Npp32f *dataIncoherentSum, int fftsize,
	Npp32f *deviceArraystd, int *arrayPos, Npp8u * pStdDeviceBuffer, int peakRange,int stdLength, Npp32f *devicearrayMean) {

	int leftPeakIndex, rightPeakIndex;
	//stdLength = (fftsize / 2) - ((peakRange) / 2)-1;
	for (int i = 0; i < numofIncoherentSums; i++) {
		
		leftPeakIndex = arrayPos[i] - peakRange/2;
		rightPeakIndex = arrayPos[i] + peakRange/2;
		
		if (rightPeakIndex >= fftsize) {//case 2
			rightPeakIndex = rightPeakIndex % fftsize;
			//stdLength = leftPeakIndex - rightPeakIndex;
			nppsMeanStdDev_32f(&dataIncoherentSum[i*fftsize+ rightPeakIndex], stdLength,&devicearrayMean[i],&deviceArraystd[i], pStdDeviceBuffer);
		}
		else if (leftPeakIndex < 0) {//case 3
			leftPeakIndex = fftsize + leftPeakIndex;
			//stdLength = leftPeakIndex-rightPeakIndex ;
			nppsMeanStdDev_32f(&dataIncoherentSum[i*fftsize + rightPeakIndex], stdLength, &devicearrayMean[i], &deviceArraystd[i], pStdDeviceBuffer);
		}
		else {//case 1
			if (arrayPos[i] < fftsize / 2) {
				//stdLength = fftsize- rightPeakIndex;
				nppsMeanStdDev_32f(&dataIncoherentSum[i*fftsize + rightPeakIndex], stdLength, &devicearrayMean[i], &deviceArraystd[i], pStdDeviceBuffer);
			}
			else {
				//stdLength = leftPeakIndex;
				nppsMeanStdDev_32f(&dataIncoherentSum[i*fftsize], stdLength, &devicearrayMean[i], &deviceArraystd[i], pStdDeviceBuffer);
			}			
		}		
	}
}

