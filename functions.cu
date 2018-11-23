#include"functions.cuh"
#include "extra/TextParser.cuh"


//INPUT CONFIG PARSER FUNCTIONS

void readConfig(const char *configFileName, int numofDataLines, int *fftsize, int *numofFFts, int *overlap, int *fSampling, int *blockSize, int *peakRangeStd, int *peakSamplesToSave,
	int* quantOfAverIncoh, int *dataOffsetBeg, int *dataOffsetEnd, int *doppler, string *fileNames,string *fileRefNames, int *ddmRes, int *ddmQuant) {

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
	TextParserSafeCall(t.seek("*PEAKRANGESTD"));
	*peakRangeStd = t.getint();
	TextParserSafeCall(t.seek("*PEAKSAMPLESTOSAVE"));
	*peakSamplesToSave = t.getint();
	TextParserSafeCall(t.seek("*REFFILENAME"));
	*fileRefNames = t.getword();
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
		
	}
}

void checkInputConfig(int argc, const char **argv, int numofDataLines, int fftsize, int numofFFts, int overlap, int fSampling,  int blockSize, int peakRangeStd, int peakSamplesToSave,
	int quantOfAverIncoh,  int *dataOffsetBeg, int *dataOffsetEnd, int *doppler, string *fileNames, string fileRefNames, int ddmRes, int ddmQuant) {

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
	cout << "Peak samples for the std: " << peakRangeStd << "\n";
	cout << "Peak samples to save: " << peakSamplesToSave << "\n";
	cout << "Ref File Name: " << fileRefNames << "\n";
	cout << "DDM Res: " << ddmRes << "\n";
	cout << "DDM Quant: " << ddmQuant << "\n";


	cout << "Num of data lines: " << numofDataLines << "\n";
	cout << "Data lines: \n";
	for (int i = 0; i < numofDataLines; i++) {
		cout << fileNames[i] << "  ";
		cout << dataOffsetBeg[i] << "  ";
		cout << dataOffsetEnd[i] << "  ";
		cout << doppler[i] << "\n";

	}

}

//READ FUNCTIONS

void readdata(int length,int offsetFromBeg, cufftComplex *data, string name)
{
	ifstream myfile;
	myfile.open(name, ios::binary);
	//float num1,num2;
	
	if (myfile.is_open())
	{
		myfile.seekg(offsetFromBeg* sizeof(cufftComplex));
		
		myfile.read((char*)data, length*sizeof(cufftComplex));
		/*int k = 0;
		while (k < length)
		{
			
			myfile.read((char*)&num1, sizeof(num1));
			myfile.read((char*)&num2, sizeof(num2));
			data[k].x = num1;
			data[k].y =  num2;
			k++;
		}*/
		myfile.close();
	}
	else cout << "Unable to open file of floats for reading " << name << "\n";
}

void readRealData(int length, int offsetFromBeg, int bytesToRead,char *data, string name)
{
	if (length > bytesToRead) {
		cout << "Error: iput length bigger than bytesToRead\n";
		exit(0);
	}

	ifstream myfile;
	myfile.open(name, ios::binary);
	if (myfile.is_open())
	{
		myfile.seekg(offsetFromBeg*sizeof(char));
		myfile.read(data, length);
				
		myfile.close();
		if(length< bytesToRead){
			cout << "Warning: length smaller than bytesToRead, " << bytesToRead - length<<" Bytes filled with 0 \n Last/s incoherents will be incomplete \n";
			memset(&data[length], 0, bytesToRead - length);
		}
	}
	else cout << "Unable to open file of Real Data for reading " << name << "\n";
}

//WRITE FUNCTIONS

void writeIncoh(int N, cuComplex *data1, string name) {

	ofstream myfile;
	myfile.open(name, ios::binary);
	if (myfile.is_open())
	{
		for (int ii = 0; ii < N/2; ii++)
		{
			

			myfile.write((char*)&data1[ii].x, sizeof(float));
			myfile.write((char*)&data1[ii].y, sizeof(float));
			

		}
		myfile.close();
	}

	else cout << "Unable to open file of incoh for writting " << name << "\n";
}

void openMaxsFile(string name) {

	ofstream myfile;
	myfile.open(name);
	if (myfile.is_open())
	{
		
		myfile.close();
	}

	else cout << "Unable to open file of Maxs " << name << "\n";
}

void writeMaxs(int N, Npp32f *dataMaxValue, int *dataMaxPos, Npp32f *hostarrayStd,int doppler,string name) {

	ofstream myfile;
	myfile.open(name, ios::app);
	if (myfile.is_open())
	{
		for (int ii = 0; ii < N; ii++)
		{
			myfile <<"Pos: "<< dataMaxPos[ii]<<" Value: " << dataMaxValue[ii] << " STD: " << hostarrayStd[ii]<< " Doppler: " << doppler<< "\n";

		}
		myfile << "\n --------------------------------------------------------------------- \n\n";
		myfile.close();
	}

	else cout << "Unable to open file of Maxs " << name << "\n";
}

void writedata(int N, cufftComplex *data1, string name) {

	ofstream myfile;
	myfile.open(name, ios::binary);
	if (myfile.is_open())
	{


		myfile.write((char*)data1, N*sizeof(cufftComplex));
		
		/*for (int ii = 0; ii < N; ii++)
		{
			
			myfile.write((char*)&data1[ii].x, sizeof(float));
			myfile.write((char*)&data1[ii].y, sizeof(float));
		}*/
		myfile.close();
	}

	else cout << "Unable to open file of data to write " << name << "\n";
}

void writetime(int N, string name, long long *readtime, long long *writetime, long long *looptime
	, long long *mask_elapsed_secs,long long *extenddop_elapsed_secs, long long *doppler_elapsed_secs,
	long long *fft_elapsed_secs, long long *mult_elapsed_secs, long long *ifft_elapsed_secs,
	long long *incho_elapsed_secs, long long *max_elapsed_secs, long long *savep_elapsed_secs, long long *std_elapsed_secs) {

	ofstream myfile;
	myfile.open(name);
	if (myfile.is_open())
	{
		myfile << "Atempt\t\tReadT.\t\tMask\t\tDoppler\t\tFFT\t\tMul\t\tIFFT\t\tScale\t\tIncoh\t\tMax\t\tSaveP.\t\tSTD\t\tWriteT.\t\tLoopT." << "\n";
		for (int ii = 0; ii < N; ii++)
		{
			myfile << ii << "\t\t" << readtime[ii] << "\t\t" << mask_elapsed_secs[ii] << "\t\t"
				<<extenddop_elapsed_secs[ii] <<"\t\t" << doppler_elapsed_secs[ii] << "\t\t"
				 << fft_elapsed_secs[ii] << "\t\t" << mult_elapsed_secs[ii] << "\t\t"
				 << ifft_elapsed_secs[ii] << "\t\t" 
				 << incho_elapsed_secs[ii] << "\t\t" << max_elapsed_secs[ii] << "\t\t"
				 << savep_elapsed_secs[ii] << "\t\t" << std_elapsed_secs[ii] << "\t\t"
				<< writetime[ii] << "\t\t" << looptime[ii] << "\n";

		}
		myfile.close();
	}

	else cout << "Unable to open file of times "<<name<<"\n";
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
	 int *deviceArrayPos, Npp8u * pDeviceBuffer) {

	for (int i = 0; i < numofIncoherentSums; i++) {

		nppsMaxIndx_32f(&deviceDataIncoherentSum[i*fftsize], fftsize, &deviceArrayMaxs[i], &deviceArrayPos[i], pDeviceBuffer);
	}
}

void stdCompute(int numofIncoherentSums, Npp32f *dataIncoherentSum, int fftsize,
	Npp32f *deviceArraystd, int *arrayPos, Npp8u * pDeviceBuffer, int peakRange) {

	int leftPeakIndex, rightPeakIndex, stdLength;
	stdLength = (fftsize / 2) - ((peakRange) / 2)-1;
	for (int i = 0; i < numofIncoherentSums; i++) {
		
		leftPeakIndex = arrayPos[i] - peakRange/2;
		rightPeakIndex = arrayPos[i] + peakRange/2;
		
		if (rightPeakIndex >= fftsize) {//case 2
			rightPeakIndex = rightPeakIndex % fftsize;
			//stdLength = leftPeakIndex - rightPeakIndex;
			nppsStdDev_32f(&dataIncoherentSum[i*fftsize+ rightPeakIndex], stdLength, &deviceArraystd[i], pDeviceBuffer);
		}
		else if (leftPeakIndex < 0) {//case 3
			leftPeakIndex = fftsize + leftPeakIndex;
			//stdLength = leftPeakIndex-rightPeakIndex ;
			nppsStdDev_32f(&dataIncoherentSum[i*fftsize + rightPeakIndex], stdLength, &deviceArraystd[i], pDeviceBuffer);
		}
		else {//case 1
			if (arrayPos[i] < fftsize / 2) {
				//stdLength = fftsize- rightPeakIndex;
				nppsStdDev_32f(&dataIncoherentSum[i*fftsize + rightPeakIndex], stdLength, &deviceArraystd[i], pDeviceBuffer);
			}
			else {
				//stdLength = leftPeakIndex;
				nppsStdDev_32f(&dataIncoherentSum[i*fftsize], stdLength, &deviceArraystd[i], pDeviceBuffer);
			}			
		}		
	}
}

//GLOBAL FUNCTIONS

__global__ void multip(int samples, cufftComplex *data1, cufftComplex *data2, int refsize)
{
	cufftComplex aux;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samples; i += stride) {
		int k = i % refsize;
		
		//(a+bi)*(c+di)=(ac−bd)+(ad+bc)i

		aux.x = data1[i].x* data2[k].x - data1[i].y*(-data2[k].y);
		aux.y = data1[i].x*(-data2[k].y) + data1[i].y*data2[k].x;
	
		
		data1[i].x = aux.x;
		data1[i].y = aux.y;
	}
}

__global__ void extendRefSignal(int samples, cufftComplex *data, int refsize) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samples; i += stride) {
		
		if (i >= refsize) {
			data[i] = data[i%refsize];

		}
	}
}

__global__ void applyDoppler(int samples, cufftComplex *data, float freqDoppler, float fs, unsigned long long samplePhaseMantain,
	int origSamples, int ddmQuant, int ddmRes, int fftsize)
{
	cufftComplex aux, aux2;
	float angle, freq, phasemantain;

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samples; i += stride) {
		phasemantain = ((i % (origSamples)) );//origninal samples signal
		freq = freqDoppler-(ddmRes * (ddmQuant/2)) +((i / (origSamples))*(ddmRes));
		angle = 2.0*PI*(phasemantain+samplePhaseMantain)*(freq /fs);
		aux2.x = cos(angle);
		aux2.y = sin(angle);

		//(a+bi)*(c+di)=(ac−bd)+(ad+bc)i
		aux.x = data[i].x*aux2.x - data[i].y*aux2.y;
		aux.y = data[i].x*aux2.y + data[i].y*aux2.x;

		data[i].x = aux.x;
		data[i].y = aux.y;
	}
}

__global__ void selectMaxs(int numOfFFT,int quantOfIncohSumAve, int ddmQuant, int *arrayPos, Npp32f *deviceArrayMaxs) {

	int step = numOfFFT / quantOfIncohSumAve;

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < step; i += stride) {
		//if (i < step) {
			for (int j = i; j < step*ddmQuant - 1; j = j + step) {
				if (deviceArrayMaxs[j + step] > deviceArrayMaxs[i]) {
					deviceArrayMaxs[i] = deviceArrayMaxs[j + step];
					arrayPos[i] = arrayPos[j + step];

				}
			}
		//}
	}
}

__global__ void savePeak(int numOfFFT, cufftComplex *dataFromIFFT, cufftComplex *dataToSave, int peakSamplesToSave,
	int quantOfIncohSumAve,int fftsize, int *arrayPos,int ddmQuant) {

	int samplesToSave = numOfFFT * peakSamplesToSave*ddmQuant;
	int posOnIFFT,fftOfThePeak,indexOfArrayPos, leftPosMaxOnOneIFFT, posOnOneIFFT;//rightPosMaxOnOneIFFT;
	
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samplesToSave; i += stride) {

		fftOfThePeak = i / peakSamplesToSave;//num of FFT in dataFromInv
		indexOfArrayPos = fftOfThePeak / quantOfIncohSumAve;//number of index in arrayPos
		//rightPosMaxOnOneIFFT = (arrayPos[indexOfPos] + peakSamplesToSave / 2);
		leftPosMaxOnOneIFFT = arrayPos[indexOfArrayPos % (numOfFFT / quantOfIncohSumAve)] - (peakSamplesToSave / 2);//begining of data to save
		posOnOneIFFT = leftPosMaxOnOneIFFT + (i%peakSamplesToSave);// sample of i in one fft

		if (posOnOneIFFT >= fftsize) {
			posOnIFFT = fftOfThePeak * fftsize + posOnOneIFFT%fftsize; //sample in the data from IFFT
			dataToSave[i] = dataFromIFFT[posOnIFFT];
			
			//case 2
		}
		else if (posOnOneIFFT < 0) {
			posOnIFFT = fftOfThePeak * fftsize + (fftsize+posOnOneIFFT);//sample in the data from IFFT
			dataToSave[i] = dataFromIFFT[posOnIFFT];
			//case 3
		}
		else {
			posOnIFFT = fftOfThePeak * fftsize + posOnOneIFFT;
			dataToSave[i] = dataFromIFFT[posOnIFFT];//sample in the data from IFFT
			//case 1
		}
	}
}

__global__ void inchoerentSum(int samplesInchoerentSum, cufftComplex *dataFromInv, Npp32f *dataStorageInocherentSum,
	int quantofAverageIncoherent, int fftsize)
{
	
	int indexofInv, numofSumM;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samplesInchoerentSum; i += stride) {
		dataStorageInocherentSum[i] = 0;
		numofSumM = i / fftsize;
		for (int k = 0; k < quantofAverageIncoherent; k++) {
			indexofInv = numofSumM*quantofAverageIncoherent*fftsize + k*fftsize+ i%fftsize;
			dataStorageInocherentSum[i] += dataFromInv[indexofInv].x*dataFromInv[indexofInv].x + dataFromInv[indexofInv].y*dataFromInv[indexofInv].y;
		}

	}
}


__global__ void scale(int samples, cufftComplex *data, int fftsize)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samples; i += stride) {
		data[i].x = data[i].x / float(fftsize);
		data[i].y = data[i].y / float(fftsize);

	}
}

__global__ void maskAndShift(char *devicedata, cuComplex *Dcomplexdata, int totalBytes)
{
	unsigned char k, aux;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < totalBytes; i += stride) {
		k = (unsigned char)(devicedata[i]);
		
		aux = k & ((unsigned) 1);
		aux = aux >> 0;
		Dcomplexdata[i * 4 + 0].x = float(2 * (aux) - 1);

		aux = k & ((unsigned)(1<<1));
		aux = aux >> 1;
		
		Dcomplexdata[i * 4 + 0].y = float(2 * (aux)-1);
		aux = k & ((unsigned)(1 << 2));
		aux = aux >> 2;
		Dcomplexdata[i * 4 + 1].x = float(2 * (aux)-1);
		aux = k & ((unsigned)(1 << 3));
		aux = aux >> 3;
		Dcomplexdata[i * 4 + 1].y = float(2 * (aux)-1);
		aux = k & ((unsigned)(1 << 4));
		aux = aux >> 4;
		Dcomplexdata[i * 4 + 2].x = float(2 * (aux)-1);
		aux = k & ((unsigned)(1 << 5));
		aux = aux >> 5;
		Dcomplexdata[i * 4 + 2].y = float(2 * (aux)-1);
		aux = k & ((unsigned)(1 << 6));
		aux = aux >> 6;
		Dcomplexdata[i * 4 + 3].x = float(2 * (aux)-1);
		aux = k & ((unsigned)(1 << 7));
		aux = aux >> 7;
		Dcomplexdata[i * 4 + 3].y = float(2 * (aux)-1);

	}
}
