#include"functions.cuh"
#include "extra/TextParser.cuh"

void readConfig(const char *configFileName, int numofDataLines, int *fftsize, int *numofFFts, int *overlap, int *fSampling, int *blockSize, int *peakRangeStd, int *peakSamplesToSave,
	int* quantOfAverIncoh, int *dataOffsetBeg, int *dataOffsetEnd, int *doppler, string *fileNames,string *fileRefNames) {

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
		fileRefNames[i] = t.getword();
	}
}

void checkInputConfig(int argc, const char **argv, int numofDataLines, int fftsize, int numofFFts, int overlap, int fSampling,  int blockSize, int peakRangeStd, int peakSamplesToSave,
	int quantOfAverIncoh,  int *dataOffsetBeg, int *dataOffsetEnd, int *doppler, string *fileNames, string *fileRefNames) {

	if (argc != 3) {
		cout << "Error: Wrong number of arguments\n"; 
		exit(0);
	}

	cout << "\n" << "Quant of args: " << argc << "\n";
	cout << "First: " << argv[0] << "\n";
	cout << "Second: " << argv[1] << "\n";
	cout << "Third: " << argv[2] << "\n\n";

	cout << "FFT Size: " << fftsize << "\n";
	cout << "Num. of FFT: " << numofFFts << "\n";
	cout << "Overlap: " << overlap << "\n";
	cout << "FSampling: " << fSampling << "\n";
	cout << "Quant of averg Inch.: " << quantOfAverIncoh << "\n";
	cout << "Blok Size: " << blockSize << "\n";
	cout << "Peak samples for the std: " << peakRangeStd << "\n";
	cout << "Peak samples to save: " << peakSamplesToSave << "\n";

	cout << "Num of data lines: " << numofDataLines << "\n";
	cout << "Data lines: \n";
	for (int i = 0; i < numofDataLines; i++) {
		cout << fileNames[i] << "  ";
		cout << dataOffsetBeg[i] << "  ";
		cout << dataOffsetEnd[i] << "  ";
		cout << doppler[i] << " ";
		cout << fileRefNames[i] << "\n";

	}

}

void readdata(int length,int offsetFromBeg, cufftComplex *data, string name)
{
	ifstream myfile;
	myfile.open(name, ios::binary);
	float num1,num2;
	
	if (myfile.is_open())
	{
		myfile.seekg(offsetFromBeg*2 * sizeof(float));
		int k = 0;
		while (k < length)
		{
			
			myfile.read((char*)&num1, sizeof(num1));
			myfile.read((char*)&num2, sizeof(num2));
			data[k].x = num1;
			data[k].y =  num2;
			k++;
		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

void readRealData(int length, int offsetFromBeg, int bytesToRead,char *data, string name)
{
	if (length > bytesToRead) {
		cout << "Error: iput length bigger than bytesToRead";
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
			cout << bytesToRead- length << "Warning: length smaller than bytesToRead, Bytes filled with 0 \n";
			memset(&data[length], 0, bytesToRead - length);
		}
	}
	else cout << "Unable to open file";
}


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

	else cout << "Unable to open file\n";
}

void writeMaxs(int N, Npp32f *dataMaxValue, int *dataMaxPos, Npp32f *hostarrayStd,string name) {

	ofstream myfile;
	myfile.open(name, ios::app);
	if (myfile.is_open())
	{
		for (int ii = 0; ii < N; ii++)
		{
			myfile <<"Pos: "<< dataMaxPos[ii]<<" Value: " << dataMaxValue[ii] << " STD: " << hostarrayStd[ii] << "\n";

		}
		myfile << "\n --------------------------------------------------------------------- \n\n";
		myfile.close();
	}

	else cout << "Unable to open file\n";
}

void writedata(int N, cufftComplex *data1, string name) {

	ofstream myfile;
	myfile.open(name, ios::binary);
	if (myfile.is_open())
	{
		for (int ii = 0; ii < N; ii++)
		{
			
			myfile.write((char*)&data1[ii].x, sizeof(float));
			myfile.write((char*)&data1[ii].y, sizeof(float));
		}
		myfile.close();
	}

	else cout << "Unable to open file\n";
}

void writetime(int N, string name, long long *readtime, long long *writetime, long long *looptime) {

	ofstream myfile;
	myfile.open(name);
	if (myfile.is_open())
	{
		myfile << "Atempt\t\tReadT.\t\tWriteT.\t\tLoopT." << "\n";
		for (int ii = 0; ii < N; ii++)
		{
			myfile << ii << "\t\t" << readtime[ii] << "\t\t" << writetime[ii] << "\t\t" << looptime[ii] << "\n";

		}
		myfile.close();
	}

	else cout << "Unable to open file";
}

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

void maxCompute(int numofIncoherentSums, Npp32f *deviceDataIncoherentSum, int fftsize, Npp32f *deviceArrayMaxs,
	 int *deviceArrayPos, Npp8u * pDeviceBuffer) {

	for (int i = 0; i < numofIncoherentSums; i++) {

		nppsMaxIndx_32f(&deviceDataIncoherentSum[i*fftsize], fftsize, &deviceArrayMaxs[i], &deviceArrayPos[i], pDeviceBuffer);
	}
}


void stdCompute(int numofIncoherentSums, Npp32f *dataIncoherentSum, int fftsize,
	Npp32f *deviceArraystd, int *arrayPos, Npp8u * pDeviceBuffer, int peakRange) {

	int leftPeakIndex, rightPeakIndex, stdLength;

	for (int i = 0; i < numofIncoherentSums; i++) {
		
		leftPeakIndex = arrayPos[i] - peakRange/2;
		rightPeakIndex = arrayPos[i] + peakRange/2;
		
		if (rightPeakIndex >= fftsize) {//case 2
			rightPeakIndex = rightPeakIndex % fftsize;
			stdLength = leftPeakIndex - rightPeakIndex;
			nppsStdDev_32f(&dataIncoherentSum[i*fftsize+ rightPeakIndex], stdLength, &deviceArraystd[i], pDeviceBuffer);
		}
		else if (leftPeakIndex < 0) {//case 3
			leftPeakIndex = fftsize + leftPeakIndex;
			stdLength = leftPeakIndex-rightPeakIndex ;
			nppsStdDev_32f(&dataIncoherentSum[i*fftsize + rightPeakIndex], stdLength, &deviceArraystd[i], pDeviceBuffer);
		}
		else {//case 1
			if (arrayPos[i] < fftsize / 2) {
				stdLength = fftsize- rightPeakIndex;
				nppsStdDev_32f(&dataIncoherentSum[i*fftsize + rightPeakIndex], stdLength, &deviceArraystd[i], pDeviceBuffer);
			}
			else {
				stdLength = leftPeakIndex;
				nppsStdDev_32f(&dataIncoherentSum[i*fftsize], stdLength, &deviceArraystd[i], pDeviceBuffer);
			}			
		}		
	}
}

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


__global__ void applyDoppler(int samples, cufftComplex *data, float freqDoppler, float fs, unsigned long long samplePhaseMantain)
{
	cufftComplex aux, aux2;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samples; i += stride) {
		aux2.x=cos(2.0*PI*(i+ samplePhaseMantain)*(float(freqDoppler)/float(fs)));
		aux2.y= sin(2.0*PI*(i+ samplePhaseMantain)*(float(freqDoppler) / float(fs)));

		//(a+bi)*(c+di)=(ac−bd)+(ad+bc)i
		aux.x = data[i].x*aux2.x - data[i].y*aux2.y;
		aux.y = data[i].x*aux2.y + data[i].y*aux2.x;

		data[i].x= aux.x;
		data[i].y= aux.y;
		 
	
	}
}

__global__ void savePeak(int numOfFFT, cufftComplex *dataFromIFFT, cufftComplex *dataToSave, int peakSamplesToSave,
	int quantOfIncohSumAve,int fftsize, int *arrayPos) {

	int samplesToSave = numOfFFT * peakSamplesToSave;
	int posOnIFFT,fftOfThePeak,indexOfArrayPos, leftPosMaxOnOneIFFT, posOnOneIFFT;//rightPosMaxOnOneIFFT;
	
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samplesToSave; i += stride) {

		fftOfThePeak = i / peakSamplesToSave;//num of FFT in dataFromInv
		indexOfArrayPos = fftOfThePeak / quantOfIncohSumAve;//number of index in arrayPos
		//rightPosMaxOnOneIFFT = (arrayPos[indexOfPos] + peakSamplesToSave / 2);
		leftPosMaxOnOneIFFT = arrayPos[indexOfArrayPos] - (peakSamplesToSave / 2);//begining of data to save
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
	unsigned k;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < totalBytes; i += stride) {

		k = (unsigned)(devicedata[i]);
		Dcomplexdata[i * 4 + 3].y = (k % 2);
		k = k / 2;
		Dcomplexdata[i * 4 + 3].x = (k % 2);
		k = k / 2;
		Dcomplexdata[i * 4 + 2].y = (k % 2);
		k = k / 2;
		Dcomplexdata[i * 4 + 2].x = (k % 2);
		k = k / 2;
		Dcomplexdata[i * 4 + 1].y = (k % 2);
		k = k / 2;
		Dcomplexdata[i * 4 + 1].x = (k % 2);
		k = k / 2;
		Dcomplexdata[i * 4 + 0].y = (k % 2);
		k = k / 2;
		Dcomplexdata[i * 4 + 0].x = (k % 2);


	}
}
