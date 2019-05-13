#include "GlobalFunc.cuh"

//GLOBAL FUNCTIONS

__global__ void multip(int samples, cufftComplex *data1, cufftComplex *data2, int refsize, bool interferometric)
{
	cufftComplex aux;
	int k;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samples; i += stride) {
		if (interferometric == true) {
			k = i;
		}
		else {
			k = i % refsize;
		}
		

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

__global__ void applyDoppler(int samples, cufftComplex *data, float freqDoppler, float fs, unsigned long long int samplePhaseMantain,
	int origSamples, int ddmQuant, int ddmRes, int fftsize)
{
	cufftComplex aux, aux2;
	double angle, freq, phasemantain;

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samples; i += stride) {
		
		if (ddmQuant > 1) {
			phasemantain = ((i % (origSamples)));//origninal samples signal
			freq = freqDoppler - (ddmRes * (ddmQuant / 2)) + ((i / (origSamples))*(ddmRes));
		}
		else{ 
			freq = freqDoppler;
			phasemantain = i;//origninal samples signal
		}
		angle = 2.0*PI*double((phasemantain) + samplePhaseMantain)*(freq / double(fs));
		//angle = 2.0*PI*(phasemantain + 0)*(freq / fs);
		aux2.x = cos(angle);
		aux2.y = sin(angle);

		//(a+bi)*(c+di)=(ac−bd)+(ad+bc)i
		aux.x = data[i].x*aux2.x - data[i].y*aux2.y;
		aux.y = data[i].x*aux2.y + data[i].y*aux2.x;

		data[i].x = aux.x;
		data[i].y = aux.y;
	}
}

__global__ void selectMaxs(int numOfFFT, int quantOfIncohSumAve, int ddmQuant, int *arrayPos, Npp32f *deviceArrayMaxs) {

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
	int quantOfIncohSumAve, int fftsize, int *arrayPos, int ddmQuant) {

	int samplesToSave = numOfFFT * peakSamplesToSave*ddmQuant;
	int posOnIFFT, fftOfThePeak, indexOfArrayPos, leftPosMaxOnOneIFFT, posOnOneIFFT;//rightPosMaxOnOneIFFT;

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < samplesToSave; i += stride) {

		fftOfThePeak = i / peakSamplesToSave;//num of FFT in dataFromInv
		indexOfArrayPos = fftOfThePeak / quantOfIncohSumAve;//number of index in arrayPos
		//rightPosMaxOnOneIFFT = (arrayPos[indexOfPos] + peakSamplesToSave / 2);
		leftPosMaxOnOneIFFT = arrayPos[indexOfArrayPos % (numOfFFT / quantOfIncohSumAve)] - (peakSamplesToSave / 2);//begining of data to save
		posOnOneIFFT = leftPosMaxOnOneIFFT + (i%peakSamplesToSave);// sample of i in one fft

		if (posOnOneIFFT >= fftsize) {
			posOnIFFT = fftOfThePeak * fftsize + posOnOneIFFT % fftsize; //sample in the data from IFFT
			dataToSave[i] = dataFromIFFT[posOnIFFT];

			//case 2
		}
		else if (posOnOneIFFT < 0) {
			posOnIFFT = fftOfThePeak * fftsize + (fftsize + posOnOneIFFT);//sample in the data from IFFT
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
			indexofInv = numofSumM * quantofAverageIncoherent*fftsize + k * fftsize + i % fftsize;
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

		aux = k & ((unsigned)1);
		aux = aux >> 0;
		Dcomplexdata[i * 4 + 0].x = float(2 * (aux)-1);

		aux = k & ((unsigned)(1 << 1));
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


__global__ void copyInt2Float(__int16 *deviceIntData, cuComplex *deviceFloatData, int length) {

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < length/2; i += stride) {
		deviceFloatData[i].x = float(deviceIntData[2*i]);
		deviceFloatData[i].y = float(deviceIntData[2*i+1]);
	}


}