#include"functions.cuh"


void readdata(int N, cufftComplex *data, string name, bool readbinary) {
	if (readbinary == true) {
		readdatabinary(N, data, name);
	}
	else {
		readdatatxt(N, data, name);
	}
}

void writedata(int N, cufftComplex *data, string name, bool writebinary) {
	if (writebinary == true) {
		writedatabinary(N, data, name);
	}
	else {
		writedatatxt(N, data, name);
	}
}


void readdatabinary(int N, cufftComplex *data, string name)
{
	ifstream myfile;
	myfile.open(name, ios::binary);
	float num1,num2;
	
	if (myfile.is_open())
	{
		int k = 0;
		while (k < N)
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


void readdatatxt(int N, cufftComplex *data, string name)
{
	ifstream myfile;
	myfile.open(name, ios::binary);
	string line;
	

	if (myfile.is_open())
	{
		int k = 0;
		while (k < N)
		{
			getline(myfile, line, '\n');
			data[k].x = stof(line);
			data[k].y = 0;

			k++;
		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

void writedatatxt(int N, cufftComplex *data1, string name) {

	ofstream myfile;
	myfile.open(name);
	if (myfile.is_open())
	{
		for (int ii = 0; ii < N; ii++)
		{
			 myfile<< ii << " " << (data1[ii].x)<< " " << data1[ii].y <<"\n"; 
	
		}
		myfile.close();
	}

	else cout << "Unable to open file\n";
}

void writedatabinary(int N, cufftComplex *data1, string name) {

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

void writetime(int N, string name, long long *readtime, long long *shifttime, long long *ffttime,
	long long *multime, long long *ifftime, long long *writetime, long long *looptime) {

	ofstream myfile;
	myfile.open(name);
	if (myfile.is_open())
	{
		myfile << "Atempt\t\tReadT.\t\tShiftT.\t\tFFTT.\t\tMulT.\t\tIFFTT.\t\tWriteT.\t\tLoopT." << "\n";
		for (int ii = 0; ii < N; ii++)
		{
			myfile << ii << "\t\t" << readtime[ii] << "\t\t" << shifttime[ii] << "\t\t" << ffttime[ii] << "\t\t" << multime[ii];
			myfile << "\t\t" << ifftime[ii] << "\t\t" << writetime[ii] << "\t\t" << looptime[ii] << "\n";

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

/*Makes the complex conjugate of data2 and multiply point by point data1 and data2
**data1 and data2 should be on device memory!!
N: length of data
data1: cufftComplex data set 1
data1: cufftComplex data set 2*/
__global__ void multip(int n, cufftComplex *data1, cufftComplex *data2, int refsize)
{
	cufftComplex aux;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < n; i += stride) {
		int k = i % refsize;
		
		//(a+bi)*(c+di)=(ac−bd)+(ad+bc)i

		aux.x = data1[i].x* data2[k].x - data1[i].y*(-data2[k].y);
		aux.y = data1[i].x*(-data2[k].y) + data1[i].y*data2[k].x;
	
		
		data1[i].x = aux.x;
		data1[i].y = aux.y;
	}
}

__global__ void extendRefSignal(int n, cufftComplex *data, int refsize) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < n; i += stride) {
		
		if (i >= refsize) {
			data[i] = data[i%refsize];

		}
	}
}



__global__
void maskandshift(char *devicedata, cuComplex *Dcomplexdata, int N)
{
	unsigned k;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;
	for (int i = index; i < N; i += stride) {

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

void readrealdata(int N, char *data, string name)
{

	ifstream myfile;
	myfile.open(name, ios::binary);
	if (myfile.is_open())
	{
		myfile.read(data, N);

		/*int k = 0;
		while (k < N)
		{
			k++
		}*/
		myfile.close();
	}
	else cout << "Unable to open file";
}

