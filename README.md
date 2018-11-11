# GNSS-REFLECTOMETRY-PROCESSING

This program is used to process GNSS-R signals using the GPU to accelerate the computation time. 

All the contents were developed for the [Passive Remote Sensing Laboratory (RSLab)](http://www.tsc.upc.edu/rslab/Passive%20Remote%20Sensing/welcome) of the [Universitat Politècnica de Catalunya (UPC)](http://www.upc.edu/?set_language=en).

New versions of this program may be found at [GitHub](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/) 

## Contents
Project is still under development, further information and modifications will be added.

## Installation of CUDA

Here are a little instructions on how to run (in Windows) the program if you do not have CUDA.

Requirements of all you need a pc with a Nvidia GPU able to run CUDA.

Insturctions: 

1. Install [VisualStudio](https://visualstudio.microsoft.com/vs/), required by CUDA. Ensure to add C++ Packages, specially SDKs.
2. If you got old GPU drivers or you want to update them I recomend running [DDU](https://www.guru3d.com/files-details/display-driver-uninstaller-download.html) to clean the GPU before installing new drivers.
3. Download [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads).
4. I recomend Personalized installation and just install CUDA option which is the CUDA itself and the Driver component if you want to update it. If you want to update drivers don't forget step 2.
5. Ensure to have on Path the following ones:
c:\program files\nvidia gpu computing toolkit\cuda\v10.0\libnvvp;
c:\program files\gnuplot\bin;c:\program files\nvidia gpu computing toolkit\cuda\v10.0\bin;
c:\program files (x86)\microsoft visual studio\2017\community\vc\tools\msvc\14.15.26726\bin\hostx64\x64;

Now CUDA should be working fine.

To build the program go to the Build and Run section.

## Build and Run

To build the program open the directory on the Powershell.

You can build the program with the following comand:

    nvcc extra/TextParser.cu functions.cu kernel.cu -lnpps -lcufft -o kernel

To run it, we will need to pass 2 arguments:
1. input.ASE (or whatever name is your input configuration file with that format)
2. number of datalines to process on that run
Also on the PowerShell: (example: 1: input.ASE 2: 4 datalines)

    .\kernel input.ASE 4


## Configuration File (input.ASE)
In this file you will need to fill each variable with the argument desired.

    *FFTSIZE 32768  <---Size of the FFT
    *NUMOFFFTS 20  <---# of FFT to do simultaneously in one dataline
    *QUANTINCOHAVER 4   <---# of coherent IFFT to averag for the incoherent
    *OVERLAP 32  <---Samples of overlap
    *FSAMPLING 32736000  <---Sampling Freq. [hz]
    *BLOCKSIZE 1024  <---Threads per blok(to know max thread/blok of your GPU go to other directory)
    *PEAKRANGESTD 64  <---Samples to avoid arround the peak for the STD computation (Peak in the mid sample)
    *PEAKSAMPLESTOSAVE 311  <---Samples to save arround the peak (Peak in the mid sample)
    *QUANTDATALINES 4  <--- Datalines to compute should be de same as argument 2

Datalines should have this format: DataFileName BeginingOfData EndOfData DopplerFreq RefSignalFileName

    *DATALINE datafiles/prn_L1CA_32_100.bin 0 654752 0 datafiles/prn_L1CA_32_100.bin



## Outputs
The results will be stored on results directory.

The outputs are:

·Maximums File:
  Here are maximum value and position, the std and doppler of each dataline.
·PeakIterationX.bin:
  Peaks and surrounding of the correlation of each coherent FFT.
  One file per dataline. X = num of dataline
  Information is stored in binary and ar complex floats (32bits real + 32 bits iamginary) stored as:
  Re1 Co1 Re2 Co2 .... Re311 Co311 (311 becouse it was decleared at PEAKSAMPLESTOSAVE, then the following peak) Re1 Co1 ...
  Find in others how to read/plot it.
·Times:
  Times of comutation in microseconds per dataline.

## Licence
You may find it in a specific licence file.

## Contact
Project developed by Oriol Cervelló.
