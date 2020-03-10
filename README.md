# GNSS-REFLECTOMETRY-PROCESSING

The main goal of the project is to process GNSS signals to obtain Earth measurements, originally for the MIR instrument. The aim is to process real GNSS-R data 
in conventional or interferometric mode
and able to compute DDMs at any coherent/incoherent correlation times. 
This program is using the GPU to accelerate the computation time.

All the contents were developed for the [Passive Remote Sensing Laboratory (RSLab)](http://www.tsc.upc.edu/rslab/Passive%20Remote%20Sensing/welcome) of the [Universitat Politècnica de Catalunya (UPC)](http://www.upc.edu/?set_language=en).

New versions of this program may be found at [GitHub](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/) 

Cite: Cervelló i Nogués, Oriol. Advanced GNSS-R signals processing with GPU. 2019
[http://hdl.handle.net/10230/42432](http://hdl.handle.net/10230/42432)

## Installation of CUDA

Here are a little instructions on how to install (in Windows) the program if you do not have CUDA.

Requirements of all you need a pc with a Nvidia GPU able to run CUDA 10.

Insturctions: 

1. Install [VisualStudio](https://visualstudio.microsoft.com/vs/), required by CUDA. Ensure to add C++ Packages, specially last version of  SDKs.
2. Install newest drivers. Look for yours [drivers](https://www.nvidia.com/Download/index.aspx) and install them.
3. Download [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads).
4. I recomend Personalized installation and just install CUDA option which is the CUDA itself. Uncheck other options.
5. Ensure to have on Path the following ones:
    c:\program files\nvidia gpu computing toolkit\cuda\v10.0\libnvvp;
    c:\program files\nvidia gpu computing toolkit\cuda\v10.0\bin;
    and the path to "cl.exe" found in the visual studio folder, can change between versions of VS but usually found inside "VC/"..."bin/" directories

Now CUDA should be working fine.

For testing CUDA I recomend go to others and try to run DeviceProperties.cu, there you will find how.

To build the program go to the Build and Run section below.

## Build and Run

To build the program you can do it with the GUI Install found in [/Gui](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/tree/master/GUI). There you will have also the instructions in how to do it.
Or you can build it with the script build.cmd found in [/scripts](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/tree/master/scripts). Also you will have also the instructions in how to do it there.
Althoug the originall programs is for the MIR instrument there are different options of building the program for adaptations to read different inputs types.


To run the program you can do it with the GUI GnssProcessing found in [/Gui](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/tree/master/GUI). There you will have also the instructions in how to do it.

Or you can do it by command line from cmd or powershell opened in the repo directory. 
To run it, we will need to pass 2 arguments:
1. input.ASE (or whatever name is your input configuration file with that format)
2. number of datalines to process on that run

Also on the PowerShell: (example: 1: input.ASE 2: 4 datalines)

    .\bin\main input.ASE 4


## Configuration File (input.ASE)

In this file you will need to fill each variable with the argument desired.

    *WRITEWAVEFORM 1 <---Produce de output of the waveform or just the maximums file, 0/1
    *FFTSIZE 32768  <---Size of the FFT, each one is a complex sample
    *NUMOFFFTS 20  <---# of FFT to do simultaneously in one dataline
    *QUANTINCOHAVER 4   <---# of coherent IFFT to averag for the incoherent
    *OVERLAP 32  <---Samples of overlap
    *FSAMPLING 32736000  <---Sampling Freq. [hz]
    *BLOCKSIZE 1024  <---Threads per blok(to know max thread/blok of your GPU go to other directory)
    *INTERFEROMETIC 1 <---Flag to run the program in interferometric(1) or conventionsl(0)
    *SAMPLESAVOIDMAX 0 <--- Samples to avoid at the begining of the correlation to compute the maximum
    *PEAKRANGESTD 65  <---Samples to avoid arround the peak for the STD computation (Peak in the mid sample)
    *PEAKSAMPLESTOSAVE 311  <---Samples to save arround the peak (Peak in the mid sample)
    *REFFILENAME datafiles/ref/clean_L1CA_14.bin <---Path of the reference file in conventional mode (avoided on interferometric)
    *RESULTSDIRECTORY results/ <---Path to store the results, be carefull they can overwrite
    *DDMFREQRES 25 <--- Resolution in Hz between to different dopplers if computing DDM
    *DDMNUMQUANT 1 <--- Quantity of dopplers to compute in a DDM, 1 for single computation
    *QUANTDATALINES 4  <--- Datalines to compute should be de same as argument 2

Datalines should be following the variables above and should have this format if INTERFEROMETRIC=0 (conventional mode): DataFileName # BeginingOfData EndOfData DopplerFreq

    *DATALINE 1 datafiles/prn_L1CA_32_100.bin 0 654752 0 datafiles/prn_L1CA_32.bin

The 1 just after the DATALINE indicates that all the samples for that line are from that archive. If a dataline have half of the samples at the end of an archive and the other half at the beggining of the following to indicate that it is the same dataline a 2 must be set in both datalines and eachone with its path and samples in that file.

if INTERFEROMETRIC=1 (interferometric mode): DataFileName BeginingOfData EndOfData DopplerFreq RefFileName BeginingOfDataRef

    *DATALINE 1 datafiles/prn_L1CA_32_100.bin 0 654752 0 datafiles/prn_L1CA_32.bin 100

The length of the Reference will be the same as the Data but we can set an offset, ex: in this line 100. If interferometric the line above * REFFILENAME will be avoided.

The BeginingOfData and EndOfData numbers are bytes by default. To compute them easy starting with offset 0, BeginingOfData=0 and EndOfData=((numofFFTs*(fftsize-overlap))+overlap)* 2/8. We multiply by 2 because they are complex samples and divide by 8 as we want bytes.

The following line would be BeginingOfData=previous_EndOfData (as it reads from BeginingOfData to EndOfData-1 ) and EndOfData=BeginingOfData+((numofFFTs*(fftsize-overlap))+overlap)* 2/8.
Any line with offset would be BeginingOfData= whatever offset in complex bytes and EndOfData=BeginingOfData+((numofFFTs*(fftsize-overlap))+overlap)* 2/8.

In case that the program is built to read different type of input (float or int16) we will be reading samples directly BeginingOfData=0 and EndOfData=((numofFFTs*(fftsize-overlap))+overlap) and the following lines as BeginingOfData=previous_EndOfData and EndOfData=BeginingOfData+((numofFFTs*(fftsize-overlap))+overlap).

In datafiles/ReadMe.md you will find how the data is structured and more useful info on how to build the datalines.
In datafiles/examples_in_float you will find that is possible to pass data in float form.

You will find an example in this repo of an 'input.ASE'

## Outputs
The results will be stored on the directory set on configFile. If we set /results the output files will be there. Be careful that if you run again the program results will overwrite if the same result folder is selected.

There is a ReadMe.md on /results directory explaining the output files and how to read them.

## Licence
You may find it in a specific licence file.

## Credit
Project developed by Oriol Cervelló.
