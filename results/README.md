Here will be stored the results

### Maximums File:

Here are maximum value and position, the std and doppler of each dataline.

To Read this files you can play the script [ReadMaxs.m](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/blob/master/results/ReadMaxs.m) found in this directory.

### PeakIterationX.bin with NO DDM (DDMQUANT=1):

Peaks and surrounding of the correlation of each coherent FFT.

One file per dataline. X = num of dataline

Information is stored in binary and are complex floats (32bits real + 32 bits imaginary) stored as:

Re1 Co1 Re2 Co2 .... Re311 Co311 (311 becouse it was decleared at PEAKSAMPLESTOSAVE, if not whatever the value setted, then the following peak) Re1 Co1 ...

To read/plot this type of file you can play the [PlotSignal.m](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/blob/master/results/PlotSignal.m) found in this directory.

### PeakIterationX.bin with DDM (DDMQUANT>1):

Peaks and surrounding of the correlation of each coherent FFT repeated DDMQUANT times to compute the different dopplers.
To build a DDM we need to read 1 coherent of each group of dopplers to have the same coherent DDMQUANT times and be able to construct the Delay Doppler Map.

Information is stored in binary and are complex floats (32bits real + 32 bits imaginary) stored as:

Re1 Co1 Re2 Co2 .... Re311 Co311 (311 becouse it was decleared at PEAKSAMPLESTOSAVE, if not whatever the value setted, then the following peak) Re1 Co1 ...

To read/plot this type of file you can play the [PlotDDM.m](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/blob/master/results/PlotDDM.m) found in this directory.

### Times:

Times of comutation in microseconds per dataline.
