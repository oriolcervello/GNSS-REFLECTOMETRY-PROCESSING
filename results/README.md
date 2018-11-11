Here will be stored the results

### Maximums File:

Here are maximum value and position, the std and doppler of each dataline.

### PeakIterationX.bin:

Peaks and surrounding of the correlation of each coherent FFT.

One file per dataline. X = num of dataline

Information is stored in binary and are complex floats (32bits real + 32 bits imaginary) stored as:

Re1 Co1 Re2 Co2 .... Re311 Co311 (311 becouse it was decleared at PEAKSAMPLESTOSAVE, then the following peak) Re1 Co1 ...

Find in others directory how to read/plot it.

### Times:

Times of comutation in microseconds per dataline.
