Here will be the files of data to process.

The data will be stored in bits, one bit real one bit complex with symbols 1 or -1(for value 0). One byte will be 4 complex numbers. Program will read size of bytes. Inside each byte bits should be stored in big-endian order.

If you want to perform a given numofFFTs of fftsize with whatever overlap, the length of data with offset 0 in bytes that you will be sting in your dataline = ((numofFFTs x (fftsize-overlap))+overlap)/4
For adding an offset just add it to the beggining and end of the data in the DATALINE.

The phase of the doppler of the end of the DATALINE 1 is maintained in the beggining of the DATALINE 2.

Num of FFTs has to be divisivable with quantofAverageIncoherent

In the ref/ directory will be the reference signals.

In the examples_in_float/ directory will be the examples of signals in float and instructions to set the program to read files of floats as input data.
