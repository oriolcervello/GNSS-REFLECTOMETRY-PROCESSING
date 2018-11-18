Here will be the files of data to process.

The data will be stored in bits, one bit real one bit complex with symbols 1 or -1(for value 0).

Length of data with offset 0 in bytes will be = ((numofFFTs x (fftsize-overlap))+overlap)/4
For adding an offset just add it to the beggining and end of the data in the DATALINE

Num of FFTs has to be divisivable with quantofAverageIncoherent

In the ref/ directory will be the reference signals.
In the examples_in_float/ directory will be the examples of signals in float and instructions to set the program to read files of floats as input data.
