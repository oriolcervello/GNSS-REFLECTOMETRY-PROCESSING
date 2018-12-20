Here you will find different scripts.

### build.cmd :

This script is to build the program.

You have 3 options that you will need to pass as argument.

· 1 -> To build and run the test program DeviceProperties found in /others
It is useful to run it the first time to check that cuda is working fine and it is recognizing your GPU.
Also with this program you will be able to check your max Blocksize (maxThreadPerBlock) that will be asked in the input config of the program and you can not exceed it (it changes in between GPUs).


· 2 -> To build the program as it is meant to be.


· 3 -> To build the program with slight variations to read the input data from a .bin file also but reading floats instead of bits as unit of data. You will find more information about this in [datafiles/examples_in_float](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/tree/master/datafiles/examples_in_float)


To call the script open cmd in this directory, if for example you want option 1:

    .\build.cmd 1

### Matlab Scripts 

PlotDDM.m PlotSignal.m ReadMaxs.m

This scripts are for reading/plotting the output data of the program and helping to understand how is stored. You will find more useful info about this in the [/results](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/tree/master/results) directory where each type of output file is explained.