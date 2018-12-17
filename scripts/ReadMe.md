### build.cmd :

This script is to build the program.

You have 3 options that you will need to pass as argument.

·1 -> To build and run the test program DeviceProperties found in /others
It is useful to run it the first time to check that cuda is working fine and it is recognizing your GPU.
Also with this program you will be able to check your max Blocksize (maxThreadPerBlock) that will be asked in the input config of the program and you can not exceed it (it changes in between GPUs).


·2 -> To build the program as it is meant to be.


·3 -> To build the program with slight variations to read the input data from a .bin file also but as floats instead of bits. You will find more information about this in [datafiles/examples_in_float](https://github.com/oriolcervello/GNSS-REFLECTOMETRY-PROCESSING/tree/master/datafiles/examples_in_float)


To call the script open cmd in this directory