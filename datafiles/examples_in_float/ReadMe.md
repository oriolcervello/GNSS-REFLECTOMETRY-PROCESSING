To build the program to be able to read this type of signals run the scripts/build.cmd with argument 3, or if you use the GUI build the version to read float inputs

This are files of example, you will be able to see how the data is structured:

### prn_L1CA_32_100.bin

Contains 100 times prn signal of a size of 32736 samples each. 
(Each sample is a complex number of float of 32bits for the real part and  float of 32bits for imaginary part)

### prn_L1CA_32_100_fd_1e3.bin

The same as prn_L1CA_32_100.bin but with a doppler of 1000 Hz.

### prn_L1CA_32_100_noisy.bin

The same as prn_L1CA_32_100.bin but with added random noise.

### prn_L1CA_32.bin

Contains 1 times prn signal of a size of 32736 samples each. Used for Reference signal for a DATALINE.
(Each sample is a complex number of float of 32bits for the real part and  float of 32bits for imaginary part)
