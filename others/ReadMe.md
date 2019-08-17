### DeviceProperties.cu

With this File you will be able to check you GPU properties. Specially Maximum Threads per block (maximum Blocksize of your GPU).

To build it you can run go to script/ directory and run build script with argument 1
    
	build.cmd 1 
	
or you can do it yourself:

    nvcc DeviceProperties.cu -o TestDevice

    .\TestDevice


