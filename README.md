# DFT using MPI
Construct a two dimensional Discrete Fourier Transform (DFT) of a square image using the Message Passing Interface (MPI) protocol. 

Sixteen CPUs are used for the calculation, each working on certain sections of the image only, with one CPU acting as master. Each CPU first performs a (row-wise) 1-D DFT on a section of the image and then sends its data to the master CPU via MPI, which then aggregates all the data into a complete 1-D DFT image. 

This 1-D DFT image is then transposed and sections of it are sent back to the worker CPUs for another round of 1-D DFT's - only this time, it is column-wise due to the prior transpose. Once this is done, the master CPU again aggregates all the data from the worker CPUs to form the final 2-D DFT of the original image.
