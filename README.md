# Stress Field Calculator
A program which calculated the stress field for given functions and grid. The input files for the plot properties and the physical constants are provided through the MATLAB program `bin_write.m` and the program `bin_read.m` can be sued to plot the results calculated by `stress_field_calculator.cpp`.

At the moment the programs does not sync automatically and must be manually run in both MATLAB and the preferred C++ editor. Also, the function for the stresses are coded in the C++ file and changing these, or appending new functions, requires changes in the main code. 

The program now includes the analytic element for gravity, and it produces the Cartesian stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{11}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{22}"/>  & <img src="https://latex.codecogs.com/gif.latex?\sigma_{12}"/> ) and the principal stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{1}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{2}"/>  & <img src="https://latex.codecogs.com/gif.latex?\theta_{p}"/> ). For a given resolution and coordinates.

The user will have to define the function files for the stresses. The program will the use those to calculate and export the stress field for a given grid.

A log file, `log.txt`, is created with all the variable data and time stamp.

## Instructions
To generate plot for stresses filed (here only for gravity at the moment) follow the following procedure:
1. Run `bin_write.m`
   - Generates the input data and plot data in `input_data.bin` and `plot_data.bin`.
2. Run `stress_field_calculator.cpp`
   - Calculates the stress fields and provide them in `data.bin` and the dimensions in `dim_data.bin`.
3. Run `bin_read.m`
   - Reads the data and plots the stress fields.

Created by,
Erik Toller
