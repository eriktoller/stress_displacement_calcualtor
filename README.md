# Stress Field Calculator
A program which calculated the stress field for given functions and grid. The input files for the plot properties and the physical constants are provided through the MATLAB program `bin_write.m` and the program `bin_read.m` can be sued to plot the results calculated by `stress_field_calculator.cpp`.

At the moment the programs does not sync automatically and must be manually run in both MATLAB and the preferred C++ editor. Also, the function for the stresses are coded in the C++ file and changing these, or appending new functions, requires changes in the main code. 

The program now includes the analytic element for gravity, and it produces the Cartesian stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{11}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{22}"/>  & <img src="https://latex.codecogs.com/gif.latex?\sigma_{12}"/> ) and the principal stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{1}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{2}"/>  & <img src="https://latex.codecogs.com/gif.latex?\theta_{p}"/> ). For a given resolution and coordinates.

The user will have to define the function files for the stresses. The program will the use those to calculate and export the stress field for a given grid.


Created by,
Erik Toller
