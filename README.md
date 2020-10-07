# Stress Field Calculator
A program which calculated the stress field for given functions and grid. The input files for the plot properties and the physical constants are provided through the MATLAB program `bin_write.m` and the program `bin_read.m` can be sued to plot the results calculated by `stress_field_calculator.cpp`.

At the moment the programs does not sync automatically and must be manually run in both MATLAB and the preferred C++ editor. Also, the function for the stresses are coded in the C++ file and changing these, or appending new functions, requires changes in the main code. 

The program now includes the analytic element for gravity, and it produces the Cartesian stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{11}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{22}"/>  & <img src="https://latex.codecogs.com/gif.latex?\sigma_{12}"/> ), the principal stress field (<img src="https://latex.codecogs.com/gif.latex?\sigma_{1}"/> , <img src="https://latex.codecogs.com/gif.latex?\sigma_{2}"/>  & <img src="https://latex.codecogs.com/gif.latex?\theta_{p}"/> ) and the principal stress trajectories. For a given resolution and coordinates.

A log file, `log.txt`, is created with all the variable data and time stamp.

This program has been developed using *Microsoft Visual Stuido* and only the `.cpp` file is included in the repository. 

## Instructions
To generate plot for stresses filed (here only for gravity at the moment) follow the following procedure:
1. Run `bin_write.m`
   - Generates the input data and plot data in `input_data.bin` and `plot_data.bin`.
2. Run `stress_field_calculator.cpp`
   - Calculates the stress fields and provide them in `data.bin` and the dimensions in `dim_data.bin`.
3. Run `bin_read.m`
   - Reads the data and plots the stress fields.
   
## Included
The program inludes several analytic element for linear elasticity. Currently the following anlaytic elements are included:
- gravity,
- cracks (not yet released) and
- circular tunnel (not yet released).

## Input data
This list contains all the definitions of the user input data defined in `bin_wirte.m`. These are the model properties:
- `H` the elevation for where the gravity is set to zero
- `rho` the density of the elastic medium
- `g` the Newtonian constant of gravitation
- `nu` the Poisson's ratio
- `kappa` the bulk modulus (is calculated by default with `kappa = 3-4*nu`)
- `nc` *inclued with cracks*
- `m` *inclued with cracks*
- `z1` *inclued with cracks*
- `z2` *inclued with cracks*
- `L` *inclued with cracks*
- `mu` *inclued with cracks*
- `beta` *inclued with cracks*
- `nt` *inclued with circular tunnel*
- `mt` *inclued with circular tunnel*
- `z0` *inclued with circular tunnel*
- `R` *inclued with circular tunnel*
- `a` *inclued with circular tunnel*
- `b` *inclued with circular tunnel*

These are the plotting properties:
- `xfrom` the starting value for x-axis
- `xto` the end value for the x-axis
- `yfrom` the starting value for the y-axis
- `yto` the end value for the y-axis
- `Nx` the number of grid points in the x-diraction
- `Ny` the numbe rof grid points in the y-direction
- `Ntraj` the number of steps for the stress trajectories
- `lvs_traj` the number of stress trajectories

__________________________________________________________________________________________
Created by,

Erik Toller
