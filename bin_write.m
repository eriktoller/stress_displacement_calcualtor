%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   BIN WRITE                                             %
%                   This function wirtes the input bin-file for the       %
%                   C++ plot program.                                     %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clearvars
clc

% Assigning variables
H = complex(0,-1);
rho = 2750;
g = 9.816;
nu = 0.3;
kappa = 3 - 4 * nu;

% Plot dim and resolution
xfrom = -1;
xto = 1;
yfrom = -1;
yto = 1;
Nx = 201;
Ny = 201;

% Write the bin-files for C++
A = [real(H),imag(H),rho,g,nu,kappa]; % Vector to write
input_file = fopen('input_data.bin','w');
fwrite(input_file,A,'double');
fclose(input_file);

B = [xfrom,xto,yfrom,yto,Nx,Ny]; % Vector to write
plot_file = fopen('plot_data.bin','w');
fwrite(plot_file,B,'double');
fclose(plot_file);