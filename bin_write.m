%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   BIN WRITE                                             %
%                   This function wirtes the input bin-file for the       %
%                   C++ plot program.                                     %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% Assigning variables
H = complex(0,-1);
rho = 2750;
g = 9.816;
nu = 0.3;
kappa = 3 - 4 * nu;

% Plot dim and resolutioin
xfrom = -1;
xto = 1;
yfrom = -1;
yto = 1;
Nx = 201;
Ny = 201;

% Write the bin-files fro C++
input_file = fopen('input_data.bin','w');
fwrite(input_file,[real(H),imag(H),rho,g,nu,kappa],'double');
fclose(input_file);

input_file = fopen('plot_data.bin','w');
fwrite(input_file,[xfrom,xto,yfrom,yto,Nx,Ny],'double');
fclose(input_file);