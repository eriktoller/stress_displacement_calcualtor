function [Nx, Ny, Ntraj, lvs_traj, x_vec, y_vec, grid_11, grid_22, grid_12, grid_1, grid_2, theta_p, traj_1, traj_2] = bin_read()
% BIN_READ Reads the binary files form the C++ computations
%   This program read the output data from the program
%   stress_field_calculator.cpp. It then prodcues the plots for the 
%   different stress fields.
close all
clearvars
clc

if_plot = 1; % set to 1 for plots to run

disp('Load the data from the bin-file')
data_file = fopen('data.bin', 'r');
dim_file = fopen('dim_data.bin', 'r');
[A,~] = fread(data_file,'double');
[B,~] = fread(dim_file,'double');
disp('Completed')
disp(' ')
disp('Getting the vecotrs and grids')
Nx = B(1);
Ny = B(2);
Ntraj = B(3);
lvs_traj = B(4);
x_vec = A(1:Nx);
y_vec = A((Nx+1):(Nx+Ny));

grid_11 = zeros(Nx,Ny);
grid_22 = zeros(Nx,Ny);
grid_12 = zeros(Nx,Ny);
grid_1 = zeros(Nx,Ny);
grid_2 = zeros(Nx,Ny);
theta_p = zeros(Nx,Ny);
traj_1 = zeros(lvs_traj,Ntraj);
traj_2 = zeros(lvs_traj,Ntraj);
start = Nx+Ny+1;
for ii = 1:Nx
    stop = start + Nx -1;
    grid_11(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_22(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_12(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_1(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_2(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    theta_p(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:lvs_traj
    stop = start + lvs_traj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_1(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:lvs_traj
    stop = start + lvs_traj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_2(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
disp('Completed')

if if_plot == 1
    % Ploting
    disp(' ')
    disp('Plotting:')
    lvs = 30;

    disp('figure (1/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_1,lvs,'blue');
    legend('\sigma_{1}')
    xlabel('x-direction')
    ylabel('y-direction')

    disp('figure (2/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_2,lvs,'red');
    legend('\sigma_{2}')
    xlabel('x-direction')
    ylabel('y-direction')

    disp('figure (3/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_11,lvs,'blue');
    legend('\sigma_{11}')
    xlabel('x-direction')
    ylabel('y-direction')

    disp('figure (4/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_22,lvs,'red');
    legend('\sigma_{22}')
    xlabel('x-direction')
    ylabel('y-direction')

    disp('figure (5/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_12,lvs,'blue');
    legend('\sigma_{12}')
    xlabel('x-direction')
    ylabel('y-direction')

    disp('figure (6/7)')
    figure
    hold on
    contour(x_vec, y_vec, theta_p,lvs,'red');
    legend('\theta_{p}')
    xlabel('x-direction')
    ylabel('y-direction')
    
    disp('figure (7/7)')
    figure
    hold on
    for ii = 1:lvs_traj
        p1 = plot(traj_1(ii,:),'blue');
        p2 = plot(traj_2(ii,:),'red');
    end
    legend([p1 p2], '\sigma_{1}','\sigma_{2}')
    xlabel('x-direction')
    ylabel('y-direction')

    disp('Completed')
else
    disp('Plotting has been disabled')
end
end