function [Nx, Ny, Ntraj, lvs_traj, x_vec, y_vec, grid_11, grid_22, grid_12, grid_1, grid_2, theta_p, traj_1, traj_2] = bin_read()
% BIN_READ Reads the binary files form the C++ computations
%   This program read the output data from the program
%   stress_field_calculator.cpp. It then prodcues the plots for the 
%   different stress fields.
if_plot = 1; % set to 1 for plots to run
close all

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
Nw = B(3);
Ntraj = B(4);
lvs_traj = B(5);
nc = B(6);
nt = B(7);
z1 = zeros(1,nc);
z2 = zeros(1,nc);
for ii = 1:nc
    re = 7 + ii;
    im = 7 + nc + ii;
    z1(ii) = complex(B(re),B(im));
    re = 7 + nc*2 + ii;
    im = 7 + nc*3 + ii;
    z2(ii) = complex(B(re),B(im));
end
z0 = zeros(1,nt);
R = [];
for ii = 1:nt
    re = 7 + nc*4 + ii;
    im = 7 + nc*4 + nt + ii;
    z0(ii) = complex(B(re),B(im));
    pos = 7 + nc*4 + nt*2 + ii;
    R = [R,B(pos)];
end


x_vec = A(1:Nx);
y_vec = A((Nx+1):(Nx+Ny));
x_vecw = A((Nx+Ny+1):(Nx+Ny+Nw));
y_vecw = A((Nx+Ny+Nw+1):(Nx+Ny+Nw+Nw));
grid_11 = zeros(Nx,Ny);
grid_22 = zeros(Nx,Ny);
grid_12 = zeros(Nx,Ny);
grid_1 = zeros(Nx,Ny);
grid_2 = zeros(Nx,Ny);
theta_p = zeros(Nx,Ny);
traj_1 = zeros(lvs_traj*2,Ntraj);
traj_2 = zeros(lvs_traj*2,Ntraj);
grid_w = zeros(Nw,Nw);
traj_w = zeros(Nw,Ntraj);
start = Nx+Ny+Nw+Nw+1;
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
for ii = 1:lvs_traj*2
    stop = start + Ntraj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_1(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:lvs_traj*2
    stop = start + Ntraj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_2(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:Nw
    stop = start + Nw*2 -1;
    for jj = 1:Nw
        re = start + jj - 1;
        im = start + Nw + jj - 1;
        grid_w(jj,ii) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:Nw*2
    stop = start + Ntraj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_w(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
disp('Completed')

lvs = 30;
if if_plot == 1
    % Ploting
    disp(' ')
    disp('Plotting:')
    disp('Contour levels set to:')
    disp(lvs)

    disp('figure (1/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_1,lvs,'blue');
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend('\sigma_{1}')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])

    disp('figure (2/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_2,lvs,'red');
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend('\sigma_{2}')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])

    disp('figure (3/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_11,lvs,'blue');
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend('\sigma_{11}')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])

    disp('figure (4/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_22,lvs,'red');
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend('\sigma_{22}')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])

    disp('figure (5/7)')
    figure
    hold on
    contour(x_vec, y_vec, grid_12,lvs,'green');
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend('\sigma_{12}')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])

    disp('figure (6/7)')
    figure
    hold on
    contour(x_vec, y_vec, theta_p,lvs,'red');
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend('\theta_{p}')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
    
    disp('figure (7/7)')
    figure
    hold on
    for ii = 1:lvs_traj*2
        p1 = plot(real(traj_1(ii,:)),imag(traj_1(ii,:)),'blue');
        p2 = plot(real(traj_2(ii,:)),imag(traj_2(ii,:)),'red');
    end
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend([p1 p2], '\sigma_{1}','\sigma_{2}')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
    
    disp('figure (8/7)')
    figure
    hold on
    quiver(x_vecw, y_vecw, real(grid_w), -imag(grid_w),'blue');
    for ii = 1:Nw*2
        p1 = plot(real(traj_w(ii,:)),imag(traj_w(ii,:)),'blue');
    end
    for ii = 1:nc
        Plot_line(z1(ii),z2(ii),'black')
    end
    for ii = 1:nt
        plot_circle(z0(ii),R(ii),'black')
    end
    legend('w')
    xlabel('x-direction')
    ylabel('y-direction')
    axis([x_vecw(1) x_vecw(end) y_vecw(1) y_vecw(end)])

    disp('Completed')
else
    disp('Plotting has been disabled')
end
end