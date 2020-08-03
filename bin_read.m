%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   BIN Read                                              %
%                   This function reades the bin-file form the            %
%                   C++ plot program.                                     %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

disp('Load the data from the bin-file')
data_file = fopen('data.bin', 'r');
dim_file = fopen('dim_data.bin', 'r');
[A,count_data] = fread(data_file,'double');
[B,count_dim] = fread(dim_file,'double');
disp('Completed')
disp(' ')
disp('Getting the vecotrs')
Nx = B(1);
Ny = B(2);
x_vec = A(1:Nx);
y_vec = A((Nx+1):(Nx+Ny));

disp(' ')
disp('Getting the grids')
grid_11 = zeros(Nx,Ny);
grid_22 = zeros(Nx,Ny);
grid_12 = zeros(Nx,Ny);
grid_1 = zeros(Nx,Ny);
grid_2 = zeros(Nx,Ny);
theta_p = zeros(Nx,Ny);
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

% Ploting
lvs = 30;

figure
hold on
contour(x_vec, y_vec, grid_1,lvs,'blue');
legend('\sigma_{1}')
xlabel('x-direction')
ylabel('y-direction')

figure
hold on
contour(x_vec, y_vec, grid_2,lvs,'red');
legend('\sigma_{2}')
xlabel('x-direction')
ylabel('y-direction')
