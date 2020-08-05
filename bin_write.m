function [] = bin_write(A, B)
% BIN_WRITE Wrties the binary files for C++ computations
%   This program provides the program stress_field_calculator.cpp with the 
%   necessary bin-files. It can either be called with the data collected in 
%   two vectors or it can be run by itself.
%
%   The length of the vectors are coded in and needs to be updated when
%   more elements are added. 
switch nargin
    case 2
        if length(A)== 6 && length(B) == 8
            % The A and B vectyors are given and the input may be skipped.
            % Write the bin-files for C++ program
            input_file = fopen('input_data.bin','w');
            fwrite(input_file, A, 'double');
            fclose(input_file);

            plot_file = fopen('plot_data.bin','w');
            fwrite(plot_file,  B, 'double');
            fclose(plot_file);
        
            disp('The output files has been written with user input values.')
        else
            disp('Dimension of A and B vectors not correct')
        end
    case 1
        disp('Only one argument was given, please give two och none')
    otherwise
        % Assigning variables
        H = complex(0,0);
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
        Ntraj = 50;
        lvs_traj = 50;

        % Write the bin-files for C++ program
        A = [real(H),imag(H),rho,g,nu,kappa]; % Vector to write
        input_file = fopen('input_data.bin','w');
        fwrite(input_file, A, 'double');
        fclose(input_file);

        B = [xfrom,xto,yfrom,yto,Nx,Ny,Ntraj,lvs_traj]; % Vector to write
        plot_file = fopen('plot_data.bin','w');
        fwrite(plot_file, B, 'double');
        fclose(plot_file);

        disp('The output files has been written with default values.')
end
end