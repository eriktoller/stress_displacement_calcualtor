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
        % The A and B vectyors are given and the input may be skipped.
        % Write the bin-files for C++ program
        input_file = fopen('input_data.bin','w');
        fwrite(input_file, A, 'double');
        fclose(input_file);

        plot_file = fopen('plot_data.bin','w');
        fwrite(plot_file,  B, 'double');
        fclose(plot_file);

        disp('The output files has been written with user input values.')
    case 1
        disp('Only one argument was given, please give two och none')
    otherwise
        % Assigning variables
        H = complex(0,0);
        rho = 2750;
        g = 9.816*0;
        nu = 0.3;
        kappa = 3 - 4 * nu;
        
        nc = 1;
        m = 4;
        z1 = complex(-.5,0);
        z2 = complex(.5,0);
        L = 1;
        mu = 0;
        beta = [5,0,0,0];
        
        nt = 0;
        mt = 0;
        R = 1;
        z0 = complex(0,0);
        
        beta_vec = [];
        for ii = 1:nc
           beta_vec = [beta_vec, real(beta(ii,:)), imag(beta(ii,:))];
        end
        a_vec = [];
        for ii = 1:mt
           a_vec = [a_vec, real(a(ii,:)), imag(a(ii,:))];
        end
        b_vec = [];
        for ii = 1:mt
           b_vec = [b_vec, real(b(ii,:)), imag(b(ii,:))];
        end
        
        % Plot dim and resolution
        xfrom = -1;
        xto = 1;
        yfrom = -1;
        yto = 1;
        Nx = 400;
        Ny = 400;
        Ntraj = 800;
        lvs_traj = 60;
        
        % Check if the data is correct
        disp('Checking the data dimensions')
        cnt = 0;
        if m == length(beta)
            cnt = cnt + 1;
        else
            disp(' m /= length(beta)')
            disp('m = ')
            disp(m)
            disp('length(beta) =')
            disp(length(beta))
        end
        if nc == length(z1)
            cnt = cnt + 1;
        else
            disp('nc /= length(z1)')
            disp('nc = ')
            disp(nc)
            disp('length(z1) =')
            disp(length(z1))
        end
        if nc == length(L)
            cnt = cnt + 1;
        else
            disp('nc /= length(L)')
            disp('nc = ')
            disp(nc)
            disp('length(L) =')
            disp(length(L))
        end
        if nc == length(mu)
            cnt = cnt + 1;
        else
            disp('nc /= length(mu)')
            disp('nc = ')
            disp(nc)
            disp('length(mu) =')
            disp(length(mu))
        end
        if cnt == 4
            disp('Data check OK')
            
            % Write the bin-files for C++ program
            A = [real(H),imag(H),rho,g,nu,kappa,nc,m,nt,mt,real(z1),imag(z1),real(z2),imag(z2),L,mu,beta_vec,real(z0),imag(z0),R,a_vec,b_vec]; % Vector to write
            input_file = fopen('input_data.bin','w');
            fwrite(input_file, A, 'double');
            fclose(input_file);

            B = [xfrom,xto,yfrom,yto,Nx,Ny,Ntraj,lvs_traj]; % Vector to write
            plot_file = fopen('plot_data.bin','w');
            fwrite(plot_file, B, 'double');
            fclose(plot_file);

            disp('The output files has been written with default values.')
        else
            disp('ERROR in data dimensions')
            disp('bin-files not created')
        end

        
end
end