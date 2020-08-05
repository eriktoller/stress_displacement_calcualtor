#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <tuple>
#include <fstream>
#include <ctime>
#include <chrono>

// Defining special types
typedef std::complex<double> dcomp;
typedef std::vector<double> dvec;
typedef std::vector< std::vector<double> > ddvec;
typedef std::vector< std::complex<double> > dcvec;
typedef std::vector< std::vector<std::complex<double> > > ddcvec;

/* --------------------------------------------------------------------

		FUNCTIONS

-----------------------------------------------------------------------*/

/* --------------------------------------------------------------------
		TRACTION GRAVITY
-----------------------------------------------------------------------*/
dcomp T_gravity(dcomp z, int H, double rho, double g, double nu, double kappa)
{
	// Defining the variables
	dcomp T;
	dcomp tau_11;
	dcomp tau_12;
	dcomp term1;
	dcomp frac1;
	dcomp frac2;
	dcomp term2;

	// Calculating the terms
	term1 = dcomp(0,2)*rho*g;
	frac1 = (1 - 2 * nu) / (kappa + 1);
	frac2 = 1 / (kappa + 1);
	term2 = (z - conj(z) - 4.0*pow(H,2));
	
	// Calculating the tau_11, tau_12 and T
	tau_11 = -term1 * frac1 * term2;
	tau_12 = term1 * frac2 * term2;
	T = -0.5*dcomp(0,1)*(tau_11 - tau_12);
	
	return T;
}

/* --------------------------------------------------------------------
		STRESS GRAVITY
-----------------------------------------------------------------------*/
std::tuple<double, double, double> sigma_gravity(dcomp z, dcomp H, double rho, double g, double nu, double kappa)
{
	// Defining the variables
	dcomp sigma_a;
	dcomp sigma_b;
	dcomp tau_11;
	dcomp tau_12;
	dcomp term1;
	dcomp frac1;
	dcomp frac2;
	dcomp term2;
	double sigma_11;
	double sigma_22;
	double sigma_12;

	// Calculating the terms
	term1 = dcomp(0, 2)*rho*g;
	frac1 = (1 - 2 * nu) / (kappa + 1);
	frac2 = 1 / (kappa + 1);
	term2 = (z - conj(z)) - 2.0*H;

	// Calculating the tau_11, tau_12 and sigma expressions
	tau_11 = -term1 * frac1 * term2;
	tau_12 = term1 * frac2 * term2;
	sigma_a = 0.5*(tau_11 + tau_12);
	sigma_b = 0.5*(tau_12 - tau_11);

	// Calculating the sigma_11, sigma_22 and sigma_12
	sigma_11 = real(sigma_a);
	sigma_22 = real(sigma_b);
	sigma_12 = imag(sigma_a);

	return {sigma_11, sigma_22, sigma_12 };
}

/* --------------------------------------------------------------------
		STRESS FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddvec, ddvec, ddvec> stress_field(double xfrom, double xto, double yfrom, double yto, int Nx, int Ny, dcomp H, double rho, double g, double nu, double kappa)
{
	// Defining the variables
	double sigma_11;
	double sigma_22;
	double sigma_12;
	ddvec grid_11(Nx, dvec(Ny));
	ddvec grid_22(Nx, dvec(Ny));
	ddvec grid_12(Nx, dvec(Ny));
	double dx;
	double dy;
	dvec x_vec(Nx);
	dvec y_vec(Ny);
	dcomp z_grid;
	int time_Nx;

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nx-1);
	dy = (yto - yfrom) / (Ny-1);
	std::cout << "Progress of the stress field calculation\n";
	for (int ii = 0; ii < Nx; ii++)
	{
		for (int jj = 0; jj < Ny; jj++)
		{
			x_vec[ii] = xfrom + ii * dx;
			y_vec[jj] = yfrom + jj * dy;
			z_grid = dcomp(x_vec[ii], y_vec[jj]);
			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z_grid, H, rho, g, nu, kappa);
			grid_11[ii][jj] = sigma_11;
			grid_22[ii][jj] = sigma_22;
			grid_12[ii][jj] = sigma_12;
		}
		time_Nx = ii*100 / (Nx - 1);
		std::cout << "\r " << time_Nx << "%";
	}
	std::cout << "\nComplete\n";
	return {x_vec, y_vec, grid_11, grid_22, grid_12};
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESSES
-----------------------------------------------------------------------*/
std::tuple<double, double, double> principal_sigma(double sigma_11, double sigma_22, double sigma_12)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	double frac1, frac2, sqrt1;

	// Calculating the terms
	frac1 = (sigma_11 + sigma_22) / 2.0;
	frac2 = ((sigma_11 - sigma_22) / 2.0);
	sqrt1 = sqrt( pow( frac2 ,2) + pow(sigma_12,2) );

	// Calculating the principal stresses and the angel of sigma_1
	sigma_1 = frac1 + sqrt1;
	sigma_2 = frac1 - sqrt1;

	// Calcuating the angel theta_p
	if (sigma_11 == sigma_22)
	{
		theta_p = 0;
	}
	else
	{
		theta_p = 0.5*atan((2 * sigma_12) / (sigma_11 - sigma_22));
	}

	return { sigma_1, sigma_2, theta_p };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESS FIELDS
-----------------------------------------------------------------------*/
std::tuple<ddvec, ddvec, ddvec> principal_stress_field(ddvec grid_11, ddvec grid_22, ddvec grid_12)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	ddvec grid_1(grid_11.size(), dvec(grid_11[0].size()));
	ddvec grid_2(grid_11.size(), dvec(grid_11[0].size()));
	ddvec grid_theta(grid_11.size(), dvec(grid_11[0].size()));
	int time_Nx;

	std::cout << "Progress of the principal stress field calculation\n";
	// Retriving the terms from the sigma function for z
	for (size_t ii = 0; ii < grid_11.size(); ii++)
	{
		for (size_t jj = 0; jj < grid_11[0].size(); jj++)
		{
			std::tie(sigma_1, sigma_2, theta_p) = principal_sigma(grid_11[ii][jj], grid_22[ii][jj], grid_12[ii][jj]);
			grid_1[ii][jj] = sigma_1;
			grid_2[ii][jj] = sigma_2;
			grid_theta[ii][jj] = theta_p;
		}
		time_Nx = ii * 100 / grid_11.size() + 1;
		std::cout << "\r " << time_Nx << "%";
	}

	std::cout << "\nComplete\n";
	return { grid_1, grid_2, grid_theta };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESS TRAJECTORIES
-----------------------------------------------------------------------*/
std::tuple<ddcvec, ddcvec> principal_stress_trajectories(double xfrom, double xto, double yfrom, double yto, int Ntraj, int lvs_traj, dcomp H, double rho, double g, double nu, double kappa)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	double sigma_11, sigma_22, sigma_12;
	ddcvec traj_1(lvs_traj, dcvec(Ntraj));
	ddcvec traj_2(lvs_traj, dcvec(Ntraj));
	int time_lvs_traj;
	double dx;
	double dy;
	dcomp z, z1, z2, z3, z4, sigma, sigma1, sigma2, sigma3;

	std::cout << "Progress of the principal stress trajectories calculation\n";
	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj - 1);
	dy = (yto - yfrom) / (Ntraj - 1);
	for (int ii = 0; ii < lvs_traj; ii++)
	{
		z = dcomp(xfrom, yfrom + ii * dy);
		traj_1[ii][0] = z;

		for (int jj = 1; jj < Ntraj; jj++)
		{
			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z, H, rho, g, nu, kappa);
			std::tie(sigma_1, sigma_2, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
			sigma = abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z1 = z + (sigma) / abs(sigma)*dx;

			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z1, H, rho, g, nu, kappa);
			sigma1 = abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z2 = z + ((sigma)+(sigma1)) / abs(sigma + sigma1)*dx;

			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z2, H, rho, g, nu, kappa);
			sigma2 = abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z3 = z + ((sigma)+(sigma2)) / abs(sigma + sigma2)*dx;

			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z3, H, rho, g, nu, kappa);
			sigma3 = abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z4 = z + ((sigma)+(sigma3)) / abs(sigma + sigma3)*dx;

			traj_1[ii][jj] = z4;
			z = z4;
		}
		z = dcomp(xfrom + ii * dx, yfrom);
		traj_2[ii][0] = z;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z, H, rho, g, nu, kappa);
			std::tie(sigma_1, sigma_2, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
			sigma = - abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z1 = z + (sigma) / abs(sigma)*dy;

			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z1, H, rho, g, nu, kappa);
			sigma1 = - abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z2 = z + ((sigma)+(sigma1)) / abs(sigma + sigma1)*dy;

			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z2, H, rho, g, nu, kappa);
			sigma2 = - abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z3 = z + ((sigma)+(sigma2)) / abs(sigma + sigma2)*dy;

			std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z3, H, rho, g, nu, kappa);
			sigma3 = - abs(sigma_1)*exp(dcomp(0, -1)*((theta_p)));
			z4 = z + ((sigma)+(sigma3)) / abs(sigma + sigma3)*dy;

			traj_2[ii][jj] = z4;
			z = z4;
		}
		time_lvs_traj = ii * 100 / (lvs_traj-1);
		std::cout << "\r " << time_lvs_traj << "%";
	}

	std::cout << "\nComplete\n";
	return { traj_1, traj_2 };
}

/* --------------------------------------------------------------------

		MAIN SCRIPT

-----------------------------------------------------------------------*/
int main()
{
	/* --------------------------------------------------------------------
			Defining variables and imprting data
	-----------------------------------------------------------------------*/

	auto start = std::chrono::system_clock::now(); // Start the clock

	std::cout << "The documentation for this program can be found on: https://github.com/eriktoller/stress_field_calcualtor \n";
	std::cout << "Written by: Erik Toller, erik.toller@geo.uu.se.\n\n";

	// Setting the data types
	dcomp z;
	dcomp T;
	dcomp H;
	double rho, g, nu, kappa;
	double xfrom, xto, yfrom, yto;
	int Nx, Ny, Ntraj, lvs_traj;
	ddvec grid_11, grid_22, grid_12;
	ddvec grid_1, grid_2, theta_p;
	ddcvec traj_1, traj_2;
	dvec x_vec, y_vec;

	/*
	// Assigning varibales
	H = dcomp(0,0);
	rho = 2750.0;
	g = 9.816;
	nu = 0.3;
	kappa = 3 - 4 * nu;
	

	// Plot dim and resolutioin
	xfrom = -1;
	xto = 1;
	yfrom = -1;
	yto = 1;
	Nx = 201;
	Ny = 201;
	*/

	// Read the input data from binary file
	double fin[6];
	std::ifstream input_file("input_data.bin", std::ios::in | std::ios::binary);
	input_file.read((char *)&fin, sizeof fin);
	H = dcomp(fin[0], fin[1]);
	rho = fin[2];
	g = fin[3];
	nu = fin[4];
	kappa = fin[5];

	// Read the plot data from binary file
	double fplot[8];
	std::ifstream plot_file("plot_data.bin", std::ios::in | std::ios::binary);
	plot_file.read((char *)&fplot, sizeof fplot);
	xfrom = fplot[0];
	xto = fplot[1];
	yfrom = fplot[2];
	yto = fplot[3];
	Nx = (int) fplot[4];
	Ny = (int) fplot[5];
	Ntraj = (int) fplot[6];
	lvs_traj = (int) fplot[7];


	// Displying the plot data in the command window
	std::cout << "This is the retrived plot data:\n";
	std::cout << "x from: " << xfrom << " to " << xto << std::endl;
	std::cout << "y from: " << yfrom << " to " << yto << std::endl;
	std::cout << "x resolution: " << Nx << std::endl;
	std::cout << "y resolution: " << Ny << std::endl;
	std::cout << "Total number of points: " << Nx * Ny << std::endl;
	std::cout << "Number of trajectory levels: " << lvs_traj << std::endl;
	std::cout << "Number of steps in trajectories: " << Ntraj;
	std::cout <<std::endl << std::endl;

	/* --------------------------------------------------------------------
			Calcuating the stress field
	-----------------------------------------------------------------------*/

	// Get the Cartisian and principal stress field
	std::tie(x_vec, y_vec, grid_11, grid_22, grid_12) = stress_field(xfrom, xto, yfrom, yto, Nx, Ny, H, rho, g, nu, kappa);
	std::tie(grid_1, grid_2, theta_p) = principal_stress_field(grid_11, grid_22, grid_12);
	std::tie(traj_1, traj_2) = principal_stress_trajectories(xfrom, xto, yfrom, yto, Ntraj, lvs_traj, H, rho, g, nu, kappa);

	std::cout << std::endl << std::endl;

	/* -------------------------------------------------------------------- 
			Saving the data as binary files
	-----------------------------------------------------------------------*/

	std::cout << "Saving the output data \n";

	// Create/open the two output files
	std::ofstream outfile("data.bin", std::ios::out | std::ios::binary);
	std::ofstream outfiledim("dim_data.bin", std::ios::out | std::ios::binary);

	// Save the x and y vectors
	dvec fx = x_vec;
	const char* pointerx = reinterpret_cast<const char*>(&fx[0]);
	std::size_t bytesx = fx.size() * sizeof(fx[0]);
	outfile.write(pointerx, bytesx);
	
	dvec fy = y_vec;
	const char* pointery = reinterpret_cast<const char*>(&fy[0]);
	std::size_t bytesy = fy.size() * sizeof(fy[0]);
	outfile.write(pointery, bytesy);

	// Save the grids
	for (size_t ii = 0; ii < grid_11.size(); ii++)
		{
			dvec fg11 = grid_11[ii];
			const char* pointerg11 = reinterpret_cast<const char*>(&fg11[0]);
			std::size_t bytesg11 = fg11.size() * sizeof(fg11[0]);
			outfile.write(pointerg11, bytesg11);
		}
	for (size_t ii = 0; ii < grid_22.size(); ii++)
		{
			dvec fg22 = grid_22[ii];
			const char* pointerg22 = reinterpret_cast<const char*>(&fg22[0]);
			std::size_t bytesg22 = fg22.size() * sizeof(fg22[0]);
			outfile.write(pointerg22, bytesg22);
		}
	for (size_t ii = 0; ii < grid_12.size(); ii++)
		{
			dvec fg12 = grid_12[ii];
			const char* pointerg12 = reinterpret_cast<const char*>(&fg12[0]);
			std::size_t bytesg12 = fg12.size() * sizeof(fg12[0]);
			outfile.write(pointerg12, bytesg12);
		}
	for (size_t ii = 0; ii < grid_1.size(); ii++)
		{
			dvec fg1 = grid_1[ii];
			const char* pointerg1 = reinterpret_cast<const char*>(&fg1[0]);
			std::size_t bytesg1 = fg1.size() * sizeof(fg1[0]);
			outfile.write(pointerg1, bytesg1);
		}
	for (size_t ii = 0; ii < grid_2.size(); ii++)
		{
			dvec fg2 = grid_2[ii];
			const char* pointerg2 = reinterpret_cast<const char*>(&fg2[0]);
			std::size_t bytesg2 = fg2.size() * sizeof(fg2[0]);
			outfile.write(pointerg2, bytesg2);
		}
	for (size_t ii = 0; ii < theta_p.size(); ii++)
	{
		dvec fgtp = theta_p[ii];
		const char* pointergtp = reinterpret_cast<const char*>(&fgtp[0]);
		std::size_t bytesgtp = fgtp.size() * sizeof(fgtp[0]);
		outfile.write(pointergtp, bytesgtp);
	}
	for (size_t ii = 0; ii < traj_1.size(); ii++)
	{
		dvec fgt1_re(Ntraj);
		dvec fgt1_im(Ntraj);
		
		for (size_t jj = 0; jj < traj_1[0].size(); jj++)
		{
			fgt1_re[jj] = real(traj_1[ii][jj]);
			fgt1_im[jj] = imag(traj_1[ii][jj]);
		}
		const char* pointergt1_re = reinterpret_cast<const char*>(&fgt1_re[0]);
		std::size_t bytesgt1_re = fgt1_re.size() * sizeof(fgt1_re[0]);
		outfile.write(pointergt1_re, bytesgt1_re);
		const char* pointergt1_im = reinterpret_cast<const char*>(&fgt1_im[0]);
		std::size_t bytesgt1_im = fgt1_im.size() * sizeof(fgt1_im[0]);
		outfile.write(pointergt1_im, bytesgt1_im);
	}
	for (size_t ii = 0; ii < traj_2.size(); ii++)
	{
		dvec fgt2_re(Ntraj);
		dvec fgt2_im(Ntraj);
		for (size_t jj = 0; jj < traj_2[0].size(); jj++)
		{
			fgt2_re[jj] = real(traj_2[ii][jj]);
			fgt2_im[jj] = imag(traj_2[ii][jj]);
		}
		const char* pointergt2_re = reinterpret_cast<const char*>(&fgt2_re[0]);
		std::size_t bytesgt2_re = fgt2_re.size() * sizeof(fgt2_re[0]);
		outfile.write(pointergt2_re, bytesgt2_re);
		const char* pointergt2_im = reinterpret_cast<const char*>(&fgt2_im[0]);
		std::size_t bytesgt2_im = fgt2_im.size() * sizeof(fgt2_im[0]);
		outfile.write(pointergt2_im, bytesgt2_im);
	}

	// Save the dimensions of Nx and Ny vectors
	dvec dim = { 1.0 * Nx, 1.0 * Ny, 1.0*Ntraj, 1.0*lvs_traj };
	const char* pointerdim = reinterpret_cast<const char*>(&dim[0]);
	std::size_t bytesdim = dim.size() * sizeof(dim[0]);
	outfiledim.write(pointerdim, bytesdim);

	// Close the output files
	outfile.close();
	outfiledim.close();

	std::cout << "Saving completed" << std::endl << std::endl;

	// Create the log file
	auto end = std::chrono::system_clock::now(); // Get the stop time
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	//const char * dt = std::ctime(&end_time);
	char str_time[26];
	ctime_s(str_time, sizeof str_time, &end_time);
	

	std::ofstream logfile;
	logfile.open("log.txt");
	logfile << "This is the log file for computation completed on: " << str_time << std::endl;
	logfile << "The computation took: " << elapsed_seconds.count() << "s";
	logfile << std::endl << std::endl;
	logfile << "These are the values of the constants:\n";
	logfile << "    H = (" << H.real() << "," << H.imag() << ")" << std::endl;
	logfile << "  rho = " << rho << std::endl;
	logfile << "    g = " << g << std::endl;
	logfile << "   nu = " << nu << std::endl;
	logfile << "kappa = " << kappa << std::endl;
	logfile << std::endl << std::endl;
	logfile << "This is the used plot data:\n";
	logfile << "x from: " << xfrom << " to " << xto << std::endl;
	logfile << "y from: " << yfrom << " to " << yto << std::endl;
	logfile << "x resolution: " << Nx << std::endl;
	logfile << "y resolution: " << Ny << std::endl;
	logfile << "Total number of points: " << Nx * Ny << std::endl;
	logfile << "Number of trajectory levels: " << lvs_traj << std::endl;
	logfile << "Number of steps in trajectories: " << Ntraj;
	logfile << std::endl << std::endl;
	logfile << "The data has been saved scussefully to:\n";
	logfile << "'data.bin'     [x_vec, y_vec, grid_11, grid_22, grid_12, grid_1, grid_2, theta_p, (traj_1.real, traj_1.imag), (traj_2.real, traj_2.imag)]\n";
	logfile << "'dim_data.bin' [Nx, Ny, Ntraj, lvs_traj]\n";
	logfile.close();

	std::cout << "Program finnished after " << elapsed_seconds.count() << " s and output data saved to binary files\n";

	return 0;
}