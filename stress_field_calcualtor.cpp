#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <tuple>
#include <fstream>

// Defining special types
typedef std::complex<double> dcomp;
typedef std::vector<double> dvec;
typedef std::vector<std::vector<double>> ddvec;

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
		time_Nx = ii*100 / Nx + 1;
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

		MAIN SCRIPT

-----------------------------------------------------------------------*/
int main()
{
	/* --------------------------------------------------------------------
			Defining variables and imprting data
	-----------------------------------------------------------------------*/

	std::cout << "This is a program under construction and is not fully developed yet.\n";
	std::cout << "The first item that will be included is the analytic element for gravity (since its published).\n\n";
	std::cout << "Reference:\nStrack, O. D. (2020). Applications of Vector Analysis and Complex Variables in Engineering. Springer International Publishing.\n\n";
	// Setting the data types
	dcomp z;
	dcomp T;
	dcomp H;
	double rho, g, nu, kappa;
	double xfrom, xto, yfrom, yto;
	int Nx, Ny;
	ddvec grid_11, grid_22, grid_12;
	ddvec grid_1, grid_2, theta_p;
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

	// Read the input and plot data
	double fin[6];
	std::ifstream input_file("input_data.bin", std::ios::in | std::ios::binary);
	input_file.read((char *)&fin, sizeof fin);
	H = dcomp(fin[0], fin[1]);
	rho = fin[2];
	g = fin[3];
	nu = fin[4];
	kappa = fin[5];

	double fplot[6];
	std::ifstream plot_file("plot_data.bin", std::ios::in | std::ios::binary);
	plot_file.read((char *)&fplot, sizeof fplot);
	xfrom = fplot[0];
	xto = fplot[1];
	yfrom = fplot[2];
	yto = fplot[3];
	Nx = (int) fplot[4];
	Ny = (int) fplot[5];

	// Displying the plot data
	std::cout << "This is the preset plot data:\n";
	std::cout << "x from: " << xfrom << " to " << xto << std::endl;
	std::cout << "y from: " << yfrom << " to " << yto << std::endl;
	std::cout << "x resolution: " << Nx << " and y resolution: " << Ny << std::endl << std::endl;

	/* --------------------------------------------------------------------
			Calcuating the stress field
	-----------------------------------------------------------------------*/

	// Get the Cartisian and principal stress field
	std::tie(x_vec, y_vec, grid_11, grid_22, grid_12) = stress_field(xfrom, xto, yfrom, yto, Nx, Ny, H, rho, g, nu, kappa);
	std::tie(grid_1, grid_2, theta_p) = principal_stress_field(grid_11, grid_22, grid_12);

	// Display the data (Example purposes)
	std::cout << "The stress feilds has been calculated \n";

	/* -------------------------------------------------------------------- 
			Wirte a binary files from data
	-----------------------------------------------------------------------*/

	std::ofstream outfile("data.bin", std::ios::out | std::ios::binary);
	std::ofstream outfiledim("dim_data.bin", std::ios::out | std::ios::binary);

	// The x and y vectors
	dvec fx = x_vec;
	const char* pointerx = reinterpret_cast<const char*>(&fx[0]);
	std::size_t bytesx = fx.size() * sizeof(fx[0]);
	outfile.write(pointerx, bytesx);

	dvec fy = y_vec;
	const char* pointery = reinterpret_cast<const char*>(&fy[0]);
	std::size_t bytesy = fy.size() * sizeof(fy[0]);
	outfile.write(pointery, bytesy);

	// The grids
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


	outfile.close();

	// The dimensions of Nx and Ny vectors
	dvec dim = { 1.0 * Nx, 1.0 * Ny };
	const char* pointerdim = reinterpret_cast<const char*>(&dim[0]);
	std::size_t bytesdim = dim.size() * sizeof(dim[0]);
	outfiledim.write(pointerdim, bytesdim);
	outfiledim.close();

	return 0;
}