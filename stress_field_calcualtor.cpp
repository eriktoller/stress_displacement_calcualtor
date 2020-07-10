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

// FUNCTIONS
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

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nx-1);
	dy = (yto - yfrom) / (Ny-1);
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
	}

	return {x_vec, y_vec, grid_11, grid_22, grid_12};
}

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
	theta_p = 0.5*atan((2 * sigma_12) / (sigma_11 - sigma_22));

	return { sigma_1, sigma_2, theta_p };
}

std::tuple<ddvec, ddvec, ddvec> principal_stress_field(ddvec grid_11, ddvec grid_22, ddvec grid_12)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	ddvec grid_1(grid_11.size(), dvec(grid_11[0].size()));
	ddvec grid_2(grid_11.size(), dvec(grid_11[0].size()));
	ddvec grid_theta(grid_11.size(), dvec(grid_11[0].size()));

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
	}

	return { grid_1, grid_2, grid_theta };
}

// MAIN
int main()
{
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
	Nx = 11;
	Ny = 11;

	std::tie(x_vec, y_vec, grid_11, grid_22, grid_12) = stress_field(xfrom, xto, yfrom, yto, Nx, Ny, H, rho, g, nu, kappa);
	std::tie(grid_1, grid_2, theta_p) = principal_stress_field(grid_11, grid_22, grid_12);

	std::cout << "The chartisian stress sigma_22\n";
	for (size_t i = 0; i < x_vec.size(); i++)
	{
		for (size_t j = 0; j < y_vec.size(); j++)
		{
			std::cout << "[" << x_vec[i] << "]" << "[" << y_vec[j] << "]= " << grid_22[i][j] << std::endl;
		}

	}

	std::cout << "The principal stress sigma_2\n";
	for (size_t i = 0; i < x_vec.size(); i++)
	{
		for (size_t j = 0; j < y_vec.size(); j++)
		{
			std::cout << "[" << x_vec[i] << "]" << "[" << y_vec[j] << "]= " << grid_2[i][j] << std::endl;
		}

	}

	std::ofstream myfile;
	myfile.open("data.txt");
	for (size_t ii = 0; ii < grid_11.size(); ii++)
	{
		for (size_t jj = 0; jj < grid_11[0].size(); jj++)
		{
			myfile << " " << grid_1[ii][jj];
		}
	}
	myfile.close();

	return 0;
}