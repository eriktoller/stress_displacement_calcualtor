#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <tuple>
#include <fstream>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <stdio.h>

// Defining special types
typedef std::complex<double> dcomp;
typedef std::vector<double> dvec;
typedef std::vector< std::vector<double> > ddvec;
typedef std::vector< std::complex<double> > dcvec;
typedef std::vector< std::vector<std::complex<double> > > ddcvec;

/* --------------------------------------------------------------------

		FUNCTIONS

-----------------------------------------------------------------------*/

inline const double pi()
{
	const double pi = std::atan(1) * 4;
	return pi;
}

/* --------------------------------------------------------------------
		CHI FROM Z
-----------------------------------------------------------------------*/
inline dcomp chi_from_z(dcomp z, dcomp z1, dcomp z2, double L, double mu)
{
	// Defining the variables
	dcomp z0, Z, chi;

	// Calculating the chi from a z value
	z0 = 0.5*(z1 + z2);
	Z = exp(dcomp(0,-1) * mu) * 2.0 * (z - z0) / L;
	chi = Z + sqrt(Z - 1.0)*sqrt(Z + 1.0);

	return chi;
}

/* --------------------------------------------------------------------
		STRESS GRAVITY
-----------------------------------------------------------------------*/
inline std::tuple<double, double, double> sigma_gravity(dcomp z, dcomp H, double rho, double g, double nu, double kappa)
{
	// Defining the variables
	dcomp sigma_a, sigma_b;
	dcomp tau_11, tau_12;
	dcomp term1, frac1, frac2, term2;
	double sigma_11,sigma_22, sigma_12;

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
		STRESS CRACK
-----------------------------------------------------------------------*/
inline std::tuple<double, double, double> sigma_crack(dcomp z, dcomp z1, dcomp z2, double L, double mu, int m, dcvec beta)
{
	// Defining the variables
	dcomp chi, chi_bar, Z, chi_pow;
	dcomp dphi, dphi_bar, ddphi, dpsi;
	dcomp tau_11, tau_12, S1, L_frac;
	double sigma_11, sigma_22, sigma_12, n;

	// Getting the chi - and Z - coordinates
	chi = chi_from_z(z, z1, z2, L, mu);
	chi_bar = conj(chi);
	Z = exp(dcomp(0,-1)*mu) * 2.0 * (z - 0.5*(z1 + z2)) / L;

	// Calculating the series
	dphi = 0;
	dphi_bar = 0;
	ddphi = 0;
	dpsi = 0;
	n = 0;
	chi_pow = chi*chi - 1.0;
	for (int ii = 0; ii < m; ii++ )
	{
		dcomp beta_n;
		n += 1; 
		beta_n = beta[ii] * n;
		dphi += conj(beta_n) * pow(chi, (1.0 - n)) / chi_pow;
		dphi_bar += beta_n * pow(chi_bar, (1.0 - n)) / (chi_bar * chi_bar - 1.0);
		ddphi += conj(beta_n) *(((1 - n) * pow(chi, (2.0 - n))) / (chi_pow * chi_pow) - (2.0 * pow(chi, (4.0 - n))) / (chi_pow * chi_pow * chi_pow));
		dpsi -= beta_n * (pow(chi, (1.0 - n))) / chi_pow;
	}

	// Multiplying the constants
	L_frac = (4.0 / L) * exp(dcomp(0, -1)*mu);
	dphi *= L_frac;
	dphi_bar *= conj(L_frac);
	ddphi *= (16.0 / (L*L))*exp(dcomp(0,-2)*mu);
	dpsi *= L_frac;

	// Calcualting simga
	tau_11 = -0.5*L*(Z - conj(Z))*ddphi - exp(dcomp(0,-1)*mu)*(dphi + dpsi);
	tau_12 = -exp(dcomp(0,1)*mu)*dphi - exp(dcomp(0,-1)*mu)*dphi_bar;

	S1 = .5*(tau_11 + tau_12);
	sigma_11 = real(S1);
	sigma_22 = real(.5*(-tau_11 + tau_12));
	sigma_12 = -imag(S1);

	return { sigma_11, sigma_22, sigma_12 };
}

/* --------------------------------------------------------------------
		STRESS CIRCULAR TUNNEL
-----------------------------------------------------------------------*/
inline std::tuple<double, double, double> sigma_circ_tunnel(dcomp z, dcomp z0, double R, int mt, dcvec a, dcvec b)
{
	// Defining the variables
	dcomp Z, Z_bar, ZZ;
	dcomp dPhi, dPhi_bar, ddPhi, dPsi;
	dcomp tau_11, tau_12, S1;
	double sigma_11, sigma_22, sigma_12, n;

	Z = (z - z0) / R;
	Z_bar = conj(Z);
	ZZ = Z * Z_bar;

	if (ZZ.real() < 0.999999)
	{
		tau_11 = dcomp(NAN, NAN);
		tau_12 = dcomp(NAN, NAN);
	} 
	else
	{
		dPsi = 0;
		dPhi = 0;
		dPhi_bar = 0;
		ddPhi = 0;
		

		for (int ii = 0; ii < mt; ii++)
		{
			n = ii;
			dcomp c = b[ii] + a[ii];
			dPsi += a[ii]*pow(Z,(-n - 2));
			dPhi += c*pow(Z,-n);
			dPhi_bar += conj(c)*pow(Z_bar,-n);
			ddPhi += n*c*pow(Z,(-n - 1));
		}
		
		dPsi = -dPsi * 2.0;
		ddPhi = -ddPhi / R;
		tau_11 = (R*Z_bar - R / Z)*ddPhi - dPsi;
		tau_12 = -dPhi - dPhi_bar;
	}

	S1 = .5*(tau_11 + tau_12);
	sigma_11 = real(S1);
	sigma_22 = real(.5*(-tau_11 + tau_12));
	sigma_12 = -imag(S1);

	return { sigma_11, sigma_22, sigma_12 };
}

/* --------------------------------------------------------------------
		STRESS TOTAL
-----------------------------------------------------------------------*/
inline std::tuple<double, double, double> sigma_total(dcomp z, dcomp H, double rho, double g, double nu, double kappa, dcvec z1, dcvec z2, dvec L, dvec mu, int m, ddcvec beta, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	dcomp sigma_a, sigma_b;
	dcomp tau_11, tau_12;
	dcomp term1, frac1, frac2, term2;
	double sigma_11, sigma_22, sigma_12;

	// Calculating the terms
	std::tie(sigma_11, sigma_22, sigma_12) = sigma_gravity(z, H, rho, g, nu, kappa);
	if (m>0)
	{
		for (size_t ii = z1.size(); ii--;)
		{
			double sigma_11c, sigma_22c, sigma_12c;
			std::tie(sigma_11c, sigma_22c, sigma_12c) = sigma_crack(z, z1[ii], z2[ii], L[ii], mu[ii], m, beta[ii]);
			sigma_11 += sigma_11c;
			sigma_22 += sigma_22c;
			sigma_12 += sigma_12c;
		}
	}
	if (mt > 0)
	{
		for (size_t ii = z0.size(); ii--;)
		{
			double sigma_11t, sigma_22t, sigma_12t;
			std::tie(sigma_11t, sigma_22t, sigma_12t) = sigma_circ_tunnel(z, z0[ii], R[ii], mt, a[ii], b[ii]);
			sigma_11 += sigma_11t;
			sigma_22 += sigma_22t;
			sigma_12 += sigma_12t;
		}
	}
	

	return { sigma_11, sigma_22, sigma_12 };
}

/* --------------------------------------------------------------------
		STRESS FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddvec, ddvec, ddvec> stress_field(double xfrom, double xto, double yfrom, double yto, int Nx, int Ny, dcomp H, double rho, double g, double nu, double kappa, dcvec z1, dcvec z2, dvec L, dvec mu, int m, ddcvec beta, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	ddvec grid_11(Nx, dvec(Ny));
	ddvec grid_22(Nx, dvec(Ny));
	ddvec grid_12(Nx, dvec(Ny));
	double dx;
	double dy;
	dvec x_vec(Nx);
	dvec y_vec(Ny);

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nx-1);
	dy = (yto - yfrom) / (Ny-1);
	#pragma omp parallel for default(none) shared(grid_11, grid_22, grid_12, x_vec, y_vec)
	for (int ii = 0; ii<Nx; ii++)
	{
		for (int jj = Ny; jj--;)
		{
			x_vec[ii] = xfrom + ii * dx;
			y_vec[jj] = yfrom + jj * dy;
			std::tie(grid_11[ii][jj], grid_22[ii][jj], grid_12[ii][jj]) = sigma_total(dcomp(x_vec[ii], y_vec[jj]), H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
		}
	}
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
	frac2 = (sigma_11 - sigma_22) / 2.0;
	sqrt1 = sqrt( frac2 * frac2 + sigma_12 * sigma_12 );

	// Calculating the principal stresses and the angel of sigma_1
	sigma_1 = frac1 + sqrt1;
	sigma_2 = frac1 - sqrt1;

	// Calcuating the angel theta_p
	theta_p = 0.5*atan( (2 * sigma_12) / abs(sigma_11 - sigma_22) ); // <- abs() of stresses

	return { sigma_1, sigma_2, theta_p };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESS FIELDS
-----------------------------------------------------------------------*/
std::tuple<ddvec, ddvec, ddvec> principal_stress_field(ddvec grid_11, ddvec grid_22, ddvec grid_12)
{
	// Defining the variables
	ddvec grid_1(grid_11.size(), dvec(grid_11[0].size()));
	ddvec grid_2(grid_11.size(), dvec(grid_11[0].size()));
	ddvec grid_theta(grid_11.size(), dvec(grid_11[0].size()));
	int Nx = (int)grid_11.size();

	// Retriving the terms from the sigma function for z
	#pragma omp parallel for default(none) shared(grid_1, grid_2, grid_theta)
	for (int ii = 0; ii < Nx; ii++)
	{
		for (size_t jj = grid_11[0].size(); jj--;)
		{
			std::tie(grid_1[ii][jj], grid_2[ii][jj], grid_theta[ii][jj]) = principal_sigma(grid_11[ii][jj], grid_22[ii][jj], grid_12[ii][jj]);
		}
	}
	return { grid_1, grid_2, grid_theta };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESS TRAJECTORIES
-----------------------------------------------------------------------*/
std::tuple<ddcvec, ddcvec> principal_stress_trajectories(double xfrom, double xto, double yfrom, double yto, int Ntraj, int lvs_traj, dcomp H, double rho, double g, double nu, double kappa, dcvec z1, dcvec z2, dvec L, dvec mu, int m, ddcvec beta, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	ddcvec traj_1(lvs_traj, dcvec(Ntraj));
	ddcvec traj_2(lvs_traj, dcvec(Ntraj));
	double dx, dy, dx_lvs, dy_lvs, pi_val;
	double xfromy, yfromx;
	double cond;
	int NIT;

	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj*0.8);
	dy = (yto - yfrom) / (Ntraj*0.8);
	dx_lvs = (xto - xfrom) / (lvs_traj/2 - 1);
	dy_lvs = (yto - yfrom) / (lvs_traj/2 - 1);
	xfromy = xfrom;
	yfromx = yto;
	pi_val = 0.5*pi();
	cond = 1e-6;
	NIT = 10;
	#pragma omp parallel for default(none) shared(traj_1, traj_2)
	for (int ii = 0; ii < lvs_traj / 2; ii++)
	{
		traj_1[ii][Ntraj - 1] = dcomp(xfromy, yfrom + ii * dy_lvs);
		for (int jj = Ntraj - 1; jj--;)
		{
			dcomp zt, z_old, sigma, sigmat;
			double sigma_1, theta_p;
			double sigma_11, sigma_22, sigma_12;
			double ee;
			std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(traj_1[ii][jj+1], H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
			sigma = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
			zt = traj_1[ii][jj + 1] + conj(sigma) / abs(sigma)*dx;
			
			ee = 1;
			for (int kk = NIT; kk--;)
			{
				z_old = zt;
				std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(zt, H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
				sigmat = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
				zt = traj_1[ii][jj + 1] + conj(sigma + sigmat) / abs(sigma + sigmat)*dx;
				ee = std::norm(z_old - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_1[ii][jj] = zt;
		}
		traj_2[ii][Ntraj - 1] = dcomp(xfrom + ii * dx_lvs, yfromx);
		for (int jj = Ntraj - 1; jj--;)
		{
			dcomp zt, z_old, sigma, sigmat;
			double sigma_2, theta_p;
			double sigma_11, sigma_22, sigma_12;
			double ee;
			std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(traj_2[ii][jj+1], H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
			std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
			sigma = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
			zt = traj_2[ii][jj + 1] + conj(sigma) / abs(sigma)*dy;

			ee = 1;
			for (int kk = NIT; kk--;)
			{
				z_old = zt;
				std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(zt, H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
				sigmat = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
				zt = traj_2[ii][jj + 1] + conj(sigma + sigmat) / abs(sigma + sigmat)*dy;
				ee = std::norm(z_old - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_2[ii][jj] = zt;
		}
	}
	dx = -abs(dx);
	dy = -abs(dy);
	xfromy = xto;
	yfromx = yfrom;
	#pragma omp parallel for default(none) shared(traj_1, traj_2)
	for (int ii = lvs_traj / 2; ii < lvs_traj; ii++)
	{
		traj_1[ii][Ntraj - 1] = dcomp(xfromy, yfrom + (ii - lvs_traj / 2) * dy_lvs);
		for (int jj = Ntraj - 1; jj--;)
		{
			dcomp zt, z_old, sigma, sigmat;
			double sigma_1, theta_p;
			double sigma_11, sigma_22, sigma_12;
			double ee;
			std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(traj_1[ii][jj + 1], H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
			sigma = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
			zt = traj_1[ii][jj + 1] + conj(sigma) / abs(sigma)*dx;

			ee = 1;
			for (int kk = NIT; kk--;)
			{
				z_old = zt;
				std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(zt, H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
				sigmat = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
				zt = traj_1[ii][jj + 1] + conj(sigma + sigmat) / abs(sigma + sigmat)*dx;
				ee = std::norm(z_old - zt);
				if (ee < cond)
				{
					break;
				}
			}
			traj_1[ii][jj] = zt;
		}
		traj_2[ii][Ntraj - 1] = dcomp(xfrom + (ii - lvs_traj / 2) * dx_lvs, yfromx);
		for (int jj = Ntraj - 1; jj--;)
		{
			dcomp zt, z_old, sigma, sigmat;
			double sigma_2, theta_p;
			double sigma_11, sigma_22, sigma_12;
			double ee;
			std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(traj_2[ii][jj + 1], H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
			std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
			sigma = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
			zt = traj_2[ii][jj + 1] + conj(sigma) / abs(sigma)*dy;

			ee = 1;
			for (int kk = NIT; kk--;)
			{
				z_old = zt;
				std::tie(sigma_11, sigma_22, sigma_12) = sigma_total(zt, H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(sigma_11, sigma_22, sigma_12);
				sigmat = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
				zt = traj_2[ii][jj + 1] + conj(sigma + sigmat) / abs(sigma + sigmat)*dy;
				ee = std::norm(z_old - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_2[ii][jj] = zt;
		}
	}
	return { traj_1, traj_2 };
}

/* --------------------------------------------------------------------

		MAIN SCRIPT

-----------------------------------------------------------------------*/
int main()
{
	/* --------------------------------------------------------------------
			Defining variables and importing data
	-----------------------------------------------------------------------*/

	// Header in console window
	auto start = std::chrono::high_resolution_clock::now(); // Start the clock

	std::cout << "The documentation for this program can be found on: https://github.com/eriktoller/stress_field_calcualtor \n";
	std::cout << "Written by: Erik Toller, erik.toller@geo.uu.se.\n\n";

	auto date_time1 = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(date_time1);
	char str_time1[26];
	ctime_s(str_time1, sizeof str_time1, &start_time);
	std::cout << "Program started: " << str_time1 << std::endl;

	// Setting the data types
	dcomp z, H;
	double rho, g, nu, kappa;
	double xfrom, xto, yfrom, yto;
	int Nx, Ny, Ntraj, lvs_traj;
	ddvec grid_11, grid_22, grid_12, grid_1, grid_2, theta_p;
	ddcvec traj_1, traj_2;
	dvec x_vec, y_vec;
	int nc, m, nt, mt;

	// Read the input data from binary file PART 1
	std::cout << "Loading input data" << std::endl;
	double fin[6000];
	std::ifstream input_file("input_data.bin", std::ios::in | std::ios::binary);
	input_file.read((char *)&fin, sizeof fin);
	H = dcomp(fin[0], fin[1]);
	rho = fin[2];
	g = fin[3];
	nu = fin[4];
	kappa = fin[5];
	nc = (int)fin[6];
	m = (int)fin[7];
	nt = (int)fin[8];
	mt = (int)fin[9];

	// Declaring the vecotrs
	dcvec z1(nc), z2(nc);
	dvec L(nc), mu(nc);
	ddcvec beta(nc, dcvec(m));
	dcvec z0(nt);
	dvec R(nt);
	ddcvec a(nt, dcvec(mt)), b(nt, dcvec(mt));

	int pos = 9 + 1;
	if (nc > 0)
	{
		for (int ii = 0; ii < nc; ii++)
		{
			int re = pos + ii;
			int im = pos + nc + ii;
			z1[ii] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * nc;
		for (int ii = 0; ii < nc; ii++)
		{
			int re = pos + ii;
			int im = pos + nc + ii;
			z2[ii] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * nc;
		for (int ii = 0; ii < nc; ii++)
		{
			L[ii] = fin[pos+ii];
		}
		pos += nc;
		for (int ii = 0; ii < nc; ii++)
		{
			mu[ii] = fin[pos + ii];
		}
		pos += nc;
		for (int ii = 0; ii < nc; ii++)
	{
		for (int jj = 0; jj < m; jj++)
		{
			int re = pos + jj;
			int im = pos + m + jj;
			beta[ii][jj] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * m;
	}
	}
	else
	{
		z1 = { dcomp(0, 0) };
		z2 = { dcomp(0, 0) };
		L = { 0 };
		mu = { 0 };
		beta = { { dcomp(0,0) } };
	}
	if (nt > 0)
	{
		for (int ii = 0; ii < nt; ii++)
		{
			int re = pos + ii;
			int im = pos + nt + ii;
			z0[ii] = dcomp(fin[re], fin[im]);
		}
		pos += 2 * nt;
		for (int ii = 0; ii < nt; ii++)
		{
			R[ii] = fin[pos + ii];
		}
		pos += nt;
		for (int ii = 0; ii < nt; ii++)
		{
			for (int jj = 0; jj < mt; jj++)
			{
				int re = pos + jj;
				int im = pos + mt + jj;
				a[ii][jj] = dcomp(fin[re], fin[im]);
			}
			pos += 2 * mt;
		}
		for (int ii = 0; ii < nt; ii++)
		{
			for (int jj = 0; jj < mt; jj++)
			{
				int re = pos + jj;
				int im = pos + mt + jj;
				b[ii][jj] = dcomp(fin[re], fin[im]);
			}
			pos += 2 * mt;
		}
	}
	else
	{
		z0 = { dcomp(0, 0) };
		R = { 0 };
		a = { { dcomp(0,0) } };
		b = { { dcomp(0,0) } };
	}
	std::cout << "Complete" << std::endl;

	double s11, s12, s13, st1, st2, st3;
	std::tie(st1, st2, st3) = sigma_crack(dcomp(0, 1), z1[0], z2[0], L[0], mu[0], m, beta[0]);
	std::tie( s11,s12,s13 ) = sigma_total(dcomp(0, 1), H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
	std::cout << " st1 = " << st1 << "\n st2 = " << st2 << "\n st3 = " << st3 << std::endl;
	std::cout << " s11 = " << s11 << "\n s12 = " << s12 << "\n s13 = " << s13 << std::endl;

	// Read the plot data from binary file PART 2
	std::cout << "Loading plot data" << std::endl;
	double fplot[80];
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
	std::cout << "Complete\n" << std::endl;

	// make sure lvs_traj is an even number
	if (lvs_traj % 2 != 0)
	{
		lvs_traj += 1;
	}


	// Displying the plot data in the command window
	std::cout << "This is the retrived plot data:\n";
	std::cout << "x from: " << xfrom << " to " << xto << std::endl;
	std::cout << "y from: " << yfrom << " to " << yto << std::endl;
	std::cout << "x resolution: " << Nx << std::endl;
	std::cout << "y resolution: " << Ny << std::endl;
	std::cout << "Total number of points: " << Nx * Ny << std::endl;
	std::cout << "Number of steps in trajectories: " << Ntraj << std::endl;
	std::cout << "Number of trajectory levels: " << lvs_traj << std::endl;
	std::cout << "Total number of points: " << Ntraj * lvs_traj * 4 << std::endl;

	// Estimate the time of the program
	int Ne = 5000 / (20+m);
	if (Ne < 2)
	{
		Ne = 2;
	}
	if (Ne % 2 != 0)
	{
		Ne += 1;
	}
	auto start0 = std::chrono::high_resolution_clock::now();
	principal_stress_trajectories(xfrom, xto, yfrom, yto, Ne, Ne, H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
	auto stop0 = std::chrono::high_resolution_clock::now();
	auto duration0 = std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0);
	long long ms0 = duration0.count();
	long long calcs = (Nx*Ny/(4*2) + Ntraj*lvs_traj);
	ms0 = ms0/(Ne *Ne) * calcs;
	long long s0 = ms0 / 1000000;
	ms0 = ms0 % 1000000;
	long long m0 = s0 / 60;
	s0 = s0 % 60;
	long long h0 = m0 / 60;
	m0 = m0 % 60;
	std::cout << std::endl;
	if (h0 < 10 && m0 < 10 && s0 < 10)
	{
		std::cout << "Estimated calculation time:  0" << h0 << ":0" << m0 << ":0" << s0 << ":" << ms0;
	}
	else if (h0 < 10 && m0 < 10)
	{
		std::cout << "Estimated calculation time:  0" << h0 << ":0" << m0<< ":" << s0 << ":" << ms0;
	}
	else if (h0 < 10)
	{
		std::cout << "Estimated calculation time:  0" << h0 << ":" << m0 << ":" << s0 << ":" << ms0;
	}
	else
	{
		std::cout << "Estimated calculation time: " << h0 << ":" << m0 << ":" << s0 << ":" << ms0;
	}
	std::cout << std::endl;

	/* --------------------------------------------------------------------
			Calcuating the stress field
	-----------------------------------------------------------------------*/

	// Get the Cartisian stress field
	std::cout << "Initiating the principal stress field calculation\n";
	auto start1 = std::chrono::high_resolution_clock::now();
	std::tie(x_vec, y_vec, grid_11, grid_22, grid_12) = stress_field(xfrom, xto, yfrom, yto, Nx, Ny, H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
	auto stop1 = std::chrono::high_resolution_clock::now();
	auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
	std::cout << "Completed, time taken by function: stress_field = " << duration1.count() << " microseconds" << std::endl << std::endl;

	// Get the principal stress field
	std::cout << "Initiating the stress field calculation\n";
	auto start2 = std::chrono::high_resolution_clock::now();
	std::tie(grid_1, grid_2, theta_p) = principal_stress_field(grid_11, grid_22, grid_12);
	auto stop2 = std::chrono::high_resolution_clock::now();
	auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
	std::cout << "Completed, time taken by function: principal_stress_field = " << duration2.count() << " microseconds" << std::endl << std::endl;

	// Get the stress trajectories
	std::cout << "Initiating the principal stress trajectories calculation\n";
	auto start3 = std::chrono::high_resolution_clock::now();
	std::tie(traj_1, traj_2) = principal_stress_trajectories(xfrom, xto, yfrom, yto, Ntraj, lvs_traj, H, rho, g, nu, kappa, z1, z2, L, mu, m, beta, z0, R, mt, a, b);
	auto stop3 = std::chrono::high_resolution_clock::now();
	auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
	std::cout << "Completed, time taken by function: principal_stress_trajectories = " << duration3.count() << " microseconds" << std::endl << std::endl;

	// Displaying the computation time
	long long ms1 = duration1.count() + duration2.count() + duration3.count();
	long long s1 = ms1 / 1000000;
	ms1 = ms1 % 1000000;
	long long m1 = s1 / 60;
	s1 = s1 % 60;
	long long h1 = m1 / 60;
	m1 = m1 % 60;
	if (h1 < 10 && m1 < 10 && s1 < 10)
	{
		std::cout << "Total calculation time: 0" << h1 << ":0" << m1 << ":0" << s1 << ":" << ms1;
	}
	else if (h1 < 10 && m1 < 10)
	{
		std::cout << "Total calculation time: 0" << h1 << ":0" << m1 << ":" << s1 << ":" << ms1;
	}
	else if (h1 < 10)
	{
		std::cout << "Total calculation time: 0" << h1 << ":" << m1 << ":" << s1 << ":" << ms1;
	}
	else
	{
		std::cout << "Total calculation time: " << h1 << ":" << m1 << ":" << s1 << ":" << ms1;
	}
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

	// Save the plot properties
	dvec dim = { 1.0 * Nx, 1.0 * Ny, 1.0*Ntraj, 1.0*lvs_traj, 1.0*nc, 1.0*nt};
	const char* pointerdim = reinterpret_cast<const char*>(&dim[0]);
	std::size_t bytesdim = dim.size() * sizeof(dim[0]);
	outfiledim.write(pointerdim, bytesdim);
	
	// saving the coordinates of the cracks
	dvec fz1_re(nc);
	dvec fz1_im(nc);
	for (int jj = 0; jj < nc; jj++)
	{
		fz1_re[jj] = real(z1[jj]);
		fz1_im[jj] = imag(z1[jj]);
	}
	const char* pointerz1_re = reinterpret_cast<const char*>(&fz1_re[0]);
	std::size_t bytesz1_re = fz1_re.size() * sizeof(fz1_re[0]);
	outfiledim.write(pointerz1_re, bytesz1_re);
	const char* pointerz1_im = reinterpret_cast<const char*>(&fz1_im[0]);
	std::size_t bytesz1_im = fz1_im.size() * sizeof(fz1_im[0]);
	outfiledim.write(pointerz1_im, bytesz1_im);

	dvec fz2_re(nc);
	dvec fz2_im(nc);
	for (int jj = 0; jj < nc; jj++)
	{
		fz2_re[jj] = real(z2[jj]);
		fz2_im[jj] = imag(z2[jj]);
	}
	const char* pointerz2_re = reinterpret_cast<const char*>(&fz2_re[0]);
	std::size_t bytesz2_re = fz2_re.size() * sizeof(fz2_re[0]);
	outfiledim.write(pointerz2_re, bytesz2_re);
	const char* pointerz2_im = reinterpret_cast<const char*>(&fz2_im[0]);
	std::size_t bytesz2_im = fz2_im.size() * sizeof(fz2_im[0]);
	outfiledim.write(pointerz2_im, bytesz2_im);

	// saving the coordinates of the circular tunnels
	dvec fz0_re(nt);
	dvec fz0_im(nt);
	dvec fR(nt);
	for (int jj = 0; jj < nt; jj++)
	{
		fz0_re[jj] = real(z0[jj]);
		fz0_im[jj] = imag(z0[jj]);
		fR[jj] = R[jj];
	}
	const char* pointerz0_re = reinterpret_cast<const char*>(&fz0_re[0]);
	std::size_t bytesz0_re = fz0_re.size() * sizeof(fz0_re[0]);
	outfiledim.write(pointerz0_re, bytesz0_re);
	const char* pointerz0_im = reinterpret_cast<const char*>(&fz0_im[0]);
	std::size_t bytesz0_im = fz0_im.size() * sizeof(fz0_im[0]);
	outfiledim.write(pointerz0_im, bytesz0_im);
	const char* pointerR = reinterpret_cast<const char*>(&fR[0]);
	std::size_t bytesR = fR.size() * sizeof(fR[0]);
	outfiledim.write(pointerR, bytesR);

	// Close the output files
	outfile.close();
	outfiledim.close();

	std::cout << "Saving completed" << std::endl << std::endl;

	// Get the date and execution time
	auto end = std::chrono::high_resolution_clock::now();
	auto date_time = std::chrono::system_clock::now();
	auto elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	int seconds = (int)elapsed_seconds.count() / 1000000;
	int mseconds = (int)elapsed_seconds.count() % 1000000;
	int minutes = seconds / 60;
	seconds = seconds % 60;
	int hours = minutes / 60;
	minutes = minutes % 60;
	std::time_t end_time = std::chrono::system_clock::to_time_t(date_time);
	char str_time[26];
	ctime_s(str_time, sizeof str_time, &end_time);
	
	// Create the log file
	std::ofstream logfile;
	logfile.open("log.txt");
	logfile << "--------------------------------------------------------------" << std::endl << std::endl;
	logfile << "    STRESS FILED CALCULATOR - LOG FILE                        " << std::endl << std::endl;
	logfile << "--------------------------------------------------------------" << std::endl;
	logfile << "This is the log file for computation completed on: " << str_time << std::endl;
	logfile << "The computation took: ";
	if (hours < 10 && minutes < 10 && seconds < 10)
	{
		logfile << "0" << hours << ":0" << minutes << ":0" << seconds << ":" << mseconds << std::endl;
	}
	else if (hours < 10 && minutes < 10)
	{
		logfile << "0" << hours << ":0" << minutes << ":" << seconds << ":" << mseconds << std::endl;
	}
	else if (hours < 10)
	{
		logfile << "0" << hours << ":" << minutes << ":" << seconds << ":" << mseconds << std::endl;
	}
	else
	{
		logfile << "" << hours << ":" << minutes << ":" << seconds << ":" << mseconds << std::endl;
	}
	logfile << std::endl << std::endl;
	logfile << "--------------------------------------------------------------" << std::endl;
	logfile << "    DATA FILES                                                " << std::endl;
	logfile << "--------------------------------------------------------------" << std::endl;
	logfile << "'data.bin'     [x_vec, y_vec, grid_11, grid_22, grid_12, grid_1, grid_2, theta_p, (traj_1.real, traj_1.imag), (traj_2.real, traj_2.imag)]\n";
	logfile << "'dim_data.bin' [Nx, Ny, Ntraj, lvs_traj, nc, z1, z2]\n";
	logfile << std::endl << std::endl;
	logfile << "--------------------------------------------------------------" << std::endl;
	logfile << "    PLOT DATA                                                 " << std::endl;
	logfile << "--------------------------------------------------------------" << std::endl;
	logfile << "x from: " << xfrom << " to " << xto << std::endl;
	logfile << "y from: " << yfrom << " to " << yto << std::endl;
	logfile << "x resolution: " << Nx << std::endl;
	logfile << "y resolution: " << Ny << std::endl;
	logfile << "Total number of points: " << Nx * Ny << std::endl;
	logfile << "Number of trajectory levels: " << lvs_traj << std::endl;
	logfile << "Number of steps in trajectories: " << Ntraj;
	logfile << std::endl << std::endl;
	logfile << "--------------------------------------------------------------" << std::endl;
	logfile << "    VARIABLES                                                 " << std::endl;
	logfile << "--------------------------------------------------------------" << std::endl;
	logfile << "Elasticity:\n";
	logfile << "     nu = " << nu << std::endl;
	logfile << "  kappa = " << kappa << std::endl;
	logfile << std::endl;
	logfile << "Gravity:\n";
	logfile << "      H = (" << H.real() << "," << H.imag() << ")" << std::endl;
	logfile << "    rho = " << rho << std::endl;
	logfile << "      g = " << g << std::endl;
	logfile << std::endl;
	if (nc>0)
	{
		logfile << "Cracks:\n";
		logfile << "     nc = " << nc << std::endl;
		logfile << "      m = " << m << std::endl;
		logfile << "     z1 = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << z1[ii];
			}
			else
			{
				logfile << z1[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		logfile << "     z2 = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << z2[ii];
			}
			else
			{
				logfile << z2[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		logfile << "      L = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << L[ii];
			}
			else
			{
				logfile << L[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		logfile << "     mu = [";
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii == (nc - 1))
			{
				logfile << mu[ii];
			}
			else
			{
				logfile << mu[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		
		for (int ii = 0; ii < nc; ii++)
		{
			logfile << "beta[" << ii <<"] = [";
			for (int jj = 0; jj < m; jj++)
			{
				if (jj == (m - 1))
				{
					logfile << beta[ii][jj];
				}
				else
				{
					logfile << beta[ii][jj] << ", ";
				}
			}
			logfile << "]" << std::endl;
		}
		logfile << std::endl << std::endl;
	}
	if (nt > 0)
	{
		logfile << "Circular Tunnel:\n";
		logfile << "     nt = " << nt << std::endl;
		logfile << "     mt = " << mt << std::endl;
		logfile << "     z0 = [";
		for (int ii = 0; ii < nt; ii++)
		{
			if (ii == (nt - 1))
			{
				logfile << z0[ii];
			}
			else
			{
				logfile << z0[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		logfile << "      R = [";
		for (int ii = 0; ii < nt; ii++)
		{
			if (ii == (nt - 1))
			{
				logfile << R[ii];
			}
			else
			{
				logfile << R[ii] << ", ";
			}
		}
		logfile << "]" << std::endl;
		for (int ii = 0; ii < nt; ii++)
		{
			logfile << "a[" << ii << "] = [";
			for (int jj = 0; jj < mt; jj++)
			{
				if (jj == (m - 1))
				{
					logfile << a[ii][jj];
				}
				else
				{
					logfile << a[ii][jj] << ", ";
				}
			}
			logfile << "]" << std::endl;
		}
		for (int ii = 0; ii < nt; ii++)
		{
			logfile << "b[" << ii << "] = [";
			for (int jj = 0; jj < mt; jj++)
			{
				if (jj == (m - 1))
				{
					logfile << b[ii][jj];
				}
				else
				{
					logfile << b[ii][jj] << ", ";
				}
			}
			logfile << "]" << std::endl;
		}
		logfile << std::endl << std::endl;
	}

	logfile.close();

	if (hours < 10 && minutes < 10 && seconds < 10)
	{
		std::cout << "Program finnished after 0" << hours << ":0" << minutes << ":0" << seconds << ":" << mseconds << " and output data saved to binary files\n";
	}
	else if (hours < 10 && minutes < 10)
	{
		std::cout << "Program finnished after 0" << hours << ":0" << minutes << ":" << seconds << ":" << mseconds << " and output data saved to binary files\n";
	}
	else if (hours < 10)
	{
		std::cout << "Program finnished after 0" << hours << ":" << minutes << ":" << seconds << ":" << mseconds << " and output data saved to binary files\n";
	}
	else
	{
		std::cout << "Program finnished after " << hours << ":" << minutes << ":" << seconds << ":" << mseconds << " and output data saved to binary files\n";
	}
	std::cout << "For log data see log.txt\n";

	return 0;
}