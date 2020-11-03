#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <numeric>
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
		FROM MICROSECONDS TO HH:MM:SS:MS
-----------------------------------------------------------------------*/
std::tuple<long long, long long, long long, long long> ms_to_time(long long ms)
{
	// Defining the variables
	long long s, m, h;

	// Calculating the time in hour:minutes:seconds:microseconds
	s = ms / 1000000;
	ms = ms % 1000000;
	m = s / 60;
	s = s % 60;
	h = m / 60;
	m = m % 60;

	return {ms,s,m,h};
}

/* --------------------------------------------------------------------
		PRINT in HH:MM:SS:MS FORMAT
-----------------------------------------------------------------------*/
void ms_to_time(long long ms, long long s, long long m, long long h)
{

	// Print the time
	if (h < 10 && m < 10 && s < 10)
	{
		std::cout << "0" << h << ":0" << m << ":0" << s << ":" << ms;
	}
	else if (h < 10 && m < 10)
	{
		std::cout << "0" << h << ":0" << m << ":" << s << ":" << ms;
	}
	else if (h < 10)
	{
		std::cout << "0" << h << ":" << m << ":" << s << ":" << ms;
	}
	else
	{
		std::cout << "" << h << ":" << m << ":" << s << ":" << ms;
	}

	return;
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
		Z FROM CHI
-----------------------------------------------------------------------*/
inline dcomp z_from_chi(dcomp chi, dcomp z1, dcomp z2, double L, double mu)
{
	// Defining the variables
	dcomp z, z0, Z;

	// Calculating the z from a chi value
	z0 = 0.5*(z1 + z2);
	Z = 0.5*(chi + 1.0 / chi);
	z = 0.5*L*Z*exp(dcomp(0, 1)*mu) + z0;

	return z;
}

/* --------------------------------------------------------------------
		ANGEL CHANGE
-----------------------------------------------------------------------*/
inline double angel_change(dcomp z, dcomp z1, dcomp z2)
{
	// Defining the variables
	double u, v, eta;

	// Calculating the angels
	u = arg(z - z1);
	v = arg(z1 - z2);

	// Correct to angel of 0 to 2pi
	if (u < 0)
	{u += 2 * pi();}
	if (v < 0)
	{v += 2 * pi();}

	// Calculate angel between the vectors
	eta = abs(u - v);

	// Correct to angel of 0 to pi
	if (eta > pi())
	{
		eta -= 2*pi();
	}

	eta = abs(eta);

	return eta;
}

/* --------------------------------------------------------------------
		W UNIFORM STRESS
-----------------------------------------------------------------------*/
inline dcomp  w_uni(dcomp z, double kappa, double G, double sigma_11inf)
{
	// Defining the variables
	dcomp phi_bar, dphi, psi, w;

	// calculating the veriables
	phi_bar = -0.5*sigma_11inf*conj(z);
	dphi = -0.5*sigma_11inf;
	psi = -0.5*sigma_11inf*z;

	// Calculating w
	w = 1 / (4 * G)*((z - conj(z))*dphi + kappa * phi_bar + psi);

	return { w };
}

/* --------------------------------------------------------------------
		W GRIVITY
-----------------------------------------------------------------------*/
inline dcomp w_gravity(dcomp z, double kappa, double G, dcomp H, double rho, double g, double nu)
{
	// Defining the variables
	dcomp w;
	dcomp term1, frac1, frac2;

	// Calculating the terms
	term1 = dcomp(0,1)*rho*g;
	frac1 = (1 - 2 * nu) / (kappa + 1);
	frac2 = 1 / (4 * G);

	// Calcualting the w
	w = term1 * frac1*frac2*((z - conj(z))*(z - conj(z)) - 4.0 * H * H);


	return { w };
}

/* --------------------------------------------------------------------
		TAU CRACK
-----------------------------------------------------------------------*/
inline dcomp w_crack(dcomp z, double kappa, double G, dcomp z1, dcomp z2, double L, double mu, int m, dcvec beta)
{
	// Defining the variables
	dcomp chi, chi_bar, Z, chi_pow;
	dcomp dphi, phi_bar, psi;
	dcomp w, L_frac;
	double n;

	// Getting the chi - and Z - coordinates
	chi = chi_from_z(z, z1, z2, L, mu);
	chi_bar = conj(chi);
	Z = exp(dcomp(0, -1)*mu) * 2.0 * (z - 0.5*(z1 + z2)) / L;

	// Calculating the series
	phi_bar = 0;
	dphi = 0;
	psi = 0;
	n = 0;
	chi_pow = chi * chi - 1.0;
	for (int ii = 0; ii < m; ii++)
	{
		dcomp beta_n;
		n += 1;
		beta_n = beta[ii] * n;
		dphi += conj(beta_n) * pow(chi, (1.0 - n)) / chi_pow;
		phi_bar -= beta[ii] * pow(chi_bar, -n);
		psi += beta[ii] * pow(chi, -n);
	}

	// Multiplying the constants
	L_frac = (4.0 / L) * exp(dcomp(0, -1)*mu);
	dphi *= L_frac;

	// Calcualting w
	w = 1 / (4 * G)*(0.5*L*(Z - conj(Z))*dphi + exp(dcomp(0, -1)*mu)*kappa*phi_bar + exp(dcomp(0, -1)*mu)*psi );

	return { w };
}

/* --------------------------------------------------------------------
		W TUNNEL
-----------------------------------------------------------------------*/
inline dcomp w_circ_tunnel(dcomp z, double kappa, double G, dcomp z0, double R, int mt, dcvec a, dcvec b)
{
	// Defining the variables
	dcomp Z, Z_bar, ZZ;
	dcomp dPhi, Phi_bar, psi, xi;
	dcomp w;
	double n;

	Z = (z - z0) / R;
	Z_bar = conj(Z);
	ZZ = Z * Z_bar;

	if (ZZ.real() < 0.999999)
	{
		w = dcomp(NAN, NAN);
	}
	else
	{
		psi = 0;
		dPhi = 0;
		Phi_bar = 0;
		xi = 0;

		for (int ii = 0; ii < mt; ii++)
		{
			n = ii;
			dcomp c = b[ii] + a[ii];
			dPhi = dPhi + c*pow(Z, -n);
			xi = xi + (n) / (n + 1)*c*pow(Z,(-n - 1));
			psi = psi + a[ii] / (n + 1)*pow(Z, (-n - 1));
			if (n > 1) 
			{
				Phi_bar = Phi_bar + conj(c) / (1 - n)*pow(Z_bar, (1 - n));
			}
		}

		psi = psi * 2.0 * R;
		Phi_bar = Phi_bar * R;

		// Getting the constant from integration
		dcomp w0;
		dcomp c = b[1] + a[1];
		if (arg(Z) >= 0)
		{
			w0 = dcomp(0,1)*R*(1.0 / (4.0 * G))* kappa*conj(c) * pi() / 2.0;
		}
		else
		{
			w0 = dcomp(0, -1)*R*(1.0 / (4.0 * G))* kappa*conj(c) * pi() / 2.0;
		}

		// Calculating the w
		double alpha = arg(Z) + .5*pi();
		w = 1.0 / (4.0 * G)*(-(R*Z_bar)*dPhi + kappa * Phi_bar + R * xi + psi) + w0;
		w = w * exp(dcomp(0,1)*alpha);
	}


	return { w };
}

/* --------------------------------------------------------------------
		W TOTAL
-----------------------------------------------------------------------*/
inline dcomp  w_total(dcomp z, double kappa, double G, dcomp H, double rho, double g, double nu, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	dcomp w, wg;

	w = w_uni(z, kappa, G, sigma_11inf);
	wg = w_gravity(z, kappa, G, H, rho, g, nu);
	w += wg;
	if (m > 0)
	{
		for (int ii = 0; ii < nc; ii++)
		{
			dcomp wc;
			wc = w_crack(z, kappa, G, z1[ii], z2[ii], L[ii], mu[ii], m, beta[ii]);
			w += wc;
		}
	}
	if (mt > 0)
	{
		for (size_t ii = z0.size(); ii--;)
		{
			dcomp wt;
			wt = w_circ_tunnel(z, kappa, G, z0[ii], R[ii], mt, a[ii], b[ii]);
			w += wt;
		}
	}

	return { w };
}

/* --------------------------------------------------------------------
		TAU UNIFORM STRESS
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp>  tau_uni(double sigma_11inf)
{
	// Defining the variables
	dcomp tau_11, tau_12;

	// calculating the tau
	tau_11 = -0.5*sigma_11inf;
	tau_12 = -0.5*sigma_11inf;

	return { tau_11, tau_12 };
}

/* --------------------------------------------------------------------
		TAU CRACK
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp> tau_crack(dcomp z, dcomp z1, dcomp z2, double L, double mu, int m, dcvec beta)
{
	// Defining the variables
	dcomp chi, chi_bar, Z, chi_pow;
	dcomp dphi, dphi_bar, ddphi, dpsi;
	dcomp tau_11, tau_12, S1, L_frac;
	double n;

	// Getting the chi - and Z - coordinates
	chi = chi_from_z(z, z1, z2, L, mu);
	chi_bar = conj(chi);
	Z = exp(dcomp(0, -1)*mu) * 2.0 * (z - 0.5*(z1 + z2)) / L;

	// Calculating the series
	dphi = 0;
	dphi_bar = 0;
	ddphi = 0;
	dpsi = 0;
	n = 0;
	chi_pow = chi * chi - 1.0;
	for (int ii = 0; ii < m; ii++)
	{
		dcomp beta_n;
		n += 1;
		beta_n = beta[ii] * n;
		dphi += conj(beta_n) * pow(chi, (1.0 - n)) / chi_pow;
		dphi_bar += beta_n * pow(chi_bar, (1.0 - n)) / (chi_bar * chi_bar - 1.0);
		ddphi -= conj(beta_n) * pow(chi, (2.0 - n)) / (chi_pow*chi_pow*chi_pow)*((n + 1.0)*chi*chi - n + 1.0);
		dpsi -= beta_n * (pow(chi, (1.0 - n))) / chi_pow;
	}

	// Multiplying the constants
	L_frac = (4.0 / L) * exp(dcomp(0, -1)*mu);
	dphi *= L_frac;
	dphi_bar *= conj(L_frac);
	ddphi *= (16.0 / (L*L))*exp(dcomp(0, -2)*mu);
	dpsi *= L_frac;

	// Calcualting tau
	tau_11 = -0.5*L*(Z - conj(Z))*ddphi - exp(dcomp(0, -1)*mu)*(dphi + dpsi);
	tau_12 = -exp(dcomp(0, 1)*mu)*dphi - exp(dcomp(0, -1)*mu)*dphi_bar;

	return { tau_11, tau_12 };
}

/* --------------------------------------------------------------------
		TAU GRAVITY
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp> tau_gravity(dcomp z, dcomp H, double rho, double g, double nu, double kappa)
{
	// Defining the variables
	dcomp sigma_a, sigma_b;
	dcomp tau_11, tau_12;
	dcomp term1, frac1, frac2, term2;

	// Calculating the terms
	term1 = dcomp(0, 2)*rho*g;
	frac1 = (1 - 2 * nu) / (kappa + 1);
	frac2 = 1 / (kappa + 1);
	term2 = (z - conj(z)) - 2.0*H;

	// Calculating the tau_11, tau_12
	tau_11 = -term1 * frac1 * term2;
	tau_12 = term1 * frac2 * term2;


	return { tau_11, tau_12 };
}


/* --------------------------------------------------------------------
		TAU CIRCULAR TUNNEL
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp> tau_circ_tunnel(dcomp z, dcomp z0, double R, int mt, dcvec a, dcvec b)
{
	// Defining the variables
	dcomp Z, Z_bar, ZZ;
	dcomp dPhi, dPhi_bar, ddPhi, dPsi;
	dcomp tau_11, tau_12;
	double n;

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
			dPsi += a[ii] * pow(Z, (-n - 2));
			dPhi += c * pow(Z, -n);
			dPhi_bar += conj(c)*pow(Z_bar, -n);
			ddPhi += n * c*pow(Z, (-n - 1));
		}

		dPsi = -dPsi * 2.0;
		ddPhi = -ddPhi / R;
		tau_11 = (R*Z_bar - R / Z)*ddPhi - dPsi;
		tau_12 = -dPhi - dPhi_bar;
	}

	return { tau_11, tau_12 };
}

/* --------------------------------------------------------------------
		TAU TOTAL
-----------------------------------------------------------------------*/
inline std::tuple<dcomp, dcomp>  tau_total(dcomp z, dcomp H, double rho, double g, double nu, double kappa, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	dcomp tau_11, tau_12;
	dcomp tau_11g, tau_12g;

	std::tie(tau_11, tau_12) = tau_uni(sigma_11inf);
	std::tie(tau_11g, tau_12g) = tau_gravity(z, H, rho, g, nu, kappa);
	tau_11 += tau_11g;
	tau_12 += tau_12g;
	if (m > 0)
	{
		for (int ii = 0; ii < nc; ii++)
		{
			if (ii != m_not)
			{
				dcomp tau_11c, tau_12c;
				std::tie(tau_11c, tau_12c) = tau_crack(z, z1[ii], z2[ii], L[ii], mu[ii], m, beta[ii]);
				tau_11 += tau_11c;
				tau_12 += tau_12c;
			}
		}
	}
	if (mt > 0)
	{
		for (size_t ii = z0.size(); ii--;)
		{
			dcomp tau_11t, tau_12t;
			std::tie(tau_11t, tau_12t) = tau_circ_tunnel(z, z0[ii], R[ii], mt, a[ii], b[ii]);
			tau_11 += tau_11t;
			tau_12 += tau_12t;
		}
	}

	return { tau_11, tau_12 };
}


/* --------------------------------------------------------------------
		STRESS TOTAL
-----------------------------------------------------------------------*/
inline std::tuple<double, double, double> sigma_total(dcomp z, dcomp H, double rho, double g, double nu, double kappa, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	dcomp S1, S2, tau_11, tau_12;
	double sigma_11, sigma_22, sigma_12;

	// Calculating the tau
	std::tie(tau_11, tau_12) = tau_total(z, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
	
	// Calculate the sigmas
	S1 = .5*(tau_11 + tau_12);
	S2 = .5*(-tau_11 + tau_12);
	sigma_11 = real(S1);
	sigma_22 = real(S2);
	sigma_12 = -imag(S1);

	return { sigma_11, sigma_22, sigma_12 };
}

/* --------------------------------------------------------------------
		STRESS FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddvec, ddvec, ddvec> stress_field(double xfrom, double xto, double yfrom, double yto, int Nx, int Ny, dcomp H, double rho, double g, double nu, double kappa, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
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
			std::tie(grid_11[ii][jj], grid_22[ii][jj], grid_12[ii][jj]) = sigma_total(dcomp(x_vec[ii], y_vec[jj]), H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
		}
	}
	return {x_vec, y_vec, grid_11, grid_22, grid_12};
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESSES
-----------------------------------------------------------------------*/
std::tuple<double, double, double> principal_sigma(dcomp z, dcomp H, double rho, double g, double nu, double kappa, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	double sigma_1;
	double sigma_2;
	double theta_p;
	double frac1, frac2, sqrt1;
	dcomp S1, S2, tau_11, tau_12;
	double sigma_11, sigma_22, sigma_12;

	// Calculating the tau
	std::tie(tau_11, tau_12) = tau_total(z, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);

	// Calculate the sigmas
	S1 = .5*(tau_11 + tau_12);
	S2 = .5*(-tau_11 + tau_12);
	sigma_11 = real(S1);
	sigma_22 = real(S2);
	sigma_12 = -imag(S1);

	// Calculating the terms
	frac1 = (sigma_11 + sigma_22) / 2.0;
	frac2 = (sigma_11 - sigma_22) / 2.0;
	sqrt1 = sqrt( frac2 * frac2 + sigma_12 * sigma_12 );

	// Calculating the principal stresses and the angel of sigma_1
	sigma_1 = frac1 + sqrt1;
	sigma_2 = frac1 - sqrt1;

	// Calcuating the angel theta_p
	theta_p = -0.5*imag(log(tau_11));

	return { sigma_1, sigma_2, theta_p };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESS FIELDS
-----------------------------------------------------------------------*/
std::tuple<ddvec, ddvec, ddvec> principal_stress_field(double xfrom, double xto, double yfrom, double yto, int Nx, int Ny, dcomp H, double rho, double g, double nu, double kappa, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	ddvec grid_1(Nx, dvec(Ny));
	ddvec grid_2(Nx, dvec(Ny));
	ddvec grid_tp(Nx, dvec(Ny));
	double dx;
	double dy;
	dvec x_vec(Nx);
	dvec y_vec(Ny);

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nx - 1);
	dy = (yto - yfrom) / (Ny - 1);
	#pragma omp parallel for default(none) shared(grid_1, grid_2, grid_tp, x_vec, y_vec)
	for (int ii = 0; ii < Nx; ii++)
	{
		for (int jj = Ny; jj--;)
		{
			x_vec[ii] = xfrom + ii * dx;
			y_vec[jj] = yfrom + jj * dy;
			std::tie(grid_1[ii][jj], grid_2[ii][jj], grid_tp[ii][jj]) = principal_sigma(dcomp(x_vec[ii], y_vec[jj]), H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
		}
	}
	return { grid_1, grid_2, grid_tp };
}

/* --------------------------------------------------------------------
		PRINCIPAL STRESS TRAJECTORIES
-----------------------------------------------------------------------*/
std::tuple<ddcvec, ddcvec> principal_stress_trajectories(double xfrom, double xto, double yfrom, double yto, dcvec xtraj, dcvec ytraj, int Ntraj, int lvs_traj, dcomp H, double rho, double g, double nu, double kappa, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, int m_not, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	ddcvec traj_1(lvs_traj*2, dcvec(Ntraj));
	ddcvec traj_2(lvs_traj*2, dcvec(Ntraj));
	double dx, dy, dx_lvsre, dx_lvsim, dy_lvsre, dy_lvsim, pi_val;
	double cond;
	int NIT;

	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj);
	dy = (yto - yfrom) / (Ntraj);
	dx_lvsre = (real(xtraj[1]) - real(xtraj[0])) / (lvs_traj - 1);
	dx_lvsim = (imag(xtraj[1]) - imag(xtraj[0])) / (lvs_traj - 1);
	dy_lvsre = (real(ytraj[1]) - real(ytraj[0])) / (lvs_traj - 1);
	dy_lvsim = (imag(ytraj[1]) - imag(ytraj[0])) / (lvs_traj - 1);
	pi_val = 0.5*pi();
	cond = 1e-6;
	NIT = 10;

	// SIGMA 1
	#pragma omp parallel for default(none) shared(traj_1)
	for (int ii = 0; ii < lvs_traj; ii++)
	{
		traj_1[ii][0] = dcomp(real(xtraj[0]) + ii * dx_lvsre, imag(xtraj[0]) + ii * dx_lvsim);
		dcomp z = traj_1[ii][0];
		dcomp z_old = z;
		double dx1 = dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_1, theta_p;
			double ee, eta;
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(z, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
			sigma = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
			z11 = z + sigma / abs(sigma)*dx1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dx1 = -dx1;
			}
			zt = z + sigma / abs(sigma)*dx1;
			
			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(zt, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
				sigmat = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dx1 = -dx1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_1[ii][jj] = zt;
			z_old = z;
			z = zt;
		}

		int kk = ii + lvs_traj;
		traj_1[kk][0] = traj_1[ii][0];
		z = traj_1[kk][0];
		z_old = z;
		dx1 = -dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_1, theta_p;
			double ee, eta;
			std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(z, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
			sigma = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
			z11 = z + sigma / abs(sigma)*dx1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dx1 = -dx1;
			}
			zt = z + sigma / abs(sigma)*dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(sigma_1, std::ignore, theta_p) = principal_sigma(zt, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
				sigmat = abs(sigma_1)*exp(dcomp(0, 1)*theta_p);
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dx1 = -dx1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_1[kk][jj] = zt;
			z_old = z;
			z = zt;
		}
	}

	// SIGMA 2
	//#pragma omp parallel for default(none) shared(traj_2)
	for (int ii = 0; ii < lvs_traj; ii++)
	{
		traj_2[ii][0] = dcomp(real(ytraj[0]) + ii * dy_lvsre, imag(ytraj[0]) + ii * dy_lvsim);
		double dy1 = dy;
		dcomp z = traj_2[ii][0];
		dcomp z_old = z;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_2, theta_p;
			double ee, eta;
			std::tie(std::ignore, sigma_2,  theta_p) = principal_sigma(z, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
			sigma = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
			z11 = z + sigma / abs(sigma)*dy1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dy1 = -dy1;
			}
			zt = z + sigma / abs(sigma)*dy1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(zt, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
				sigmat = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
				z11 = z + (sigma + sigmat) / abs(sigma + sigmat)*dy1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dy1 = -dy1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dy1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_2[ii][jj] = zt;
			z_old = traj_2[ii][jj-1];
			z = zt;
		}

		int kk = ii + lvs_traj;
		traj_2[kk][0] = traj_2[ii][0];
		dy1 = -dy;
		z = traj_2[kk][0];
		z_old = z;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z11, z_oldt, sigma, sigmat;
			double sigma_2, theta_p;
			double ee, eta;
			std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(z, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
			sigma = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
			z11 = z + sigma / abs(sigma)*dy1;
			eta = angel_change(z11, z, z_old);
			if (eta > pi_val && jj > 1) {
				dy1 = -dy1;
			}
			zt = z + sigma / abs(sigma)*dy1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				std::tie(std::ignore, sigma_2, theta_p) = principal_sigma(zt, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, m_not, z0, R, mt, a, b);
				sigmat = abs(sigma_2)*exp(dcomp(0, 1)*(theta_p + pi_val));
				z11 = z + (sigma + sigmat) / abs(sigma)*dy1;
				eta = angel_change(z11, z, z_old);
				if (eta > pi_val && jj > 1) {
					dy1 = -dy1;
				}
				zt = z + (sigma + sigmat) / abs(sigma + sigmat)*dy1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_2[kk][jj] = zt;
			z_old = z;
			z = zt;
		}
	}
	return { traj_1, traj_2 };
}

/* --------------------------------------------------------------------
		DISPLACEMENT FIELD
-----------------------------------------------------------------------*/
std::tuple<dvec, dvec, ddcvec> w_field(double xfrom, double xto, double yfrom, double yto, int Nw, double kappa, double G, dcomp H, double rho, double g, double nu, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)

{
	// Defining the variables
	ddcvec grid_w(Nw, dcvec(Nw));
	dvec x_vecw(Nw), y_vecw(Nw);
	double dx;
	double dy;

	// Retriving the terms from the sigma function for z
	dx = (xto - xfrom) / (Nw - 1);
	dy = (yto - yfrom) / (Nw - 1);
	//#pragma omp parallel for default(none) shared(grid_w, x_vecw, y_vecw)
	for (int ii = 0; ii < Nw; ii++)
	{
		for (int jj = 0; jj < Nw; jj++)
		{
			
			x_vecw[ii] = xfrom + ii * dx;
			y_vecw[jj] = yfrom + jj * dy;
			grid_w[ii][jj] = w_total(dcomp(x_vecw[ii], y_vecw[jj]), kappa, G, H, rho, g, nu, sigma_11inf, z1, z2, L, mu, m, nc, beta, z0, R, mt, a, b);
		}
	}
	return { x_vecw, y_vecw, grid_w };
}

/* --------------------------------------------------------------------
		DISPLACEMENT TRAJECTORIES
-----------------------------------------------------------------------*/
ddcvec w_trajectories(double xfrom, double xto, double yfrom, double yto, int Ntraj, int Nw, double kappa, double G, dcomp H, double rho, double g, double nu, double sigma_11inf, dcvec z1, dcvec z2, dvec L, dvec mu, int m, int nc, ddcvec beta, dcvec z0, dvec R, int mt, ddcvec a, ddcvec b)
{
	// Defining the variables
	ddcvec traj_w(Nw * 2, dcvec(Ntraj));
	double dx, dy_lvs;
	double cond, pi_val;
	int NIT;

	// Getting the starting points
	dx = (xto - xfrom) / (Ntraj);
	dy_lvs = (yto - yfrom) / (Nw-1);
	pi_val = 0.5*pi();
	cond = 1e-6;
	NIT = 10;

	// w trajectories
	#pragma omp parallel for default(none) shared(traj_w)
	for (int ii = 0; ii < Nw; ii++)
	{
		traj_w[ii][0] = dcomp(xfrom, yfrom + ii * dy_lvs);
		dcomp z = traj_w[ii][0];
		dcomp z_old = z;
		double dx1 = dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z_oldt, w, w1;
			double ee;
			w = w_total(z, kappa, G, H, rho, g, nu, sigma_11inf, z1, z2, L, mu, m, nc, beta, z0, R, mt, a, b);
			zt = z + conj(w) / abs(w)*dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				w1 = w_total(zt, kappa, G, H, rho, g, nu, sigma_11inf, z1, z2, L, mu, m, nc, beta, z0, R, mt, a, b);
				zt = z + conj(w + w1) / abs(w + w1)*dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}

			traj_w[ii][jj] = zt;
			z_old = z;
			z = zt;
		}

		int kk = ii + Nw;
		traj_w[kk][0] = dcomp(xto, yfrom + ii * dy_lvs);
		z = traj_w[kk][0];
		z_old = z;
		dx1 = dx;
		for (int jj = 1; jj < Ntraj; jj++)
		{
			dcomp zt, z_oldt, w, w1;
			double ee;
			w = w_total(z, kappa, G, H, rho, g, nu, sigma_11inf, z1, z2, L, mu, m, nc, beta, z0, R, mt, a, b);
			zt = z + conj(w) / abs(w)*dx1;

			ee = 1;
			for (int rr = NIT; rr--;)
			{
				z_oldt = zt;
				w1 = w_total(zt, kappa, G, H, rho, g, nu, sigma_11inf, z1, z2, L, mu, m, nc, beta, z0, R, mt, a, b);
				zt = z + conj(w + w1) / abs(w + w1)*dx1;
				ee = std::norm(z_oldt - zt);
				if (ee < cond)
				{
					break;
				}
			}
				
			traj_w[kk][jj] = zt;
			z_old = z;
			z = zt;

		}
	}
	return { traj_w };
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

	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		STRESS AND DISPLACEMENT CALCULATOR	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "The documentation for this program can be found on:\nhttps://github.com/eriktoller/stress_field_calcualtor \n";
	std::cout << "Written by: Erik Toller, erik.toller@geo.uu.se.\n\n";
	

	auto date_time1 = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(date_time1);
	char str_time1[26];
	ctime_s(str_time1, sizeof str_time1, &start_time);
	std::cout << "Program started: " << str_time1 << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		LODING THE DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	

	// Setting the data types
	dcomp z, H;
	double rho, g, nu, kappa, sigma_11inf, G;
	double xfrom, xto, yfrom, yto;
	int Nx, Ny, Nw, Ntraj, lvs_traj;
	ddvec grid_11, grid_22, grid_12, grid_1, grid_2, theta_p;
	ddcvec traj_1, traj_2, grid_w, traj_w;
	dvec x_vec, y_vec, x_vecw, y_vecw;
	dcvec xtraj, ytraj;
	int nc, m, nt, mt;

	// Read the input data from binary file PART 1
	std::cout << "Loading input data" << std::endl;
	double fin[6000];
	std::ifstream input_file("input_data.bin", std::ios::in | std::ios::binary);
	input_file.read((char *)&fin, sizeof fin);
	H = dcomp(fin[0], fin[1]);
	rho = fin[2];
	g = fin[3];
	sigma_11inf = fin[4];
	nu = fin[5];
	kappa = fin[6];
	G = fin[7];
	nc = (int)fin[8];
	m = (int)fin[9];
	nt = (int)fin[10];
	mt = (int)fin[11];

	// Declaring the vecotrs
	dcvec z1(nc), z2(nc);
	dvec L(nc), mu(nc);
	ddcvec beta(nc, dcvec(m));
	dcvec z0(nt);
	dvec R(nt);
	ddcvec a(nt, dcvec(mt)), b(nt, dcvec(mt));

	int pos = 11 + 1;
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
	std::cout << " -> Complete" << std::endl;

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
	Nw = (int)fplot[6];
	Ntraj = (int) fplot[7];
	lvs_traj = (int) fplot[8];
	xtraj = { fplot[9] + dcomp(0,1)*fplot[10], fplot[11] + dcomp(0,1)*fplot[12] };
	ytraj = { fplot[13] + dcomp(0,1)*fplot[14], fplot[15] + dcomp(0,1)*fplot[16] };
	std::cout << " -> Complete\n" << std::endl;

	// Displying the plot data in the command window
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		THE READ PLOT DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "This is the retrived plot data:\n";
	std::cout << "x from: " << xfrom << " to " << xto << std::endl;
	std::cout << "y from: " << yfrom << " to " << yto << std::endl;
	std::cout << "x resolution: " << Nx << std::endl;
	std::cout << "y resolution: " << Ny << std::endl;
	std::cout << "Total number of points: " << Nx * Ny << std::endl;
	std::cout << "Number of steps in trajectories: " << Ntraj << std::endl;
	std::cout << "Number of trajectory levels: " << lvs_traj << std::endl;
	std::cout << "Total number of trajectory points: " << Ntraj * lvs_traj * 4 << std::endl;
	std::cout << "sigma_1 starting line from: " << xtraj[0] << " to " << xtraj[1] << std::endl;
	std::cout << "sigma_2 starting line from: " << ytraj[0] << " to " << ytraj[1] << std::endl;
	std::cout << "x and y quiver resolution: " << Nw << std::endl << std::endl;

	// Estimate the time of the program
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		ESTIMATED CALCULATION TIME	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
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
	principal_stress_field(xfrom, xto, yfrom, yto, Ne, Ne, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1, z0, R, mt, a, b);
	auto stop0 = std::chrono::high_resolution_clock::now();
	auto start01 = std::chrono::high_resolution_clock::now();
	principal_stress_trajectories(xfrom, xto, yfrom, yto, xtraj, ytraj, Ne, Ne, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1, z0, R, mt, a, b);
	auto stop01 = std::chrono::high_resolution_clock::now();
	auto duration0 = std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0);
	long long ms0 = duration0.count();
	auto duration01 = std::chrono::duration_cast<std::chrono::microseconds>(stop01 - start01);
	long long ms01 = duration01.count();
	long long calcs0 = (Nx*Ny + Nx*Ny + Nw*Nw);
	long long calcs01 = (Ntraj * lvs_traj + Ntraj * Nw);
	ms0 = ms0/(Ne*Ne) * calcs0 + ms01 / (Ne*Ne*2) * calcs01;
	long long s0, m0, h0;
	std::tie(ms0, s0, m0, h0) = ms_to_time(ms0);
	std::cout << "Estimated calculation time:  ";
	ms_to_time(ms0, s0, m0, h0);
	std::cout << std::endl << std::endl;

	/* --------------------------------------------------------------------
			Calcuating the stress field
	-----------------------------------------------------------------------*/
	
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		COMPUTING THE PLOTS	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;

	// Get the Cartisian stress field
	std::cout << "Initiating the stress field calculation\n";
	auto start1 = std::chrono::high_resolution_clock::now();
	std::tie(x_vec, y_vec, grid_11, grid_22, grid_12) = stress_field(xfrom, xto, yfrom, yto, Nx, Ny, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1, z0, R, mt, a, b);
	auto stop1 = std::chrono::high_resolution_clock::now();
	auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
	std::cout << "Completed, time taken by function: stress_field = ";
	long long ms1, s1, m1, h1;
	ms1 = duration1.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	ms_to_time(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	// Get the principal stress field
	std::cout << "Initiating the principal stress field calculation\n";
	auto start2 = std::chrono::high_resolution_clock::now();
	std::tie(grid_1, grid_2, theta_p) = principal_stress_field(xfrom, xto, yfrom, yto, Nx, Ny, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1, z0, R, mt, a, b);
	auto stop2 = std::chrono::high_resolution_clock::now();
	auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
	std::cout << "Completed, time taken by function: principal_stress_field = ";
	ms1 = duration2.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	ms_to_time(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	// Get the stress trajectories
	std::cout << "Initiating the principal stress trajectories calculation\n";
	auto start3 = std::chrono::high_resolution_clock::now();
	std::tie(traj_1, traj_2) = principal_stress_trajectories(xfrom, xto, yfrom, yto, xtraj, ytraj, Ntraj, lvs_traj, H, rho, g, nu, kappa, sigma_11inf, z1, z2, L, mu, m, nc, beta, -1, z0, R, mt, a, b);
	auto stop3 = std::chrono::high_resolution_clock::now();
	auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
	std::cout << "Completed, time taken by function: principal_stress_trajectories = ";
	ms1 = duration3.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	ms_to_time(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	// Get the displacement field
	std::cout << "Initiating the displacement field calculation\n";
	auto start4 = std::chrono::high_resolution_clock::now();
	std::tie(x_vecw, y_vecw, grid_w) = w_field(xfrom, xto, yfrom, yto, Nw, kappa, G, H, rho, g, nu, sigma_11inf, z1, z2, L, mu, m, nc, beta, z0, R, mt, a, b);
	auto stop4 = std::chrono::high_resolution_clock::now();
	auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
	std::cout << "Completed, time taken by function: w_field = ";
	ms1 = duration4.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	ms_to_time(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;

	
	// Get the displacement trajectories
	std::cout << "Initiating the displacement trajectories calculation\n";
	auto start5 = std::chrono::high_resolution_clock::now();
	traj_w = w_trajectories(xfrom, xto, yfrom, yto, Ntraj, Nw, kappa, G, H, rho, g, nu, sigma_11inf, z1, z2, L, mu, m, nc, beta, z0, R, mt, a, b);
	auto stop5 = std::chrono::high_resolution_clock::now();
	auto duration5 = std::chrono::duration_cast<std::chrono::microseconds>(stop5 - start5);
	std::cout << "Completed, time taken by function: w_trajectories = ";
	ms1 = duration5.count();
	std::tie(ms1, s1, m1, h1) = ms_to_time(ms1);
	ms_to_time(ms1, s1, m1, h1);
	std::cout << std::endl << std::endl;


	// Displaying the computation time
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		TOTAL CALCULATION TIME	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;
	long long mstime = duration1.count() + duration2.count() + duration3.count() + duration4.count() + duration5.count();
	long long stime, mtime, htime;
	std::tie(mstime, stime, mtime, htime) = ms_to_time(mstime);
	std::cout << "Total calculation time: ";
	ms_to_time(mstime, stime, mtime, htime);
	std::cout << std::endl << std::endl;

	/* -------------------------------------------------------------------- 
			Saving the data as binary files
	-----------------------------------------------------------------------*/
	std::cout << "=================================================================" << std::endl << std::endl;
	std::cout << "		SAVING THE OUTPUT DATA	" << std::endl << std::endl;
	std::cout << "=================================================================" << std::endl << std::endl;

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

	dvec fxw = x_vecw;
	const char* pointerxw = reinterpret_cast<const char*>(&fxw[0]);
	std::size_t bytesxw = fxw.size() * sizeof(fxw[0]);
	outfile.write(pointerxw, bytesxw);

	dvec fyw = y_vecw;
	const char* pointeryw = reinterpret_cast<const char*>(&fyw[0]);
	std::size_t bytesyw = fyw.size() * sizeof(fyw[0]);
	outfile.write(pointeryw, bytesyw);

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
	for (size_t ii = 0; ii < grid_w.size(); ii++)
	{
		dvec fgw_re(Nw);
		dvec fgw_im(Nw);
		for (size_t jj = 0; jj < grid_w[0].size(); jj++)
		{
			fgw_re[jj] = real(grid_w[ii][jj]);
			fgw_im[jj] = imag(grid_w[ii][jj]);
		}
		const char* pointergw_re = reinterpret_cast<const char*>(&fgw_re[0]);
		std::size_t bytesgw_re = fgw_re.size() * sizeof(fgw_re[0]);
		outfile.write(pointergw_re, bytesgw_re);
		const char* pointergw_im = reinterpret_cast<const char*>(&fgw_im[0]);
		std::size_t bytesgw_im = fgw_im.size() * sizeof(fgw_im[0]);
		outfile.write(pointergw_im, bytesgw_im);
	}
	for (size_t ii = 0; ii < traj_w.size(); ii++)
	{
		dvec fgtw_re(Ntraj);
		dvec fgtw_im(Ntraj);
		for (size_t jj = 0; jj < traj_w[0].size(); jj++)
		{
			fgtw_re[jj] = real(traj_w[ii][jj]);
			fgtw_im[jj] = imag(traj_w[ii][jj]);
		}
		const char* pointergtw_re = reinterpret_cast<const char*>(&fgtw_re[0]);
		std::size_t bytesgtw_re = fgtw_re.size() * sizeof(fgtw_re[0]);
		outfile.write(pointergtw_re, bytesgtw_re);
		const char* pointergtw_im = reinterpret_cast<const char*>(&fgtw_im[0]);
		std::size_t bytesgtw_im = fgtw_im.size() * sizeof(fgtw_im[0]);
		outfile.write(pointergtw_im, bytesgtw_im);
	}

	// Save the plot properties
	dvec dim = { 1.0 * Nx, 1.0 * Ny, 1.0*Nw, 1.0*Ntraj, 1.0*lvs_traj, 1.0*nc, 1.0*nt};
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

	std::cout << " -> Complete" << std::endl << std::endl;

	// Get the date and execution time
	auto end = std::chrono::high_resolution_clock::now();
	auto date_time = std::chrono::system_clock::now();
	auto elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	long long mseconds = elapsed_seconds.count();
	long long seconds, hours, minutes;
	std::tie(mseconds, seconds, minutes, hours) = ms_to_time(mseconds);
	std::time_t end_time = std::chrono::system_clock::to_time_t(date_time);
	char str_time[26];
	ctime_s(str_time, sizeof str_time, &end_time);
	
	// Create the log file
	std::ofstream logfile;
	logfile.open("log.txt");
	logfile << "=================================================================" << std::endl << std::endl;
	logfile << "    STRESS AND DISPLACEMENT CALCULATOR - LOG FILE                 " << std::endl << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "This is the log file for computation completed on: " << str_time << std::endl;
	logfile << "The program took: ";
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
	logfile << "=================================================================" << std::endl;
	logfile << "    DATA FILES								           " << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "'data.bin'     [x_vec, y_vec, x_vecw, y_vecw, grid_11, grid_22, grid_12, grid_1, grid_2, theta_p, (traj_1.real, traj_1.imag), (traj_2.real, traj_2.imag), (grid_w.real, grid_w.imag), (traj_w.real, traj_w.imag)]\n";
	logfile << "'dim_data.bin' [Nx, Ny, Nw, Ntraj, lvs_traj, nc, z1, z2]\n";
	logfile << std::endl << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "    PLOT DATA                                          " << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "x from: " << xfrom << " to " << xto << std::endl;
	logfile << "y from: " << yfrom << " to " << yto << std::endl;
	logfile << "x resolution: " << Nx << std::endl;
	logfile << "y resolution: " << Ny << std::endl;
	logfile << "Total number of points: " << Nx * Ny << std::endl;
	logfile << "Number of trajectory levels: " << lvs_traj << std::endl;
	logfile << "Number of steps in trajectories: " << Ntraj;
	logfile << "sigma_1 starting line from: " << xtraj[0] << " to " << xtraj[1] << std::endl;
	logfile << "sigma_2 starting line from: " << ytraj[0] << " to " << ytraj[1] << std::endl;
	logfile << "x and y quiver resolution: " << Nw;
	logfile << std::endl << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "    VARIABLES										   " << std::endl;
	logfile << "=================================================================" << std::endl;
	logfile << "Elasticity:\n";
	logfile << "         nu = " << nu << std::endl;
	logfile << "      kappa = " << kappa << std::endl;
	logfile << std::endl;
	logfile << "Unifrom stress state:\n";
	logfile << "sigma_11inf = " << sigma_11inf << std::endl;
	logfile << std::endl;
	logfile << "Gravity:\n";
	logfile << "          H = (" << H.real() << "," << H.imag() << ")" << std::endl;
	logfile << "        rho = " << rho << std::endl;
	logfile << "          g = " << g << std::endl;
	logfile << std::endl;
	if (nc>0)
	{
		logfile << "Cracks:\n";
		logfile << "         nc = " << nc << std::endl;
		logfile << "          m = " << m << std::endl;
		logfile << "         z1 = [";
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
		logfile << "         z2 = [";
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
		logfile << "          L = [";
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
		logfile << "         mu = [";
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
			logfile << "    beta[" << ii <<"] = [";
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
		logfile << "         nt = " << nt << std::endl;
		logfile << "         mt = " << mt << std::endl;
		logfile << "         z0 = [";
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
		logfile << "          R = [";
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
			logfile << "    a[" << ii << "] = [";
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
			logfile << "    b[" << ii << "] = [";
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
	std::cout << "Program finnished after ";
	ms_to_time(mseconds, seconds, minutes, hours);
	std::cout << std::endl;
	std::cout << "Uutput data saved to binary files: data.bin and dim_data.bin" << std::endl;
	std::cout << "For log data see log.txt" << std::endl;

	return 0;
}