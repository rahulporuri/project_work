/* dammit! since the rk4_stepper is returning more than one value,
 we need to define an array that will hold its values or a struct! */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double A(double N, double a0);

double an(double n, double n0, double a0);
double aN_nve(double N, double n0, double a0, double p);
double aN_pve(double N, double n0, double a0, double p);

double h(double n, double n0, double a0);
double H_nve(double N, double n0, double a0, double p);
double H_pve(double N, double n0, double a0, double p);

double fn(double n, double n0, double a0);
double fN_nve(double N, double n0, double a0, double p);
double fN_pve(double N, double n0, double a0, double p);

double n_nve(double N, double n0, double a0, double p);
double n_pve(double N, double n0, double a0, double p);

double DH_nve(double N);
double DH_pve(double N);

double DDhk_nve(double k, double N, double hk, double Dhk, double n0, double a0, double p);
double DDhk_pve(double k, double N, double hk, double Dhk, double n0, double a0, double p);

void initialize_hk(double k, double Nics, double *hk, double a0);
void initialize_Dhk(double k, double Nics, double *Dhk, double a0, double n0, double p);

void rk4_stepper_nve(double k, double N, double *hk, double *Dhk, double *update_hk, double *update_Dhk, double n0, double a0, double p, double step);
void rk4_stepper_pve(double k, double N, double *hk, double *Dhk, double *update_hk, double *update_Dhk, double n0, double a0, double p, double step);

int
main(void)
{
	/* define constants */
	double n0, a0, p;

	double k0;
	double N, Nics, Nshss;
	
	double tps;

	int npts;
	double step;

	double hk[2], Dhk[2];
	double increment_hk[2], increment_Dhk[2];

	a0 = pow(10,5);
	n0 = 1;
	p = 1;

	npts = 100000;

	k0 = 1e-30;

	Nics = -17.20805585;
	Nshss = 7.52650892;

	printf("%le, %le, %le, %le\n", hk[0], hk[1], Dhk[0], Dhk[1]);
	printf("%le, %le, %le, %le\n", hk[0], hk[1], Dhk[0], Dhk[1]);
	step = (Nshss -Nics)/(npts);

	while (k0 < 1e+2*k0)
	{
		initialize_hk(k0, Nics, hk, a0);
		initialize_Dhk(k0, Nics, Dhk, a0, n0, p);

		N = Nics;
		while (N < 0 -5*step)
		{
			rk4_stepper_nve(k0, N, hk, Dhk, increment_hk, increment_Dhk, n0, a0, p, step);

			hk[0] += increment_hk[0];
			hk[1] += increment_hk[1];
			Dhk[0] += increment_Dhk[0];
			Dhk[1] += increment_Dhk[1];

			N += step;
		}

		printf("%le, %le, %le, %le\n", hk[0], hk[1], Dhk[0], Dhk[1]);

		N = 0 +5*step;
		while (N < Nshss +step)
		{
			rk4_stepper_pve(k0, N, hk, Dhk, increment_hk, increment_Dhk, n0, a0, p, step);

			hk[0] += increment_hk[0];
			hk[1] += increment_hk[1];
			Dhk[0] += increment_Dhk[0];
			Dhk[1] += increment_Dhk[1];
	
			N += step;
		}

		printf("%le, %le, %le, %le\n", hk[0], hk[1], Dhk[0], Dhk[1]);

		tps = 8*pow(k0, 3)/(2*pow(M_PI,2))*(hk[0]*hk[0] +hk[1]*hk[1]);
		printf("%le \n", tps);

		k0 = 10*k0;
	}
	
	return (0);
}


double A(double N, double a0)
{
	return a0*exp(N*N/2);
}

double an(double n, double n0, double a0)
{
	return a0*(1+(n/n0)*(n/n0));
}

double aN_nve(double N, double n0, double a0, double p)
{
	return an(n_nve(N, n0, a0, p), n0, a0);
}

double aN_pve(double N, double n0, double a0, double p)
{
	return an(n_pve(N, n0, a0, p), n0, a0);
}

double h(double n, double n0, double a0)
{
	return (2*a0*n/(n0*n0))*(1/(an(n, n0, a0)*an(n, n0, a0)));
}

double H_nve(double N, double n0, double a0, double p)
{
	return h(n_nve(N, n0, a0, p), n0, a0);
}

double H_pve(double N, double n0, double a0, double p)
{
	return h(n_pve(N, n0, a0, p), n0, a0);
}

double fn(double n, double n0, double a0)
{
	return (2*a0/(n0*n0))*(1/an(n, n0, a0));
}

double fN_nve(double N, double n0, double a0, double p)
{
	return fn(n_nve(N, n0, a0, p), n0, a0);
}

double fN_pve(double N, double n0, double a0, double p)
{
	return fn(n_pve(N, n0, a0, p), n0, a0);
}

double n_nve(double N, double n0, double a0, double p)
{
	return -n0*sqrt(pow(A(N, a0)/a0, 1/p)-1);
}

double n_pve(double N, double n0, double a0, double p)
{
	return n0*sqrt(pow(A(N, a0)/a0, 1/p)-1);
}

double DH_nve(double N)
{
	return -100000.0*N/sqrt(exp(N*N/2) -1.0)*exp(N*N/2)/pow((100000.0*(exp(N*N/2) -1) +100000.0),2) +40000000000.0*N*sqrt(exp(N*N/2) -1.0)*exp(N*N/2)/pow((100000.0*(exp(N*N/2) -1) +100000.0),3);
}

double DH_pve(double N)
{
	return +100000.0*N/sqrt(exp(N*N/2) -1.0)*exp(N*N/2)/pow((100000.0*(exp(N*N/2) -1) +100000.0),2) -40000000000.0*N*sqrt(exp(N*N/2) -1.0)*exp(N*N/2)/pow((100000.0*(exp(N*N/2) -1) +100000.0),3);
}

/*initialize hk and Dhk values*/

void initialize_hk(double k, double Nics, double *hk, double a0)
{
	hk[0] = 1/(sqrt(2*k))/(A(Nics, a0));
	hk[1] = 0;
	return;
}

void initialize_Dhk(double k, double Nics, double *Dhk, double a0, double n0, double p)
{
	Dhk[0] = -Nics/sqrt(2*k)/A(Nics, a0);
	Dhk[1] = -Nics*(sqrt(k/2))/(A(Nics, a0)*A(Nics, a0)*H_nve(Nics, n0, a0, p));
	return;
}

/* evolve hk values from Nics to Nshss */

double DDhk_nve(double k, double N, double hk, double Dhk, double n0, double a0, double p)
{
	return -(3*N -(1/N) +(DH_nve(N)/H_nve(N, n0, a0, p)))*Dhk -pow(k/(A(N, a0)*H_nve(N, n0, a0, p)),2)*hk;
}

double DDhk_pve(double k, double N, double hk, double Dhk, double n0, double a0, double p)
{
	return -(3*N -(1/N) +(DH_pve(N)/H_pve(N, n0, a0, p)))*Dhk -pow(k/(A(N, a0)*H_pve(N, n0, a0, p)),2)*hk;
}

void rk4_stepper_nve(double k, double N, double *hk, double *Dhk, double *update_hk, double *update_Dhk, double n0, double a0, double p, double step)
{
	double F1_real, f1_real, F2_real, f2_real, F3_real, f3_real, F4_real, f4_real;
	double F1_imag, f1_imag, F2_imag, f2_imag, F3_imag, f3_imag, F4_imag, f4_imag;

	F1_real = Dhk[0];
	f1_real = DDhk_nve(k, N, hk[0], Dhk[0], n0, a0, p);
	F2_real = Dhk[0] +f1_real*step/2.;
	f2_real = DDhk_nve(k, N +step/2., hk[0] +F1_real*step/2., Dhk[0] +f1_real*step/2., n0, a0, p);
	F3_real = Dhk[0] +f2_real*step/2.;
	f3_real = DDhk_nve(k, N +step/2., hk[0] +F2_real*step/2., Dhk[0] +f2_real*step/2., n0, a0, p);
	F4_real = Dhk[0] +f3_real*step;
	f4_real = DDhk_nve(k, N +step, hk[0] +F3_real*step, Dhk[0] +f3_real*step, n0, a0, p);

	F1_imag = Dhk[1];
	f1_imag = DDhk_nve(k, N, hk[1], Dhk[1], n0, a0, p);
	F2_imag = Dhk[1] +f1_imag*step/2.;
	f2_imag = DDhk_nve(k, N +step/2., hk[1] +F1_imag*step/2., Dhk[1] +f1_imag*step/2., n0, a0, p);
	F3_imag = Dhk[1] +f2_imag*step/2.;
	f3_imag = DDhk_nve(k, N +step/2., hk[1] +F2_imag*step/2., Dhk[1] +f2_imag*step/2., n0, a0, p);
	F4_imag = Dhk[1] +f3_imag*step;
	f4_imag = DDhk_nve(k, N +step, hk[1] +F3_imag*step, Dhk[1] +f3_imag*step, n0, a0, p);

	update_hk[0] = (F1_real +2*F2_real +2*F3_real +F4_real)*step/6.;
	update_Dhk[0] = (f1_real +2*f2_real +2*f3_real +f4_real)*step/6.;
	update_hk[1] = (F1_imag +2*F2_imag +2*F3_imag +F4_imag)*step/6.;
	update_Dhk[1] = (f1_imag +2*f2_imag +2*f3_imag +f4_imag)*step/6.;

	return;
}

void rk4_stepper_pve(double k, double N, double *hk, double *Dhk, double *update_hk, double *update_Dhk, double n0, double a0, double p, double step)
{
	double F1_real, f1_real, F2_real, f2_real, F3_real, f3_real, F4_real, f4_real;
	double F1_imag, f1_imag, F2_imag, f2_imag, F3_imag, f3_imag, F4_imag, f4_imag;

	F1_real = Dhk[0];
	f1_real = DDhk_pve(k, N, hk[0], Dhk[0], n0, a0, p);
	F2_real = Dhk[0] +f1_real*step/2.;
	f2_real = DDhk_pve(k, N +step/2., hk[0] +F1_real*step/2., Dhk[0] +f1_real*step/2., n0, a0, p);
	F3_real = Dhk[0] +f2_real*step/2.;
	f3_real = DDhk_pve(k, N +step/2., hk[0] +F2_real*step/2., Dhk[0] +f2_real*step/2., n0, a0, p);
	F4_real = Dhk[0] +f3_real*step;
	f4_real = DDhk_pve(k, N +step, hk[0] +F3_real*step, Dhk[0] +f3_real*step, n0, a0, p);

	F1_imag = Dhk[1];
	f1_imag = DDhk_pve(k, N, hk[1], Dhk[1], n0, a0, p);
	F2_imag = Dhk[1] +f1_imag*step/2.;
	f2_imag = DDhk_pve(k, N +step/2., hk[1] +F1_imag*step/2., Dhk[1] +f1_imag*step/2., n0, a0, p);
	F3_imag = Dhk[1] +f2_imag*step/2.;
	f3_imag = DDhk_pve(k, N +step/2., hk[1] +F2_imag*step/2., Dhk[1] +f2_imag*step/2., n0, a0, p);
	F4_imag = Dhk[1] +f3_imag*step;
	f4_imag = DDhk_pve(k, N +step, hk[1] +F3_imag*step, Dhk[1] +f3_imag*step, n0, a0, p);

	update_hk[0] = (F1_real +2*F2_real +2*F3_real +F4_real)*step/6.;
	update_Dhk[0] = (f1_real +2*f2_real +2*f3_real +f4_real)*step/6.;
	update_hk[1] = (F1_imag +2*F2_imag +2*F3_imag +F4_imag)*step/6.;
	update_Dhk[1] = (f1_imag +2*f2_imag +2*f3_imag +f4_imag)*step/6.;

	return;
}
