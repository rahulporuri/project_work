/* dammit! since the rk4_stepper is returning more than one value,
 we need to define an array that will hold its values or a struct! */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double A(double N, double a0);

double an(double n, double n0, double a0)
double aN_nve(double N, double n0, double a0, double p)
double aN_pve(double N, double n0, double a0, double p)

double h(double n, double n0, double a0)
double H_nve(double N, double n0, double a0, double p)
double H_pve(double N, double n0, double a0, double p)

double fn(double n, double n0, double a0)
double fN_nve(double N, double n0, double double a0)
double fN_pve(double N, double n0, double double a0)

double n_nve(double N, double n0, double a0, double p)
double n_pve(double N, double n0, double a0, double p)

double DH_nve(double N)
double DH_pve(double N)

int
main(void)
{
	/* define constants */

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

double fN_nve(double N, double n0, double double a0)
{
	return fn(n_nve(N, n0, a0, p), n0, a0);
}

double fN_pve(double N, double n0, double double a0)
{
	return fn(n_pve(N, n0, a0, p), n0, a0);
}

double n_nve(double N, double n0, double a0, double p)
{
	return -n0*sqrt(pow(A(N)/a0, 1/p)-1);
}

double n_pve(double N, double n0, double a0, double p)
{
	return n0*sqrt(pow(A(N)/a0, 1/p)-1);
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

/* evolve hk values from Nics to Nshss */

double DDhk_nve(double k, double N, double hk, double Dhk, double n0, double a0, double p)
{
	return -(3*N -(1/N) +(DH_nve(N)/H_nve(N, n0, a0, p)))*Dhk -pow(k/(A(N, a0)*H_nve(N, n0, a0, p)),2)*hk;
}

double DDhk_pve(double k, double N, double hk, double Dhk, double n0, double a0, double p)
{
	return -(3*N -(1/N) +(DH_pve(N)/H_pve(N, n0, a0, p)))*Dhk -pow(k/(A(N, a0)*H_pve(N, n0, a0, p)),2)*hk;
}
