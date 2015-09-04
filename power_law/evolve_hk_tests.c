/* dammit! since the rk4_stepper is returning more than one value,
 we need to define an array that will hold its values or a struct! */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double V(double phi, double V0, double q, double phi0);
double dV(double phi, double V0, double q, double phi0);
double DDphi(double N, double phi, double Dphi, double V0, double q, double phi0);
void rk4_stepper(double N, double phi, double Dphi, double step, double V0, double q, double phi0, double *update);

int
main(void)
{
	/* define constants */

	double q, V0, t0;
	double phi0, dphi0;
	double phi, Dphi;
	double Ni, Nf;
	double N;
	double H0, Dphi0;

	int npts;
	double step;

	double increment[2];

	q = 51.0;
	V0 = (204/100)*pow(10,-8);
	printf("%le, %le \n", 1e-8, (204/100)*1e-8);
	t0 = sqrt(q*(3*q -1)/V0);

	phi0 = 1.0;
	dphi0 = sqrt(2*q)/t0;

	Ni = 0.0;
	Nf = 70.0;

//	H0 = sqrt((1./3)*((dphi0*dphi0)/2 +V(phi0, V0, q, phi0)));
	H0 = sqrt((1./3)*((dphi0*dphi0)/2 +V0));
	Dphi0 = dphi0/H0;

	/* solve ODE and obtain phi values */

	npts = 1000000;
	step = (Nf-Ni)/(npts);

	phi = phi0;
	Dphi = Dphi0;

/*	phi_array = [];
	Dphi_array = [];
	N_array = []; */
	
	printf("%lf, %lf, %lf, %lf, %lf \n", Ni, phi, Dphi, dphi0, H0);

	N = Ni;

//	printf("%llf, %llf, %llf, %llf \n", V(phi0, V0, q, phi0), V0*exp(-sqrt(2/q)*0), V0, V0*exp(0));

	while (N < Nf)
	{

/*		N_array = numpy.append(N_array, N);
		phi_array = numpy.append(phi_array, phi_);
		Dphi_array = numpy.append(Dphi_array, Dphi_); */
		
/* how to load the arrays so they can be saved for later use? */

		rk4_stepper(N, phi, Dphi, step, V0, q, phi0, increment);
		phi = phi +increment[0];
		Dphi = Dphi +increment[1];

		N += step;
	}

	printf("%lf, %lf, %lf, %lf \n", N, step, phi, Dphi);
	
	return (0);
}

double V(double phi, double V0, double q, double phi0)
{
	return V0*exp(-sqrt(2./q)*(phi -phi0));
}

double dV(double phi, double V0, double q, double phi0)
{
	return -sqrt(2/q)*V0*exp(-sqrt(2/q)*(phi -phi0));
}

double DDphi(double N, double phi, double Dphi, double V0, double q, double phi0)
{
	return -(3 -(Dphi*Dphi)/2.)*Dphi -(dV(phi, V0, q, phi0)/(2*V(phi, V0, q, phi0)))*(6 -(Dphi*Dphi));
}

void rk4_stepper(double N, double phi, double Dphi, double step, double V0, double q, double phi0, double *update)
{
	double F1, f1, F2, f2, F3, f3, F4, f4;

	F1 = Dphi;
	f1 = DDphi(N, phi, Dphi, V0, q, phi0);
	F2 = Dphi +f1*step/2.;
	f2 = DDphi(N +step/2., phi +F1*step/2., Dphi +f1*step/2., V0, q, phi0);
	F3 = Dphi +f2*step/2.;
	f3 = DDphi(N +step/2., phi +F2*step/2., Dphi +f2*step/2., V0, q, phi0);
	F4 = Dphi +f3*step;
	f4 = DDphi(N +step, phi +F3*step, Dphi +f3*step, V0, q, phi0);

	update[0] = (F1 +2*F2 +2*F3 +F4)*step/6.;
	update[1] = (f1 +2*f2 +2*f3 +f4)*step/6.;

	return; 
	/* (phi, Dphi) update */
}

/* root finding to estimate Nics and Nshss*/

/*initialize hk and Dhk values*/

/* evolve hk values from Nics to Nshss */


