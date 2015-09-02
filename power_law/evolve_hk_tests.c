/* dammit! since the rk4_stepper is returning more than one value,
 we need to define an array that will hold its values or a struct! */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double V(double phi, double V0, double q, double phi0);
double dV(double phi, double V0, double q, double phi0);
double DDphi(double N, double phi, double Dphi, double V0, double q, double phi0);
double rk4_stepper(double N, double phi, double Dphi, double step, double V0, double q, double phi0);

int
main(void)
{
	/* define constants */

	double q, V0, t0;
	double ph0, dphi0;
	double Ni, Nf;
	double H0, Dphi0;

	int npts;
	double step;


	q = 51.;
	V0 = (204./100.)*1e-08;
	t0 = (q*(3.*q -1.)/V0)**(1./2);

	phi0 = 1.;
	dphi0 = (2.*q)**(1./2)/t0;

	Ni = 0.;
	Nf = 70.;

	H0 = ((1./3)*(dphi0**2/2. +V(phi0, V0, q, phi0)))**(1./2.);
	Dphi0 = dphi0/H0;

	/* solve ODE and obtain phi values */

	npts = 250000;
	step = (Nf-Ni)/(npts);

	phi = phi0;
	Dphi = Dphi0;

	phi_array = [];
	Dphi_array = [];
	N_array = [];
	
	inc_array = [];

	N = Ni;
	while (N < Nf)
	{
/*		N_array = numpy.append(N_array, N);
		phi_array = numpy.append(phi_array, phi_);
		Dphi_array = numpy.append(Dphi_array, Dphi_); */
		
/* how to load the arrays so they can be saved for later use? */

		rk4_step(N, phi_, Dphi_, step);
		phi_ = phi_ +inc_array[0];
		Dphi_ = Dphi_ +inc_array[1];

		N += step;
	}
	
	return (0);
}

double V(double phi, double V0, double q, double phi0)
{
	return V0*numpy.exp(-(2./q)**(1./2)*(phi -phi0));
}

double dV(double phi, double q, double V0, double phi0)
{
	return -(2./q)**(1./2)*V0*numpy.exp(-(2./q)**(1./2)*(phi -phi0));
}

double DDphi(double N, double phi, double Dphi, double V0, double q, double phi0)
{
	return -(3 -Dphi**2/2.)*Dphi -(dV(phi, V0, q, phi0)/(2*V(phi, V0, q, phi0)))*(6 -Dphi**2);
}

double rk4_stepper(double N, double phi, double Dphi, double step, double V0, double q, double phi0)
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

	return (F1 +2*F2 +2*F3 +F4)*step/6., (f1 +2*f2 +2*f3 +f4)*step/6.;
	/* (phi, Dphi) update */
}

/* root finding to estimate Nics and Nshss*/

/*initialize hk and Dhk values*/

/* evolve hk values from Nics to Nshss */


