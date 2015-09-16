/* dammit! since the rk4_stepper is returning more than one value,
 we need to define an array that will hold its values or a struct! */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double V(double phi, double V0, double p, double mu);
double dV(double phi, double V0, double p, double mu);
double DDphi(double N, double phi, double Dphi, double V0, double p, double mu);
void rk4_stepper_phi(double N, double phi, double Dphi, double step, double V0, double p, double mu, double *update);

double phi_function(double *phi_array, double N, double Ni, double step);
double Dphi_function(double *Dphi_array, double N, double Ni, double step);
double H(double N, double V0, double q, double phi0, double *phi_array, double *Dphi_array, double Ni, double step);
double DH(double N, double V0, double q, double phi0, double *phi_array, double *Dphi_array, double Ni, double step);

double A(double N, double ai);

double find_Nics(double k, double *N_array, int npts, double ai, double V0, double p, double mu, double *phi_array, double *Dphi_array, double N, double Ni, double step);
double find_Nshss(double k, double *N_array, int npts, double ai, double V0, double p, double mu, double *phi_array, double *Dphi_array, double N, double Ni, double step);

void initialize_hk(double k, double Nics, double ai, double *hk);
void initialize_Dhk(double k, double Nics, double ai, double *Dhk, double N, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step);

double DDhk(double k, double N, double hk, double Dhk, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step, double ai);
void rk4_stepper_hk(double k, double N, double *hk, double *Dhk, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step, double *update_hk, double *update_Dhk, double ai);

int
main(void)
{
	/* define constants */

	FILE *phi_ptr;
	FILE *H_ptr;
	FILE *DH_ptr;
	FILE *eps1_ptr;
	FILE *tps_data_ptr;

	double p, V0;
	double mu, phi0, dphi0;
	double phi, Dphi;
	double Ni, Nf;
	double N;
	double H0, Dphi0;

	double ai;
	double k;

	int npts;
	double tps;

	int i, j;
	double step;

	double increment_phi[2];
	double increment_hk[2];
	double increment_Dhk[2];

	double N_array[100000 +1];
	double phi_array[100000 +1];
	double Dphi_array[100000 +1];

	double H_array[100000 +1];
	double DH_array[100000 +1];
	double eps1_array[100000 +1];

	double hk[2];
	double Dhk[2];
	double DDhk[2];

	double Nics, Nshss;

	phi_ptr = fopen("phi_vs_N_c.txt", "w");
	H_ptr = fopen("H_vs_N_c.txt", "w");
	DH_ptr = fopen("DH_vs_N_c.txt", "w");
	eps1_ptr = fopen("eps_vs_N_c.txt","w");

	tps_data_ptr = fopen("tps_c.txt","w");

	p = 4.0;
	mu = 15.0;
	V0 = (5.55702)*pow(10,-10);

	phi0 = 7.3;
	dphi0 = 8*pow(10, -7);

	Ni = 0.0;
	Nf = 70.0;

	ai = pow(10, -5);
	k = pow(10, -6);

//	H0 = sqrt((1./3)*((dphi0*dphi0)/2 +V(phi0, V0, q, phi0)));
	H0 = sqrt((1.0/3.0)*((dphi0*dphi0)/2.0 +V0));
	Dphi0 = dphi0/H0;

	/* solve ODE and obtain phi values */
	npts = 100000;
	step = (Nf-Ni)/(npts);

	phi = phi0;
	Dphi = Dphi0;

	printf("%lf, %lf, %lf, %lf, %lf \n", Ni, phi, Dphi, dphi0, H0);

	N = Ni;
	j = 0;

	while (N < Nf +step)
	{
		N_array[j] = N;
		phi_array[j] = phi;
		Dphi_array[j] = Dphi;
		
/* how to load the arrays so they can be saved for later use? */

		rk4_stepper_phi(N, phi, Dphi, step, V0, p, mu, increment_phi);
		phi = phi +increment_phi[0];
		Dphi = Dphi +increment_phi[1];

		N += step;
		j += 1;
	}

/*	printf("%lf, %lf, %lf, %lf \n", N, step, phi, Dphi);
	printf("%lf, %lf, %lf, %lf, %lf, %lf \n", N_array[0], N_array[npts], phi_array[0], phi_array[npts], Dphi_array[0], Dphi_array[npts]);
	printf("%lf, %lf, %lf, %lf \n", phi_function(phi_array, Ni, Ni, step), Dphi_function(Dphi_array, Ni, Ni, step), phi_function(phi_array, Nf, Ni, step), Dphi_function(Dphi_array, Nf, Ni, step));
	printf("%d \n", j); */

	for (i=0; i<npts+1; i++)
	{
		H_array[i] = sqrt(V(phi_array[i], V0, p, mu)/(3.0 -Dphi_array[i]*Dphi_array[i]/2.0));
		DH_array[i] = (-1.0/2.0)*sqrt(V(phi_array[i], V0, p, mu)/(3.0 -Dphi_array[i]*Dphi_array[i]/2.0))*Dphi_array[i]*Dphi_array[i];
		eps1_array[i] = Dphi_array[i]*Dphi_array[i]/2.0;
	}

	for (i=0; i<npts+1; i++)
	{
		fprintf(phi_ptr, "%lf , %lf \n", N_array[i], phi_array[i]);
		fprintf(H_ptr, "%lf , %lf \n", N_array[i], H_array[i]);
		fprintf(DH_ptr, "%lf , %lf \n", N_array[i], DH_array[i]);
		fprintf(eps1_ptr, "%lf , %lf \n", N_array[i], eps1_array[i]);
	}

	while (k < pow(10,0))
	{
		Nics = find_Nics(k, N_array, npts, ai, V0, p, mu, phi_array, Dphi_array, N, Ni, step);
		Nshss = find_Nshss(k, N_array, npts, ai, V0, p, mu, phi_array, Dphi_array, N, Ni, step);

		printf("===================================");
		printf("%le, %lf, %lf \n", k, Nics, Nshss);

		initialize_hk(k, Nics, ai, hk);
		initialize_Dhk(k, Nics, ai, Dhk, Nics, V0, p, mu, phi_array, Dphi_array, Ni, step);

		printf("%lf, %lf \n", hk[0], hk[1]);
		printf("%lf, %lf \n", Dhk[0], Dhk[1]);

/*		printf("%lf \n", Nics);
		printf("%lf \n", Nshss);
		printf("%lf \n", find_Nics(pow(10,-5), N_array, npts, ai, V0, q, phi0, phi_array, Dphi_array, N, Ni, step));
		printf("%lf \n", find_Nshss(pow(10,-5), N_array, npts, ai, V0, q, phi0, phi_array, Dphi_array, N, Ni, step));
		printf("%lf \n", find_Nics(pow(10,-4), N_array, npts, ai, V0, q, phi0, phi_array, Dphi_array, N, Ni, step));
		printf("%lf \n", find_Nshss(pow(10,-4), N_array, npts, ai, V0, q, phi0, phi_array, Dphi_array, N, Ni, step)); */

		N = Nics;
		while (N < Nshss +step)
		{
			rk4_stepper_hk(k, N, hk, Dhk, V0, p, mu, phi_array, Dphi_array, Ni, step, increment_hk, increment_Dhk, ai);
			hk[0] += increment_hk[0];
			hk[1] += increment_hk[1];
			Dhk[0] += increment_Dhk[0];
			Dhk[1] += increment_Dhk[1];

			N += step;
		}

		printf("%lf, %lf, %le, %le \n", hk[0], hk[1], Dhk[0], Dhk[1]);

		tps = 2*pow(k,3)/(2*pow(M_PI,2))*(hk[0]*hk[0] +hk[1]*hk[1]);
		printf("%le \n", tps);

		fprintf(tps_data_ptr, "%le, %lf, %lf, %le \n", k, Nics, Nshss, tps);

		k = pow(10,1./2)*k;
		printf("\n");
	}

	fclose(phi_ptr);
	fclose(H_ptr);
	fclose(DH_ptr);
	fclose(eps1_ptr);

	fclose(tps_data_ptr);

	return (0);
}

double V(double phi, double V0, double p, double mu)
{
	return V0*(1.0 -pow(phi/mu, p));
}

double dV(double phi, double V0, double p, double mu)
{
	return -V0*p*pow(phi/mu, p-1)/mu;
}

double DDphi(double N, double phi, double Dphi, double V0, double p, double mu)
{
	return -(3. -(Dphi*Dphi)/2.)*Dphi -(dV(phi, V0, p, mu)/(2.*V(phi, V0, p, mu)))*(6. -(Dphi*Dphi));
}

void rk4_stepper_phi(double N, double phi, double Dphi, double step, double V0, double p, double mu, double *update_phi)
{
	double F1, f1, F2, f2, F3, f3, F4, f4;

	F1 = Dphi;
	f1 = DDphi(N, phi, Dphi, V0, p, mu);
	F2 = Dphi +f1*step/2.0;
	f2 = DDphi(N +step/2.0, phi +F1*step/2.0, Dphi +f1*step/2.0, V0, p, mu);
	F3 = Dphi +f2*step/2.0;
	f3 = DDphi(N +step/2.0, phi +F2*step/2.0, Dphi +f2*step/2.0, V0, p, mu);
	F4 = Dphi +f3*step;
	f4 = DDphi(N +step, phi +F3*step, Dphi +f3*step, V0, p, mu);

	update_phi[0] = (F1 +2*F2 +2*F3 +F4)*step/6.0;
	update_phi[1] = (f1 +2*f2 +2*f3 +f4)*step/6.0;

	return; 
	/* (phi, Dphi) update */
}

double phi_function(double *phi_array, double N, double Ni, double step)
{
	int index;
	index = floor((N-Ni)/step);
	return phi_array[index];
}

double Dphi_function(double *Dphi_array, double N, double Ni, double step)
{
	int index;
	index = floor((N-Ni)/step);
	return Dphi_array[index];
}

double A(double N, double ai)
{
	return ai*exp(N);
}

double H(double N, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step)
{
	return sqrt(V(phi_function(phi_array, N, Ni, step), V0, p, mu)/(3 -Dphi_function(Dphi_array, N, Ni, step)*Dphi_function(Dphi_array, N, Ni, step)/2));
}

double DH(double N, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step)
{
	return (-1.0/2)*H(N, V0, p, mu, phi_array, Dphi_array, Ni, step)*Dphi_function(Dphi_array, N, Ni, step)*Dphi_function(Dphi_array, N, Ni, step);
}

/* root finding to estimate Nics and Nshss */

double find_Nics(double k, double *N_array, int npts, double ai, double V0, double p, double mu, double *phi_array, double *Dphi_array, double N, double Ni, double step)
{
	double test_array[npts+1];
	int i, j;
	double min;

	for (i=0; i<npts+2; i++)
	{
		test_array[i] = fabs(k -pow(10,2)*A(N_array[i], ai)*H(N_array[i], V0, p, mu, phi_array, Dphi_array, Ni, step));
	}

//	printf("%lf, %lf\n", test_array[0], test_array[npts+1]);

	min = test_array[0];
	j = 0;
	for (i=1; i<npts+2; i++)
	{
		if (test_array[i] < min)
		{
			min = test_array[i];
			j = i;
		}
	}

//	printf("%lf, %lf \n", min, test_array[j]);

	return N_array[j];
}


double find_Nshss(double k, double *N_array, int npts, double ai, double V0, double p, double mu, double *phi_array, double *Dphi_array, double N, double Ni, double step)
{
	double test_array[npts+1];
	int i, j;
	double min;

	for (i=0; i<npts+2; i++)
	{
		test_array[i] = fabs(k -pow(10,-5)*A(N_array[i], ai)*H(N_array[i], V0, p, mu, phi_array, Dphi_array, Ni, step));
	}

//	printf("%lf, %lf\n", test_array[0], test_array[npts+1]);

	min = test_array[0];
	j = 0;
	for (i=1; i<npts+2; i++)
	{
		if (test_array[i] < min)
		{
			min = test_array[i];
			j = i;
		}
	}

//	printf("%lf, %lf \n", min, test_array[j]);

	return N_array[j];
}

/*initialize hk and Dhk values*/

void initialize_hk(double k, double Nics, double ai, double *hk)
{
	hk[0] = 1/sqrt(2*k)/A(Nics, ai);
	hk[1] = 0;
	return;	
}

void initialize_Dhk(double k, double Nics, double ai, double *Dhk, double N, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step)
{
	Dhk[0] = -1/sqrt(2*k)/A(Nics, ai);
	Dhk[1] = -sqrt(k/2)/(A(Nics, ai)*A(Nics, ai)*H(Nics, V0, p, mu, phi_array, Dphi_array, Ni, step));
	return;
}

/* evolve hk values from Nics to Nshss */

double DDhk(double k, double N, double hk, double Dhk, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step, double ai)
{
	double DDhk_val;
	DDhk_val = -(3 +(DH(N, V0, p, mu, phi_array, Dphi_array, Ni, step)/H(N, V0, p, mu, phi_array, Dphi_array, Ni, step)))*(Dhk)
			-pow(k/(A(N, ai)*H(N, V0, p, mu, phi_array, Dphi_array, Ni, step)),2)*hk;
	return DDhk_val;
}

void rk4_stepper_hk(double k, double N, double *hk, double *Dhk, double V0, double p, double mu, double *phi_array, double *Dphi_array, double Ni, double step, double *update_hk, double *update_Dhk, double ai)
{
	double F1_real, f1_real, F2_real, f2_real, F3_real, f3_real, F4_real, f4_real;
	double F1_imag, f1_imag, F2_imag, f2_imag, F3_imag, f3_imag, F4_imag, f4_imag;

	F1_real = Dhk[0];
	f1_real = DDhk(k, N, hk[0], Dhk[0], V0, p, mu, phi_array, Dphi_array,  Ni, step, ai);
	F2_real = Dhk[0] +f1_real*step/2.;
	f2_real = DDhk(k, N +step/2., hk[0] +F1_real*step/2., Dhk[0] +f1_real*step/2., V0, p, mu, phi_array, Dphi_array, Ni, step, ai);
	F3_real = Dhk[0] +f2_real*step/2.;
	f3_real = DDhk(k, N +step/2., hk[0] +F2_real*step/2., Dhk[0] +f2_real*step/2., V0, p, mu, phi_array, Dphi_array, Ni, step, ai);
	F4_real = Dhk[0] +f3_real*step;
	f4_real = DDhk(k, N +step, hk[0] +F3_real*step, Dhk[0] +f3_real*step, V0, p, mu, phi_array, Dphi_array, Ni, step, ai);

	F1_imag = Dhk[1];
	f1_imag = DDhk(k, N, hk[1], Dhk[1], V0, p, mu, phi_array, Dphi_array, Ni, step, ai);
	F2_imag = Dhk[1] +f1_imag*step/2.;
	f2_imag = DDhk(k, N +step/2., hk[1] +F1_imag*step/2., Dhk[1] +f1_imag*step/2., V0, p, mu, phi_array, Dphi_array, Ni, step, ai);
	F3_imag = Dhk[1] +f2_imag*step/2.;
	f3_imag = DDhk(k, N +step/2., hk[1] +F2_imag*step/2., Dhk[1] +f2_imag*step/2., V0, p, mu, phi_array, Dphi_array, Ni, step, ai);
	F4_imag = Dhk[1] +f3_imag*step;
	f4_imag = DDhk(k, N +step, hk[1] +F3_imag*step, Dhk[1] +f3_imag*step, V0, p, mu, phi_array, Dphi_array, Ni, step, ai);

	update_hk[0] = (F1_real +2*F2_real +2*F3_real +F4_real)*step/6.;
	update_Dhk[0] = (f1_real +2*f2_real +2*f3_real +f4_real)*step/6.;
	update_hk[1] = (F1_imag +2*F2_imag +2*F3_imag +F4_imag)*step/6.;
	update_Dhk[1] = (f1_imag +2*f2_imag +2*f3_imag +f4_imag)*step/6.;

	return;
}
