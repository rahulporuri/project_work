/* dammit! since the rk4_stepper is returning more than one value,
 we need to define an array that will hold its values or a struct! */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double V(double phi, double V0, double q, double phi0);
double dV(double phi, double V0, double q, double phi0);
double DDphi(double N, double phi, double Dphi, double V0, double q, double phi0);
void rk4_stepper_phi(double N, double phi, double Dphi, double step, double V0, double q, double phi0, double *update);

double phi_function(double *phi_array, double N, double Ni, double step);
double Dphi_function(double *Dphi_array, double N, double Ni, double step);
double H(double N, double Ni, double step, double *H_array);
double DH(double N, double Ni, double step, double *DH_array);

double A(double N, double ai);

double find_Nics(double k, double *N_array, int npts, double ai, double Ni, double step, double *H_array);
double find_Nshss(double k, double *N_array, int npts, double ai, double Ni, double step, double *H_array);

void initialize_hk(double k, double Nics, double ai, double *hk);
void initialize_Dhk(double k, double Nics, double ai, double *Dhk, double Ni, double step, double *H_array);

double DDhk(double k, double N, double hk, double Dhk, double Ni, double step, double ai, double *H_array, double *DH_array);
void rk4_stepper_hk(double k, double N, double *hk, double *Dhk, double Ni, double step, double *update_hk, double *update_Dhk, double ai, double *H_array, double *DH_array);

void evolve_hk(double k, double Nics, double Nshss, double ai, double Ni, double step, double *H_array, double *DH_array, double **hk_array);

//void calG(double *CalG);
//void calG_cc(double *CalG_cc);

int
main(void)
{
	/* define constants */

	FILE *tps_data_ptr;

	double q, V0, t0;
	double phi0, dphi0;
	double phi, Dphi;
	double Ni, Nf;
	double N;
	double H0, Dphi0;

	double ai;
	double k;

 	double k1, k2, k3;

	int npts;
	double tps;
	double tps_k1, tps_k2, tps_k3;

	int i, j;
	double step;

	double increment_phi[2];



	double *N_array = malloc(10000001*sizeof(double));
	double *phi_array = malloc(10000001*sizeof(double));
	double *Dphi_array = malloc(10000001*sizeof(double));

	double *H_array = malloc(10000001*sizeof(double));
	double *DH_array = malloc(10000001*sizeof(double));
	double *eps1_array = malloc(10000001*sizeof(double));



	double DDhk[2];

	double Nics, Nshss;

	double CalG[2];
	double CalG_cc[2];

	double **hk_k1_array = malloc(10000001*sizeof(double));
	for (i=0; i<10000001; i++)
	{
		hk_k1_array[i] = (double*) malloc(2*sizeof(double));
	}

	double **hk_k2_array = malloc(10000001*sizeof(double));
	for (i=0; i<10000001; i++)
	{
		hk_k2_array[i] = (double*) malloc(2*sizeof(double));
	}

	double **hk_k3_array = malloc(10000001*sizeof(double));
	for (i=0; i<10000001; i++)
	{
		hk_k3_array[i] = (double*) malloc(2*sizeof(double));
	}

	double *N_int_range = malloc(10000001*sizeof(double));

	int size_hk_array;

//	tps_data_ptr = fopen("data_files/tps_c.txt","w");

	q = 51.0;
	V0 = (204.0/100.0)*pow(10,-8);
	printf("%le, %le \n", 1e-8, (204/100)*1e-8);
	t0 = sqrt(q*(3.0*q -1.0)/V0);

	phi0 = 1.0;
	dphi0 = sqrt(2.0*q)/t0;

	Ni = 0.0;
	Nf = 70.0;

	ai = pow(10, -5);
	k = 5*pow(10, -4);

	H0 = sqrt((1.0/3.0)*((dphi0*dphi0)/2.0 +V0));
	Dphi0 = dphi0/H0;

	/* solve ODE and obtain phi values */
	npts = 10000000;
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
		
/*		how to load the arrays so they can be saved for later use? */

		rk4_stepper_phi(N, phi, Dphi, step, V0, q, phi0, increment_phi);
		phi = phi +increment_phi[0];
		Dphi = Dphi +increment_phi[1];

		N += step;
		j += 1;
	}

	for (i=0; i<npts+1; i++)
	{
		H_array[i] = sqrt(V(phi_array[i], V0, q, phi0)/(3.0 -Dphi_array[i]*Dphi_array[i]/2.0));
		DH_array[i] = (-1.0/2.0)*sqrt(V(phi_array[i], V0, q, phi0)/(3.0 -Dphi_array[i]*Dphi_array[i]/2.0))*Dphi_array[i]*Dphi_array[i];
		eps1_array[i] = Dphi_array[i]*Dphi_array[i]/2.0;
	}

/*	Nics = find_Nics(5*pow(10,-4), N_array, npts, ai, Ni, step);
	Nshss = find_Nshss(5*pow(10,-2), N_array, npts, ai, Ni, step); */

	Nics = 8.880899999991;
	Nshss = 30.018590000083;

	size_hk_array = floor((Nshss-Nics)/step);

	while (k < 5*pow(10,-1))
	{
		printf("=================================== \t");
		printf("%le \n", k);
		k1 = k; k2 = k; k3 = k;

		evolve_hk(k1, Nics, Nshss, ai, Ni, step, H_array, DH_array, hk_k1_array);
		evolve_hk(k2, Nics, Nshss, ai, Ni, step, H_array, DH_array, hk_k2_array);
		evolve_hk(k3, Nics, Nshss, ai, Ni, step, H_array, DH_array, hk_k3_array);

		tps = 2*pow(k1,3)/(2*pow(M_PI,2))*(hk_k1_array[size_hk_array][0]*hk_k1_array[size_hk_array][0] +hk_k1_array[size_hk_array][1]*hk_k1_array[size_hk_array][1]);
		tps = 2*pow(k2,3)/(2*pow(M_PI,2))*(hk_k2_array[size_hk_array][0]*hk_k2_array[size_hk_array][0] +hk_k2_array[size_hk_array][1]*hk_k2_array[size_hk_array][1]);
		tps = 2*pow(k3,3)/(2*pow(M_PI,2))*(hk_k3_array[size_hk_array][0]*hk_k3_array[size_hk_array][0] +hk_k3_array[size_hk_array][1]*hk_k3_array[size_hk_array][1]);

		printf("%le, %le \n", k1, tps_k1);
		printf("%le, %le \n", k2, tps_k2);
		printf("%le, %le \n", k3, tps_k3);

//		fprintf(tps_data_ptr, "%le, %lf, %lf, %le \n", k, Nics, Nshss, tps);

		k = pow(10,1./2)*k;
	}

//	fclose(tps_data_ptr);

	free(N_array);
	free(phi_array);
	free(Dphi_array);

	free(H_array);
	free(DH_array);
	free(eps1_array);

/*	for(i=0; i<10000001; i++)
	{
		free(hk_k1_array[i]);
	}
	free(hk_k1_array);

	for(i=0; i<10000001; i++)
	{
		free(hk_k1_array[i]);
	}
	free(H_array);

	for(i=0; i<10000001; i++)
	{
		free(hk_k1_array[i]);
	}
	free(H_array);*/

	free(N_int_range);

	return (0);
}

double V(double phi, double V0, double q, double phi0)
{
	return V0*exp(-sqrt(2.0/q)*(phi -phi0));
}

double dV(double phi, double V0, double q, double phi0)
{
	return -sqrt(2/q)*V0*exp(-sqrt(2/q)*(phi -phi0));
}

double DDphi(double N, double phi, double Dphi, double V0, double q, double phi0)
{
	return -(3 -(Dphi*Dphi)/2.)*Dphi -(dV(phi, V0, q, phi0)/(2*V(phi, V0, q, phi0)))*(6 -(Dphi*Dphi));
}

void rk4_stepper_phi(double N, double phi, double Dphi, double step, double V0, double q, double phi0, double *update_phi)
{
	double F1, f1, F2, f2, F3, f3, F4, f4;

	F1 = Dphi;
	f1 = DDphi(N, phi, Dphi, V0, q, phi0);
	F2 = Dphi +f1*step/2.0;
	f2 = DDphi(N +step/2.0, phi +F1*step/2.0, Dphi +f1*step/2.0, V0, q, phi0);
	F3 = Dphi +f2*step/2.0;
	f3 = DDphi(N +step/2.0, phi +F2*step/2.0, Dphi +f2*step/2.0, V0, q, phi0);
	F4 = Dphi +f3*step;
	f4 = DDphi(N +step, phi +F3*step, Dphi +f3*step, V0, q, phi0);

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

double H(double N, double Ni, double step, double *H_array)
{
	int index;
	index = floor((N-Ni)/step);
	return H_array[index];
//	return sqrt(V(phi_function(phi_array, N, Ni, step), V0, q, phi0)/(3 -Dphi_function(Dphi_array, N, Ni, step)*Dphi_function(Dphi_array, N, Ni, step)/2));
}

double DH(double N, double Ni, double step, double *DH_array)
{
	int index;
	index = floor((N-Ni)/step);
	return DH_array[index];
//	return (-1.0/2)*H(N, V0, q, phi0, phi_array, Dphi_array, Ni, step)*Dphi_function(Dphi_array, N, Ni, step)*Dphi_function(Dphi_array, N, Ni, step);
}

/* root finding to estimate Nics and Nshss */

double find_Nics(double k, double *N_array, int npts, double ai, double Ni, double step, double *H_array)
{
	double test_array[npts+1];
	int i, j;
	double min;

	for (i=0; i<npts+2; i++)
	{
		test_array[i] = fabs(k -pow(10,2)*A(N_array[i], ai)*H(N_array[i], Ni, step, H_array));
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


double find_Nshss(double k, double *N_array, int npts, double ai, double Ni, double step, double *H_array)
{
	double test_array[npts+1];
	int i, j;
	double min;

	for (i=0; i<npts+2; i++)
	{
		test_array[i] = fabs(k -pow(10,-5)*A(N_array[i], ai)*H(N_array[i], Ni, step, H_array));
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

void initialize_Dhk(double k, double Nics, double ai, double *Dhk, double Ni, double step, double *H_array)
{
	Dhk[0] = -1/sqrt(2*k)/A(Nics, ai);
	Dhk[1] = -sqrt(k/2)/(A(Nics, ai)*A(Nics, ai)*H(Nics, Ni, step, H_array));
	return;
}

/* evolve hk values from Nics to Nshss */

double DDhk(double k, double N, double hk, double Dhk, double Ni, double step, double ai, double *H_array, double *DH_array)
{
	double DDhk_val;
	DDhk_val = -(3 +(DH(N, Ni, step, DH_array)/H(N, Ni, step, H_array)))*(Dhk)
			-pow(k/(A(N, ai)*H(N, Ni, step, H_array)),2)*hk;
	return DDhk_val;
}

void rk4_stepper_hk(double k, double N, double *hk, double *Dhk, double Ni, double step, double *update_hk, double *update_Dhk, double ai, double *H_array, double *DH_array)
{
	double F1_real, f1_real, F2_real, f2_real, F3_real, f3_real, F4_real, f4_real;
	double F1_imag, f1_imag, F2_imag, f2_imag, F3_imag, f3_imag, F4_imag, f4_imag;

	F1_real = Dhk[0];
	f1_real = DDhk(k, N, hk[0], Dhk[0], Ni, step, ai, H_array, DH_array);
	F2_real = Dhk[0] +f1_real*step/2.;
	f2_real = DDhk(k, N +step/2., hk[0] +F1_real*step/2., Dhk[0] +f1_real*step/2., Ni, step, ai, H_array, DH_array);
	F3_real = Dhk[0] +f2_real*step/2.;
	f3_real = DDhk(k, N +step/2., hk[0] +F2_real*step/2., Dhk[0] +f2_real*step/2., Ni, step, ai, H_array, DH_array);
	F4_real = Dhk[0] +f3_real*step;
	f4_real = DDhk(k, N +step, hk[0] +F3_real*step, Dhk[0] +f3_real*step, Ni, step, ai, H_array, DH_array);

	F1_imag = Dhk[1];
	f1_imag = DDhk(k, N, hk[1], Dhk[1], Ni, step, ai, H_array, DH_array);
	F2_imag = Dhk[1] +f1_imag*step/2.;
	f2_imag = DDhk(k, N +step/2., hk[1] +F1_imag*step/2., Dhk[1] +f1_imag*step/2., Ni, step, ai, H_array, DH_array);
	F3_imag = Dhk[1] +f2_imag*step/2.;
	f3_imag = DDhk(k, N +step/2., hk[1] +F2_imag*step/2., Dhk[1] +f2_imag*step/2., Ni, step, ai, H_array, DH_array);
	F4_imag = Dhk[1] +f3_imag*step;
	f4_imag = DDhk(k, N +step, hk[1] +F3_imag*step, Dhk[1] +f3_imag*step, Ni, step, ai, H_array, DH_array);

	update_hk[0] = (F1_real +2*F2_real +2*F3_real +F4_real)*step/6.;
	update_Dhk[0] = (f1_real +2*f2_real +2*f3_real +f4_real)*step/6.;
	update_hk[1] = (F1_imag +2*F2_imag +2*F3_imag +F4_imag)*step/6.;
	update_Dhk[1] = (f1_imag +2*f2_imag +2*f3_imag +f4_imag)*step/6.;

	return;
}

void evolve_hk(double k, double Nics, double Nshss, double ai, double Ni, double step, double *H_array, double *DH_array, double **hk_array)
{
	double hk[2];
	double Dhk[2];

	double increment_hk[2];
	double increment_Dhk[2];

	double tps;

	double N;
	int i;

//	Nics = find_Nics(k, N_array, npts, ai, V0, q, phi0, phi_array, Dphi_array, N, Ni, step);
//	Nshss = find_Nshss(k, N_array, npts, ai, V0, q, phi0, phi_array, Dphi_array, N, Ni, step);

	initialize_hk(k, Nics, ai, hk);
	initialize_Dhk(k, Nics, ai, Dhk, Ni, step, H_array);

	i = 0;
	N = Nics;
	while (N < Nshss +step)
	{
		hk_array[i][0] = hk[0];
		hk_array[i][1] = hk[1];

		rk4_stepper_hk(k, N, hk, Dhk, Ni, step, increment_hk, increment_Dhk, ai, H_array, DH_array);
		hk[0] += increment_hk[0];
		hk[1] += increment_hk[1];
		Dhk[0] += increment_Dhk[0];
		Dhk[1] += increment_Dhk[1];

		N += step;
		i += 1;
	}

/*	tps = 2*pow(k,3)/(2*pow(M_PI,2))*(hk[0]*hk[0] +hk[1]*hk[1]);
	printf("=================== \n");
	printf("%le, %le \n", k, tps);
	printf("size : %d", i-1); */

	return;
}

/*

void calG(double k1, double k2, double k3, double *hk_k1_array, double *hk_k2_array, double *hk_k3_array, double *N_int_range, double *CalG, double V0, double q, double phi0, double *phi_array, double *Dphi_array, double Ni, double step, double ai)
{
	double *func_int_real = malloc(10000001*sizeof(double));
	double *func_int_imag = malloc(10000001*sizeof(double));

	double e;
	e = 1/50;

	for(i=0;i<N;i++)
	{
		func_int_real[i] = (A(N_int_range[i], ai)/H(N_int_range[i, V0, q, phi0, phi_array, Dphi_array, Ni, step]))
					*hk_k1_array[i][0]*hk_k2_array[i][0]*hk_k3_array[i][0]
					*exp(e*(k1+k2+k3)/(3*A(N_int_range[i], ai)*H(N_int_range[i], V0, q, phi0, phi_array, Dphi_array, Ni, step)));

		func_int_imag[i] = (A(N_int_range[i], ai)/H(N_int_range[i], V0, q, phi0, phi_array, Dphi_array, Ni, step))
					*(-hk_k1_array[i][1])*(-hk_k2_array[i][1])*(-hk_k3_array[i][1])
					*exp(e*(k1+k2+k3)/(3*A(N_int_range[i], ai)*H(N_int_range[i], V0, q, phi0, phi_array, Dphi_array, Ni, step)));
	}

	return;
}

void calG_cc(double *hk_array, double *N_array, double *CalG_cc)
{
	return;
}

*/
