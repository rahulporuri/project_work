/*
Note that in this code, I use the prefix 'd' to represent derivative with respect to time 
(except for the case of dV where the derivative is with respect to phi) and the prefix 'D' 
to represent derivative with respect to e-fold N. Also, the suffix '0' is used to represent 
the initial conditions in various cases. Also, as can be seen here, 
we evaluate the scalar field in the e-fold N range Ni to Nf.
*/

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

void evolve_hk(double k, int npts, double ai, double Ni, double step, double *N_array, double *H_array, double *DH_array, double **hk_array, int size_hk_array, double Nics, double Nshss);

void calG(double k1, double k2, double k3, int size_hk_array, double Nics, double Nshss, double **hk_k1_array, double **hk_k2_array, double **hk_k3_array, double *CalG, double Ni, double step, double ai, double *H_array);
void calG_cc(double k1, double k2, double k3, int size_hk_array, double Nics, double Nshss, double **hk_k1_array, double **hk_k2_array, double **hk_k3_array, double *CalG, double Ni, double step, double ai, double *H_array);

void G_func(double *G, double a, double Ib, double c, double Id, double e, double If, double *CalG, double *CalG_cc);

int
main(void)
{
	double q, V0, t0;
	double phi0, dphi0;
	double phi, Dphi;
	double Ni, Nf;
	double N;
	double H0, Dphi0;

	double ai;
 	double k; /* the triangular set of modes*/
	int npts; /* number of points of integration need to be > 10^8! */
	/* otherwise, larger modes will die off faster and becomes NaN at Nshss */
	double tps; /* tensor power spectrum values for the three modes k1, k2 and k3*/

	int i, j;
	double step; /* step size for integration*/

	double increment_phi[2]; /* increment_phi[0] is the phi increment and increment_phi[1] is the Dphi increment */

	/* hold the values for N, phi and Dphi which are needed to estimate */
	double *N_array = malloc(10000001*sizeof(double));
	double *phi_array = malloc(10000001*sizeof(double));
	double *Dphi_array = malloc(10000001*sizeof(double));
	/* the following arrays containing H and DH corresponding to the above N */
	double *H_array = malloc(10000001*sizeof(double));
	double *DH_array = malloc(10000001*sizeof(double));

	/* limits of hk integration*/
	double Nics;
	double Nshss;

	/* stores the real and imaginary values of the second derivative of hk during integration */
	double DDhk[2];

	/* stores the value of script G and it's complex conjugate which are used to estimate */
	double CalG[2];
	double CalG_cc[2];
	/* the tensor bi-spectrum G and the non gaussianity parameter h_NL*/
	double G[2];
	double h_NL;

	/* the following arrays hold the real and imaginary values of hk
	corresponding to the three modes k1, k2, k3 - which are used to estimate calG and calG_cc.
	dynamic memory allocation is done because static arrays don't support 10^8 data points! */
	double **hk_array = malloc(10000001*sizeof(double));
	for (i=0; i<10000001; i++)
	{
		hk_array[i] = (double*) malloc(2*sizeof(double));
	}

	int size_hk_array;

	/* define initial conditions */
	q = 51.0;
	V0 = (204.0/100.0)*pow(10,-8);
	printf("%le, %le \n", 1e-8, (204/100)*1e-8);
	t0 = sqrt(q*(3.0*q -1.0)/V0);

	phi0 = 1.0;
	dphi0 = sqrt(2.0*q)/t0;

	Ni = 0.0;
	Nf = 70.0;

	ai = pow(10, -5);

	H0 = sqrt((1.0/3.0)*((dphi0*dphi0)/2.0 +V0));
	Dphi0 = dphi0/H0;

//	printf("done setting initial params \n");

	/* solve ODE and obtain phi values */
	npts = 10000000;
	step = (Nf-Ni)/(npts);

	phi = phi0;
	Dphi = Dphi0;

//	printf("%lf, %lf, %lf, %lf, %lf \n", Ni, phi, Dphi, dphi0, H0);

	N = Ni;
	j = 0;

	while (N < Nf +step)
	{
		N_array[j] = N;
		phi_array[j] = phi;
		Dphi_array[j] = Dphi;

		rk4_stepper_phi(N, phi, Dphi, step, V0, q, phi0, increment_phi);
		phi = phi +increment_phi[0];
		Dphi = Dphi +increment_phi[1];

		N += step;
		j += 1;
	}

//	printf("done with phi evolution\n");

	for (i=0; i<npts+1; i++)
	{
		H_array[i] = sqrt(V(phi_array[i], V0, q, phi0)/(3.0 -Dphi_array[i]*Dphi_array[i]/2.0));
		DH_array[i] = (-1.0/2.0)*sqrt(V(phi_array[i], V0, q, phi0)/(3.0 -Dphi_array[i]*Dphi_array[i]/2.0))*Dphi_array[i]*Dphi_array[i];
	}

//	printf("done with loading H_arrays \n");

//	Nics = 2.5;
//	Nshss = 32;

//	size_hk_array = floor((Nshss-Nics)/step);

/*	k = pow(10,-6);
	while(k < pow(10,0))
	{
		Nics = find_Nics(k, N_array, npts, ai, Ni, step, H_array);
		Nshss = find_Nshss(k, N_array, npts, ai, Ni, step, H_array);
		printf("%lf, %lf, %d \n", Nics, Nshss, size_hk_array);
		k = k*pow(10,1);
	}
*/	

	k = pow(10,-6);
	/* generate values k2 and k3 and evolve the hk for the corresponding modes */
	while (k < pow(10,-2))
	{
		evolve_hk(k, npts, ai, Ni, step, N_array, H_array, DH_array, hk_array, size_hk_array, Nics, Nshss);
		printf("now to get tps %lf, %lf, %d \n", Nics, Nshss, size_hk_array);
		tps = 2*pow(k,3)/(2*pow(M_PI,2))*(hk_array[size_hk_array-1][0]*hk_array[size_hk_array-1][0] +hk_array[size_hk_array-1][1]*hk_array[size_hk_array-1][1]);

/*		* now that we have hk corresponding to k1, k2 and k3; we estimate the value of script G and it's complex conjugate *
		calG(k, k, k, size_hk_array, Nics, Nshss, hk_array, hk_array, hk_array, CalG, Ni, step, ai, H_array);
		calG_cc(k, k, k, size_hk_array, Nics, Nshss, hk_array, hk_array, hk_array, CalG_cc, Ni, step, ai, H_array);

		* and using the above, calculate the tensor bi-spectrum and *
		G_func(G, hk_array[size_hk_array-1][0], hk_array[size_hk_array-1][1],
				hk_array[size_hk_array-1][0], hk_array[size_hk_array-1][1],
				hk_array[size_hk_array-1][0], hk_array[size_hk_array-1][1], CalG, CalG_cc);

		* the non-gaussianity parameter h_NL *
		h_NL = (-(4/(2*M_PI*M_PI))*(4/(2*M_PI*M_PI))*(k*k*k*k*k*k*k*k*k)*G[0]/
			(2*k*k*k*tps*tps +2*k*k*k*tps*tps +2*k*k*k*tps*tps));

//		printf("%lf, %lf \n", hk_array[size_hk_array][0], hk_array[size_hk_array][1]);
		printf("%le, %le, %le, %le, %le, %lf \n", k, tps, CalG[0]*CalG[0] +CalG[1]*CalG[1], G[0], (k*k*k*k*k*k)*G[0], h_NL);*/
		printf("%le, %lf, %lf, %lf, %le \n", k, Nics, Nics +size_hk_array*step, Nshss, tps);
		k = k*pow(10,1);
	}

	free(N_array);
	free(phi_array);
	free(Dphi_array);

	free(H_array);
	free(DH_array);

	for(i=0; i<10000001; i++)
	{
		free(hk_array[i]);
	}
	free(hk_array);

	/* DONT FORGET TO FREE THE MEMORY USED USING MALLOC */
	/* IF YOU DECIDE TO KILL THE PROGRAM HALF-WAY, IT'LL TAKE */
	/* THE SYSTEM A VERY LONG TIME TO CLEANUP ALL OF THE GARBAGE! */

	return (0);
}

double V(double phi, double V0, double q, double phi0)
/* The value of the potential V as a function of the scalar field phi */
{
	return V0*exp(-sqrt(2.0/q)*(phi -phi0));
}

double dV(double phi, double V0, double q, double phi0)
/* the value of the derivative of the potential V with respect to phi, 
as a function of the scalar field phi*/
{
	return -sqrt(2/q)*V0*exp(-sqrt(2/q)*(phi -phi0));
}

double DDphi(double N, double phi, double Dphi, double V0, double q, double phi0)
/* the value of the second derivative of the scalar field phi with respect to e-fold N */
/* it is used when we are evaluating and integrating phi */
{
	return -(3 -(Dphi*Dphi)/2.)*Dphi -(dV(phi, V0, q, phi0)/(2*V(phi, V0, q, phi0)))*(6 -(Dphi*Dphi));
}

void rk4_stepper_phi(double N, double phi, double Dphi, double step, double V0, double q, double phi0, double *update_phi)
/* the stepper function used to evaluate phi as a function of e-fold N */
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
/* given the array phi_array, this function gives us access to the values 
as a function of N instead of using an array index */
{
	int index;
	index = floor((N-Ni)/step);
	return phi_array[index];
}

double Dphi_function(double *Dphi_array, double N, double Ni, double step)
/* given the array Dphi_array, this function gives us access to the values 
 as a function of N instead of using an array index */
{
	int index;
	index = floor((N-Ni)/step);
	return Dphi_array[index];
}

double A(double N, double ai)
/* the value of the scale factor as a function of e-fold N */
{
	return ai*exp(N);
}

double H(double N, double Ni, double step, double *H_array)
/* given the array H_array, this function gives us access to the values 
as a function of N instead of using an array index */
{
	int index;
	index = floor((N-Ni)/step);
	return H_array[index];
}

double DH(double N, double Ni, double step, double *DH_array)
/* given the array DH_array, this function gives us access to the values 
as a function of N instead of using an array index */
{
	int index;
	index = floor((N-Ni)/step);
	return DH_array[index];
}

/* root finding to estimate Nics and Nshss */
double find_Nics(double k, double *N_array, int npts, double ai, double Ni, double step, double *H_array)
{
	double *test_array = malloc(10000001*sizeof(double));
	int i, j;
	double min;

//	printf("evaluating Nics \n");

	/* creates a test array containing the absolute values of k -10^2*A*H fr various N */
	for (i=0; i<npts+2; i++)
	{
		test_array[i] = fabs(k -pow(10,2)*A(N_array[i], ai)*H(N_array[i], Ni, step, H_array));
	}

	/* look for the minimum value of the test array, the array index and corresponding value of N */
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

	free(test_array);

	return N_array[j];
}


double find_Nshss(double k, double *N_array, int npts, double ai, double Ni, double step, double *H_array)
{
	double *test_array = malloc(10000001*sizeof(double));
	int i, j;
	double min;

//	printf("evaluating Nshss \n");

	/* creates a test array containing the absolute values of k -10^(-5)*A*H fr various N */
	for (i=0; i<npts+2; i++)
	{
		test_array[i] = fabs(k -pow(10,-5)*A(N_array[i], ai)*H(N_array[i], Ni, step, H_array));
	}

	/* look for the minimum value of the test array, the array index and corresponding value of N */
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

	free(test_array);

	return N_array[j];
}

/* initialize hk and Dhk values */
/* according to the Bunch-Davies initial conditions */
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

void evolve_hk(double k, int npts, double ai, double Ni, double step, double *N_array, double *H_array, double *DH_array, double **hk_array, int size_hk_array, double Nics, double Nshss)
{
	double hk[2];
	double Dhk[2];

	double increment_hk[2];
	double increment_Dhk[2];

	double tps;

	double N;
	int i;

	Nics = find_Nics(k, N_array, npts, ai, Ni, step, H_array);
	Nshss = find_Nshss(k, N_array, npts, ai,  Ni, step, H_array);

	size_hk_array = floor((Nshss-Nics)/step);

	printf("got Nics and Nshss %lf, %lf, %d\n", Nics, Nshss, size_hk_array);

	initialize_hk(k, Nics, ai, hk);
	initialize_Dhk(k, Nics, ai, Dhk, Ni, step, H_array);

//	printf("done with hk and Dhk initialization \n");

	i = 0;
	N = Nics;

	hk_array[i][0] = hk[0];
	hk_array[i][1] = hk[1];

//	printf("%lf, %lf \t", N, Nics);
	while (N < Nshss)
	{
		/* stores the real and imaginary values of hk in hk_kx_arrays */
	//	printf("c");
		rk4_stepper_hk(k, N, hk, Dhk, Ni, step, increment_hk, increment_Dhk, ai, H_array, DH_array);
		hk[0] += increment_hk[0];
		hk[1] += increment_hk[1];
		Dhk[0] += increment_Dhk[0];
		Dhk[1] += increment_Dhk[1];

		N += step;
		i += 1;

		hk_array[i][0] = hk[0];
		hk_array[i][1] = hk[1];
	}

//	printf("%lf, %lf \n", N, Nshss);
	printf("done with hk evolution \n");

	return;
}

void calG(double k1, double k2, double k3, int size_hk_array, double Nics, double Nshss, double **hk_k1_array, double **hk_k2_array, double **hk_k3_array, double *CalG, double Ni, double step, double ai, double *H_array)
/* function to estimate of script G */
{
	double *func_int_real = malloc(size_hk_array*sizeof(double));
	double *func_int_imag = malloc(size_hk_array*sizeof(double));

	double int_real;
	double int_imag;

	double N;
	N = Nics;

	int i;

	int_real = 0;
	int_imag = 0;

	double e;
	e = 1/10;

	/* evaluates and stores the values of (A/H)*hk(k1)*hk(k2)*hk(k3)*exp(-e*k/(A*H)) */
	/* which is to be integrated over Nics to Nshss.*/
	for(i=0; i<size_hk_array; i++)
	{
		func_int_real[i] = (A(N, ai)/H(N, Ni, step, H_array))
					*(hk_k1_array[i][0]*hk_k2_array[i][0]*hk_k3_array[i][0] -hk_k1_array[i][1]*hk_k2_array[i][1]*hk_k3_array[i][0]
						-(hk_k1_array[i][0]*hk_k2_array[i][1] +hk_k1_array[i][1]*hk_k2_array[i][0])*hk_k3_array[i][1])
					*exp(e*(k1+k2+k3)/(3*A(N, ai)*H(N, Ni, step, H_array)));

		func_int_imag[i] = (A(N, ai)/H(N, Ni, step, H_array))
					*(hk_k1_array[i][0]*hk_k2_array[i][0]*hk_k3_array[i][1] -hk_k1_array[i][1]*hk_k2_array[i][1]*hk_k3_array[i][1]
						+(hk_k1_array[i][0]*hk_k2_array[i][1] +hk_k1_array[i][1]*hk_k2_array[i][0])*hk_k3_array[i][0])
					*exp(e*(k1+k2+k3)/(3*A(N, ai)*H(N, Ni, step, H_array)));
		N += step;
	}

	/* note that the exponential factor is being put in by hand.*/
	int_real = func_int_real[0] +func_int_real[size_hk_array];
	int_imag = func_int_imag[0] +func_int_imag[size_hk_array];

	/* simpson's rule of integration */
	for(i=0;i<size_hk_array/2;i++)
	{
		int_real += 2*func_int_real[2*i];
		int_real += 4*func_int_real[2*i-1];
		int_imag += 2*func_int_imag[2*i];
		int_imag += 4*func_int_imag[2*i-1];
	}

	CalG[0] = ((k1*k1 +k2*k2 +k3*k3)/4)*int_imag*step/3;
	CalG[1] = -((k1*k1 +k2*k2 +k3*k3)/4)*int_real*step/3;

	free(func_int_real);
	free(func_int_imag);
	
	return;
}

void calG_cc(double k1, double k2, double k3, int size_hk_array, double Nics, double Nshss, double **hk_k1_array, double **hk_k2_array, double **hk_k3_array, double *CalG_cc, double Ni, double step, double ai, double *H_array)
/* function to estimate the complex conjugate of script G */
{
	double *func_int_real = malloc(size_hk_array*sizeof(double));
	double *func_int_imag = malloc(size_hk_array*sizeof(double));

	double int_real;
	double int_imag;

	double N;
	N = Nics;

	int i;

	int_real = 0;
	int_imag = 0;

	double e;
	e = 1/10;

	/* evaluates and stores the values of (A/H)*hk(k1)*hk(k2)*hk(k3)*exp(-e*k/(A*H)) */
	/* which is to be integrated over Nics to Nshss.*/
	for(i=0; i<size_hk_array; i++)
	{
		func_int_real[i] = (A(N, ai)/H(N, Ni, step, H_array))
					*(hk_k1_array[i][0]*hk_k2_array[i][0]*hk_k3_array[i][0] -hk_k1_array[i][1]*hk_k2_array[i][1]*hk_k3_array[i][0]
						-(hk_k1_array[i][0]*hk_k2_array[i][1] +hk_k1_array[i][1]*hk_k2_array[i][0])*hk_k3_array[i][1])
					*exp(e*(k1+k2+k3)/(3*A(N, ai)*H(N, Ni, step, H_array)));

		func_int_imag[i] = (A(N, ai)/H(N, Ni, step, H_array))
					*(hk_k1_array[i][0]*hk_k2_array[i][0]*hk_k3_array[i][1] -hk_k1_array[i][1]*hk_k2_array[i][1]*hk_k3_array[i][1]
						+(hk_k1_array[i][0]*hk_k2_array[i][1] +hk_k1_array[i][1]*hk_k2_array[i][0])*hk_k3_array[i][0])
					*exp(e*(k1+k2+k3)/(3*A(N, ai)*H(N, Ni, step, H_array)));
		N += step;
	}

	/* note that the exponential factor is being put in by hand.*/
	int_real = func_int_real[0] +func_int_real[size_hk_array];
	int_imag = func_int_imag[0] +func_int_imag[size_hk_array];

	/* simpson's rule of integration */
	for(i=0;i<size_hk_array/2;i++)
	{
		int_real += 2*func_int_real[2*i];
		int_real += 4*func_int_real[2*i-1];
		int_imag += 2*func_int_imag[2*i];
		int_imag += 4*func_int_imag[2*i-1];
	}

	CalG_cc[0] = ((k1*k1 +k2*k2 +k3*k3)/4)*int_imag*step/3;
	CalG_cc[1] = ((k1*k1 +k2*k2 +k3*k3)/4)*int_real*step/3;

	free(func_int_real);
	free(func_int_imag);
	
	return;
}

void G_func(double *G, double a, double Ib, double c, double Id, double e, double If, double *CalG, double *CalG_cc)
/* function to estimate the value of the tensor bi-spectrum G using hk_kx[Nshss] and CalG, CalG_cc */
{
	G[0] = (a*c*e -Ib*Id*e -(a*Id +Ib*c)*If)*CalG[0] -(a*c*If -Ib*Id*If +(a*Id +Ib*c)*e)*CalG[1]
		+(a*c*e -Ib*Id*e -(a*Id +Ib*c)*If)*CalG_cc[0] -(a*c*If -Ib*Id*If +(a*Id +Ib*c)*e)*CalG_cc[1];

	G[1] = (a*c*e -Ib*Id*e -(a*Id +Ib*c)*If)*CalG[1] +(a*c*If -Ib*Id*If +(a*Id +Ib*c)*e)*CalG[0]
		+(a*c*e -Ib*Id*e -(a*Id +Ib*c)*If)*CalG_cc[1] -(a*c*If -Ib*Id*If +(a*Id +Ib*c)*e)*CalG_cc[0];
	
//	printf("%le, %le, %le, %le \n", (a*c*e -Ib*Id*e -(a*Id +Ib*c)*If)*CalG[0], -(a*c*If -Ib*Id*If +(a*Id +Ib*c)*e)*CalG[1], +(a*c*e -Ib*Id*e - (a*Id +Ib*c)*If)*CalG_cc[0], -(a*c*If -Ib*Id*If +(a*Id +Ib*c)*e)*CalG_cc[1]);

	return;
}
