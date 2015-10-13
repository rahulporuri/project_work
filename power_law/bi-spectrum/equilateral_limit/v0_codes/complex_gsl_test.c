#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

int
main(void)
{

	gsl_complex z;
	GSL_SET_COMPLEX(&z, 2, 2);
//	printf("%lf, \n", gsl_complex_abs(z));
	printf("%lf, %lf \n", z, z);
	return (0);
}
