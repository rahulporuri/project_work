#include <stdio.h>
#include <stdlib.h>

void set_val(double *var, double cnst);

int
main(void)
{
	double x;
	double y;

	set_val(&x, 3);
	set_val(&y, 4);

	printf("in main : %lf, %lf \n", x, y);

	return;
}

void set_val(double *var, double c)
{
	*var = c;
	printf("inside the function : %lf %lf\n", *var, c);
	return;
}
