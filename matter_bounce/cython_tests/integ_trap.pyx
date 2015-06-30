# python setup.py build_ext --inplace
# call the functions from the shared library

cimport numpy
import numpy




def fn(double x):
	return numpy.exp(x)

def integ_trap(double a, double b, double h):
	cdef double output
	output = fn(a)/2.
	while a<b:
		output += fn(a+h)
		a += h
	return (output + fn(b)/2.)*h

def integ_simp(double a, double b, double h):
	cdef double output
	output = 0
	output = fn(a)/3.
	while a<b:
		output += (2./3.)*fn(a+2*h) + (4./3.)*fn(a+h)
		a += 2*h
	return (output + fn(b)/3.)*h

cdef double exact
exact = fn(1) - fn(0)

def main():
	for i in range(1,10):
		h = 1e-01**i
		print integ_trap(0,1,h) - exact, integ_simp(0,1,h) - exact, h
