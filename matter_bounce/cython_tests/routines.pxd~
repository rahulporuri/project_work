# python setup.py build_ext --inplace
# http://docs.cython.org/src/userguide/fusedtypes.html
# http://docs.cython.org/src/userguide/sharing_declarations.html
# http://docs.cython.org/src/tutorial/pxd_files.html and
# http://stackoverflow.com/questions/7295638/how-do-you-get-cimport-to-work-in-cython

cimport numpy
cimport cython

import numpy

cdef DDhk(double k0, double N, cython.numeric hk0, cython.numeric Dhk0):
    return -((3.*N -(1./N) +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(A(N)*H(N)))**2)*hk0)

cdef rk4_step(double k0, double N, cython.numeric hk0, cython.numeric Dhk0, double step):
    cdef complex F1, f1, F2, f2, F3, f3, F4, f4
    F1 = Dhk0
    f1 = DDhk(k0, N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(k0, N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(k0, N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(k0, N +step, hk0 +F3*step, Dhk0 +f3*step)    

    return [numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex)] # [Dhk, hk] update
    
A = lambda double N : a0*numpy.exp(N**2/2.)
an = lambda double n : a0*((1.+(n/n0)**2)**p)
aN = lambda double N : an(n(N))
h = lambda double n : (2.*a0*n/n0**2)*(1./an(n)**2)
H = lambda double N : h(n(N))

fn = lambda double n : (2.*a0/n0**2)*(1./an(n))
fN = lambda double N : fn(n(N))    

Heavi = lambda double N : (1./2)*(numpy.sign(N)+1.)

nm = lambda double N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
np = lambda double N : n0*((A(N)/a0)**(1./p)-1)**(1./2)

DHm = lambda double N : -1.3887943864964e-16*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**(-0.5)*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**2 + 5.55517754598561e-21*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**0.5*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**3

DHp = lambda double N : +1.3887943864964e-16*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**(-0.5)*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**2 - 5.55517754598561e-21*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**0.5*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**3

n = lambda double N : Heavi(N)*np(N) + Heavi(-N)*nm(N)
DH = lambda double N : Heavi(N)*DHp(N) + Heavi(-N)*DHm(N)

cdef double a0, n0, p

a0 = 1e+05
n0 = 1.0
p = 1
