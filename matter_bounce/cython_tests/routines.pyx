# python setup.py build_ext --inplace
# http://docs.cython.org/src/userguide/fusedtypes.html
# http://docs.cython.org/src/userguide/sharing_declarations.html
# http://docs.cython.org/src/tutorial/pxd_files.html and
# http://stackoverflow.com/questions/7295638/how-do-you-get-cimport-to-work-in-cython

cimport numpy
cimport cython

import numpy

def DDhk(double k0, double N, cython.numeric hk0, cython.numeric Dhk0):
    return -((3.*N -(1./N) +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(A(N)*H(N)))**2)*hk0)

def rk4_step(double k0, double N, double hk0_real, double hk0_imag, double Dhk0_real, double Dhk0_imag, double step):
    cdef double F1, f1, F2, f2, F3, f3, F4, f4
    cdef double iF1, if1, iF2, if2, iF3, if3, iF4, if4
    F1 = Dhk0_real
    iF1 = Dhk0_imag
    f1 = DDhk(k0, N, hk0_real, Dhk0_real)
    if1 = DDhk(k0, N, hk0_imag, Dhk0_imag)
    F2 = Dhk0_real +f1*step/2.
    iF2 = Dhk0_imag +if1*step/2.
    f2 = DDhk(k0, N +step/2., hk0_real +F1*step/2., Dhk0_real +f1*step/2.)
    if2 = DDhk(k0, N +step/2., hk0_imag +iF1*step/2., Dhk0_imag +if1*step/2.)
    F3 = Dhk0_real +f2*step/2.
    iF3 = Dhk0_imag +if2*step/2.
    f3 = DDhk(k0, N +step/2., hk0_real +F2*step/2., Dhk0_real +f2*step/2.)
    if3 = DDhk(k0, N +step/2., hk0_imag +iF2*step/2., Dhk0_imag +if2*step/2.)
    F4 = Dhk0_real +f3*step
    iF4 = Dhk0_imag +if3*step
    f4 = DDhk(k0, N +step, hk0_real +F3*step, Dhk0_real +f3*step)    
    if4 = DDhk(k0, N +step, hk0_imag +iF3*step, Dhk0_imag +if3*step)    

    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6., (if1 +2*if2 +2*if3 +if4)*step/6., (F1 +2*F2 +2*F3 +F4)*step/6., (iF1 +2*iF2 +2*iF3 +iF4)*step/6.], dtype =numpy.double)
    # [Dhk_real, Dhk_imag, hk_real, hk_imag] update
    
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
