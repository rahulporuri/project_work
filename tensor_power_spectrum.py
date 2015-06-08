import numpy as numpy
import matplotlib.pyplot as plt
import sympy as sym
from sympy.functions import exp
import scipy.optimize as opt

x = sym.symbols('x')

a0 = 1e-05
n0 = numpy.exp(25)
p = 1

def DDhk(k0, N, hk0, Dhk0):
    return -((3.*N -(1./N) +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(A(N)*H(N)))**2)*hk0)

def rk4_step(k0, N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(k0, N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(k0, N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(k0, N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(k0, N +step, hk0 +F3*step, Dhk0 +f3*step)   

    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)

np = lambda N : +n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : np(N)

aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))
DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))    

Nshss = numpy.array([4])
for i in range(200):
	Nshss = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),Nshss)

#Nshss = 3.02287778
print Nshss

print 'lift off!'

k_min = 1e-25
k_max = 1e-10

k_vs_hk = numpy.zeros(1,dtype=complex)

k0 = k_min
while k0 < k_max:
    print 'k0 = ', k0

    nm = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
    n = lambda N : nm(N)
    
    Nics = numpy.array([-4])
    for i in range(200):
        Nics = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),Nics)
    
    hk0 = numpy.zeros(1,dtype=complex)
    hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)
    
    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -(sym.diff(A(x),x).subs(x,Nics).evalf()/(A(Nics)**2))*((2*k0)**(-1./2))
    Dhk0.imag = -Nics*((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))

    print 'got Nics, hk0 and Dhk0'

    npts = 7500
    step = (0-Nics)/(npts-1)

    print 'starting from Nics'
    
    N = Nics
    while N < 0 -2*step:
        #array = euler_step(N, hk0, Dhk0, step)
        array = rk4_step(k0, N, hk0, Dhk0, step)
        hk0 = hk0 + array[1]
        Dhk0 = Dhk0 + array[0]
        N += step

    print N-step, 0, hk0, Dhk0
    print 'half way there!'

    A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
    an = lambda n : a0*((1+(n/n0)**2)**p)

    np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)
    n = lambda N : np(N)
    
    aN = lambda N : an(n(N))
    h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
    H = lambda N : h(n(N))
    DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()

    fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
    fN = lambda N : fn(n(N))    

    step = (Nshss-0)/(npts-1)

    N = 0 +step
    while N < Nshss + step:
        #array = euler_step()
        array = rk4_step(k0, N, hk0, Dhk0, step)
        hk0 = hk0 + array[1]
        Dhk0 = Dhk0 + array[0]
        N += step
        
    k_vs_hk = numpy.append(k_vs_hk, hk0)
    print N-step, Nshss, hk0, Dhk0, Nics
    print '\n'
    k0 = 10*k0


k_list = numpy.array([10**(-25 + i) for i in range(5)])
print len(k_list), len(k_vs_hkhk)
TPS = [8*(k_list[i])**3/(2*numpy.pi**2)*(numpy.absolute(k_vs_hk[i+1]))**2 for i in range(len(k_list))]
#print k_list, TPS

plt.loglog(k_list, TPS)
