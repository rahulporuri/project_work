import numpy as numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt
from sympy.functions import exp
import sympy as sym

x = sym.symbols('x')

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)
aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))
fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))
DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()

def DDhk(N,hk0,Dhk0):
    return -((3*N -1./N +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(aN(N)*H(N)))**2)*hk0)

def rk4_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(N +step, hk0 +F3*step, Dhk0 +f3*step)

    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

nn = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : nn(N)

a0 = 1e-05
n0 = numpy.exp(25)
k0 = 1e-10*numpy.exp(-25)
p = 1

npts = 5000

Nics = numpy.array([-10])
for i in range(200):
    Nics = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),Nics)

print 'Nics =', Nics
print 'we have Nics'

lolol = Nics

step = (0-Nics)/(npts-1)

hk0 = numpy.zeros(1,dtype=complex)
hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)

Dhk0 = numpy.zeros(1,dtype=complex)
Dhk0.real = -(sym.diff(A(x),x).subs(x,Nics).evalf()/(A(Nics)**2))*((2*k0)**(-1./2))
Dhk0.imag = -Nics*((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))

hkm = numpy.array(hk0, dtype = complex)

print 'starting Nics to 0'

N = Nics
while N < 0 - 2*step:
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hkm = numpy.append(hkm, hk0)
    N += step

print Nics, N-step, hk0, Dhk0

np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : np(N)
    
Nshss = numpy.array([10])
for i in range(200):
    Nshss = opt.newton_krylov(lambda N : (k0)**2 - 1e+06*fN(N),Nshss)

print 'Nshss =', Nshss
print 'we have Nshss'

step = (Nshss-0)/(npts-1)
hkp = numpy.array(hk0, dtype = complex)

print 'starting 0 to Nshss'

N = 0 + 2*step
while N < Nshss:
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hkp = numpy.append(hkp, hk0)
    N += step

print N-step, Nshss, hk0, Dhk0

Nics = -10.57877019

plt.cla()
plt.hold(True)
plt.semilogy(numpy.linspace(Nics,0,len(hkm)), numpy.absolute(hkm))
plt.semilogy(numpy.linspace(0,Nshss,len(hkp[:4500])), numpy.absolute(hkp[:4500]))
plt.savefig('abs_hk_vs_N_k3.png')

plt.cla()
plt.hold(True)
plt.plot(numpy.linspace(Nics,0,len(hkm)), hkm.imag)
plt.plot(numpy.linspace(0, Nshss,len(hkp[:4500])), hkp.imag[:4500])
plt.savefig('imag_hk_vs_N_k3.png')
