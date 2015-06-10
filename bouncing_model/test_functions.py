import numpy as numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt
from sympy.functions import exp
import sympy as sym

x = sym.symbols('x')

a0 = 1e-05
n0 = numpy.exp(25)
k0 = numpy.exp(-25)
p = 1

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)

plt.cla()
plt.plot(numpy.arange(-1e06, 1e06,100),an(numpy.arange(-1e06, 1e06,100)))
plt.savefig('scale_factor_vs_conformal_time.png')

aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))
DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))

#################################################

Nm = numpy.linspace(-5,0,51)

nm = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : nm(N)

plt.cla()
plt.hold(True)
plt.plot(Nm,[H(i) for i in Nm])

Np = numpy.linspace(0,5,51)

np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : np(N)

plt.plot(Np,[H(i) for i in Np])
plt.savefig('H_vs_N.png')

##################################################

Nm = numpy.linspace(-5,0,51)

nm = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : nm(N)

plt.cla()
plt.hold(True)
plt.plot(Nm,[DH(i) for i in Nm])

Np = numpy.linspace(0,5,51)

np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : np(N)

plt.plot(Np,[DH(i) for i in Np])
plt.savefig('DH_vs_N.png')

#################################################

plt.cla()
plt.plot(numpy.linspace(-5, 5,101), [fN(i) for i in numpy.linspace(-5, 5,101)])
plt.savefig('blah.png')
