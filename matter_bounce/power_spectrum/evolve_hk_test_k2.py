import numpy as numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt

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

A = lambda N : a0*numpy.exp(N**2/2.)
an = lambda n : a0*((1.+(n/n0)**2)**p)
aN = lambda N : an(n(N))
h = lambda n : (2.*a0*n/n0**2)*(1./an(n)**2)
H = lambda N : h(n(N))

fn = lambda n : (2.*a0/n0**2)*(1./an(n))
fN = lambda N : fn(n(N))    

Heavi = lambda N : (1./2)*(numpy.sign(N)+1.)

nm = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)

DHm = lambda N : -1.3887943864964e-16*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**(-0.5)*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**2 + 5.55517754598561e-21*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**0.5*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**3

DHp = lambda N : +1.3887943864964e-16*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**(-0.5)*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**2 - 5.55517754598561e-21*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**0.5*numpy.exp(N**2/2)**1.0/(1.0e-5*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 1.0e-5)**3

n = lambda N : Heavi(N)*np(N) + Heavi(-N)*nm(N)
DH = lambda N : Heavi(N)*DHp(N) + Heavi(-N)*DHm(N)

f = open('evolve_hk_test_2.dat','w')

a0 = 1e-05
n0 = numpy.exp(25)
p = 1

k0 = 10**(-5)*numpy.exp(-25)

Nics = numpy.array([-4])
for i in range(200):
    Nics = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),Nics)

print 'Nics =', Nics
print 'we have Nics'

Nshss = numpy.array([5])
for i in range(200):
    Nshss = opt.newton_krylov(lambda N : (k0)**2 - 1e+06*fN(N),Nshss)

print 'Nshss =', Nshss
print 'we have Nshss'

hk0 = numpy.zeros(1,dtype=complex)
hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)

Dhk0 = numpy.zeros(1,dtype=complex)
Dhk0.real = -(Nics/A(Nics))*((2*k0)**(-1./2))
Dhk0.imag = -Nics*((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))

print 'starting Nics to Nshss'

hk_array = numpy.array(hk0, dtype = complex)

f.write(str(Nics)+"\t"+str(hk0)+"\n")

npts = 5000
step = (Nshss-Nics)/(npts)
N = Nics
while N < Nshss:
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hk_array = numpy.append(hk_array, hk0)
    N += step
    f.write(str(Nics)+"\t"+str(hk0.real)+"\t"+str(hk0.imag)+"\n")    

print Nics, N-step, Nshss, hk0, Dhk0

f.close()

Nics = -8.11533591

plt.subplot(121)

plt.cla()
plt.hold(True)
plt.semilogy(numpy.linspace(Nics,Nshss,len(hk_array[:4700])), numpy.absolute(hk_array[:4700]))
#plt.savefig('abs_hk_vs_N_k2.png')

plt.subplot(122)

plt.cla()
plt.hold(True)
plt.plot(numpy.linspace(Nics,Nshss,len(hk_array[:4700])), hk_array.imag[:4700])
#plt.savefig('imag_hk_vs_N_k2.png')

plt.savefig('k2.png')
