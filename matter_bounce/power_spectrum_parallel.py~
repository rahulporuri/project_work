import numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt

import time
import multiprocessing as mp
import random
import string

parallel_output = mp.Queue()

tps_file = open('power_spectrum_bouncing_model_functional.dat','w')

A = lambda N : a0*numpy.exp(N**2/2.)
an = lambda n : a0*(1.+(n/n0)**2)
aN = lambda N : an(n(N))
h = lambda n : (2.*a0*n/n0**2)*(1./an(n)**2)
H = lambda N : h(n(N))

fn = lambda n : (2.*a0/n0**2)*(1./an(n))
fN = lambda N : fn(n(N))    

Heavi = lambda N : (1./2)*(numpy.sign(N)+1.)

nm = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)

DHm = lambda N : -100000.0*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**(-0.5)*numpy.exp(N**2/2)**1.0/(100000.0*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 100000.0)**2 + 40000000000.0*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**0.5*numpy.exp(N**2/2)**1.0/(100000.0*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 100000.0)**3
DHp = lambda N : 100000.0*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**(-0.5)*numpy.exp(N**2/2)**1.0/(100000.0*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 100000.0)**2 - 40000000000.0*N*(1.0*numpy.exp(N**2/2)**1.0 - 1.0)**0.5*numpy.exp(N**2/2)**1.0/(100000.0*(1.0*numpy.exp(N**2/2)**1.0 - 1)**1.0 + 100000.0)**3

n = lambda N : Heavi(N)*np(N) + Heavi(-N)*nm(N)
DH = lambda N : Heavi(N)*DHp(N) + Heavi(-N)*DHm(N)

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

    return [numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex)] # [Dhk, hk] update

def solve_Nics(k0):
    Nics = numpy.array([-5])
    for i in range(200):
        Nics = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),Nics)

    return Nics

def initialize(Nics, k0):
    hk0 = numpy.zeros(1,dtype=complex)
    hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)

    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -(Nics/A(Nics))*((2*k0)**(-1./2))
    Dhk0.imag = -Nics*((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))
    
    return [hk0, Dhk0]

def evolve_hk(Nics, Nshss, k0, hk0, Dhk0, step):
    N = Nics
    while N < Nshss:
        array = rk4_step(k0, N, hk0, Dhk0, step)
        hk0 = hk0 + array[1]
        Dhk0 = Dhk0 + array[0]
        N += step
        
    return hk0

a0 = 1e+05
n0 = 1.0
p = 1.0

k_min = 1e-30
k_max = 1e-4
npts = 30000

Nshss = numpy.array([7])
for i in range(200):
    Nshss = opt.newton_krylov(lambda N : (k_max)**2 - 1e+04*fN(N),Nshss)

print 'Nshss =', Nshss

def main(k0, Nshss):
    
    Nics = solve_Nics(k0)
    
    temp_array = initialize(Nics, k0)
    hk0 = temp_array[0]
    Dhk0 = temp_array[1]

    step = (Nshss-Nics)/(npts)
    
    hk0 = evolve_hk(Nics, Nshss, k0, hk0, Dhk0, step)
    
    tps_temp = 8*(k0)**3/(2*numpy.pi**2)*(numpy.absolute(hk0))**2
    
    print k0, numpy.absolute(hk0), str(tps_temp).strip('[]')
    
    return str(tps_temp).strip('[]')

k_list = numpy.array([10**(-10 + i) for i in range(4)])

pool = mp.Pool(processes = 4)
temp_results = [pool.apply_async(main, args = (k,Nshss,)) for k in k_list]
results = []

for i in range(len(temp_results)):
    results.append(temp_results[i].get())

results = numpy.asarray(results, dtype=numpy.float)
print results

'''
plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal{P}}_{\rm T}(k)$')
plt.title(r'${\mathcal{P}}_{\rm T}(k)$ as a function of $k$')
numerics, = plt.loglog(k_list, TPS)
plt.legend([numerics],['numerical results'])
plt.savefig('power_spectrum_bouncing_model.png')

data = [line.split() for line in open('fort.30')]
data = numpy.asarray(data, dtype=numpy.float)

plt.cla()
plt.hold(True)
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal{P}}_{\rm T}(k)$')
plt.title(r'${\mathcal{P}}_{\rm T}(k)$ as a function of $k$')
py_numerics, = plt.loglog(k_list, TPS)
f_numerics, = plt.loglog(data[:,0], data[:,1])
plt.legend([py_numerics, f_numerics],['python', 'fortran'], loc='lower left')
plt.savefig('power_spectrum_matter_bounce.png')

plt.cla()
plt.hold(True)
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal{P}}_{\rm T}(k)$')
plt.title(r'${\mathcal{P}}_{\rm T}(k)$ as a function of $k$')
py_numerics, = plt.loglog(k_list[:24], TPS[:24])
f_numerics, = plt.loglog(data[:19,0], data[:19,1])
plt.legend([py_numerics, f_numerics],['python', 'fortran'], loc='lower left')
plt.savefig('power_spectrum_matter_bounce_2.png')
'''
