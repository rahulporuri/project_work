import numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.integrate import romb, simps

import time
import multiprocessing as mp
import random
import string

parallel_output = mp.Queue()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

q = 51.
V0 = (204./100.)*1e-08
t0 = (q*(3.*q -1.)/V0)**(1./2)

phi0 = 1.
dphi0 = (2.*q)**(1./2)/t0

Ni = 0.
Nf = 70. 

V = lambda phi : V0*numpy.exp(-(2./q)**(1./2)*(phi -phi0))
dV = lambda phi : -(2./q)**(1./2)*V0*numpy.exp(-(2./q)**(1./2)*(phi -phi0))

H0 = ((1./3)*(dphi0**2/2. +V(phi0)))**(1./2.)
Dphi0 = dphi0/H0

def DDphi(N, phi0, Dphi0):
    return -(3 -Dphi0**2/2.)*Dphi0 -(dV(phi0)/(2*V(phi0)))*(6 -Dphi0**2)

def rk4_step(N, phi0, Dphi0, step):
    F1 = Dphi0
    f1 = DDphi(N, phi0, Dphi0)
    F2 = Dphi0 +f1*step/2.
    f2 = DDphi(N +step/2., phi0 +F1*step/2., Dphi0 +f1*step/2.)
    F3 = Dphi0 +f2*step/2.
    f3 = DDphi(N +step/2., phi0 +F2*step/2., Dphi0 +f2*step/2.)
    F4 = Dphi0 +f3*step
    f4 = DDphi(N +step, phi0 +F3*step, Dphi0 +f3*step)  

    return [(f1 +2*f2 +2*f3 +f4)*step/6., (F1 +2*F2 +2*F3 +F4)*step/6.] # [Dhk, hk] update

npts = 100000
step = (Nf-Ni)/(npts)

phi_ = phi0
Dphi_ = Dphi0

phi_array = numpy.empty(0)
Dphi_array = numpy.empty(0)
N_array = numpy.empty(0) 

N = Ni
while N < Nf:
    phi_array = numpy.append(phi_array, phi_)
    Dphi_array = numpy.append(Dphi_array, Dphi_)
    N_array = numpy.append(N_array, N)

    array = rk4_step(N, phi_, Dphi_, step)
    phi_ = phi_ + array[1]
    Dphi_ = Dphi_ + array[0]
    N += step

phi = lambda N : phi_array[int((N-Ni)/step)]
Dphi = lambda N : Dphi_array[int((N-Ni)/step)]

H = lambda N : (V(phi(N))/(3 -Dphi(N)**2/2))**(1./2)
DH = lambda N : -(1./2)*H(N)*Dphi(N)**2

ai = 1e-05
A = lambda N : ai*numpy.exp(N)
k0 = numpy.empty(0)

def DDhk(k0, N, hk0, Dhk0):
    return -((3. +(DH(N)/H(N)))*Dhk0 +((k0/(A(N)*H(N)))**2)*hk0)

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

def solve_Nics(k0, N_array):
    step = N_array[1] -N_array[0]
    Nics_temp = numpy.asarray([k0 - 1e+02*A(N)*H(N) for N in N_array])   
    nics_test = numpy.where(Nics_temp > 0)
    return Ni + nics_test[0][-1]*step

def solve_Nshss(k0, N_array):
    step = N_array[1] -N_array[0]
    Nshss_temp = numpy.asarray([k0 - 1e-05*A(N)*H(N) for N in N_array])
    nshss_test = numpy.where(Nshss_temp > 0)
    return Ni + nshss_test[0][-1]*step

def initialize_hk(k0, Nics):
    hk0 = numpy.zeros(1,dtype=complex)             
    hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)
    return hk0

def initialize_Dhk(k0, Nics):
    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -(1/A(Nics))*((2*k0)**(-1./2))
    Dhk0.imag = -((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))
    return Dhk0

def evolve_hk(k0, hk0, Dhk0, Nics, Nshss, step):
    hk_array = numpy.empty(0, dtype=complex)
    N = Nics
    while N < Nshss:
        hk_array = numpy.append(hk_array, hk0)
        array = rk4_step(k0, N, hk0, Dhk0, step)
        hk0 = hk0 + array[1]
        Dhk0 = Dhk0 + array[0]
        N += step

    return hk_array

e = 10**(-1.)

def calG(hk_array, k0, Nics, Nshss):
    N_range = numpy.linspace(Nics, Nshss, len(hk_array))
    func_int = (A(N_range)/numpy.asarray([H(N) for N in N_range]))*numpy.conj(hk_array)**3*numpy.exp(-e*k0/(A(N_range)*numpy.asarray([H(N) for N in N_range])))
    #result = romb(func_int, N_Array[1]-N_Array[0])
    result = simps(func_int, N_range)
    return (-1/4.)*result*numpy.array([0.+1.j], dtype=complex)

def calG_cc(hk_array, k0, Nics, Nshss):
    N_range = numpy.linspace(Nics, Nshss, len(hk_array))
    func_int = (A(N_range)/numpy.asarray([H(N) for N in N_range]))*(hk_array)**3*numpy.exp(-e*k0/(A(N_range)*numpy.asarray([H(N) for N in N_range])))
    #result = romb(func_int, N_Array[1]-N_Array[0])
    result = simps(func_int, N_range)
    return (+1/4.)*result*numpy.array([0.+1.j], dtype=complex)

def main(k0, N_array):
    Nics = solve_Nics(k0, N_array)
    Nshss = solve_Nshss(k0, N_array)
    step = N_array[1] -N_array[0]

    hk0 = numpy.empty(0,dtype=complex)             
    Dhk0 = numpy.empty(0,dtype=complex)

    hk0 = initialize_hk(k0, Nics)
    Dhk0 = initialize_Dhk(k0, Nics)

    hk_array = numpy.empty(0, dtype=complex)
    hk_array = evolve_hk(k0, hk0, Dhk0, Nics, Nshss, step)
    tps= 2.*(k0)**3/(2.*numpy.pi**2)*(numpy.absolute(hk_array[-1]))**2

    CalG = calG(hk_array, k0, Nics, Nshss)
    CalG_cc = calG_cc(hk_array, k0, Nics, Nshss)

    G = (3.*(k0)**2)*((hk_array[-1]**3)*CalG + (numpy.conj(hk_array[-1])**3)*CalG_cc)
    h_NL = (-1./6)*((4./(2.*numpy.pi**2))**2)*(k0**6)*G/(tps**2)

    print k0, str(tps).strip('[]'), str((k0**(3./2))*numpy.absolute(CalG)).strip('[]'), str(G.real).strip('[]'), str((k0**6)*G.real).strip('[]'), str(h_NL.real).strip('[]')
    return None

k_list = numpy.array([10**((-12 + i)/2.) for i in range(13)])

pool = mp.Pool(processes = 4)
temp_results = [pool.apply_async(main, args = (k0, N_array,)) for k0 in k_list]
results = []

for i in range(len(k_list)):
    results.append(temp_results[i].get())
#    k_vs_hk[i] = temp_results[i].get()

#results = numpy.asarray(results, dtype=numpy.float)
#print k_list, results
#print k_vs_hk
#print '\n'
