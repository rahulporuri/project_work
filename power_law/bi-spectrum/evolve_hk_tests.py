
# coding: utf-8

# In[1]:

import numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.integrate import romb, simps

import time
import multiprocessing as mp
import random
import string


# In[2]:

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# In[3]:

q = 51.
V0 = (204./100.)*1e-08
t0 = (q*(3.*q -1.)/V0)**(1./2)

phi0 = 1.
dphi0 = (2.*q)**(1./2)/t0

Ni = 0.
Nf = 70.

phi_ptr = open("../mathematica_codes/phi_vs_N_python.txt", "w")
H_ptr = open("../mathematica_codes/H_vs_N_python.txt", "w")
DH_ptr = open("../mathematica_codes/DH_vs_N_python.txt", "w")
eps_ptr = open("../mathematica_codes/eps_vs_N_python.txt","w")
# In[4]:

V = lambda _phi : V0*numpy.exp(-(2./q)**(1./2)*(_phi -phi0))
dV = lambda _phi : -(2./q)**(1./2)*V0*numpy.exp(-(2./q)**(1./2)*(_phi -phi0)) 

H0 = ((1./3)*(dphi0**2/2. +V(phi0)))**(1./2.)
Dphi0 = dphi0/H0

def DDphi(_N, _phi, _Dphi):
    return -(3 -_Dphi**2/2.)*_Dphi -(dV(_phi)/(2*V(_phi)))*(6 -_Dphi**2)

def rk4_step(_N, _phi, _Dphi, _step):
    F1 = _Dphi
    f1 = DDphi(_N, _phi, _Dphi)
    F2 = _Dphi +f1*_step/2.
    f2 = DDphi(_N +_step/2., _phi +F1*_step/2., _Dphi +f1*_step/2.)
    F3 = _Dphi +f2*step/2.
    f3 = DDphi(_N +_step/2., _phi +F2*_step/2., _Dphi +f2*_step/2.)
    F4 = _Dphi +f3*step
    f4 = DDphi(_N +_step, _phi +F3*_step, _Dphi +f3*_step)  

    return (F1 +2*F2 +2*F3 +F4)*_step/6., (f1 +2*f2 +2*f3 +f4)*_step/6. # (phi, Dphi) update


# In[5]:

npts = 100000
step = (Nf-Ni)/(npts)

phi_ = phi0
Dphi_ = Dphi0

phi_array = numpy.empty(0)
Dphi_array = numpy.empty(0)
N_array = numpy.empty(0)
H_array = numpy.empty(0)
DH_array = numpy.empty(0)
eps_array = numpy.empty(0)

N = Ni

print N, phi_, Dphi_, dphi0, H0
while N < Nf+step:
    N_array = numpy.append(N_array, N)
    phi_array = numpy.append(phi_array, phi_)
    Dphi_array = numpy.append(Dphi_array, Dphi_)

    phi_inc, Dphi_inc = rk4_step(N, phi_, Dphi_, step)
    phi_ = phi_ +phi_inc
    Dphi_ = Dphi_ +Dphi_inc

    N += step


print N, phi_, Dphi_
print len(N_array)

# In[6]:

phi = lambda _N : phi_array[int((_N-Ni)/step)]
Dphi = lambda _N : Dphi_array[int((_N-Ni)/step)]

H = lambda _N : (V(phi(_N))/(3. -Dphi(_N)**2/2.))**(1./2)
DH = lambda _N : H(_N)*Dphi(_N)

for i in range(len(N_array)):
	H_array = numpy.append(H_array, (V(phi_array[i])/(3. -Dphi_array[i]**2/2.))**(1./2))
	DH_array = numpy.append(DH_array, -(1./2)*(V(phi_array[i])/(3. -Dphi_array[i]**2/2.))**(1./2)*Dphi_array[i]**2)
	eps_array = numpy.append(eps_array, Dphi_array[i]**2/2.0)

print H_array[0], DH_array[0], H_array[-1], DH_array[-1]

for i in range(len(N_array)):
	phi_ptr.write(str(N_array[i])+" , "+str(phi_array[i])+"\n")
	H_ptr.write(str(N_array[i])+" , "+str(H_array[i])+"\n")
	DH_ptr.write(str(N_array[i])+" , "+str(DH_array[i])+"\n")
	eps_ptr.write(str(N_array[i])+" , "+str(eps_array[i])+"\n")

phi_ptr.close()
H_ptr.close()
DH_ptr.close()
eps_ptr.close()

'''

ai = 1e-05
A = lambda _N : ai*numpy.exp(_N)

k0 = numpy.empty(0)

def DDhk(_k, _N, _hk, _Dhk):
    return -(3. +(DH(_N)/H(_N)))*_Dhk -((_k/(A(_N)*H(_N)))**2)*_hk

def rk4_step(_k, _N, _hk, _Dhk, _step):
    F1 = _Dhk
    f1 = DDhk(_k, _N, _hk, _Dhk)
    F2 = _Dhk +f1*_step/2.
    f2 = DDhk(_k, _N +_step/2., _hk +F1*_step/2., _Dhk +f1*_step/2.)
    F3 = _Dhk +f2*_step/2.
    f3 = DDhk(_k, _N +_step/2., _hk +F2*_step/2., _Dhk +f2*_step/2.)
    F4 = _Dhk +f3*_step
    f4 = DDhk(_k, _N +_step, _hk +F3*_step, _Dhk +f3*_step)

    return (numpy.array([(F1 +2*F2 +2*F3 +F4)*_step/6.], dtype=complex),
            numpy.array([(f1 +2*f2 +2*f3 +f4)*_step/6.], dtype=complex))
            # [hk, Dhk] update


# In[7]:

def solve_Nics(k, eN_array):
    Ni = eN_array[0]
    step = eN_array[1] -eN_array[0]
    Nics_temp = numpy.asarray([k -1e+02*A(N)*H(N) for N in eN_array])
    nics_test = numpy.where(Nics_temp > 0)
    return Ni + nics_test[0][-1]*step

def solve_Nshss(k, eN_array):
    Ni = eN_array[0]
    step = eN_array[1] -eN_array[0]
    Nshss_temp = numpy.asarray([k -1e-03*A(N)*H(N) for N in eN_array])
    nshss_test = numpy.where(Nshss_temp > 0)
    return Ni + nshss_test[0][-1]*step

def initialize_hk(k, _Nics):
    hk0 = numpy.zeros(1,dtype=complex)             
    hk0.real = 1./((2.*k)**(1./2.))/A(_Nics)
    return hk0

def initialize_Dhk(k, _Nics):
    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -1./((2.*k)**(1./2.))/A(_Nics)
    Dhk0.imag = -((k/2.)**(1./2.))/(A(_Nics)*A(_Nics)*H(_Nics))
    return Dhk0 

def evolve_hk(k, _Nics, _Nshss, _step):    
    hk = numpy.empty(0, dtype=complex)
    Dhk = numpy.empty(0, dtype=complex)

    hk = initialize_hk(k, _Nics)
    Dhk = initialize_Dhk(k, _Nics)

    #print _Nics, str(hk).strip('[]'), str(Dhk).strip('[]')

    hk_array = numpy.empty(0, dtype=complex)
    N = _Nics

    while N < _Nshss:
        hk_array = numpy.append(hk_array, hk)
        hk_inc, Dhk_inc = rk4_step(k, N, hk, Dhk, _step)
        hk = hk + hk_inc
        Dhk = Dhk + Dhk_inc
        N += _step

    #print N, _Nshss, str(hk).strip('[]'), str(Dhk).strip('[]'), '\n'
    return hk_array


# In[8]:

k0 = 1e-06
while k0 < 1e-01:
    Nics = solve_Nics(k0, N_array)
    print k0, Nics, initialize_hk(k0, Nics), initialize_Dhk(k0, Nics)
    k0 = 10*k0


# In[9]:

k0 = 1e-06
Nics = solve_Nics(k0, N_array)
while k0 < 1e-01:
    print k0, Nics, initialize_hk(k0, Nics), initialize_Dhk(k0, Nics)
    k0 = 10*k0


# In[8]:

k0 = 1e-06
while k0 < 1e-01:
    Nics = solve_Nics(k0, N_array)
    Nshss = solve_Nshss(k0, N_array)

    hk_k0_array = numpy.empty(0, dtype=complex)
    hk_k0_array = evolve_hk(k0, Nics, Nshss, step)
    tps_k0 = 8.*(k0)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k0_array[-1]))**2.

    print k0, Nics, Nshss, str(hk_k0_array[-1]).strip('[]'), str(tps_k0).strip('[]')
    k0 = 10*k0


# In[11]:

k0 = 1e-06
while k0 < 1e-01:
    Nics = solve_Nics(k0, N_array)
    Nshss = solve_Nshss(k0, N_array)

    hk_k0_array = numpy.empty(0, dtype=complex)
    hk_k0_array = evolve_hk(k0, Nics-1, Nshss, step)
    tps_k0 = 8.*(k0)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k0_array[-1]))**2.

    print k0, Nics-1, Nshss, str(hk_k0_array[-1]).strip('[]'), str(tps_k0).strip('[]')
    k0 = 10*k0


# In[12]:

k0 = 1e-06
while k0 < 1e-01:
    Nics = solve_Nics(k0, N_array)
    Nshss = solve_Nshss(k0, N_array)

    hk_k0_array = numpy.empty(0, dtype=complex)
    hk_k0_array = evolve_hk(k0, Nics, Nshss+1, step)
    tps_k0 = 8.*(k0)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k0_array[-1]))**2.

    print k0, Nics, Nshss+1, str(hk_k0_array[-1]).strip('[]'), str(tps_k0).strip('[]')
    k0 = 10*k0


# In[15]:

k0 = 1e-06
while k0 < 1e-01:
    Nics = solve_Nics(k0, N_array)
    Nshss = solve_Nshss(k0, N_array)

    hk_k0_array = numpy.empty(0, dtype=complex)
    hk_k0_array = evolve_hk(k0, Nics, Nshss, step)
    tps_k0 = 8.*(k0)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k0_array[-1]))**2.

    print k0, Nics, Nshss, str(hk_k0_array[0]).strip('[]'), str(hk_k0_array[0]*A(Nics)).strip('[]')
    print str(hk_k0_array[-1]).strip('[]'), str(tps_k0).strip('[]')
    k0 = 10*k0


# In[9]:

k0 = 1e-06
Nics = solve_Nics(k0, N_array)
while k0 < 1e-01:
    Nshss = solve_Nshss(k0, N_array)

    hk_k0_array = numpy.empty(0, dtype=complex)
    hk_k0_array = evolve_hk(k0, Nics, Nshss, step)
    tps_k0 = 8.*(k0)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k0_array[-1]))**2.

    print k0, Nics, Nshss, str(hk_k0_array[-1]).strip('[]'), str(tps_k0).strip('[]')
    k0 = 10*k0


# In[14]:

k0 = 1e-06
Nics = solve_Nics(k0, N_array)
while k0 < 1e-01:
    Nshss = solve_Nshss(k0, N_array)

    hk_k0_array = numpy.empty(0, dtype=complex)
    hk_k0_array = evolve_hk(k0, Nics, Nshss, step)
    tps_k0 = 8.*(k0)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k0_array[-1]))**2.

    print k0, Nics, Nshss, str(hk_k0_array[0]).strip('[]'), str(hk_k0_array[0]*A(Nics)).strip('[]')
    print str(hk_k0_array[-1]).strip('[]'), str(tps_k0).strip('[]')
    k0 = 10*k0


# In[15]:




# In[ ]:



'''
