
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

get_ipython().magic(u'matplotlib inline')


# In[3]:

parallel_output = mp.Queue()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# In[4]:

q = 51.
V0 = (204./100.)*1e-08
t0 = (q*(3.*q -1.)/V0)**(1./2)

phi0 = 1.
dphi0 = (2.*q)**(1./2)/t0

Ni = 0.
Nf = 70.


# In[5]:

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

    return [(f1 +2*f2 +2*f3 +f4)*_step/6., (F1 +2*F2 +2*F3 +F4)*_step/6.] # [Dhk, hk] update


# In[6]:

npts = 100000
step = (Nf-Ni)/(npts)

phi_ = phi0
Dphi_ = Dphi0

phi_array = numpy.empty(0)
Dphi_array = numpy.empty(0)
N_array = numpy.empty(0)

N = Ni
while N < Nf +step:
    phi_array = numpy.append(phi_array, phi_)
    Dphi_array = numpy.append(Dphi_array, Dphi_)
    N_array = numpy.append(N_array, N)
    
    array = rk4_step(N, phi_, Dphi_, step)
    phi_ = phi_ + array[1]
    Dphi_ = Dphi_ + array[0]
    
    N += step


# In[7]:

phi = lambda _N : phi_array[int((_N-Ni)/step)]
Dphi = lambda _N : Dphi_array[int((_N-Ni)/step)]

H = lambda _N : (V(phi(_N))/(3 -Dphi(_N)**2/2))**(1./2)
DH = lambda _N : -(1./2)*H(_N)*Dphi(_N)**2.

ai = 1e-05
A = lambda _N : ai*numpy.exp(_N)

k0 = numpy.empty(0)

def DDhk(_k, _N, _hk, _Dhk):
    return -((3. +(DH(_N)/H(_N)))*_Dhk +((_k/(A(_N)*H(_N)))**2)*_hk)


def rk4_step(_k, _N, _hk, _Dhk, _step):
    F1 = _Dhk
    f1 = DDhk(_k, _N, _hk, _Dhk)
    F2 = _Dhk +f1*_step/2.
    f2 = DDhk(_k, _N +_step/2., _hk +F1*_step/2., _Dhk +f1*_step/2.)
    F3 = _Dhk +f2*_step/2.
    f3 = DDhk(_k, _N +_step/2., _hk +F2*_step/2., _Dhk +f2*_step/2.)
    F4 = _Dhk +f3*_step
    f4 = DDhk(_k, _N +_step, _hk +F3*_step, _Dhk +f3*_step)

    return (numpy.array([(f1 +2*f2 +2*f3 +f4)*_step/6.], dtype=complex),
            numpy.array([(F1 +2*F2 +2*F3 +F4)*_step/6.], dtype=complex))
            # [Dhk, hk] update


# In[8]:

def solve_Nics(k, eN_array):
    Ni = eN_array[0]
    step = eN_array[1] -eN_array[0]
    Nics_temp = numpy.asarray([k - 1e+02*A(N)*H(N) for N in eN_array])
    nics_test = numpy.where(Nics_temp > 0)
    return Ni + nics_test[0][-1]*step

def solve_Nshss(k, eN_array):
    Ni = eN_array[0]
    step = eN_array[1] -eN_array[0]
    Nshss_temp = numpy.asarray([k - 1e-05*A(N)*H(N) for N in eN_array])
    nshss_test = numpy.where(Nshss_temp > 0)
    return Ni + nshss_test[0][-1]*step

def initialize_hk(k, _Nics):
    hk0 = numpy.zeros(1,dtype=complex)             
    hk0.real = (((2.*k)**(1./2))*A(_Nics))**(-1.)
    return hk0

def initialize_Dhk(k, _Nics):
    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -(1/A(_Nics))*((2*k)**(-1./2))
    Dhk0.imag = -((k/2)**(1./2))/(A(_Nics)*A(_Nics)*H(_Nics))
    return Dhk0

def evolve_hk(k, _Nics, _Nshss, _step):    
    hk = numpy.empty(0, dtype=complex)
    Dhk = numpy.empty(0, dtype=complex)
    
    hk = initialize_hk(k, _Nics)
    Dhk = initialize_Dhk(k, _Nics)
    
    hk_array = numpy.empty(0, dtype=complex)
    
    N = _Nics
    while N < _Nshss:
        hk_array = numpy.append(hk_array, hk)

        array = rk4_step(k, N, hk, Dhk, _step)
        hk = hk + array[1]
        Dhk = Dhk + array[0]

        N += _step

    return hk_array


# In[9]:

e = 10**(-1)

def calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, _Nics, _Nshss):
    N_range = numpy.linspace(_Nics, _Nshss, len(hk_k1_array))
    func_int = ((A(N_range)/numpy.asarray([H(N) for N in N_range]))*
                (numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))*
                (numpy.exp(-e*(k1 +k2 +k3)/(3.*A(N_range)*numpy.asarray([H(N) for N in N_range])))))
    result = simps(func_int, N_range)

    return -(k1**2. +k2**2. +k3**2)/4.*result*numpy.array([0.+1.j], dtype=complex)


def calG_cc(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, _Nics, _Nshss):
    N_range = numpy.linspace(_Nics, _Nshss, len(hk_k1_array))
    func_int = ((A(N_range)/numpy.asarray([H(N) for N in N_range]))*
                (hk_k1_array*hk_k2_array*hk_k3_array)*
                (numpy.exp(-e*(k1 +k2 +k3)/(3.*A(N_range)*numpy.asarray([H(N) for N in N_range])))))
    result = simps(func_int, N_range)

    return (k1**2. +k2**2. +k3**2)/4.*result*numpy.array([0.+1.j], dtype=complex)


# In[10]:

#k1 = 1e-06
k1 = 5e-02
k3 = numpy.arange(0.,1.,.1)*k1

k_array = []

for i in range(len(k3)):
    if k3[i]/k1 < 0.5:
        k2 = numpy.linspace(1. -k3[i]/k1, 1., 2 +int(k3[i]/k1/0.1))*k1
        [k_array.append([kx, k3[i]]) for kx in k2]
    else :
        k2 = numpy.linspace(k3[i]/k1, 1., 2 +int((1. -k3[i]/k1)/0.1))*k1
        [k_array.append([kx, k3[i]]) for kx in k2]


# In[11]:

k_array


# In[12]:

for k_set in k_array[2:]:
    
    k2, k3 = k_set
    
    if k2 > k3:
        Nics = solve_Nics(k3, N_array)
    else :
        Nics = solve_Nics(k2, N_array)
    
    Nshss = solve_Nshss(k1, N_array)

    hk_k1_array = numpy.empty(0, dtype=complex)
    hk_k2_array = numpy.empty(0, dtype=complex)
    hk_k3_array = numpy.empty(0, dtype=complex)

    hk_k1_array = evolve_hk(k1, Nics, Nshss, step)
    hk_k2_array = evolve_hk(k2, Nics, Nshss, step)
    hk_k3_array = evolve_hk(k3, Nics, Nshss, step)
    
    tps_k1 = 2.*(k1)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k1_array[-1]))**2.
    tps_k2 = 2.*(k2)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k2_array[-1]))**2.
    tps_k3 = 2.*(k3)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k3_array[-1]))**2.

    CalG = calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)
    CalG_cc = calG_cc(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)

    G = (hk_k1_array[-1]*hk_k2_array[-1]*hk_k3_array[-1])*CalG +(numpy.conj(hk_k1_array[-1])*numpy.conj(hk_k2_array[-1])*numpy.conj(hk_k3_array[-1]))*CalG_cc
    h_NL = -(4./(2.*numpy.pi**2.))**2.*(k1**3.*k2**3.*k3**3*G)/(k1**3.*tps_k2*tps_k3 +k2**3.*tps_k3*tps_k1 +k3**3.*tps_k1*tps_k2)

    print Nics, Nshss
    print k1, k2, k3, str(tps_k1).strip('[]'), str(tps_k2).strip('[]'), str(tps_k3).strip('[]')
    print str(CalG).strip('[]'), str(G).strip('[]'), str(h_NL).strip('[]')


# In[ ]:

'''    N_range = numpy.linspace(Nics, Nshss, len(hk_k1_array))
    func_int = ((A(N_range)/numpy.asarray([H(N) for N in N_range]))*
                (numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))*
                (numpy.exp(-e*k1/(A(N_range)*numpy.asarray([H(N) for N in N_range])))))
    CalG = -(k1**2. +k2**2. +k3**2)/4.*simps(func_int, N_range)*numpy.array([0.+1.j], dtype=complex)'''


# In[37]:

temp = A(N_range)/numpy.asarray([H(N) for N in N_range])
test = temp*(numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))
final = test*(numpy.exp(-e*k1/(A(N_range)*numpy.asarray([H(N) for N in N_range]))))
print len(test), len(temp), len(final)


# In[36]:

numpy.exp(-e*k1/(A(N_range)*numpy.asarray([H(N) for N in N_range])))


# In[35]:

-e*k1/(A(N_range)*numpy.asarray([H(N) for N in N_range]))


# In[ ]:



