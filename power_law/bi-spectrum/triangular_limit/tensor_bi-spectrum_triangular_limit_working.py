import numpy
from scipy.integrate import simps

import time
import multiprocessing as mp
import random
import string

parallel_output = mp.Queue()

q = 51.
V0 = (204./100.)*1e-08
t0 = (q*(3.*q -1.)/V0)**(1./2)

phi0 = 1.
dphi0 = (2.*q)**(1./2)/t0

Ni = 0.
Nf = 70.

# Note that in this code, I use the prefix 'd' to represent derivative with respect to time (except for the case of dV where the derivative is with respect to phi) and the prefix 'D' to represent derivative with respect to e-fold N. Also, the suffix '0' is used to represent the initial conditions in various cases. Also, as can be seen here, we evaluate the scalar field in the e-fold N range Ni to Nf.

V = lambda _phi : V0*numpy.exp(-(2./q)**(1./2)*(_phi -phi0))
dV = lambda _phi : -(2./q)**(1./2)*V0*numpy.exp(-(2./q)**(1./2)*(_phi -phi0))

''' Functions to evaluate the values of the potential function V(phi)
and the derivative of V with respect to phi.
Note that functions can be defined using the lambda notation, as shown 
above or using the usual def and return statements, as shown below.'''

H0 = ((1./3)*(dphi0**(2.)/2. +V(phi0)))**(1./2.)
Dphi0 = dphi0/H0

def DDphi(_N, _phi, _Dphi):
    ''' Returns the value of the second derivative of 
    phi with respect to e-fold N.'''
    return -(3. -_Dphi**(2.)/2.)*_Dphi -(dV(_phi)/(2.*V(_phi)))*(6. -_Dphi**(2.))

def phi_rk4_step(_N, _phi, _Dphi, _step):
    ''' Returns 2 values, the first of the two is the value by which phi 
    needs to be updated and the second of the two is the value by which the 
    first derivative of phi with respect to e-fold N needs to be updated.'''
    F1 = _Dphi
    f1 = DDphi(_N, _phi, _Dphi)
    F2 = _Dphi +f1*_step/2.
    f2 = DDphi(_N +_step/2., _phi +F1*_step/2., _Dphi +f1*_step/2.)
    F3 = _Dphi +f2*step/2.
    f3 = DDphi(_N +_step/2., _phi +F2*_step/2., _Dphi +f2*_step/2.)
    F4 = _Dphi +f3*step
    f4 = DDphi(_N +_step, _phi +F3*_step, _Dphi +f3*_step)  

    return (F1 +2.*F2 +2.*F3 +F4)*_step/6., (f1 +2.*f2 +2.*f3 +f4)*_step/6.

'''We evolve the scalar field phi for e-fold N ranging from Ni to Nf.'''

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
    
    phi_update, Dphi_update = phi_rk4_step(N, phi_, Dphi_, step)
    phi_ = phi_ +phi_update
    Dphi_ = Dphi_ +Dphi_update
    
    N += step

#2000001
#2000000
N_new = numpy.linspace(Ni,Nf,2000001)
phi_array_new = numpy.interp(N_new, N_array, phi_array)
Dphi_array_new = numpy.interp(N_new, N_array, Dphi_array)

phi_array = phi_array_new
Dphi_array = Dphi_array_new
N_array = N_new
step = (Nf-Ni)/(2000000)

phi = lambda _N : phi_array[int((_N-Ni)/step)]
Dphi = lambda _N : Dphi_array[int((_N-Ni)/step)]

H = lambda _N : (V(phi(_N))/(3. -Dphi(_N)**(2.)/2.))**(1./2)
DH = lambda _N : -(1./2)*H(_N)*Dphi(_N)**2.

'''The above functions let us access the values of H(N) and DH(N) 
when we try to evaluate the tensor perturbations h_k. We have obtained 
these values from the phi and Dphi values earlier.'''

ai = 1e-05
A = lambda _N : ai*numpy.exp(_N)
'''The scale factor in terms of e-fold N.'''

def DDhk(_k, _N, _hk, _Dhk):
    '''Returns the value of the second derivative of the tensor perturbatons
    h_k th respec to e-fold N. We need this value when we are trying to 
    evluate h_k'''
    return -((3. +(DH(_N)/H(_N)))*_Dhk +((_k/(A(_N)*H(_N)))**(2.))*_hk)


def hk_rk4_step(_k, _N, _hk, _Dhk, _step):
    '''a runge-kutta 4 stepper function that returns the value by which
    h_k nd Dh_k need to be updated.'''
    F1 = _Dhk
    f1 = DDhk(_k, _N, _hk, _Dhk)
    F2 = _Dhk +f1*_step/2.
    f2 = DDhk(_k, _N +_step/2., _hk +F1*_step/2., _Dhk +f1*_step/2.)
    F3 = _Dhk +f2*_step/2.
    f3 = DDhk(_k, _N +_step/2., _hk +F2*_step/2., _Dhk +f2*_step/2.)
    F4 = _Dhk +f3*_step
    f4 = DDhk(_k, _N +_step, _hk +F3*_step, _Dhk +f3*_step)

#    print f1, f2, f3, f4, F1, F2, F3, F4, _step

    return (f1 +2.*f2 +2.*f3 +f4)*_step/6., (F1 +2.*F2 +2.*F3 +F4)*_step/6.
           # [Dhk, hk] update

def solve_Nics(k, eN_array):
    '''Returns the value of e-fold N when the mode is
    in the sub-Hubble domain, which we define as k/(A*H) =10^2.'''
    Ni = eN_array[0]
    step = eN_array[1] -eN_array[0]
    Nics_temp = numpy.asarray([k - 1e+02*A(N)*H(N) for N in eN_array])
    nics_test = numpy.where(Nics_temp > 0)
    return Ni + nics_test[0][-1]*step

def solve_Nshss(k, eN_array):
    '''Returns the value of e-fold N when the mode is
    in the super-Hubble domain, which we define as k/(A*H) =10^(-5).'''
    Ni = eN_array[0]
    step = eN_array[1] -eN_array[0]
    Nshss_temp = numpy.asarray([k - 1e-05*A(N)*H(N) for N in eN_array])
    nshss_test = numpy.where(Nshss_temp > 0)
    return Ni + nshss_test[0][-1]*step

def initialize_hk(k, _Nics):
    '''Returns the value of h_k for the mode k at e-fold N of _Nics.
    We obtain his value by imposing the Bunch-Davies initial conditions'''
    hk0 = numpy.zeros(1,dtype=complex)
    hk0.real = (((k)**(1./2))*A(_Nics))**(-1.)
    return hk0

def initialize_Dhk(k, _Nics):
    '''Returns the value of h_k for the mode k at e-fold N of _Nshss.
    We obtain his value by imposing the Bunch-Davies initial conditions'''
    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -(1./A(_Nics))*((k)**(-1./2))
    Dhk0.imag = -((k)**(1./2))/(A(_Nics)*A(_Nics)*H(_Nics))
    return Dhk0

def evolve_hk(k, _Nics, _Nshss, _step):
    '''Returns the values of h_k for the mode k for e-fold N ranging from
    _Nics to _Nshss. We use the h_k values later on to estimate calG.'''
    hk = numpy.empty(0, dtype=complex)
    Dhk = numpy.empty(0, dtype=complex)
    
    hk = initialize_hk(k, _Nics)
    Dhk = initialize_Dhk(k, _Nics)
    
    hk_array = numpy.empty(0, dtype=complex)
    
    N = _Nics
    while N < _Nshss:
        hk_array = numpy.append(hk_array, hk)

        array = hk_rk4_step(k, N, hk, Dhk, _step)
        hk = hk + array[1]
        Dhk = Dhk + array[0]

        N += _step

    return hk_array

def calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, _Nics, _Nshss):
    '''Returns the value of \mathcal{G} which is in turn used to estimate G, the 
    tensor bi-spectrum. The integral is evaluated for e-fold N ranging from 
    _Nics till _Nshss. Note that the extra factor exp(-(e*k)/(A*H)) is put in by 
    hand to satisfy the consistency relation.'''
    N_range = numpy.linspace(_Nics, _Nshss, len(hk_k1_array))
    func_int = ((A(N_range)/numpy.asarray([H(N) for N in N_range]))*
                (numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))*
                (numpy.exp(-e*(k1 +k2 +k3)/(3.*A(N_range)*numpy.asarray([H(N) for N in N_range])))))
    result = simps(func_int, N_range)

    return -(k1**2. +k2**2. +k3**2)/4.*result*numpy.array([0.+1.j], dtype=complex)

def calG_cc(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, _Nics, _Nshss):
    '''Returns the value of the complex conjugate of \mathcal{G} which is 
    in turn used to estimate G, the tensor bi-spectrum. The integral is 
    evaluated for e-fold N ranging from _Nics till _Nshss. Note that the 
    extra factor exp(-(e*k)/(A*H)) is put in by hand to satisfy the consistency relation.'''
    N_range = numpy.linspace(_Nics, _Nshss, len(hk_k1_array))
    func_int = ((A(N_range)/numpy.asarray([H(N) for N in N_range]))*
                (hk_k1_array*hk_k2_array*hk_k3_array)*
                (numpy.exp(-e*(k1 +k2 +k3)/(3.*A(N_range)*numpy.asarray([H(N) for N in N_range])))))
    result = simps(func_int, N_range)

    return (k1**2. +k2**2. +k3**2)/4.*result*numpy.array([0.+1.j], dtype=complex)

k1 = 1e-02
kmin = 1e-03
kmax = 1e-02
k3 = numpy.arange(0.,1.,.05)*k1

k_array = []
'''Since k1 is fixed, k_array will contain the set of [k2, k3] values.'''
for i in range(len(k3)):
    if k3[i]/k1 < 0.5:
        k2 = numpy.linspace(1. -k3[i]/k1, 1., 2 +int(k3[i]/k1/0.05))*k1
        [k_array.append([kx, k3[i]]) for kx in k2]
    else :
        k2 = numpy.linspace(k3[i]/k1, 1., 2 +int((1. -k3[i]/k1)/0.05))*k1
        [k_array.append([kx, k3[i]]) for kx in k2]

print len(k_array)
print k_array

hk_k1_array = numpy.empty(0, dtype=complex)
hk_k2_array = numpy.empty(0, dtype=complex)
hk_k3_array = numpy.empty(0, dtype=complex)

e =1./50

Nics = solve_Nics(kmin, N_array)
Nshss = solve_Nshss(kmax, N_array)

hk_k1_array = evolve_hk(k1, Nics, Nshss, step)
tps_k1 = 4.*(k1)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k1_array[-1]))**2.   

def main(k_set):
    k2, k3 = k_set

    hk_k2_array = evolve_hk(k2, Nics, Nshss, step)
    hk_k3_array = evolve_hk(k3, Nics, Nshss, step)

    tps_k2 = 4.*(k2)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k2_array[-1]))**2.
    tps_k3 = 4.*(k3)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k3_array[-1]))**2.
    
    CalG = calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)
    CalG_cc = calG_cc(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)

    G = ((hk_k1_array[-1]*hk_k2_array[-1]*hk_k3_array[-1])*CalG 
		+(numpy.conj(hk_k1_array[-1])*numpy.conj(hk_k2_array[-1])*numpy.conj(hk_k3_array[-1]))*CalG_cc)

    h_NL = -((4./(2.*numpy.pi**2.))**2.*(k1**3.*k2**3.*k3**3*G)/
             (2.*k1**3.*tps_k2*tps_k3 +2.*k2**3.*tps_k3*tps_k1 +2.*k3**3.*tps_k1*tps_k2))
   
    print k1, k2, k3, str(tps_k1).strip('[]'), str(tps_k2).strip('[]'), str(tps_k3).strip('[]'), str(CalG.real).strip('[]'), str(CalG.imag).strip('[]'), str(G.real).strip('[]'), str(G.imag).strip('[]'), str(h_NL.real).strip('[]')

    return None

pool = mp.Pool(processes = 4)
temp_results = [pool.apply_async(main, args = (k_set, )) for k_set in k_array[2:]]
results = []

for i in range(len(temp_results)):
        results.append(temp_results[i].get())

print results

CalG = calG(hk_k1_array, hk_k1_array, hk_k1_array, k1, k1, k1, Nics, Nshss)
CalG_cc = calG_cc(hk_k1_array, hk_k1_array, hk_k1_array, k1, k1, k1, Nics, Nshss)

G = ((hk_k1_array[-1]*hk_k1_array[-1]*hk_k1_array[-1])*CalG 
	+(numpy.conj(hk_k1_array[-1])*numpy.conj(hk_k1_array[-1])*numpy.conj(hk_k1_array[-1]))*CalG_cc)

h_NL = -((4./(2.*numpy.pi**2.))**2.*(k1**3.*k1**3.*k1**3*G)/
        (2.*k1**3.*tps_k1*tps_k1 +2.*k1**3.*tps_k1*tps_k1 +2.*k1**3.*tps_k1*tps_k1))

print k1, k1, k1, str(tps_k1).strip('[]'), str(tps_k1).strip('[]'), str(tps_k1).strip('[]'), str(CalG.real).strip('[]'), str(CalG.imag).strip('[]'), str(G.real).strip('[]'), str(G.imag).strip('[]'), str(h_NL.real).strip('[]')
