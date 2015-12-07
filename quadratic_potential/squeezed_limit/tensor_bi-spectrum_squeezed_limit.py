import numpy
import matplotlib.pyplot as plt
from scipy.integrate import simps

import time
import multiprocessing as mp
import random
import string

parallel_output = mp.Queue()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

phi0 = 165./10.
Dphi0 = -10**(-5.)

m = 7147.*10**(-9.)

Ni = 0.
Nf = 69.

''' Note that in this code, I use the prefix 'd' to represent derivative with respect 
to time (except for the case of dV where the derivative is with respect to phi) and 
the prefix 'D' to represent derivative with respect to e-fold N. Also, the suffix '0' 
is used to represent the initial conditions in various cases. Also, as can be seen here, 
we evaluate the scalar field in the e-fold N range Ni to Nf.'''

V = lambda _phi : m**2.*_phi**2./2.
dV = lambda _phi : m**2.*_phi

''' Functions to evaluate the values of the potential function V(phi)
and the derivative of V with respect to phi.
Note that functions can be defined using the lambda notation, as shown 
above or using the usual def and return statements, as shown below.'''

def DDphi(_N, _phi, _Dphi):
    ''' Returns the value of the second derivative of 
    phi with respect to e-fold N.'''
    return -(3 -_Dphi**2/2.)*_Dphi -(dV(_phi)/(2*V(_phi)))*(6 -_Dphi**2)

def rk4_step(_N, _phi, _Dphi, _step):
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

    return (F1 +2*F2 +2*F3 +F4)*_step/6., (f1 +2*f2 +2*f3 +f4)*_step/6. # [phi, Dphi] update

'''We evolve the scalar field phi for e-fold N ranging from Ni to Nf.'''
npts = 500000
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
    
    phi_update, Dphi_update = rk4_step(N, phi_, Dphi_, step)
    phi_ = phi_ + phi_update
    Dphi_ = Dphi_ + Dphi_update
    
    N += step

phi = lambda _N : phi_array[int((_N-Ni)/step)]
Dphi = lambda _N : Dphi_array[int((_N-Ni)/step)]

H = lambda _N : (V(phi(_N))/(3 -Dphi(_N)**2/2))**(1./2)
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
    return -((3. +(DH(_N)/H(_N)))*_Dhk +((_k/(A(_N)*H(_N)))**2)*_hk)

def rk4_step(_k, _N, _hk, _Dhk, _step):
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

    return (numpy.array([(f1 +2*f2 +2*f3 +f4)*_step/6.], dtype=complex),
            numpy.array([(F1 +2*F2 +2*F3 +F4)*_step/6.], dtype=complex))
            # [Dhk, hk] update

def solve_Nics(k, eN_array):
    '''Returns the value of e-fold N when the mode is
    in the sub-Hubble domain, which we define as k/(A*H) =10^2.'''
    Ni = eN_array[0]
    step = eN_array[1] -eN_array[0]
    Nics_temp = numpy.asarray([k - 1e+03*A(N)*H(N) for N in eN_array])
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
    hk0.real = 1./(k**(1./2))/A(_Nics)
    return hk0

def initialize_Dhk(k, _Nics):
    '''Returns the value of h_k for the mode k at e-fold N of _Nshss.
    We obtain his value by imposing the Bunch-Davies initial conditions'''
    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -1./(k**(1./2))/A(_Nics)
    Dhk0.imag = -k**(1./2)/(A(_Nics)*A(_Nics)*H(_Nics))
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

        array = rk4_step(k, N, hk, Dhk, _step)
        hk = hk + array[1]
        Dhk = Dhk + array[0]

        N += _step

    return hk_array

e = 50**(-1)

def calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, _Nics, _Nshss):
    '''Returns the value of \mathcal{G} which is in turn used to estimate G, the 
    tensor bi-spectrum. The integral is evaluated for e-fold N ranging from 
    _Nics till _Nshss. Note that the extra factor exp(-(e*k)/(A*H)) is put in by 
    hand to satisfy the consistency relation.'''
    N_range = numpy.linspace(_Nics, _Nshss, len(hk_k1_array))
    func_int = ((A(N_range)/numpy.asarray([H(N) for N in N_range]))*
                (numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))*
                (numpy.exp(-e*k2/(A(N_range)*numpy.asarray([H(N) for N in N_range])))))
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
                (numpy.exp(-e*k2/(A(N_range)*numpy.asarray([H(N) for N in N_range])))))
    result = simps(func_int, N_range)

    return (k1**2. +k2**2. +k3**2)/4.*result*numpy.array([0.+1.j], dtype=complex)

'''We first evolve the pseudo-zero mode 'k0' from Nics 
corresponding to k0 to Nshss corresponding to the largest mode.'''
k0 = 1000*ai*numpy.exp(Ni)*H(Ni)

k_min = 1e-05
k_max = 1e-01

print 'evolving k0'

#Nics_k0 = solve_Nics(k0, N_array)
Nics_k0 = 0.
Nshss_k0 = solve_Nshss(k_max, N_array)

hk_k0_array = numpy.empty(0, dtype=complex)
hk_k0_array = evolve_hk(k0, Nics_k0, Nshss_k0, step)
tps_k0 = 4.*(k0)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_k0_array[-1]))**2.

print k0, str(tps_k0).strip('[]')

def main(ki, N_array):
    '''The main routine that calls the other functions. It takes the mode k as input and 
    estimates Nics and Nshss. Afterwards, it evaluates the h_k array from Nics to Nshss. 
    Following that, it estimates the tensor power spectrum, \mathcal{G}, G and h_NL.'''
    Nics = solve_Nics(ki, N_array)
    Nshss = solve_Nshss(ki, N_array)

    hk_ki_array = numpy.empty(0, dtype=complex)
    hk_ki_array = evolve_hk(ki, Nics, Nshss, step)
    tps_ki = 4.*(ki)**3./(2.*numpy.pi**2.)*(numpy.absolute(hk_ki_array[-1]))**2.

    test_array= hk_k0_array[int((Nics -Nics_k0)/step):int((Nics -Nics_k0)/step)+len(hk_ki_array)]

    CalG = calG(test_array, hk_ki_array, hk_ki_array, k0, ki, ki, Nics, Nshss)
    CalG_cc = calG_cc(test_array, hk_ki_array, hk_ki_array, k0, ki, ki, Nics, Nshss)

    G = (test_array[-1]*hk_ki_array[-1]**2)*CalG +(numpy.conj(test_array[-1])*numpy.conj(hk_ki_array[-1])**2)*CalG_cc
    h_NL = -(((4./(2.*numpy.pi**2.))**2.)*(k0**3.*ki**3.*ki**3.*G)/
        (2*k0**3.*tps_ki*tps_ki +2.*ki**3.*tps_k0*tps_ki +2.*ki**3.*tps_k0*tps_ki))

    return k0, ki, str(tps_ki).strip('[]'), str(numpy.absolute(CalG)).strip('[]'), str(G.real).strip('[]'), str(h_NL.real).strip('[]')

k_list = numpy.asarray([10**((-12.+i)/2.) for i in range(4)])
'''k_list contains the mode k =k1 =k2.'''

pool = mp.Pool(processes = 4)
temp_results = [pool.apply_async(main, args = (ki, N_array,)) for ki in k_list]
results = []

for i in range(len(k_list)):
    results.append(temp_results[i].get())

data = numpy.asarray(results,dtype=float)
numpy.savetxt('squeezed_limit.dat', data)

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal P}_T$')
plt.title(r'${\mathcal P}_T$ as a function of $k$')
python, = plt.loglog(data[:,1], data[:,2], '--')
plt.legend([python],['python'])
plt.savefig('tps_vs_k_sq.png')

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal G}$')
plt.title(r'${\mathcal G}$ as a function of $k$')
python, = plt.loglog(data[:,1], data[:,3], '--')
plt.legend([python],['python'])
plt.savefig('calG_vs_k_sq.png')
'''
plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'$k^{3/2}{\mathcal G}$')
plt.title(r'$k^{3/2}{\mathcal G}$ as a function of $k$')
python, = plt.semilogx(data[:,1], data[:,0]**(3/2)*data[:,3], '--')
plt.legend([python],['python'])
plt.savefig('k32calG_vs_k_sq.png')
'''
plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'G')
plt.title(r'G as a function of $k$')
python, = plt.loglog(data[:,1], -data[:,4], '--')
plt.legend([python],['python'])
plt.savefig('G_vs_k_sq.png')

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'G')
plt.title(r'$k^6$G as a function of $k$')
python, = plt.loglog(data[:,1], -data[:,0]**3.*data[:,1]**3.*data[:,4], '--')
plt.legend([python],['python'])
plt.savefig('k6G_vs_k_sq.png')

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\rm h}_{NL}$')
plt.title(r'${\rm h}_{NL}$ as a function of $k$')
python, = plt.semilogx(data[:,1], data[:,5], '--')
plt.legend([python],['python'])
plt.savefig('h_NL_vs_k_sq.png')
