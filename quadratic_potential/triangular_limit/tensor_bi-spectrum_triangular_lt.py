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

'''Note that in this code, I use the prefix 'd' to represent derivative with 
respect to time (except for the case of dV where the derivative is with respect 
to phi) and the prefix 'D' to represent derivative with respect to e-fold N. 
Also, the suffix '0' is used to represent the initial conditions in various cases. 
Also, as can be seen here, we evaluate the scalar field in the e-fold N range Ni to Nf.'''

V = lambda phi : V0*numpy.exp(-(2./q)**(1./2)*(phi -phi0))
dV = lambda phi : -(2./q)**(1./2)*V0*numpy.exp(-(2./q)**(1./2)*(phi -phi0))

H0 = ((1./3)*(dphi0**2/2. +V(phi0)))**(1./2.)
Dphi0 = dphi0/H0

''' Functions to evaluate the values of the potential function V(phi)
and the derivative of V with respect to phi.
Note that functions can be defined using the lambda notation, as shown 
above or using the usual def and return statements, as shown below.'''

def DDphi(N, phi0, Dphi0):
	''' Returns the value of the second derivative of 
	phi with respect to e-fold N.'''
	return -(3 -Dphi0**2/2.)*Dphi0 -(dV(phi0)/(2*V(phi0)))*(6 -Dphi0**2)

def rk4_step(N, phi0, Dphi0, step):
	''' Returns 2 values, the first of the two is the value by which phi 
	needs to be updated and the second of the two is the value by which the 
	first derivative of phi with respect to e-fold N needs to be updated.'''
	F1 = Dphi0
	f1 = DDphi(N, phi0, Dphi0)
	F2 = Dphi0 +f1*step/2.
	f2 = DDphi(N +step/2., phi0 +F1*step/2., Dphi0 +f1*step/2.)
	F3 = Dphi0 +f2*step/2.
	f3 = DDphi(N +step/2., phi0 +F2*step/2., Dphi0 +f2*step/2.)
	F4 = Dphi0 +f3*step
	f4 = DDphi(N +step, phi0 +F3*step, Dphi0 +f3*step)  

	return [(f1 +2*f2 +2*f3 +f4)*step/6., (F1 +2*F2 +2*F3 +F4)*step/6.] # [Dhk, hk] update

'''We evolve the scalar field phi for e-fold N ranging from Ni to Nf.'''

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

H = lambda N : (V(phi(N))/(3. -Dphi(N)**2/2.))**(1./2)
DH = lambda N : -(1./2)*H(N)*Dphi(N)**2

'''The above functions let us access the values of H(N) and DH(N) 
when we try to evaluate the tensor perturbations h_k. We have obtained 
these values from the phi and Dphi values earlier.'''

ai = 1e-05
A = lambda N : ai*numpy.exp(N)
'''The scale factor in terms of e-fold N.'''

k0 = numpy.empty(0)

def DDhk(k0, N, hk0, Dhk0):
	'''Returns the value of the second derivative of the tensor perturbatons
	h_k th respec to e-fold N. We need this value when we are trying to 
	evluate h_k'''
	return -((3. +(DH(N)/H(N)))*Dhk0 +((k0/(A(N)*H(N)))**2)*hk0)

def rk4_step(k0, N, hk0, Dhk0, step):
	'''a runge-kutta 4 stepper function that returns the value by which
	h_k nd Dh_k need to be updated.'''
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
	'''Returns the value of e-fold N when the mode is
	in the sub-Hubble domain, which we define as k/(A*H) =10^2.'''
	step = N_array[1] -N_array[0]
	Nics_temp = numpy.asarray([k0 - 1e+02*A(N)*H(N) for N in N_array])   
	nics_test = numpy.where(Nics_temp > 0)
	return Ni + nics_test[0][-1]*step

def solve_Nshss(k0, N_array):
	'''Returns the value of e-fold N when the mode is
	in the super-Hubble domain, which we define as k/(A*H) =10^(-5).'''
	step = N_array[1] -N_array[0]
	Nshss_temp = numpy.asarray([k0 - 1e-05*A(N)*H(N) for N in N_array])
	nshss_test = numpy.where(Nshss_temp > 0)
	return Ni + nshss_test[0][-1]*step

def initialize_hk(k0, Nics):
	'''Returns the value of h_k for the mode k at e-fold N of _Nics.
	We obtain his value by imposing the Bunch-Davies initial conditions'''
	hk0 = numpy.zeros(1,dtype=complex)             
	hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)
	return hk0

def initialize_Dhk(k0, Nics):
	'''Returns the value of h_k for the mode k at e-fold N of _Nshss.
	We obtain his value by imposing the Bunch-Davies initial conditions'''
	Dhk0 = numpy.zeros(1,dtype=complex)
	Dhk0.real = -(1/A(Nics))*((2*k0)**(-1./2))
	Dhk0.imag = -((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))
	return Dhk0

def evolve_hk(k0, Nics, Nshss, step):
	'''Returns the values of h_k for the mode k for e-fold N ranging from
	_Nics to _Nshss. We use the h_k values later on to estimate calG.'''
	hk = numpy.empty(0, dtype=complex)
	Dhk = numpy.empty(0, dtype=complex)

	hk = initialize_hk(k, _Nics)
	Dhk = initialize_Dhk(k, _Nics)

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

def calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss):
	'''Returns the value of \mathcal{G} which is in turn used to estimate G, the 
	tensor bi-spectrum. The integral is evaluated for e-fold N ranging from 
	_Nics till _Nshss. Note that the extra factor exp(-(e*k)/(A*H)) is put in by 
	hand to satisfy the consistency relation.'''
	N_range = numpy.linspace(Nics, Nshss, len(hk_k1_array))

	func_int = (A(N_range)/numpy.asarray([H(N) for N in N_range])*
		(numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))*
		numpy.exp(-e*(k1 +k2 +k3)/(3*A(N_range)*numpy.asarray([H(N) for N in N_range]))))

	result = simps(func_int, N_range)
	return (-1/4.)*(k1**2+k2**2+k3**2)*result*numpy.array([0.+1.j], dtype=complex)

def calG_cc(hk_kx_array, hk_ky_array, hk_kz_array, kx, ky, kz, N_ics, N_shss):
	 '''Returns the value of the complex conjugate of \mathcal{G} which is 
	 in turn used to estimate G, the tensor bi-spectrum. The integral is 
	 evaluated for e-fold N ranging from _Nics till _Nshss. Note that the 
	 extra factor exp(-(e*k)/(A*H)) is put in by hand to satisfy the consistency relation.'''
	N_range = numpy.linspace(N_ics, N_shss, len(hk_kx_array))

	AoverH = numpy.asarray([A(N) for N in N_range])/numpy.asarray([H(N) for N in N_range])
	factor = numpy.exp(-e*(k1+k2+k3)/(3*A(N_range)*numpy.asarray([H(N) for N in N_range])))
	func_int = AoverH*(hk_kx_array*hk_ky_array*hk_kz_array)*factor

	result = simps(func_int, N_range)
	return (+1/4.)*(kx**2+ky**2+kz**2)*result*numpy.array([0.+1.j], dtype=complex)

def main(k_set, N_array):
	k2, k3 = k_set

	if k2 > k3:
		Nics = solve_Nics(k3, N_array)
	else :
		Nics = solve_Nics(k2, N_array)

	Nshss = solve_Nshss(k1, N_array)
	step = N_array[1] -N_array[0]

	hk_k1_array = numpy.empty(0, dtype=complex)
	hk_k2_array = numpy.empty(0, dtype=complex)
	hk_k3_array = numpy.empty(0, dtype=complex)

	hk_k1_array = evolve_hk(k1, hk_k1, Dhk_k1, Nics, Nshss, step)
	hk_k2_array = evolve_hk(k2, hk_k2, Dhk_k2, Nics, Nshss, step)
	hk_k3_array = evolve_hk(k3, hk_k3, Dhk_k3, Nics, Nshss, step)

	tps_k1= 2.*(k1)**3/(2.*numpy.pi**2)*(numpy.absolute(hk_k1_array[-1]))**2
	tps_k2= 2.*(k2)**3/(2.*numpy.pi**2)*(numpy.absolute(hk_k2_array[-1]))**2
	tps_k3= 2.*(k3)**3/(2.*numpy.pi**2)*(numpy.absolute(hk_k3_array[-1]))**2

	CalG = calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)
	CalG_cc = calG_cc(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)

	G = (hk_k1_array[-1]*hk_k2_array[-1]*hk_k3_array[-1])*CalG +(numpy.conj(hk_k1_array[-1])*numpy.conj(hk_k2_array[-1])*numpy.conj(hk_k3_array[-1]))*CalG_cc
	h_NL = -((4./(2.*numpy.pi**2))**2)*(k1**3*k2**3*k3**3)*G/(2*k3**2*tps_k1*tps_k2 +2*k2**2*tps_k3*tps_k1 +2*k1**2*tps_k2*tps_k3)

	print k1, str(tps_k1).strip('[]'), k2, str(tps_k2).strip('[]'), k3, str(tps_k3).strip('[]'),
	print str(numpy.absolute(CalG)).strip('[]'), str(G.real).strip('[]'), str(h_NL.real).strip('[]')

	return None

#k_list = numpy.array([10**((-12 + i)/2.) for i in range(13)])

k1 = 1e-6
k3 = numpy.arange(0.,1.,0.1)*k1

k_list = []
'''Since k1 is fixed, k_list will contain the set of [k2, k3] values.'''
for i in range(len(k3)):
	if k3[i] < 0.5:
		k2 = numpy.linspace(1. -k3[i]/k1, 1.,int(k3[i]/k1/0.1))*k1
		[k_list.append([kx, k3[i]]) for kx in k2]
	else :
		k2 = numpy.linspace(k3[i]/k1, 1.,int((1 -k3[i]/k1) /0.1))*k1
		[k_list.append([kx, k3[i]]) for kx in k2]

#k_list = k_list[1]
print len(k_list)
#print k_list
pool = mp.Pool(processes =4)
#for i in range(len(k_list)):
#	k2, k3 = k_list[i]
#	print k2, k3

temp_results = [pool.apply_async(main, args = (k_set, N_array,)) for k_set in k_list[20:]]
results = []

for i in range(len(k_list)):
	results.append(temp_results[i].get())
#    k_vs_hk[i] = temp_results[i].get()

#results = numpy.asarray(results, dtype=numpy.float)
#print k_list, results
#print k_vs_hk
#print '\n'
