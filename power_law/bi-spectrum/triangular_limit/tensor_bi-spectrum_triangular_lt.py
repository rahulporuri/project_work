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

H = lambda N : (V(phi(N))/(3. -Dphi(N)**2/2.))**(1./2)
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
'''
def calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss):
	N_range = numpy.linspace(Nics, Nshss, len(hk_k1_array))
	print len(N_range), len(A(N_range)), len(numpy.asarray([H(N) for N in N_range])), len(hk_k1_array), len(hk_k2_array), len(hk_k3_array)
	func_int = A(N_range)/numpy.asarray([H(N) for N in N_range])*(numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))*numpy.exp(-e*k0/(A(N_range)*numpy.asarray([H(N) for N in N_range])))
	#result = romb(func_int, N_Array[1]-N_Array[0])
	result = simps(func_int, N_range)
	return (-1/4.)*(k1**2+k2**2+k3**2)*result*numpy.array([0.+1.j], dtype=complex)

def calG_cc(hk_kx_array, hk_ky_array, hk_kz_array, kx, ky, kz, N_ics, N_shss):
	N_range = numpy.linspace(N_ics, N_shss, len(hk_kx_array))
	print len(N_range), len(A(N_range)), len(numpy.asarray([H(N) for N in N_range])), len(hk_kx_array), len(hk_ky_array), len(hk_kz_array)
	AoverH = numpy.asarray([A(N) for N in N_range])/numpy.asarray([H(N) for N in N_range])
	factor = numpy.exp(-e*k0/(numpy.asarray([A(N) for N in N_range])*numpy.asarray([H(N) for N in N_range])))
	print len(AoverH), len(factor)
	func_int = AoverH*(hk_kx_array*hk_ky_array*hk_kz_array)*factor
	#result = romb(func_int, N_Array[1]-N_Array[0])
	result = simps(func_int, N_range)
	return (+1/4.)*(kx**2+ky**2+kz**2)*result*numpy.array([0.+1.j], dtype=complex)
'''
def main(k_set, N_array):
	k2, k3 = k_set

	if k2 > k3:
		Nics = solve_Nics(k3, N_array)
	else :
		Nics = solve_Nics(k2, N_array)

	Nshss = solve_Nshss(k1, N_array)
	step = N_array[1] -N_array[0]

	hk_k1 = numpy.empty(0,dtype=complex) 
	Dhk_k1 = numpy.empty(0,dtype=complex)

	hk_k1 = initialize_hk(k1, Nics)
	Dhk_k1 = initialize_Dhk(k1, Nics)

	hk_k2 = numpy.empty(0,dtype=complex) 
	Dhk_k2 = numpy.empty(0,dtype=complex)

	hk_k2 = initialize_hk(k2, Nics)
	Dhk_k2 = initialize_Dhk(k2, Nics)

	hk_k3 = numpy.empty(0,dtype=complex) 
	Dhk_k3 = numpy.empty(0,dtype=complex)

	hk_k3 = initialize_hk(k3, Nics)
	Dhk_k3 = initialize_Dhk(k3, Nics)

	hk_k1_array = numpy.empty(0, dtype=complex)
	hk_k2_array = numpy.empty(0, dtype=complex)
	hk_k3_array = numpy.empty(0, dtype=complex)

	hk_k1_array = evolve_hk(k1, hk_k1, Dhk_k1, Nics, Nshss, step)
	hk_k2_array = evolve_hk(k2, hk_k2, Dhk_k2, Nics, Nshss, step)
	hk_k3_array = evolve_hk(k3, hk_k3, Dhk_k3, Nics, Nshss, step)

	tps_k1= 2.*(k1)**3/(2.*numpy.pi**2)*(numpy.absolute(hk_k1_array[-1]))**2
	tps_k2= 2.*(k2)**3/(2.*numpy.pi**2)*(numpy.absolute(hk_k2_array[-1]))**2
	tps_k3= 2.*(k3)**3/(2.*numpy.pi**2)*(numpy.absolute(hk_k3_array[-1]))**2

#	CalG = calG(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)
#	CalG_cc = calG_cc(hk_k1_array, hk_k2_array, hk_k3_array, k1, k2, k3, Nics, Nshss)

	N_range = numpy.linspace(Nics, Nshss, len(hk_k1_array))
#	print len(N_range), len(A(N_range)), len(numpy.asarray([H(N) for N in N_range])), len(hk_k1_array), len(hk_k2_array), len(hk_k3_array)
	func_int = A(N_range)/numpy.asarray([H(N) for N in N_range])*(numpy.conj(hk_k1_array)*numpy.conj(hk_k2_array)*numpy.conj(hk_k3_array))*numpy.exp(-e*k0/(A(N_range)*numpy.asarray([H(N) for N in N_range])))

	CalG =  (-1/4.)*(k1**2+k2**2+k3**2)*simps(func_int, N_range)*numpy.array([0.+1.j], dtype=complex)
	CalG_cc = numpy.conj(CalG_cc)

	G = (hk_k1_array[-1]*hk_k2_array[-1]*hk_k3_array[-1])*CalG +(numpy.conj(hk_k1_array[-1])*numpy.conj(hk_k2_array[-1])*numpy.conj(hk_k3_array[-1]))*CalG_cc
	h_NL = ((4./(2.*numpy.pi**2))**2)*(k1**3*k2**3*k3**3)*G/(2*k3**2*tps_k1*tps_k2 +2*k2**2*tps_k3*tps_k1 +2*k1**2*tps_k2*tps_k3)

	print k1, str(tps_k1).strip('[]'), k2, str(tps_k2).strip('[]'), k3, str(tps_k3).strip('[]'),
	print str(numpy.absolute(CalG)).strip('[]'), str(G.real).strip('[]'), str(h_NL.real).strip('[]')

	return None

#k_list = numpy.array([10**((-12 + i)/2.) for i in range(13)])

k1 = 1e-6
k3 = numpy.arange(0.,1.,0.1)*k1

k_list = []

for i in range(len(k3)):
	if k3[i] < 0.5:
		k2 = numpy.arange(1. -k3[i]/k1, 1., 0.1)*k1
		[k_list.append([kx, k3[i]]) for kx in k2]
	else :
		k2 = numpy.arange(k3[i]/k1, 1., 0.1)*k1
		[k_list.append([kx, k3[i]]) for kx in k2]

#k_list = k_list[1]
print len(k_list)
#print k_list
pool = mp.Pool(processes = 4)
#for i in range(len(k_list)):
#	k2, k3 = k_list[i]
#	print k2, k3

temp_results = [pool.apply_async(main, args = (k_set, N_array,)) for k_set in k_list]
results = []

for i in range(len(k_list)):
	results.append(temp_results[i].get())
#    k_vs_hk[i] = temp_results[i].get()

#results = numpy.asarray(results, dtype=numpy.float)
#print k_list, results
#print k_vs_hk
#print '\n'
