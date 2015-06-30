import numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt

#cimport routines
#from routines import *
#from routines cimport rk4_step
#from routines cimport DDhk


tps_file = open('power_spectrum_bouncing_model.dat','w')

k_min = 1e-30
k_max = 1e-4

Nshss = 3.02287778 

print 'lift off!' 

k_vs_hk = numpy.zeros(1,dtype=complex)

k0 = k_min
while k0 < k_max:
    print 'k0 = ', k0

    Nics = numpy.array([-1])
    for i in range(200):
        Nics = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),Nics)

    print Nics
    
    hk0 = numpy.zeros(1,dtype=complex)
    hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)

    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -(Nics/A(Nics))*((2*k0)**(-1./2))
    Dhk0.imag = -Nics*((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))
    
    hkm = numpy.array(hk0, dtype = complex)

    print 'got Nics, hk0 and Dhk0'

    npts = 10000
    step = (Nshss-Nics)/(npts)
    print 'starting from Nics'

    N = Nics
    while N < Nshss:
        array = rk4_step(k0, N, hk0, Dhk0, step)
        hk0 = hk0 + array[1]
        Dhk0 = Dhk0 + array[0]
        N += step

    k_vs_hk = numpy.append(k_vs_hk, hk0)
    
    temp = 8*(k0)**3/(2*numpy.pi**2)*(numpy.absolute(hk0))**2
    tps_file.write(str(k0)+"\t"+str(temp).strip('[]')+"\n")  

    print N, temp
    print '\n'
    
    k0 = 10*k0

k_list = numpy.array([10**(-30 + i) for i in range(27)])
TPS = [8*(k_list[i])**3/(2*numpy.pi**2)*(numpy.absolute(k_vs_hk[i+1]))**2 for i in range(len(k_list))]
print k_list, TPS

tps_file.close()

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal{P}}_{\rm T}(k)$')
plt.title(r'${\mathcal{P}}_{\rm T}(k)$ as a function of $k$')
numerics, = plt.loglog(k_list, TPS)
plt.legend([numerics],['numerical results'])
plt.savefig('power_spectrum_bouncing_model.png')
