import numpy as numpy
import matplotlib.pyplot as plt
import scipy.optimize as opt

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

tps_file = open('power_spectrum_bouncing_model.dat','w')

a0 = 1e-05
n0 = numpy.exp(25)
p = 1

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
