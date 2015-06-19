import numpy
import matplotlib.pyplot as plt

tps_file = open('power_spectrum_power_law.dat','w')
phi_file = open('phi_vs_N_power_law.dat','w')
h_file = open('H_vs_N_power_law.dat','w')
eps_file = open('eps1_vs_N_power_law.dat','w')

q = 51.
V0 = (204./100.)*1e-08
t0 = (q*(3.*q -1.)/V0)**(1./2)

phi0 = 1.
dphi0 = (2.*q)**(1./2)/t0

Ni = 0.
Nf = 70. 

kp = 5.*1e-02
beta = -((2.*q -1.)/(q -1.))
eps1a = ((beta +2.)/(beta +1.))

#V = lambda phi : V0*numpy.exp(-(2*q)**(1./2.)*(phi-phi_i))
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

phi_array = numpy.array([phi_])
Dphi_array = numpy.array([Dphi_])
N_array = numpy.array([Ni]) 

phi_theory = lambda N : (2./q)**(1./2)*N + phi0

N = Ni
phi_file.write(str(N)+"\t"+str(phi_)+"\t"+str(Dphi_)+"\t"+str(phi_theory(N))+"\n")
while N < Nf:
    array = rk4_step(N, phi_, Dphi_, step)
    phi_ = phi_ + array[1]
    Dphi_ = Dphi_ + array[0]
    phi_array = numpy.append(phi_array,phi_)
    Dphi_array = numpy.append(Dphi_array,Dphi_)
    N += step
    N_array = numpy.append(N_array,N)
    phi_file.write(str(N)+"\t"+str(phi_)+"\t"+str(Dphi_)+"\t"+str(phi_theory(N))+"\n")

phi_file.close()
#plt.plot(numpy.linspace(0,70,npts+1), phi_array)

phi = lambda N : phi_array[int((N-Ni)/step)]
Dphi = lambda N : Dphi_array[int((N-Ni)/step)]

plt.cla()
plt.hold()
#plt.plot(N_array, [phi(N) for N in N_array],'*')
plt.xlim([Ni,Nf])
numerical, = plt.plot(N_array, phi_array, label = 'numerical results')
theory, = plt.plot(N_array, [phi_theory(N) for N in N_array],'*', label = 'theory results')
plt.legend([numerical, theory], ['numerical results', 'theoretical results'])
plt.savefig('phi_vs_N_power_law.png')
#plt.show()

eps0 = (3./2)*((dphi0**2)/(dphi0**2/2. + V(phi0)))
eps = 1./q 

#H = [((V(phi_array[i]))/(3 -Dphi_array[i]**2/2))**(1./2) for i in range(len(phi_array))]

H = lambda N : (V(phi(N))/(3 -Dphi(N)**2/2))**(1./2)
DH = lambda N : H(N)*Dphi(N)

H_theory = lambda N : H0*numpy.exp(-N/q)

for N in N_array:
	h_file.write(str(N)+"\t"+str(H(N)/H0)+"\t"+str(H_theory(N)/H0)+"\n")

plt.cla()
plt.hold()
plt.xlim([Ni,Nf])
numerical, = plt.plot(N_array, numpy.asarray([H(i) for i in N_array], dtype= numpy.float64)/H0, '*', label = 'numerical results')
theory, = plt.plot(N_array, [H_theory(N)/H0 for N in N_array], label = 'theory')
plt.legend([numerical, theory], ['numerical results', 'theoretical results'])
plt.savefig('H_vs_N_power_law.png')

ai = 1e-05

eps1 = lambda N : Dphi_array[int((N-Ni)/step)]**2/2.
eps1_theory = eps0

for N in N_array:
	eps_file.write(str(N)+"\t"+str(eps1(N))+"\t"+str(eps1_theory)+"\n")

plt.cla()
plt.xlim([Ni,Nf])
numerical, = plt.plot(N_array, [str(eps1(i)).strip('[]') for i in N_array], '*', label = 'numerical results')
#theory, = plt.axhline(y=eps0, label = 'theory')
plt.axhline(y=eps0)
plt.legend([numerical], ['numerical results'])
plt.savefig('eps1_vs_N_power_law.png')

#z = [ai*numpy.exp(N_array[i])*Dphi_array[i] for i in range(len(N_array))]
z = lambda N: ai*numpy.exp(N)*Dphi(N)
A = lambda N : ai*numpy.exp(N)

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


k_list = numpy.array([10**((-12 + i)/2.) for i in range(13)])
Nics_array = []
Nshs_array = []

for k in k_list:
    #temp = numpy.asarray([k/(A(N)*H(N)) for N in N_array])
    Nics_temp = numpy.asarray([k - 1e+02*A(N)*H(N) for N in N_array])
    Nshss_temp = numpy.asarray([k - 1e-05*A(N)*H(N) for N in N_array])
    nics_test = numpy.where(Nics_temp > 0)
    nshss_test = numpy.where(Nshss_temp > 0)
    Nics_array.append(Ni + nics_test[0][-1]*step)
    Nshs_array.append(Ni + nshss_test[0][-1]*step)

Nics_arr = numpy.asarray(Nics_array)
Nshs_arr = numpy.asarray(Nshs_array)

k_min = 1e-6
k_max = 10**(1./2)

print 'lift off!'

k_vs_hk = numpy.zeros(1,dtype=complex)

i = 0
k0 = k_min

while k0 < k_max:
    print 'k0 = ', k0

    Nics = Nics_arr[i]
    Nshss = Nshs_arr[i]

    hk0 = numpy.zeros(1,dtype=complex)
    hk0.real = (((2.*k0)**(1./2))*A(Nics))**(-1.)

    Dhk0 = numpy.zeros(1,dtype=complex)
    Dhk0.real = -(1/A(Nics))*((2*k0)**(-1./2))
    Dhk0.imag = -((k0/2)**(1./2))/(A(Nics)*A(Nics)*H(Nics))
 
    print 'got Nics, hk0 and Dhk0'
    print 'starting from Nics'

    N = Nics
    while N < Nshss:
        #array = euler_step()
        array = rk4_step(k0, N, hk0, Dhk0, step)
        hk0 = hk0 + array[1]
        Dhk0 = Dhk0 + array[0]
        N += step

    k_vs_hk = numpy.append(k_vs_hk, hk0) 
    print N-step, Nshss, hk0, Dhk0, Nics
    print '\n'
    
    temp = 8*(k0)**3/(2*numpy.pi**2)*(numpy.absolute(hk0))**2
    tps_file.write(str(k0)+"\t"+str(temp).strip('[]')+"\n")
#    tps_file.write(str(k0)+"\t"+str(hk0.real)+"\t"+str(hk0.imag)+"\t"+str(temp).strip('[]')+"\n")
    k0 = 10**(1./2)*k0
    i += 1

k_list = numpy.array([10**((-12 + i)/2.) for i in range(13)])
#print len(k_list), len(k_vs_hkhk)
TPS = [8*(k_list[i])**3/(2*numpy.pi**2)*(numpy.absolute(k_vs_hk[i+1]))**2 for i in range(len(k_list))]
print k_list, TPS

tps_file.close()

plt.cla()
plt.xlabel('k')
plt.ylabel('P(k)')
plt.title('P(k) vs k')
numerics, = plt.loglog(k_list, TPS)
plt.legend([numerics],['numerical results'])
plt.savefig('power_spectrum_power_law.png')
