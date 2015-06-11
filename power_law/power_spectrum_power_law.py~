import numpy
import matplotlib.pyplot as plt

f = open('power_spectrum_power_law.dat','w')

q = 51.
V0 = (204./100.)*1e-08
t0 = (q*(3*q -1)/V0)**(1./2)

phi0 = 1.
dphi0 = (2.*q)**(1./2)/t0

Ni = 0.
Nf = 70. 

kp = 5.*1e-02
beta = -((2*q -1)/(q -1))
eps1a = ((beta +2)/(beta +1))

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

    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.]), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.]) # [Dhk, hk] update

npts = 20000
step = (Nf-Ni)/(npts)

phi_array = numpy.array([phi0])
Dphi_array = numpy.array([Dphi0])
N_array = numpy.array([Ni]) 

N = Ni
while N < Nf:
    array = rk4_step(N, phi0, Dphi0, step)
    phi0 = phi0 + array[1]
    Dphi0 = Dphi0 + array[0]
    phi_array = numpy.append(phi_array,phi0)
    Dphi_array = numpy.append(Dphi_array,Dphi0)
    N += step
    N_array = numpy.append(N_array,N)

#plt.plot(numpy.linspace(0,70,npts+1), phi_array)
#plt.plot(N_array, phi_array)

eps0 = (3./2)*((dphi0**2)/(dphi0**2/2. + V(phi0)))
eps = 1./q 

#H = [((V(phi_array[i]))/(3 -Dphi_array[i]**2/2))**(1./2) for i in range(len(phi_array))]

phi = lambda N : phi_array[int((N-Ni)/step)]
Dphi = lambda N : Dphi_array[int((N-Ni)/step)]

H = lambda N : ((V(phi(N))/(3 -Dphi(N))**2/2))**(1./2)
DH = lambda N : H(N)*Dphi(N)

ai = 1e-05

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


Nics_array = numpy.array([2.54198038578949765072641948787, 3.71629878321646094957559512976, 4.89061718064342424842477077165, 6.06493557807038754727394641354, 7.23925397549735084612312205543, 8.41357237292431414497229769732, 9.58789077035127744382147333808, 10.7622091677782407426704414080, 11.9365275652052040415198238617, 13.1108459626321673403690002649, 14.2851643600591306392181759068, 15.4594827574860939380673515487, 16.6338011549130572369165271905, 17.8081195523400205357657028324])

Nshs_array = numpy.array([18.9824379497669838346148777580, 20.1567563471939471334640541162, 21.3310747446209104323132297581, 22.5053931420478737311624054000, 23.6797115394748370300115810419, 24.8540299369018003288607566837, 26.0283483343287636277099323256, 27.2026667317557269265591079675, 28.3769851291826902254082834314, 29.5513035266096535242574592513, 30.7256219240366168231066348932, 31.8999403214635801219558105351, 33.0742587188905434208049861770, 34.2485771163175067196541618192])

Nics_arr = Ni + numpy.array((Nics_array-Ni)/step,dtype=int)*step
Nshss_arr = Ni + numpy.array((Nshs_array-Ni)/step,dtype=int)*step

k_min = 1e-6
k_max = 10

print 'lift off!'

k_vs_hk = numpy.zeros(1,dtype=complex)

i = 0
k0 = k_min

while k0 < k_max:
    print 'k0 = ', k0

    Nics = Nics_arr[i]
    Nshss = Nshss_arr[i]

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
    f.write(str(k0)+"\t"+str(hk0.real)+"\t"+str(hk0.imag)+"\t"+str(temp)+"\n")
    
    k0 = 10**(1./2)*k0
    i += 1

k_list = numpy.array([10**((-12 + i)/2.) for i in range(13)])
#print len(k_list), len(k_vs_hkhk)
TPS = [8*(k_list[i])**3/(2*numpy.pi**2)*(numpy.absolute(k_vs_hk[i+1]))**2 for i in range(len(k_list))]
print k_list, TPS

f.close()

plt.loglog(k_list, TPS)
plt.savefig('power_spectrum_power_law.png')
