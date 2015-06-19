import numpy
import matplotlib.pyplot as plt

tps_file = open('power_spectrum_small_field.dat','w')
phi_file = open('phi_vs_N_small_field.dat','w')
h_file = open('H_vs_N_small_field.dat','w')
eps_file = open('eps1_vs_N_small_field.dat','w')

p = 4.
mu = 15.
V0 = 5.55702*1e-10
phi0 = 7.3
dphi0 = 8*1e-07

Ni = 0.
Nf = 70.

#V = lambda phi : V0*numpy.exp(-(2*q)**(1./2.)*(phi-phi_i))
V = lambda phi : V0*(1-(phi/mu)**p)
dV = lambda phi : -V0*(-p*(phi/mu)**(p-1))*(1./mu)

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

N = Ni
phi_file.write(str(N)+"\t"+str(phi_)+"\t"+str(Dphi_)+"\n")
while N < Nf:
    array = rk4_step(N, phi_, Dphi_, step)
    phi_ = phi_ + array[1]
    Dphi_ = Dphi_ + array[0]
    phi_array = numpy.append(phi_array,phi_)
    Dphi_array = numpy.append(Dphi_array,Dphi_)
    N += step
    N_array = numpy.append(N_array,N)
    phi_file.write(str(N)+"\t"+str(phi_)+"\t"+str(Dphi_)+"\n")

phi_file.close()

#plt.plot(numpy.linspace(0,70,npts+1), phi_array)
plt.cla()
plt.xlabel('e-fold N')
plt.ylabel('phi(N)')
plt.title('plot of phi(N) vs e-fold N')
numerics, = plt.plot(N_array, phi_array, '*', label = 'numerical results')
plt.legend([numerics],['numerical results'])
plt.savefig('phi_vs_N_small_field.png')

eps0 = (3./2)*((dphi0**2)/(dphi0**2/2. + V(phi0)))
#eps = 1./q 

#H = [((V(phi_array[i]))/(3 -Dphi_array[i]**2/2))**(1./2) for i in range(len(phi_array))]

phi = lambda N : phi_array[int((N-Ni)/step)]
Dphi = lambda N : Dphi_array[int((N-Ni)/step)]

H = lambda N : (V(phi(N))/(3 -Dphi(N)**2/2))**(1./2)
DH = lambda N : H(N)*Dphi(N)

for N in N_array:
	h_file.write(str(N)+"\t"+str(H(N)/H0)+"\t"+"\n")

plt.cla()
plt.xlabel('e-fold N')
plt.ylabel('H(N)')
plt.title('plot of Hubble parameter H(N) vs e-fold N')
numerics, = plt.plot(N_array, numpy.asarray([str(H(i)).strip('[]') for i in N_array], dtype= numpy.float64)/H0, '*', label = 'numerical results')
plt.legend([numerics],['numerical results'])
plt.savefig('H_vs_N_small_field.png')

ai = 1e-05

eps1 = lambda N : Dphi_array[int((N-Ni)/step)]**2/2.

for N in N_array:
	eps_file.write(str(N)+"\t"+str(eps1(N))+"\t"+str(eps0)+"\n")

plt.cla()
plt.xlabel('e-fold N')
plt.ylabel('eps(N)')
plt.title('eps(N) vs e-fold N')
numerics, = plt.plot(N_array, [str(eps1(i)).strip('[]') for i in N_array], '*', label = 'numerical results')
plt.legend([numerics],['numerical results'])
plt.savefig('eps1_vs_N_small_field.png')

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

Nics_array = numpy.asarray(Nics_array)
Nshs_array = numpy.asarray(Nshs_array)

#Nics_array = numpy.array([2.54198038578949765072641948787, 3.71629878321646094957559512976, 4.89061718064342424842477077165, 6.06493557807038754727394641354, 7.23925397549735084612312205543, 8.41357237292431414497229769732, 9.58789077035127744382147333808, 10.7622091677782407426704414080, 11.9365275652052040415198238617, 13.1108459626321673403690002649, 14.2851643600591306392181759068, 15.4594827574860939380673515487, 16.6338011549130572369165271905, 17.8081195523400205357657028324])

#Nshs_array = numpy.array([18.9824379497669838346148777580, 20.1567563471939471334640541162, 21.3310747446209104323132297581, 22.5053931420478737311624054000, 23.6797115394748370300115810419, 24.8540299369018003288607566837, 26.0283483343287636277099323256, 27.2026667317557269265591079675, 28.3769851291826902254082834314, 29.5513035266096535242574592513, 30.7256219240366168231066348932, 31.8999403214635801219558105351, 33.0742587188905434208049861770, 34.2485771163175067196541618192])

k_min = 1e-6
k_max = 10

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
    tps_file.write(str(k0)+"\t"+str(hk0.real)+"\t"+str(hk0.imag)+"\t"+str(temp)+"\n")
    
    k0 = 10**(1./2)*k0
    i += 1

k_list = numpy.array([10**((-12 + i)/2.) for i in range(13)])
#print len(k_list), len(k_vs_hkhk)
TPS = [8*(k_list[i])**3/(2*numpy.pi**2)*(numpy.absolute(k_vs_hk[i+1]))**2 for i in range(len(k_list))]
print k_list, TPS

f.close()

plt.loglog(k_list, TPS)
plt.savefig('power_spectrum_small_field.png')
