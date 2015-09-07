import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data  = [line.split() for line in open('power_spectrum_power_law.dat')]

k_list = [data[i][0] for i in range(len(data))]
TPS = [data[i][1] for i in range(len(data))]

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal{P}}_{\rm T}(k)$')
plt.title(r'${\mathcal{P}}_{\rm T}(k)$ as a function of $k$')
numerics, = plt.loglog(k_list, TPS)
plt.legend([numerics],['numerical results'])
plt.savefig('power_spectrum_power_law.png')
