import matplotlib.pyplot as plt

data  = [line.split() for line in open('power_spectrum_power_law.dat')]

k_list = [data[i][0] for i in range(len(data))]
TPS = [data[i][1] for i in range(len(data))]

plt.cla()
plt.xlabel('k')
plt.ylabel('P(k)')
plt.title('P(k) vs k')
numerics, = plt.loglog(k_list, TPS)
plt.legend([numerics],['numerical results'])
plt.savefig('power_spectrum_power_law.png')

data  = [line.split() for line in open('power_spectrum_small_field.dat')]

k_list = [data[i][0] for i in range(len(data))]
TPS = [data[i][1] for i in range(len(data))]

plt.cla()
plt.xlabel('k')
plt.ylabel('P(k)')
plt.title('P(k) vs k')
numerics, = plt.loglog(k_list, TPS)
plt.legend([numerics],['numerical results'])
plt.savefig('power_spectrum_small_field.png')
