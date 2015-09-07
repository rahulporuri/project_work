import matplotlib.pyplot as plt
import numpy

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#data  = [line.split() for line in open('power_spectrum_power_law.dat')]

c_data = [line.split() for line in open("phi_vs_N_c.txt")]
py_data = [line.split() for line in open("phi_vs_N_python.txt")]

print c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2]

N_array_c = [c_data[i][0] for i in range(100000 +1)]
N_array_py = [py_data[i][0] for i in range(100000 +1)]
phi_array_c = [c_data[i][2] for i in range(100000 +1)]
phi_array_py = [py_data[i][2] for i in range(100000 +1)]

plt.cla()
plt.xlabel(r'$N$')
plt.ylabel(r'$\phi$')
plt.title(r'$\phi$ as a function of $N$')
c, = plt.plot(N_array_c, phi_array_c, '.')
python, = plt.plot(N_array_py, phi_array_py, '--')
plt.legend([c, python],['c', 'python'])
plt.savefig('phi_vs_N.png')
