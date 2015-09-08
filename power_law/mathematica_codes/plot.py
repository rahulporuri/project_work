import matplotlib.pyplot as plt
import numpy

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#data  = [line.split() for line in open('power_spectrum_power_law.dat')]

c_data = [line.split() for line in open("phi_vs_N_c.txt")]
py_data = [line.split() for line in open("phi_vs_N_python.txt")]
m_data = [line.split(",") for line in open("phi_vs_N_mathematica.txt")]

#print c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2]

N_array_c = [c_data[i][0] for i in range(100000 +1)]
N_array_py = [py_data[i][0] for i in range(100000 +1)]
N_array_m = [m_data[i][1] for i in range(100000 +1)]
N_array_m = numpy.asarray(N_array_m, dtype=float)

phi_array_c = [c_data[i][2] for i in range(100000 +1)]
phi_array_py = [py_data[i][2] for i in range(100000 +1)]
phi_array_m = [m_data[i][2] for i in range(100000 +1)]

plt.cla()
plt.xlabel(r'$N$')
plt.ylabel(r'$\phi$')
plt.title(r'$\phi$ as a function of $N$')
c, = plt.plot(N_array_c, phi_array_c, '.')
python, = plt.plot(N_array_py, phi_array_py, '--')
m, = plt.plot(N_array_m, phi_array_m, '-')
plt.legend([c, python, m],['c', 'python', 'mathematica'])
plt.savefig('phi_vs_N.png')
