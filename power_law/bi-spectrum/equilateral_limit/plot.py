import numpy
import matplotlib.pyplot as plt
'''
data = numpy.loadtxt('test.dat')

plt.semilogx(data[:,0], data[:,8])
plt.show()
'''
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#data  = [line.split() for line in open('power_spectrum_power_law.dat')]

#h_NL_c_data = [line.split() for line in open("phi_vs_N_c.txt")]
h_NL_py_data = numpy.loadtxt('test.dat')
#h_NL_m_data = [line.split(",") for line in open("phi_vs_N_mathematica.txt")]

#print c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2]

#N_array_c = [phi_c_data[i][0] for i in range(100000 +1)]
#N_array_py = [phi_py_data[i][0] for i in range(100000 +1)]
#N_array_m = [phi_m_data[i][1] for i in range(100000 +1)]
#N_array_m = numpy.asarray(N_array_m, dtype=float)

#phi_array_c = [phi_c_data[i][2] for i in range(100000 +1)]
#phi_array_py = [phi_py_data[i][2] for i in range(100000 +1)]
#phi_array_m = [phi_m_data[i][2] for i in range(100000 +1)]

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'$h_{NL}$')
plt.title(r'$h_{NL}$ as a function of $k$')
plt.ylim([0.87,0.89])
#c, = plt.plot(h_NL_c_data[:,], h_NL_c_data[:,], '.')
python, = plt.semilogx(h_NL_py_data[:,0], h_NL_py_data[:,8], '--')
#m, = plt.plot(h_NL_m_data[:,], h_NL_m_data[:,], '.')
#plt.legend([c, python, m],['c', 'python', 'mathematica'])
plt.legend([python],['python'])
plt.savefig('h_NL_vs_k.png')
