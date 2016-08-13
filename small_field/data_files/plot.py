from __future__ import print_function

import matplotlib.pyplot as plt
import numpy

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#data  = [line.split() for line in open('power_spectrum_power_law.dat')]

#phi_c_data = [line.split() for line in open("phi_vs_N_c.txt")]
phi_py_data = [line.split() for line in open("phi_vs_N_python.txt")]
#phi_m_data = [line.split(",") for line in open("phi_vs_N_mathematica.txt")]

#print(c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2])

#N_array_c = [phi_c_data[i][0] for i in range(100000 +1)]
N_array_py = [phi_py_data[i][0] for i in range(100000 +1)]
#N_array_m = [phi_m_data[i][1] for i in range(100000 +1)]
#N_array_m = numpy.asarray(N_array_m, dtype=float)

#phi_array_c = [phi_c_data[i][2] for i in range(100000 +1)]
phi_array_py = [phi_py_data[i][1] for i in range(100000 +1)]
#phi_array_py_theory = [phi_py_data[i][4] for i in range(100000 +1)]
#phi_array_m = [phi_m_data[i][2] for i in range(100000 +1)]

plt.cla()
plt.xlabel(r'$N$')
plt.ylabel(r'$\phi$')
plt.title(r'$\phi$ as a function of $N$')
#c, = plt.plot(N_array_c, phi_array_c, '.')
python, = plt.plot(N_array_py, phi_array_py, '--')
#py_theory = plt.plot(N_array_py, phi_array_py_theory, '*')
#m, = plt.plot(N_array_py, phi_array_m, '-')
#plt.legend([c, python, m],['c', 'python', 'mathematica'])
#plt.legend([c, python],['c', 'python'], loc = 'upper left')
plt.legend([python],['numerical results'], loc = 'upper left')
plt.savefig('phi_vs_N.png')

##################################

#H_c_data = [line.split() for line in open("H_vs_N_c.txt")]
H_py_data = [line.split() for line in open("H_vs_N_python.txt")]
#H_m_data = [line.split(",") for line in open("H_vs_N_mathematica.txt")]

#print(c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2])

#N_array_c = [H_c_data[i][0] for i in range(100000 +1)]
N_array_py = [H_py_data[i][0] for i in range(100000 +1)]
#N_array_m = [H_m_data[i][1] for i in range(100000 +1)]
#N_array_m = numpy.asarray(N_array_m, dtype=float)

#H_array_c = [H_c_data[i][2] for i in range(100000 +1)]
H_array_py = [H_py_data[i][1] for i in range(100000 +1)]
#H_array_py_theory = [H_py_data[i][4] for i in range(100000 +1)]
#H_array_m = [H_m_data[i][2] for i in range(100000 +1)]

plt.cla()
plt.xlabel(r'$N$')
plt.ylabel(r'${\rm H}$')
plt.title(r'${\rm H}$ as a function of $N$')
#c, = plt.plot(N_array_c, H_array_c, '.')
python, = plt.plot(N_array_py, H_array_py, '--')
#py_theory, = plt.plot(N_array_py, H_array_py_theory, '*')
#m, = plt.plot(N_array_py, H_array_m, '-')
#plt.legend([c, python, m],['c', 'python', 'mathematica'])
plt.legend([python],['numerical results'], loc = 'upper left')
#plt.legend([c, python],['c', 'python'], loc = 'lower left')
plt.savefig('H_vs_N.png')

##################################

#DH_c_data = [line.split() for line in open("DH_vs_N_c.txt")]
DH_py_data = [line.split() for line in open("DH_vs_N_python.txt")]
#DH_m_data = [line.split(",") for line in open("DH_vs_N_mathematica.txt")]

#print(c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2])

#N_array_c = [DH_c_data[i][0] for i in range(100000 +1)]
N_array_py = [DH_py_data[i][0] for i in range(100000 +1)]
#N_array_m = [DH_m_data[i][1] for i in range(100000 +1)]
#N_array_m = numpy.asarray(N_array_m, dtype=float)

#DH_array_c = [DH_c_data[i][2] for i in range(100000 +1)]
DH_array_py = [DH_py_data[i][1] for i in range(100000 +1)]
#DH_array_py_theory = [DH_py_data[i][4] for i in range(100000 +1)]
#DH_array_m = [DH_m_data[i][2].strip(' ') for i in range(100000 +1)]

plt.cla()
plt.xlabel(r'$N$')
plt.ylabel(r'$\frac{{\rm d}H}{{\rm d}N}$')
plt.title(r'$\frac{{\rm d}H}{{\rm d}N}$ as a function of $N$')
#c, = plt.plot(N_array_c, DH_array_c, '.')
python, = plt.plot(N_array_py, DH_array_py, '--')
#py_theory, = plt.plot(N_array_py, DH_array_py_theory, '*')
#m, = plt.plot(N_array_py, numpy.asarray(DH_array_m, dtype=float)*(10**(-6)), '-')
#plt.legend([c, python, m],['c', 'python', 'mathematica'])
#plt.legend([c, python],['c', 'python'], loc = 'lower left')
plt.legend([python],['numerical results'], loc = 'upper left')
plt.savefig('DH_vs_N.png')

#################################
#eps_c_data = [line.split() for line in open("eps_vs_N_c.txt")]
eps_py_data = [line.split() for line in open("eps_vs_N_python.txt")]
#eps_m_data = [line.split(",") for line in open("eps_vs_N_mathematica.txt")]

#print(c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2])

#N_array_c = [eps_c_data[i][0] for i in range(100000 +1)]
N_array_py = [eps_py_data[i][0] for i in range(100000 +1)]
#N_array_m = [eps_m_data[i][1] for i in range(100000 +1)]
#N_array_m = numpy.asarray(N_array_m, dtype=float)

#eps_array_c = [eps_c_data[i][2] for i in range(100000 +1)]
eps_array_py = [eps_py_data[i][1] for i in range(100000 +1)]
#eps_array_py_theory = [eps_py_data[i][4] for i in range(100000 +1)]
#eps_array_m = [eps_m_data[i][2] for i in range(100000 +1)]

plt.cla()
plt.xlabel(r'$N$')
plt.ylabel(r'$\epsilon_1$')
plt.title(r'$\epsilon_1$ as a function of $N$')
#c, = plt.plot(N_array_c, eps_array_c, '.')
python, = plt.plot(N_array_py, eps_array_py, '--')
#py_theory, = plt.plot(N_array_py, eps_array_py_theory, '*')
#m, = plt.plot(N_array_py, eps_array_m, '-')
#plt.legend([c, python, m],['c', 'python', 'mathematica'], loc='upper left')
plt.legend([python],['numerical results'], loc = 'upper left')
#plt.legend([c, python],['c', 'python'])
plt.savefig('eps_vs_N.png')

#tps_c_data = [line.split() for line in open("tps_c.txt")]
tps_py_data = [line.split() for line in open("tps_python.txt")]
#tps_m_data = [line.split() for line in open("tps_mathematica.txt")]

#print(c_data[0][0], c_data[0][2], py_data[0][0], py_data[0][2])

#k_array_c = [tps_c_data[i][0] for i in range(11)]
k_array_py = [tps_py_data[i][0] for i in range(11)]
#k_array_m = [phi_m_data[i][1] for i in range(100000 +1)]
#k_array_m = numpy.asarray(N_array_m, dtype=float)

#tps_array_c = [tps_c_data[i][3] for i in range(11)]
tps_array_py = [tps_py_data[i][1] for i in range(11)]
#tps_array_m = [tps_m_data[i][3].strip('}') for i in range(11)]

#print(tps_array_py)

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal P}_T$')
plt.title(r'$k^{3/2}{\mathcal P}_T$ as a function of $k$')
#c, = plt.semilogx(numpy.asarray(k_array_py,dtype=float), tps_array_c, '.')
python, = plt.semilogx(numpy.asarray(k_array_py, dtype=float), tps_array_py, '--')
#m, = plt.semilogx(numpy.asarray(k_array_py, dtype=float), numpy.asarray(tps_array_m,dtype=float)*(10**(-10)), '-')
#plt.legend([c, python, m],['c', 'python', 'mathematica'])
#plt.legend([c, python],['c', 'python'], loc = 'center right')
plt.legend([python],['numerical results'], loc = 'upper left')
plt.show()
plt.savefig('tps.png')
