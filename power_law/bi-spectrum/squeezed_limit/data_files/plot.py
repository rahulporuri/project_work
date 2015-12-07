import numpy
import matplotlib.pyplot as plt

#h_NL_py_data = numpy.loadtxt('h_NL_py.txt')
#h_NL_py_data = [line.split() for line in open('h_NL_py.txt')]
G_py_data = numpy.loadtxt('h_NL_py.txt')
G_m_data = numpy.loadtxt('p-law-ggg-one-sl-n.txt')

#print G_py_data[:,0]**6.*G_py_data[:,1]

#print h_NL_py_data

#k_py_data = numpy.asarray([h_NL_py_data[i][0] for i in  range(9)])
#h_NL_py = numpy.asarray([h_NL_py_data[i][2] for i in  range(9)])

#plt.xlim([0.,1.])
#plt.semilogx(k_py_data, h_NL_py)
#plt.show()

#G_py = [h_NL_py_ata[i][1] for i in range(9)]

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'$-G$')
plt.title(r'$-G$ as a function of $k$')
#plt.ylim([0.87,0.89])
#c, = plt.plot(h_NL_c_data[:,], h_NL_c_data[:,], '.')
python, = plt.loglog(G_py_data[:,0], (G_py_data[:,0]**3.)*(5e-07**3.)*numpy.absolute(G_py_data[:,4]), '.')
mathematica, = plt.loglog(G_m_data[:,0], numpy.absolute(G_m_data[:,1]), '--')
#m, = plt.plot(h_NL_m_data[:,], h_NL_m_data[:,], '.')
#plt.legend([c, python, m],['c', 'python', 'mathematica'])
plt.legend([python, mathematica],['mine','sreenath'])
plt.savefig('G_vs_k.png')

 data[:,0]**(4.*(-(2.*q-1.)/(q-1.)+2.))*10**(-17)

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'$h_{NL}$')
plt.title(r'$h_{NL}$ as a function of $k$')
#plt.ylim([0.87,0.89])
#c, = plt.plot(h_NL_c_data[:,], h_NL_c_data[:,], '.')
python, = plt.semilogx(G_py_data[:,0], G_py_data[:,5], '.')
#mathematica, = plt.loglog(h_NL_m_data[:,0], h_NL_m_data[:,1], '--')
#m, = plt.plot(h_NL_m_data[:,], h_NL_m_data[:,], '.')
#plt.legend([c, python, m],['c', 'python', 'mathematica'])
plt.legend([python,],['python'])
plt.savefig('h_NL_vs_k.png')
