import numpy
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data = numpy.genfromtxt('equilateral_limit.dat')

'''

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal P}_T$')
plt.title(r'${\mathcal P}_T$ as a function of $k$')
python, = plt.loglog(data[:,0], data[:,1], '--')
plt.legend([python],['python'])
plt.savefig('tps_vs_k.png')

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\mathcal G}$')
plt.title(r'${\mathcal G}$ as a function of $k$')
python, = plt.loglog(data[:,0], data[:,2], '--')
plt.legend([python],['python'])
plt.savefig('calG_vs_k.png')

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'$k^{3/2}{\mathcal G}$')
plt.title(r'$k^{3/2}{\mathcal G}$ as a function of $k$')
python, = plt.loglog(data[:,0], data[:,0]**(3./2)*data[:,2], '--')
plt.legend([python],['python'])
plt.savefig('k32calG_vs_k.png')

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'G')
plt.title(r'G as a function of $k$')
python, = plt.loglog(data[:,0], -data[:,3], '--')
plt.legend([python],['python'])
plt.savefig('G_vs_k.png')

'''

q = 51.

print data[0][3]

plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'G')
plt.title(r'$k^6$G as a function of $k$')
python, = plt.loglog(data[:,0], -data[:,0]**(6.)*data[:,3], '--')
th, = plt.loglog(data[:,0], data[:,0]**(4.*(-(2.*q-1.)/(q-1.)+2.))*10**(-17),'x')
plt.legend([python, th],['python', 'th'])
plt.savefig('k6G_vs_k.png')
'''
plt.cla()
plt.xlabel(r'$k$')
plt.ylabel(r'${\rm h}_{NL}$')
plt.title(r'${\rm h}_{NL}$ as a function of $k$')
python, = plt.semilogx(data[:,0], data[:,4], '--')
plt.legend([python],['python'])
plt.ylim([0.47, 0.475])
plt.savefig('h_NL_vs_k.png')
'''
