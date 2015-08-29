import numpy
import matplotlib.pyplot as plt

k1 = 1
k3 = numpy.arange(0,1,10**(-2))

plt.cla()
plt.hold(True)

for i in range(len(k3)):
	if k3[i] < 0.5:
		k2 = numpy.arange(1- k3[i]/k1, 1, 10**(-2))
		x = numpy.zeros(len(k2))
		x[:] = k3[i]
		plt.plot(x, k2,'.')
	else:
		k2 = numpy.arange(k3[i]/k1, 1, 10**(-2))
		x = numpy.zeros(len(k2))
		x[:] = k3[i]
		plt.plot(x, k2,'.')

plt.show()
