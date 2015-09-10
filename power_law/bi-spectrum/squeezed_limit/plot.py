import numpy
import matplotlib.pyplot as plt

#h_NL_py_data = numpy.loadtxt('h_NL_py.txt')
h_NL_py_data = [line.split() for line in open('h_NL_py.txt')]

print h_NL_py_data

k_py_data = numpy.asarray([h_NL_py_data[i][0] for i in  range(9)])
h_NL_py = numpy.asarray([h_NL_py_data[i][2] for i in  range(9)])

plt.xlim([0.,1.])
plt.semilogx(k_py_data, h_NL_py)
plt.show()

G_py = [h_NL_py_data[i][1] for i in range(9)]
