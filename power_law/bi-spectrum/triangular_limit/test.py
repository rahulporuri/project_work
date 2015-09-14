import numpy
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')

import random

k1 = 1.
k3 = numpy.arange(0.,1.,10**(-1))

plt.cla()
plt.hold(True)

#print k3

k_array = []

for i in range(len(k3)):
	if k3[i] < 0.5:
		k2 = numpy.arange(1.- k3[i]/k1, 1., 10**(-1))
		x = numpy.zeros(len(k2))
#		print k2
		[k_array.append([kx, k3[i]]) for kx in k2]
#		plt.plot(x, k2,'.')
	else:
		k2 = numpy.arange(k3[i]/k1, 1., 10**(-1))
#		print k2
		[k_array.append([kx, k3[i]]) for kx in k2]
#		plt.plot(x, k2,'.')
#plt.show()
#print k_array
print numpy.asarray(k_array)
'''
# In[12]:

#random.random()
k1 = 1
k3 = numpy.arange(0,1,10**(-1))

plt.cla()
plt.hold(True)

for i in range(len(k3)):
	if k3[i] < 0.5:
		k2 = numpy.arange(1- k3[i]/k1, 1, 10**(-1))
		x = numpy.zeros(len(k2))
		x[:] = k3[i]
		plt.plot(x, k2,'.')
	else:
		k2 = numpy.arange(k3[i]/k1, 1, 10**(-1))
		x = numpy.zeros(len(k2))
		x[:] = k3[i]
		plt.plot(x, k2,'.')
#plt.show()


# In[27]:

k1 = 1
k3 = numpy.arange(0,1,10**(-2))

matrix = []

for i in range(len(k3)):
    if k3[i] < 0.5:
        k2 = numpy.arange(1- k3[i]/k1, 1, 10**(-2))
        color = numpy.zeros(len(k2))
        color[:] = random.random()
        temp_matrix = numpy.empty(shape=3)
        temp_matrix[:,0] = k3[i]
        temp_matrix[:,1] = k2[i]
        temp_matrix[:,2] = color[i]
        matrix.append(temp_matrix)
    else:
        k2 = numpy.arange(k3[i]/k1, 1, 10**(-2))
        color = numpy.zeros(len(k2))
        color[:] = random.random()
        temp_matrix = numpy.empty(shape=3)
        temp_matrix[:,0] = k3[i]
        temp_matrix[:,1] = k2[i]
        temp_matrix[:,2] = color[i]
        matrix.append(temp_matrix)


# In[33]:

#numpy.empty(shape=3)
temp_matrix[:,0]


# In[20]:

numpy.array([[1,2,3],[4,5,6]])[:,0]


# In[39]:

numpy.meshgrid(numpy.array([[1,2]]),numpy.array([[1,2]]))


# In[61]:

x = numpy.arange(-5, 5, 0.1)
y = numpy.arange(-5, 5, 0.1)
xx, yy = numpy.meshgrid(x, y, sparse=True)
print len(xx), len(yy)
#print xx
z = numpy.sin(xx**2 + yy**2) / (xx**2 + yy**2)
#h = plt.imshow(x,y,z)


# In[64]:

#numpy.concatenate((xx,yy,z),axis=1)
#print len(xx), len(yy), len(z)
numpy.hstack((x,yy,z))


# In[ ]:




# In[ ]:




# In[133]:

#random.random()
k1 = 1
k3 = numpy.arange(0,1,10**(-1))

matrix = numpy.empty((1,3))
#matrix = []

for i in range(len(k3)):
    if k3[i] < 0.5:
        k2 = numpy.arange(1 -k3[i]/k1, 1, 10**(-1))
#        x = numpy.zeros(len(k2))
        #x[:] = k3[i]
        color = random.random()
        for j in range(len(k2)):
            #matrix = numpy.vstack((matrix, numpy.array([k3[i], k2[j], color])))
            matrix = numpy.vstack((matrix, numpy.array([color, k3[i], k2[j]])))
            #matrix.append(numpy.array([k3[i], k2[j]]))
    else:
        k2 = numpy.arange(k3[i]/k1, 1, 10**(-1))
        #x = numpy.zeros(len(k2))
        #x[:] = k3[i]
        for j in range(len(k2)):
#            matrix = numpy.append(matrix, numpy.array([k3[i],k2[j]]))
            #matrix = numpy.vstack((matrix, numpy.array([k3[i], k2[j], color])))
            matrix = numpy.vstack((matrix, numpy.array([color, k3[i], k2[j]])))
#plt.show()


# In[134]:

#plt.plot(matrix[:,0], matrix[:,1])
matrix
#numpy.vstack((numpy.zeros((1,2)), numpy.ones((1,2))))


# In[131]:

plt.scatter(matrix[:,0], matrix[:,1])


# In[138]:

plt.imshow(matrix)


# In[144]:

x = matrix[:,1]
y = matrix[:,2]
#triangles = [[0, 1, 4], [1, 2, 5], [2, 3, 6], [1, 5, 4], [2, 6, 5], [4, 5, 7],
#             [5, 6, 8], [5, 8, 7], [7, 8, 9]]
# for an upright triangle.
# need to figure this out for an inverted triangle,
# which is what we need!

triang = mtri.Triangulation(x,y,triangles)

import matplotlib.tri as mtri

plt.cla()
plt.hold(True)
plt.tricontourf(triang, matrix[:,0])
plt.triplot(triang, 'ko-')

'''
