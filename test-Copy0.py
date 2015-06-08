
# coding: utf-8

# $$a(\mathcal{N}) = a_0*e^{\frac{\mathcal{N}^2}{2}}$$
# 
# $$a(\eta) = a_0*(1+(\frac{\eta}{\eta_0})^2)^p$$
# 
# $$\eta^+(\mathcal{N}) = \eta_0\left[\left(\frac{a(\mathcal{N})}{a_0}\right)^{\frac{1}{p}}-1\right]^{\frac{1}{2}}$$
# $$\eta^-(\mathcal{N}) = -\eta_0\left[\left(\frac{a(\mathcal{N})}{a_0}\right)^{\frac{1}{p}}-1\right]^{\frac{1}{2}}$$
# 
# <!--i could also use \big(, \Big(, \bigg(, \Bigg(-->
# 
# $$\eta(\mathcal{N}) = \theta(\mathcal{N})*\eta^+(\mathcal{N}) + \theta(\mathcal{-N})*\eta^-(\mathcal{N})$$
# 
# $$a(\mathcal{N}) = a(\eta(\mathcal{N}))$$
# 
# $$\mathcal{H}(\eta) = \frac{a^{'}(\eta)}{a(\eta)^2}$$

# function prefix of 'd_' for differentiation with respect to $\eta$ and 'D_' for differentiation with respect to $\mathcal{N}$. Lower case function names for functions in terms of $\eta$ and higher case function names for functions in terms of $\mathcal{N}$

# $$\frac{d^2h_k}{d\mathcal{N}^2} + \left(3\mathcal{N} - \frac{1}{\mathcal{N}} + \frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{d\mathcal{N}}\right)\frac{dh_k}{d\mathcal{N}} + \left(\frac{k_s\mathcal{N}}{a\mathcal{H}}\right)^2 h_k = 0$$
# 
# where $h = h(\mathcal{N})$, $a = a(\eta(\mathcal{N}))$ and $\mathcal{H} = \mathcal{H}(\eta(\mathcal{N}))$

# In[1]:

import numpy as numpy
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
import sympy as sym
x = sym.symbols('x')
# references -
# http://docs.scipy.org/doc/scipy-0.14.0/reference/optimize.html
#import scipy
import scipy.optimize as opt


# In[2]:

a0 = 1e-05
n0 = numpy.exp(25)
k0 = numpy.exp(-25)
p = 1


# In[3]:

#sym.diff(np(x),x)
#sym.diff(np(x),x)
from sympy.functions import exp


# In[4]:

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)
nm = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : nm(N)
#aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))


# In[5]:

DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()


# In[6]:

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))


# In[7]:

aN = lambda N : an(n(N))


# In[8]:

N = numpy.linspace(-5,0,51)


# the behavior of a($\eta$) for various $\eta$
# ===

# In[9]:

plt.plot(numpy.arange(-1e06, 1e06,100),an(numpy.arange(-1e06, 1e06,100)))


# the behavior of $\mathcal{H}(\mathcal{N})$ for various $\mathcal{N}$
# ===

# In[10]:

plt.plot(N,[H(i) for i in N])


# the behavior of $\mathcal{H}^{'}(\mathcal{N})$ for various $\mathcal{N}$
# ===

# In[11]:

plt.plot(N,[DH(i) for i in N])


# the behavior of $\frac{a^{''}(\mathcal{N})}{a(\mathcal{N})}$ for various $\mathcal{N}$
# ===

# In[12]:

plt.plot(N,[fN(i) for i in N])


# finding the roots and initial conditions to solve the diff. eqn
# ===
# 
# $$h(Nics) = \frac{1}{(2k_s)^\frac{1}{2}*a(Nics)}$$
# $$h^{'}(Nics) = -\frac{I*Nics*(2ks)^\frac{1}{2}}{2*a(Nics)*a(Nics)*H(Nics)} - \frac{a^{'}(Nics)}{a^2(Nics)*(2ks)^\frac{1}{2}}$$
# 
# results from mathematica -
# 
# h0 = 948715.57472053941534264
# dh0 = 4.2222613144777582640577*1e+06

# In[14]:

temp = numpy.array([-4])
for i in range(10):
    temp = opt.newton_krylov(lambda N : (1./n0)**2 - 1e+04*fN(N),temp)
print temp


# evolution of $h_k$ from Nics to 0 for $k_0 = exp(-25)$
# ====

# In[14]:

def DDhk(N, hk0, Dhk0):
    return -((3.*N -(1./N) +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(A(N)*H(N)))**2)*hk0)
'''
def euler_step(N, hk0, Dhk0, step):
    F = Dhk0
    f = DDhk(N,hk0,Dhk0)
    
    return [f*step, F*step] #[Dhk, hk] update

def rk2_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N,Dhk0,hk0)
    F2 = Dhk0 +f1*step
    f2 = DDhk(N+step/2.,Dhk0+f1*step/2.,hk0+F1*step/2.)
    
    return [(f1+f2)*step/2., (F1+F2)*step/2.] # [Dhk, hk] update
'''
def rk4_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(N +step, hk0 +F3*step, Dhk0 +f3*step)
    
    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

npts = 1000
Nics = -4.450498032115866
Nshss = 0

step = (Nshss-Nics)/(npts-1)

hk0 = numpy.array([948715.57472053941534264 + 0.j], dtype = complex)
Dhk0 = numpy.array([4222261.3144777582640577 -298566424.99481408456441j], dtype = complex)

hk = numpy.array([948715.57472053941534264 + 0.j], dtype = complex)

N = Nics
while N < Nshss-2*step:
    #array = euler_step(N, hk0, Dhk0, step)
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hk = numpy.append(hk,hk0)
    N += step


# In[15]:

print N-step, hk0, Dhk0


# In[16]:

plt.plot(numpy.absolute(hk))


# In[17]:

plt.semilogy(numpy.absolute(hk))


# In[18]:

plt.plot(hk.imag)


# In[19]:

N = numpy.linspace(0,5,51)


# In[20]:

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)
np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : np(N)
#aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))


# In[21]:

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))


# In[22]:

aN = lambda N : an(n(N))


# In[23]:

DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()


# the behavior of a($\eta$) for various $\eta$
# ===

# In[24]:

plt.plot(numpy.arange(-1e06, 1e06,100),an(numpy.arange(-1e06, 1e06,100)))


# the behavior of $\mathcal{H}(\mathcal{N})$ for various $\mathcal{N}$
# ===

# In[25]:

plt.plot(N,[H(i) for i in N])


# the behavior of $\mathcal{H}^{'}(\mathcal{N})$ for various $\mathcal{N}$
# ===

# In[26]:

plt.plot(N,[DH(i) for i in N])


# In[27]:

plt.plot(N,[fN(i) for i in N])


# finding the roots and initial conditions to solve the diff. eqn
# ===
# 
# $$h(Nics) = \frac{1}{(2k_s)^\frac{1}{2}*a(Nics)}$$
# $$h^{'}(Nics) = -\frac{I*Nics*(2ks)^\frac{1}{2}}{2*a(Nics)*a(Nics)*H(Nics)} - \frac{a^{'}(Nics)}{a^2(Nics)*(2ks)^\frac{1}{2}}$$
# 
# results from mathematica -
# 
# h0 = 948715.57472053941534264
# dh0 = 4.2222613144777582640577*1e+06

# In[28]:

temp = numpy.array([5])
for i in range(10):
    temp = opt.newton_krylov(lambda N : (1./n0)**2 - 1e+06*fN(N),temp)
    
print temp


# evolution of $h_k$ from 0 to Nshss for $k_0 = exp(-25)$
# ====

# In[29]:

def DDhk(N,hk0,Dhk0):
    return -((3*N -1./N +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(aN(N)*H(N)))**2)*hk0)
'''
def euler_step():
    F = Dhk0
    f = DDhk(N,Dhk0,hk0)
    return [f*step, F*step] #[dh, h] update
'''
def rk4_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(N +step, hk0 +F3*step, Dhk0 +f3*step)
    
    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

npts = 1000

Nics = 0
Nshss = 5.386772255310419

step = (Nshss-Nics)/(npts-1)

hk0 = numpy.array([  1.15806723e+09 -1.92903590e+10j], dtype = complex) 
Dhk0 = numpy.array([ -2.60841343e+09 -1.06553978e+10j], dtype = complex) 

hk = numpy.array([  1.15806723e+09 -1.92903590e+10j], dtype = complex)

N = Nics+step
while N < Nshss:
    #array = euler_step()
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hk = numpy.append(hk, hk0)
    N += step


# In[30]:

plt.plot(numpy.absolute(hk))


# In[32]:

plt.semilogy(numpy.absolute(hk[:850]))


# In[34]:

plt.plot(hk.imag[:850])


# In[35]:

print hk0, Dhk0


# In[ ]:




# references to help handle complex numbers using numpy
# http://docs.scipy.org/doc/numpy/reference/routines.math.html#handling-complex-numbers
# https://docs.python.org/2/library/cmath.html

# evolution of $h_k$ from Nics to Nshss for $k_0 = 10^{-5}*exp^{-25}$
# ====

# In[15]:

k0 = 1e-5*numpy.exp(-25)


# In[ ]:

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)
nn = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : nn(N)
aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))

DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()


# In[19]:

temp = numpy.array([-8])
for i in range(20):
    temp = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),temp)

print temp


# In[20]:

Nics = -8.11533591


# In[21]:

print 'hk0', (((2.*k0)**(1./2))*A(Nics))**(-1.)
print 'imag part dhk0', -Nics*((k0/2)**(1/2))/(A(Nics)*A(Nics)*H(Nics))
print 'real part dhk0', -(H(Nics))*((2*k0)**(-1/2))
print 'real part dhk0', -(sym.diff(an(x),x).subs(x,n(Nics)).evalf()/(A(Nics)**2))*((2*k0)**(-1/2))


# In[22]:

def DDhk(N, hk0, Dhk0):
    return -((3.*N -(1./N) +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(A(N)*H(N)))**2)*hk0)

'''

def euler_step(N, hk0, Dhk0, step):
    F = Dhk0
    f = DDhk(N,hk0,Dhk0)    

    return [f*step, F*step] #[Dhk, hk] update

 

def rk2_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N,Dhk0,hk0)
    F2 = Dhk0 +f1*step
    f2 = DDhk(N+step/2.,Dhk0+f1*step/2.,hk0+F1*step/2.)

    return [(f1+f2)*step/2., (F1+F2)*step/2.] # [Dhk, hk] update
'''

def rk4_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(N +step, hk0 +F3*step, Dhk0 +f3*step)

    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

npts = 1000

Nics = -8.115335912021952142866641
Nshss = 0

step = (Nshss-Nics)/(npts-1)

hk0 = numpy.array([0.030001020677925666319695 + 0.j], dtype = complex)
Dhk0 = numpy.array([0.24346836050488333224282 - 17.215812871737443813285j], dtype = complex)

hk = numpy.array([0.030001020677925666319695 + 0.j], dtype = complex)

N = Nics
while N < Nshss-2*step:
    #array = euler_step(N, hk0, Dhk0, step)
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hk = numpy.append(hk,hk0)
    N += step


# In[23]:

print N-step, hk0, Dhk0


# In[24]:

plt.plot(numpy.absolute(hk))


# In[25]:

plt.semilogy(numpy.absolute(hk))


# In[26]:

plt.plot(hk.imag)


# In[ ]:




# In[32]:

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)
np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : np(N)
aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))

DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()


# In[33]:

temp = numpy.array([8])
for i in range(20):
    temp = opt.newton_krylov(lambda N : (k0)**2 - 1e+06*fN(N),temp)

print temp


# In[29]:

Nshss = 8.66423784


# In[34]:

def DDhk(N,hk0,Dhk0):
    return -((3*N -1./N +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(aN(N)*H(N)))**2)*hk0)
'''
def euler_step():
    F = Dhk0
    f = DDhk(N,Dhk0,hk0)
    return [f*step, F*step] #[dh, h] update
'''
def rk4_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(N +step, hk0 +F3*step, Dhk0 +f3*step)
    
    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

npts = 1000

Nics = 0
Nshss = 8.664237839356059165326164

step = (Nshss-Nics)/(npts-1)

hk0 = numpy.array([-5.18773965e+24 -1.26789704e+25j], dtype = complex) 
Dhk0 = numpy.array([-4.73915369e+24 -1.15826147e+25j], dtype = complex) 

hk = numpy.array([-5.18773965e+24 -1.26789704e+25j], dtype = complex)

N = Nics+step

while N < Nshss:
    #array = euler_step()
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hk = numpy.append(hk, hk0)
    N += step


# In[44]:

plt.plot(numpy.absolute(hk[:900]))


# In[41]:

plt.semilogy(numpy.absolute(hk[:900]))


# In[40]:

plt.plot(hk.imag[:900])


# In[ ]:




# In[ ]:




# evolution of $h_k$ from Nics to Nshss for $k_0 = 10^{-10}*exp^{-25}$
# ====

# In[45]:

k0 = 1e-10*numpy.exp(-25)


# In[46]:

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)
nn = lambda N : -n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : nn(N)
aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))

DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()


# In[48]:

temp = numpy.array([-10])
for i in range(20):
    temp = opt.newton_krylov(lambda N : (k0)**2 - 1e+04*fN(N),temp)

print temp


# In[49]:

Nics = -10.57877019


# In[50]:

print 'hk0', (((2.*k0)**(1./2))*A(Nics))**(-1.)
print 'imag part dhk0', -Nics*((k0/2)**(1/2))/(A(Nics)*A(Nics)*H(Nics))
print 'real part dhk0', -(H(Nics))*((2*k0)**(-1/2))
print 'real part dhk0', -(sym.diff(an(x),x).subs(x,n(Nics)).evalf()/(A(Nics)**2))*((2*k0)**(-1/2))


# In[62]:

def DDhk(N,hk0,Dhk0):
    return -((3*N -1./N +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(aN(N)*H(N)))**2)*hk0)

'''
def euler_step():
    F = Dhk0
    f = DDhk(N,Dhk0,hk0)
    return [f*step, F*step] #[dh, h] update
'''
def rk4_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(N +step, hk0 +F3*step, Dhk0 +f3*step)

    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

npts = 1000

Nics = -10.57877019434839994065461
Nshss = 0

step = (Nshss-Nics)/(npts-1)

hk0 = numpy.array([9.487155747205394153426*1e-10 + 0.j], dtype = complex) 
Dhk0 = numpy.array([1.0036244044767754696560*1e-8 -7.096696221698383396916*1e-7j], dtype = complex)

hk = numpy.array([9.487155747205394153426*1e-10 + 0j], dtype = complex)

N = Nics
while N < Nshss - 2*step:
    #array = euler_step()
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hk = numpy.append(hk, hk0)
    N += step


# In[63]:

plt.plot(numpy.absolute(hk))


# In[64]:

plt.semilogy(numpy.absolute(hk))


# In[65]:

plt.plot(hk.imag)


# In[66]:

print N-step, hk0, Dhk0


# In[ ]:




# In[55]:

A = lambda N : a0*exp(x**2/2).subs(x,N).evalf()
an = lambda n : a0*((1+(n/n0)**2)**p)
np = lambda N : n0*((A(N)/a0)**(1./p)-1)**(1./2)
n = lambda N : np(N)
aN = lambda N : an(n(N))
h = lambda n : (sym.diff(an(x),x).subs(x,n).evalf())/(an(n)**2)
H = lambda N : h(n(N))

fn = lambda n : sym.diff(sym.diff(an(x),x),x).subs(x,n).evalf()/an(n)
fN = lambda N : fn(n(N))

DH = lambda N : sym.diff(H(x),x).subs(x,N).evalf()


# In[57]:

temp = numpy.array([10])
for i in range(20):
    temp = opt.newton_krylov(lambda N : (k0)**2 - 1e+06*fN(N),temp)
    
print temp


# In[ ]:

Nshss = 11.00548587


# In[ ]:

print 'hk0', (((2.*k0)**(1./2))*A(Nics))**(-1.)
print 'imag part dhk0', -Nics*((k0/2)**(1/2))/(A(Nics)*A(Nics)*H(Nics))
print 'real part dhk0', -(H(Nics))*((2*k0)**(-1/2))
print 'real part dhk0', -(sym.diff(an(x),x).subs(x,n(Nics)).evalf()/(A(Nics)**2))*((2*k0)**(-1/2))


# In[67]:

def DDhk(N,hk0,Dhk0):
    return -((3*N -1./N +(DH(N)/H(N)))*Dhk0 +(((k0*N)/(aN(N)*H(N)))**2)*hk0)

'''
def euler_step():
    F = Dhk0
    f = DDhk(N,Dhk0,hk0)
    return [f*step, F*step] #[dh, h] update
'''

def rk4_step(N, hk0, Dhk0, step):
    F1 = Dhk0
    f1 = DDhk(N, hk0, Dhk0)
    F2 = Dhk0 +f1*step/2.
    f2 = DDhk(N +step/2., hk0 +F1*step/2., Dhk0 +f1*step/2.)
    F3 = Dhk0 +f2*step/2.
    f3 = DDhk(N +step/2., hk0 +F2*step/2., Dhk0 +f2*step/2.)
    F4 = Dhk0 +f3*step
    f4 = DDhk(N +step, hk0 +F3*step, Dhk0 +f3*step)

    return numpy.array([(f1 +2*f2 +2*f3 +f4)*step/6.], dtype=complex), numpy.array([(F1 +2*F2 +2*F3 +F4)*step/6.], dtype=complex) # [Dhk, hk] update

npts = 1000

Nics = 0
Nshss = 11.00548586827543321741572

step = (Nshss-Nics)/(npts-1)

hk0 = numpy.array([-1.24919723e+41 -1.10490234e+41j], dtype = complex) 
Dhk0 = numpy.array([-1.14621079e+41 -1.01381187e+41j], dtype = complex) 

hk = numpy.array([-1.24919723e+41 -1.10490234e+41j], dtype = complex)

N = Nics+step
while N < Nshss:
    #array = euler_step()
    array = rk4_step(N, hk0, Dhk0, step)
    hk0 = hk0 + array[1]
    Dhk0 = Dhk0 + array[0]
    hk = numpy.append(hk, hk0)
    N += step


# In[69]:

plt.plot(numpy.absolute(hk[:950]))


# In[71]:

plt.semilogy(numpy.absolute(hk[:925]))


# In[72]:

plt.plot(hk.imag[:925])


# In[ ]:




# In[ ]:



