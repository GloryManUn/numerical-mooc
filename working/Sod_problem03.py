# -*- coding: utf-8 -*-
"""
Created on Thu May 12 23:25:09 2016

@author: Сергей
"""

print("dasda")

import numpy 
from matplotlib import pyplot
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

#Basic initial condition parameters
#defining grid size, time steps, CFL condition, etc...
nx = 1281
T = 0.01
dx = 20./(nx-1)
sigma1 = 0.0008
dt = sigma1*dx
nt = T/dt
gamma = 1.4


x = numpy.linspace(-10,10,nx)


def u_initial(nx):
    
    Rho = numpy.ones(nx)
    Rho[int((nx-1)/2.):] = 0.125
    v = numpy.zeros(nx)
    p = 100000*numpy.ones(nx)
    p[int((nx-1)/2.):] = 10000
    return numpy.array([Rho, v, p])
    
    
z = u_initial(nx)        #initial conditions
print(z)


def f(u, gamma):
    
    """Rho = u[0]
    v = u[1]
    p = u[2]
    e = p/((gamma-1)*Rho)
    eT = e+0.5*v*v
    """
    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
    
    
    return numpy.array([u2, u2*u2/u1+(gamma-1)*(u3-0.5*u2*u2/u1),
                      (u3+(gamma-1)*(u3-0.5*u2*u2/u1))*u2/u1])
                      
                      
                      
def f1(u, gamma):
    Rho = u[0]
    v = u[1]
    p = u[2]
    e = p/((gamma-1)*Rho)
    eT = e+0.5*v*v
    return numpy.array([Rho, Rho*v, Rho*eT])  


u = f1(z, gamma) #flux variables 


u_star = numpy.zeros((len(u),len(u[0])))
u_n = numpy.zeros_like(u)      
    #copy the initial u array into each row of our new array
print(u.shape,u_star.shape)


u_plus = numpy.zeros_like(u)
u_minus = numpy.zeros_like(u)
flux = numpy.zeros_like(u)

"""     Godunov
for n in range(0,int(nt)):
    #print('next step')
    
    u_n = u.copy() 
    #print(n)
    
    
   
    
    
    
    u_plus[:,:-1] = u[:,1:] # Can't do i+1/2 indices, so cell boundary
    u_minus = u.copy() # arrays at index i are at location i+1/2
    flux = 0.5 * (f(u_minus, gamma) + 
                      f(u_plus, gamma) + 
                      dx / dt * (u_minus - u_plus))
    u_n[:,1:-1] = u[:,1:-1] + dt/dx*(flux[:,:-2]-flux[:,1:-1])
    u_n[:,0] = u[:,0]
    u_n[:,-1] = u[:,-1]
    u = u_n.copy()
 """



""" Rihtmayer

for n in range(0,int(nt)):
    #print('next step')
    #print(u)
    u_n = u.copy() 
    #print(n)
    #print(u_n)
    
    F = f(u, gamma)
    #print("u is ", type(u))
    #print(F)
    u_star[:,:] = 0.5*((u[:,1:]+u[:,:-1])-dt/dx*(F[:,1:]-F[:,:-1]))
    #print(u_star.shape)
        
    #print("u_star is ", type(u_star))
    #print(u_star)
    F_star = f(u_star, gamma)
    #print('F_star = ', F_star)
    u_n[:,1:-1] = u[:,1:-1]-dt/dx*(F_star[:,1:]-F_star[:,:-1])
    u_n[:,0] = u[:,0]
    u_n[:,-1] = u[:,-1]
    #print('u_n = ', u_n)
    u = u_n.copy()
    #print(u.shape)
    #print("u is ", type(u))
 """   
    
 
def minmod(e, dx):
    """
    Compute the minmod approximation to the slope
    
    Parameters
    ----------
    e : array of float 
        input data
    dx : float 
        spacestep
    
    Returns
    -------
    sigma : array of float 
            minmod slope
    """
    
    sigma = numpy.zeros_like(e)
    de_minus = numpy.ones_like(e)
    de_plus = numpy.ones_like(e)
    
    de_minus[:,1:] = (e[:,1:] - e[:,:-1])/dx
    de_plus[:,:-1] = (e[:,1:] - e[:,:-1])/dx
    
    # The following is inefficient but easy to read
    for i in range(1, len(e[0])-1):
        if any(de_minus[:,i]*de_plus[:,i] <0.0):
            sigma[:,i] = 0.0
        elif any(numpy.abs(de_minus[:,i]) < numpy.abs(de_plus[:,i])):
            sigma[:,i] = de_minus[:,i]
        else:
            sigma[:,i] = de_plus[:,i]
            
    return sigma



#MUSCL
for t in range(0,int(nt)):
        u_n = u.copy()   
        sigma = minmod(u,dx) #calculate minmod slope

        #reconstruct values at cell boundaries
        u_left = u + sigma*dx/2.
        u_right = u - sigma*dx/2.     
        
        flux_left = f(u_left, gamma) 
        flux_right = f(u_right, gamma)
        
        #flux i = i + 1/2
        flux[:,:-1] = 0.5 * (flux_right[:,1:] + flux_left[:,:-1] - dx/dt *\
                          (u_right[:,1:] - u_left[:,:-1] ))
        
        #rk2 step 1
        u_star[:,1:-1] = u[:,1:-1] + dt/dx * (flux[:,:-2] - flux[:,1:-1])
        
        u_star[:,0] = u[:,0]
        u_star[:,-1] = u[:,-1]
        
        
        sigma = minmod(u_star,dx) #calculate minmod slope
    
        #reconstruct values at cell boundaries
        u_left = u_star + sigma*dx/2.
        u_right = u_star - sigma*dx/2.
        
        flux_left = f(u_left, gamma) 
        flux_right = f(u_right, gamma)
        
        flux[:,:-1] = 0.5 * (flux_right[:,1:] + flux_left[:,:-1] - dx/dt *\
                          (u_right[:,1:] - u_left[:,:-1] ))
        
        u_n[:,1:-1] = .5 * (u[:,1:-1] + u_star[:,1:-1] + dt/dx * (flux[:,:-2] - flux[:,1:-1]))
        
        u_n[:,0] = u[:,0]
        u_n[:,-1] = u[:,-1]
        u = u_n.copy()











   
def f3(u, gamma):
    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
   
    return numpy.array([u1, u2/u1, (u3/u1-0.5*u2*u2/u1/u1)*((gamma-1)*u1)])
    
z2 = f3(u, gamma)

#count velocity, pressure and density
print(z2[1,50])
print(z2[2,50])
print(z2[0,50])


p = numpy.arange(len(u[0]))
p = z2[2,:]
print(p.shape,x.shape)
rho = numpy.arange(len(u[0]))
rho = z2[0,:]
v = numpy.arange(len(u[0]))
v = z2[1,:]
#pyplot.plot(x, rho, color='#003366', ls='-', lw=3)
#pyplot.plot(x, v, color='blue', ls='-', lw=3)
pyplot.plot(x, p, color='green', ls='-', lw=3)
#pyplot.ylabel('pressure')
pyplot.xlabel('Distance')