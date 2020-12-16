"""
Project PHY407
@author: Genevieve Beauregard

Final Project

"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from PHY407_Functions import Gaussian, FDLaplaceEstimate, GetNextP
import scipy.ndimage
from matplotlib import cm



#%%

# Initialize Parameters
Lx= 300.0 # length of x
Ly = 300.0 # length of y
dx = 1  # resolution

v = 20 # constant wave velocity
v2 = v**2

# mesh 
Nx = int(Lx/dx) # no of points
Ny = int(Ly/dx) # no of points
x = np.linspace(0, Lx, Nx + 1)
y = np.linspace(0, Ly, Ny + 1)
X, Y = np.meshgrid(x, y)


end_time = 10 # time scale
dt = 0.005 # time grid
# create the time mesh,
time = np.arange(0, end_time+dt, dt)

# Courant Stability, we need this to be below 1
Courant  = (2* v * dt) / dx
print("Courant Stability is "+ str(round(Courant, 2))) 


# Source 


# initial conditions 


# set the initial Gaussian here, setting centre to be zero 
# amplitude 
a = 10
P0 = Gaussian(X, Y, Lx/2, Ly/2, a, Lx/100, Ly/100 )

# set initial V, wave velocity array 
V0 = np.zeros(np.shape(X))


# define updating variable for loop 
P = np.copy(P0)

# to kickstart the algorithm we need a P0^{-1} term, or a non existent 
# (-dt) term. We use our initial V condition using central difference approx
P_next = P + dt* V0 - 0.5 *(v*dt)**2 *  scipy.ndimage.laplace(P)

# create Stack (change later)
Pstack = []
Pstack.append(P)
Pstack.append(P_next)

P_prev = np.copy(P)
P = np.copy(P_next)

# now the current p is associated with n = 1

# Now loop from current time = 0 till and including end_time  
for i in range(1, len(time)):
    P_next = GetNextP(P, P_prev, dt, dx, v)
    #P_next = 2*P - P_prev + (v * dt)**2 * FDLaplaceEstimate(P, dx)
    
    # P_next =  (v * dt)**2 * scipy.ndimage.laplace(P)
    # P_next = 2*P + P_prev 
    
    P_prev = np.copy(P) # update
    P = np.copy(P_next) 
    
    Pstack.append(P) # collect
    
Pstack = np.array(Pstack)
# plot 
for i in range(0,len(time), 100):
    plt.clf() # clear the plot
    plt.imshow(Pstack[i]) 
    plt.colorbar(label = "Pressure (Pa)")
    plt.clim(-2, 2)
    plt.xlabel("x, m")
    plt.ylabel("y, m")
    plt.title("time = "+str(time[i])+"s")
    plt.draw() 
    plt.savefig("Plot_Neuman/FD_{}.pdf".format(round(time[i],2)))
    plt.pause(0.01) #pause to allow a smooth animation

# plt.imshow(P0)
# plt.colorbar()

# index = 0 
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# surf = ax.plot_surface(X, Y, Pstack[index],cmap=cm.coolwarm, linewidth=0,\
#                        antialiased=False )
# ax.set_title("$T ={}$s, for $\Delta x ={}$, $\Delta t = {}$s".format(round(time[index],2), dx, dt))
# ax.set_xticks([0, Lx])
# ax.set_yticks([0, Ly])
# ax.set_zticks([0])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('P')
# fig.colorbar(surf)
# fig.clim(-2.0, 2.0)
# plt.savefig("Plot_Neuman/FD_Initial.pdf")


# for i in range(1,10):
#     index = int(i * len(time)/10) 
#     fig = plt.figure()
#     ax = plt.axes(projection='3d')
#     surf = ax.plot_surface(X, Y, Pstack[index],cmap=cm.coolwarm,\
#                            linewidth=0, antialiased=False )
#     ax.set_title("$T ={}$s, for $\Delta x ={}$, $\Delta t = {}$s".\
#                  format(round(time[index],2), dx, dt))
#     ax.set_xticks([0, Lx])
#     ax.set_yticks([0, Ly])
#     ax.set_zticks([0])
#     ax.set_xlabel('x')
#     ax.set_ylabel('y')
#     ax.set_zlabel('P')
#     cbar = fig.colorbar(surf)
#     cbar.set_clim(-a, a)
#     plt.savefig("Plot_Neuman/FD_{}.pdf".format(round(time[index],2)))

