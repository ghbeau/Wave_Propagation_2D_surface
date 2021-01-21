"""
Project PHY407
@author: Genevieve Beauregard

Final Project
Neumann
"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from PHY407_Functions import Gaussian, FDLaplaceEstimate, GetNextP
import scipy.ndimage
from matplotlib import cm
import matplotlib.animation as animation





#%% # Parameters and Initial Values 

# Initialize Parameters in SI units
Lx= 300.0 # length of x [m]
Ly = 300.0 # length of y [m]
dx = 1  # resolution [m]

v = 340 # constant wave velocity speed of sound [ms-2]
v2 = v**2 

# mesh 
Nx = int(Lx/dx) # no of points
Ny = int(Ly/dx) # no of points
x = np.linspace(0, Lx, Nx + 1)
y = np.linspace(0, Ly, Ny + 1)
X, Y = np.meshgrid(x, y)


end_time = 1 # end time s 
dt = 0.0005 # time grid s
# create the time mesh,
time = np.arange(0, end_time+dt, dt)

# Courant Stability, we need this to be below 1
Courant  = (v * dt) / dx
print("Courant Stability is "+ str(round(Courant, 2))) 


# set the initial Gaussian here, setting centre to be zero 
# amplitude 
a = 10
P0 = Gaussian(X, Y, Lx/2, Ly/2, a, Lx/100, Ly/100 )


# set initial V, wave velocity array 
V0 = np.zeros(np.shape(X))



# %% First timestep from n = 0 to n = 1 

# define updating variable for loop 
P = np.copy(P0)

# to kickstart the algorithm we need a P0^{-1} term, or a non existent 
# (-dt) term. We use our initial V condition using central difference approx
P_next = P + dt* V0 - 0.5 *(v*dt)**2 *  FDLaplaceEstimate(P,dx)


# Declare memory to store P values
Pstack = np.zeros((len(time), Nx+1 , Ny+ 1))
Pstack[0] = np.copy(P)
Pstack[1] = np.copy(P_next)

# update
P_prev = np.copy(P)
P = np.copy(P_next)

# now the current p is associated with n = 1

#%%  Now loop from current time = 0 till len(t)

for i in range(1, len(time)):
    P_next = GetNextP(P, P_prev, dt, dx, v)
    P_prev = np.copy(P) # update
    P = np.copy(P_next)
    Pstack[i] = np.copy(P) # store
    
    
    
    
#%% Plotting and animations. 



plt.rc('text', usetex=True)             # use LaTeX for text
plt.rc('font', family='serif')          # use serif font
plt.rcParams.update({'font.size': 30})  # set font size 



for i in range(0,len(time), 100):
    plt.clf() # clear the plot
    plt.imshow(Pstack[i]) 
    plt.colorbar(label = "Pressure (Pa)")
    plt.set_cmap('seismic')
    plt.xlabel("x, m")
    plt.ylabel("y, m")
    plt.yticks([0, Ly])
    plt.xticks([0,Lx])
    plt.title("T = "+str(round(time[i],2))+"s")
    plt.draw() 
    plt.tight_layout()
    plt.savefig("Plot_Neuman/FD_{}.pdf".format(round(time[i],2)))
    plt.pause(0.01) #pause to allow a smooth animation
    

#%% Animate!

# fig = plt.figure()
# a = Pstack[0]

# im = plt.imshow(a)

# included_frames = np.arange(0, len(time), 100)  #skip every 5


# def animate_func(i):
#     im.set_array(Pstack[i])
#     return [im]

# anim = animation.FuncAnimation(fig, animate_func,\
#                                 frames=included_frames\
#                                     , blit=True)
# plt.xlabel('x [m]')
# plt.ylabel('y [m]')
# plt.title("Wave with Neumann Boundary over 1 seconds")
# plt.colorbar()
# plt.set_cmap('seismic')         
# anim.save('Wave.mp4', writer='ffmpeg')

    

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



