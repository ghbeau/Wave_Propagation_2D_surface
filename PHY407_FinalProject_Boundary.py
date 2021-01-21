"""
Project PHY407
@author: Genevieve Beauregard

Final Project
Damped Waves with Sponge Boundaries

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PHY407_Functions import Gaussian, FDLaplaceEstimate, GetNextP
import scipy.ndimage
from matplotlib import cm
import matplotlib.animation as animation



# just to avoid doing loops 
i_want_to_compute_all = True



#%% Parameters and Initial Values 

# Initialize Parameters in SI units
Lx= 300.0 # length of x [m]
Ly = 300.0 # length of y [m]
dx = 1  # resolution [m]

v = 340 # constant wave velocity speed of sound [ms-2]
v2 = v**2 

nb = 60 # no of layers for damping 


end_time = 1 # end time s 
dt = 0.0005 # time grid s
# create the time mesh,
time = np.arange(0, end_time+dt, dt)

# Courant Stability, we need this to be below 1
Courant  = (2* v * dt) / dx
print("Courant Stability is "+ str(round(Courant, 2))) 


#%% Introduce Mesh with Damping Layers
# mesh 
Nx = int(Lx/dx) # no of points in Solution Space
Ny = int(Ly/dx) # no of points in Solution Space


allNx = Nx + 2*nb+1 # Total Points inclusive of Absorbing Boundary
allNy =  Ny + 2*nb+1 # Total Points inclusive of Absorbing Boundary

#grid with boundary
x = np.linspace(-nb*dx , Lx + nb*dx, allNx)
y = np.linspace(-nb*dx , Ly + nb*dx, allNy)
X, Y = np.meshgrid(x, y)

#%% Initial Conditions

# We set the initial Gaussian here, setting centre to be zero 
# amplitude 
a = 10
P0 = Gaussian(X, Y, Lx/2, Ly/2, a, Lx/100, Ly/100 )

# set initial V, wave velocity array 
V0 = np.zeros(np.shape(X))

#%% Damping Coefficient for D_x and D_y 

# We want to create a Damping Numpy matrix, for which we take the exponent
# Define reflection coefficient reflected wave/incident wave amplitudes

Dx = np.zeros(np.shape(X))
Dy = np.zeros(np.shape(X))

nbL = nb*dx # damping layer width 


factor = 0.015 # K multiplicated factor 


# the damping grid for Dx occurs at 0 to nb-1 columns and from -nb to the end
Dx[:,:nb] = factor * ( (X[:,:nb] - X[0][nb]) / dx)**2
Dx[:,-nb:] = factor* ( (X[:,-nb:] - X[0][-nb-1]) / dx)**2

# the damping grid for Dy occurs similarly but for Rows, switching indices
Dy[:nb,:] = factor * ( (Y[:nb,:] - Y[nb][0]) /dx)**2
Dy[-nb:,:] = factor* ((  Y[-nb:,] - Y[-nb-1][0]) /dx)**2

D = Dx + Dy

# Create exponential damping term for all time steps 
expD = np.exp(-D)




# %% First timestep from n = 0 to n = 1 

# define updating variable for loop 
P = np.copy(P0)

# to kickstart the algorithm we need a P0^{-1} term, or a non existent 
# (-dt) term. We use our initial V condition using central difference approx
P_next = P + dt* V0 - 0.5 *(v*dt)**2 *  FDLaplaceEstimate(P,dx)

# Damp 
# P_next = P_next*expDstack[1]
P_next = P_next*expD


# Declare memory to store P values
Pstack = np.zeros((len(time), allNx, allNy))

# Store
Pstack[0] = np.copy(P) 
Pstack[1] = np.copy(P_next) 

# update
P_prev = np.copy(P)
P = np.copy(P_next)

# now the current p is associated with n = 1

#%% Loop from time step n = 1 till len(time)

if i_want_to_compute_all:
    for i in range(1, len(time)):
        # P_next = GetNextP(P, P_prev, dt, dx, v)*expDstack[i]
        P_next = GetNextP(P, P_prev, dt, dx, v)*expD

        P_prev = np.copy(P) # update
        P = np.copy(P_next)
        Pstack[i] = np.copy(P) # store
    np.savez('output_file', Pstack = Pstack)
        
else: 
    npzfile = np.load('output_file.npz')
    Pstack = npzfile['Pstack']
    
      
#%% Plotting and animations. 

Pstack_soln = np.copy(Pstack[:,nb:-nb,nb:-nb])



plt.rc('text', usetex=True)             # use LaTeX for text
plt.rc('font', family='serif')          # use serif font
plt.rcParams.update({'font.size': 15})  # set font size 


# The solution space arrays 
Pstack_soln = np.copy(Pstack[:,nb:-nb,nb:-nb])
X_soln = np.copy(X[-nb:nb, nb:-nb])
Y_soln = np.copy(Y[-nb:nb, nb:-nb])

# for i in range(0,len(time), 100):
#     plt.clf() # clear the plot
#     plt.imshow(Pstack_soln[i])
#     plt.colorbar(label = "P(Pa)")
#     plt.set_cmap('seismic')
#     plt.clim(-0.5,1.0)
#     plt.xlabel("x, m")
#     plt.ylabel("y, m")
#     plt.yticks([0, Ly])
#     plt.xticks([0,Lx])
#     plt.title("T = "+str(round(time[i],2))+"s")
#     plt.draw() 
#     plt.tight_layout()
#     plt.savefig("Plot_Neuman/nb={}_factor={}_DampledFD_{}.pdf".format(nb, factor, round(time[i],2)))
#     plt.pause(0.01) #pause to allow a smooth animation
    




#%% Animate! 


fig = plt.figure()
a = Pstack_soln[0]

im = plt.imshow(a)

included_frames = np.arange(0, len(time), 5)  #skip every 5


def animate_func(i):
    im.set_array(Pstack_soln[i])
    return [im]

anim = animation.FuncAnimation(fig, animate_func,\
                                frames=included_frames\
                                    , blit=True)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title("nb="+str(nb)+",factor = "+str(factor))
plt.colorbar()
plt.clim(-0.5, 1)
plt.set_cmap('seismic')         
anim.save('DampedWavefactor={}.mp4'.format(factor), writer='ffmpeg')



