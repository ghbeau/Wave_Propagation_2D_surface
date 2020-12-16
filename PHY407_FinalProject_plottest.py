"""
Project PHY407
@author: Genevieve Beauregard

This is just the script to generate the square plate
"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from PHY407_Functions import Gaussian, FDLaplaceEstimate, GetNextP
import scipy.linalg




# %% Parameters
Lx= 1.0 # length of the square
Ly = 1.0
nx = 200 # grid size
ny = 200

# nx = 3
# ny = 3
#N test
#N = 5
v = 30 # constant wave velocity



# %% creating arrays

#space arrays for functions
x = np.linspace(0, Lx,nx)
y = np.linspace(0, Ly, ny)
dx = x[1]-x[0]
X, Y = np.meshgrid(x, y)
time = 10 # time scale
dt = 0.1 # time grid
time = np.arange(0, time+dt, dt)


# P0 = np.zeros(np.shape(X)) # the P0 = 0 boundary conditions are already set
# P1 = np.zeros( np.shape(X))# the P1 = 0 boundary conditions set in here 
P1 = Gaussian(X, Y, 0.5, 0.5, 1, 0.1, 0.1) -0.1
P0 = Gaussian(X, Y, 0.5, 0.5, 1, 0.1, 0.1)

P_prev = np.copy(P0)
P = np.copy(P1) 



# D0 = np.zeros((3,3))
# D1 = np.copy(D0)
# D0[1][1] = 1
# D1[1][1] = 1 -0.1
# Stack = []
# Stack.append(D0)
# Stack.append(D1)
# for i in range(1, len(time)): 
#     print(i)
#     NextD = GetNextP(D1, D0, dt, dx, 1)
#     print(NextD)
#     Stack.append(NextD)
# # F = np.zeros((nx,ny))
# # F[sx][sy] = 60

# Pstack = np.zeros((len(time), nx, ny)) # create stack array 





#Pstack = P_prev
Pstack = np.array([np.copy(P_prev)])
Pstack = np.vstack((Pstack, np.array([np.copy(P)]))) #first index correspending to time

for i in range(1,len(time)):
    P_new = 2 * P - P_prev + (v ** 2) * (dt ** 2) * (FDLaplaceEstimate(P, dx))
    # P_next = GetNextP(Pstack[i],Pstack[i-1], dx, dt, v)
    Pstack = np.vstack((Pstack, np.array([P_new])))
    P_prev = np.copy(P)
    P = np.copy(P_new)
    
    
# plt.figure()
# plt.title("First")
# plt.imshow(Pstack[1])
# plt.colorbar()

# plt.figure()
# plt.title("Second")
# plt.imshow(Pstack[2])
# plt.colorbar()


# plt.figure()
# plt.title("3")
# plt.imshow(Pstack[3])
# plt.colorbar()

# plt.figure()
# plt.title("4")
# plt.imshow(Pstack[4])
# plt.colorbar()

# plt.figure()
# plt.title("5")
# plt.imshow(Pstack[5])
# plt.colorbar()

# plt.figure()
# plt.title("first")
# plt.imshow(Pstack[1])
# plt.colorbar()

# # %%
# # # plt.rc('text', usetex=True)             # use LaTeX for text
# # # plt.rc('font', family='serif')          # use serif font
# # # plt.rcParams.update({'font.size': 10})  # set font size 
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_wireframe(X, Y, Pstack[0])
# ax.set_title("First")
# #ax.set_xticks([0, 1])
# #ax.set_yticks([0, 1])
# ax.set_zticks([0])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# #ax.set_zlabel('z')
# #ax.view_init(90, 0)


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_surface(X, Y, Pstack[5])
# ax.set_title("Second")
# #ax.set_xticks([0, 1])
# #ax.set_yticks([0, 1])
# #ax.set_zticks([0])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_surface(X, Y, Pstack[30])
# ax.set_title("Second")
# #ax.set_xticks([0, 1])
# #ax.set_yticks([0, 1])
# ax.set_zticks([0])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')



# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_surface(X, Y, Pstack[50])
# ax.set_title("Second")
# #ax.set_xticks([0, 1])
# #ax.set_yticks([0, 1])
# ax.set_zticks([0])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
