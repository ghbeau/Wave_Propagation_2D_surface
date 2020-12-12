"""
Project PHY407
@author: Genevieve Beauregard

This is just the script to generate the square plate
"""

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from PHY407_Functions import *
import scipy.linalg



# %% Parameters
L = 1.0 # length of the square
N = 20 # grid size
#N test
#N = 5
a = 1 # amplitude initial conditions
b = 1 
m = 2 # modes
n = 2 # modes
time = 10 # time scale
dt = 0.5  # time grid
time = np.arange(0, time+dt, dt)
T = 1 
sigma = 1 


# %% Initial Value Conditions
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
dx = x[1]-x[0]
X, Y = np.meshgrid(x, y)

# Initial Condition (Make sure Boundary is set)
Z0 = getAnalyticEigenmode(m, n, a, b, X, Y)
D = getDynamicMatrixSquare(N)

# %% Computing Eigenmodes (This step may be skipped, this takes awhile)

I_want_to_compute_everything = True 

if I_want_to_compute_everything:

    # Find eigenmodes
    eigenval, eigenvec = scipy.linalg.eig(D)
    print(eigenval)
    np.savez('output_file', eigenval=eigenval, eigenvec=eigenvec, D=D)

else: 
    npzfile = np.load('output_file.npz') 
    eigenval = npzfile['eigenval']
    eigenvec = npzfile['eigenvec']    

# %% Now we  process them and reorder these modes: 
#omega_square = omega_square/((sigma * (dx)**2 )/T )
#omega_square = np.abs(omega_square)

    

omega = np.real(eigenval)/ np.sqrt(T/(sigma*(dx)**2))
# %%
# # plt.rc('text', usetex=True)             # use LaTeX for text
# # plt.rc('font', family='serif')          # use serif font
# # plt.rcParams.update({'font.size': 10})  # set font size 
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_wireframe(X, Y, )
# ax.set_title("Side View")
# #ax.set_xticks([0, 1])
# #ax.set_yticks([0, 1])
# ax.set_zticks([0])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# #ax.view_init(90, 0)
# plt.savefig('Neutral_Square.pdf')
y = np.arange(1, (N**2)+1)
plt.figure
plt.plot(omega)

