"""
Project PHY407
@author: Genevieve Beauregard

These are just some of the functions.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy.fft as fft 

def Gaussian(XX, YY, X0, Y0, A, sigmax, sigmay):
    """
    This is a gaussian bell function we use for an initial parameter.
    
    Input:
        XX: x coordinate mesh, array  
        YY: y coordinate mesh, array
        X0: centre of gaussian, scalar
        Y0: centre of Gaussian, scalar
        A: Amplitude, scalar
        sigmax: parameter wrt to x, scalar
        sigmay: parameter wrt to y, scalar
    Output:
        Gaussian curve: Array
    
    """
    
    return A*np.exp((-0.5)*(((XX - X0)/sigmax)**2 +((YY -Y0)*(1/sigmay))**2))

def FDLaplaceEstimate(P, dx):
    """This will calculate the finite difference version of laplace for 
    neumann boundary conditions. 
    Input:
        P_now: Current pressure array, array NxN
        dx: grid step (we assume same dx = dy), scalar
        
    Output: 
        Second Derivative Matrix of Laplace, (N-1) x (N-1), with zeros
        at boundary
    """
    # calculate second derivative, Had to reference 
    # https://hplgit.github.io/fdm-book/doc/pub/book/sphinx/._book008.html
    # and Igel's coursera course on wave equations for Geophysics
    # as I kept on messing up the array handling
    # P_xx is a vector of lower dimensionality as the edges P_0 are zero
    shape = np.shape(P)
    
    #define second derivative matrix 
    d2PdX = np.zeros(shape)
    d2PdY = np.zeros(shape)
    # for loop over the rows of P_now to calculate d2Pdy
    for i in range(1,shape[0]-1):
        d2PdY[i,:] = (P[i-1,:] - 2* P[i,:] + P[i+1,:])*(1/dx)**2
        
    # again with this but with the rows
    for j in range(1, shape[1] - 1):
        d2PdX[:, j] = (P[:, j-1] - 2*P[:,j] + P[:,j+1]) * (1/ dx)**2
    
    
    return d2PdX + d2PdY

def GetNextP(P_now, P_prev, dt, dx, v,  F_now=0):
    """ This will go in the forward timestep
    Input:
        P_now: Current pressure, array
        P_old: prev Pressure, array
        dt: time step, scalar
        dx: grid size, scalar
        Fnow: Current Source value, array or scalar, set to zero
        v: wave speed, same units as dx, scalar
    Output: 
        Next P: Next Pressure Array, array
    """
    
    # Get laplacian estimate
    P_nowLaplace = FDLaplaceEstimate(P_now, dx)
    P_next = (v**2) * (dt**2 )* P_nowLaplace + (dt**2) * F_now +\
        2* P_now - P_prev # forward time
    return P_next


    
