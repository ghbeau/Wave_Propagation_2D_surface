"""
Project PHY407
@author: Genevieve Beauregard

This is just the script to generate the square plate
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



#def getNextP(P_prev, P_now, F_now, dx, dt, v, boundary=True,): 
    
    # # This is a non-vectorized scheme I wrote. Its simple but extremely slow
    # # and does not work very well
    # # Initiazing arrays
    # shape = np.shape(P_prev)
    # P_next = np.zeros(shape)
    
    # factor = ((dt/dx)**2) * (v**2) # the multiplicative fact

    # # Populate array and sum
    # for i in range(1, shape[0]-1) :
    #     for j in range(1,shape[1]-1):
            
    #         # gross sum is the equation used to calculate P_next
    #         gross_sum = (P_now[i+1][j] + P_now[i-1][j] + P_now[i][j-1] +\
    #             P_now[i][j])* factor
    
        
            
    #         gross_sum += (dt**2) * F_now[i][j] + 2 * P_now[i,j] - P_prev[i][j]
    
    #         P_next[i,j] = gross_sum
    
    # if boundary:
    #     P_next[0,:] = 0
    #     P_next[-1,:] = 0
    #     P_next[:,0] = 0
    #     P_next[:,-1] = 0
    # print(np.shape(P_next))
    
    # return P_next
    
    

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
    """This will calculate the finite difference version of laplace.
    This may be replaced by a built in function scipy.sparse.csgraph.laplacian¶

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
        2* P_now - P_prev
    # P_next[:,0] = 0
    # P_next[:,-1] = 0
    # P_next[0,:] = 0
    # P_next[-1,:] = 0
    return P_next

    
