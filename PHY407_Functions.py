"""
Project PHY407
@author: Genevieve Beauregard

This is just the script to generate the square plate
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def getAnalyticEigenmode(m, n, a, b, X, Y):
    """
    Returns Analytical Eigenmode Solution of a resonant plate for Z

    Parameters
    ----------
    m : int
        mode pertaining to x
    n : int
        mode pertaining to y.
    a : float
        amplitude.
    b : float
        amplitude.
    X : array_like
        X grid N x N
    Y : array_like
        Y grid N x N

    Returns
    -------
    array_like

    """
    
    return a*np.sin(m * np.pi * X) * np.sin(n * np.pi * Y)

def getDynamicMatrixSquare(N):
    """
    Returns the dynamical matrix of a square membrane, without correction
    for edges. 

    Parameters
    ----------
    N : Int
        Length of Square grid.

    Returns
    -------
    D : array_like
        Dynamical Matrix.

    """
    
    Dlen = N**2
    D = np.zeros((Dlen, Dlen))
    
    for i in range(Dlen): # iterate of k  
        D[i][i] = 4 #populate the centre 
        

        if i+1 < Dlen: 
            D[i][i+1] = -1
        
        if i-1 <= Dlen and (i-1) >= 0: 
            D[i][i-1] = -1
            
        if i + N < Dlen: 
            D[i][i + N] = -1
        
        if i - N < Dlen and i-N >=0: 
            D[i][i -N] = -1
    return D

def getzvals(A, omega, t):
    
    return A* np.exp(-1j * omega * t)
    

    
    
