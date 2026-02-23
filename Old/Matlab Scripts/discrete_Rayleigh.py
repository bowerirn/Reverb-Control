# This function computes the sound pressure from a panel of size (Lx,Ly) at
# a point in space (x,y,z), where z is the distance outward from the panel
# center.

# f is the frequency at which the calculation takes place (Hz)
# A is the panel displacement matrix 
# (for simplicity, use linspace so panel is 100x100 matrix)
import numpy as np

def discrete_Rayleigh(A, Lx, Ly, f, x, y, z):
    A = np.array(A) # sanity check

    # convert to radians
    w = 2 * np.pi * f

    # air density, frequency etc
    Constant = 1.224 * (w ** 2) / (2 * np.pi)

    # number of boxes (should be 100 x 100 if panel is made using linspace)
    SzX = A.shape[0]
    SzY = A.shape[1]

    # divide panel into small boxes with area LxLy/SzXSxY
    Dx = np.linspace(-Lx/2, Lx/2, SzX)
    Dy = np.linspace(-Ly/2, Ly/2, SzY)

    # Panel area
    dS = (Lx * Ly) / (SzX * SzY)

    # wavenumber
    k = w / 343

    out = 0

    # perform rayleigh integral (discrete)
    for i in range(SzX):
        
        for h in range(SzY):
            
            R = np.sqrt((Dx(i) - x)^2 + (Dy(h) - y)^2 + z^2)
            
            out = out + dS * Constant * A[i, h] * np.exp(-1j * k*  R) / R
    
    return out
            