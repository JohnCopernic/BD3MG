#from tkinter import Y
import warnings
import multiprocessing as mp
import time
import os, sys
import numpy as np
from scipy.io import loadmat
from scipy.io import savemat
warnings.filterwarnings("ignore")
sys.path.append(os.getcwd())


from BP3MG.PAR3MG_MPI import PAR3MG_MPI
from PSF_tools.gaussian_kernel_3D import gaussian_kernel_3D
from BD3MG.blur_alt_z import blur_alt_z
from BD3MG.adjblur_alt_z import adjblur_alt_z



print("imports ok")
# Loading image
y = loadmat(os.getcwd() + "/Images/ImageSmall.mat")["y"]
Nx, Ny, Nz = y.shape
###ATTENTION : Il faut Nz=nb cpu
# Slicing the layers along the z-axis
y = y[:, :, 25:45]
Nz = 20

print("image loaded")

# loading h_full
h_full = loadmat(os.getcwd() + "/Images/PSF.mat")["h_full"]
print("h_full loaded")

Nh = np.array(h_full.shape).astype(int)

h = lambda z: h_full[:, :, :, z]

print("size kernel: Nx = {}, Ny = {}, Nz = {}".format(*h(0).shape))

print("calculating the inputs...")
# variables to initialize : hty, h1
from scipy.signal import convolve

H1Z = [convolve(np.ones((Nx, Ny, Nz)), h(z), "same") for z in range(Nz)]

from PSF_tools.applyPSFadjvar3Dz import applyPSFadjvar3Dz_2

Hty = [applyPSFadjvar3Dz_2(y, z, h) for z in range(Nz)] 

Hty = np.dstack(Hty)
H1 = np.dstack([H1Z[z][:, :, z] for z in range(len(H1Z))])

# I is the estimation. Because we don't have any estimation, let's define :
from copy import deepcopy
I = deepcopy(y)


print("Initialization done")
print("Starting restoration...")


# Initialisation
Timemax = 24*3600
NbIt = 5000  # Max iterations number
# Regularization parameters:
lambda_ = 1
delta = 2
phi = 4
# Bounds of the constrained domain:
xmin = 0
xmax = 1
# Elastic net parameter:
tau = 1e-3
# Weight of the quadratic distance function :Gradz
eta = 0.1
# Initialization
x0 = np.zeros((Nx, Ny, Nz))

cores=20
print(cores)
# Optimization loop
Times = {}
Crits = {}
ratio = {}
SNR = {}

BMMD = PAR3MG_MPI(y, h, Hty, H1, eta, tau, lambda_, delta, xmin, xmax, phi, x0, I, Nx, Ny, Nz, NbIt, cores, Timemax)

start=time.time()
BMMD.optimize()
timer=time.time()-start
print("timer is :", timer)
"""
    :param y: observed data
    :param h: Gaussian blur operator
    :param Hty: adjoint of the gaussian blur operator applied to y
    :param eta: regularization parameter on the distance to the hypercube [xmin, xmax]
    :param lambda_: regularization parameter on the horizontal total variation norm
    :param delta: second regularization parameter on the horizontal total variation norm
    :param kappa: regularization parameter on the vertical total variation norm
    :param phi: choice of regularization function
    :param I: estimation
    :param xstar: ground truth
    :param xbar: minimal estimation
    :param xmin: lower bound on image pixel value
    :param xmax: upper bound on image pixel value
    :param NbIt: nb of iteration for the algorithm to reach
    :param timemax: maximal time of computation

"""
# Collecting and saving statistics
Crits[cores] = BMMD.Crit
Times[cores] = np.cumsum(BMMD.Time)
SNR[cores] = BMMD.SNR

layer = 20
y = BMMD.y
x = BMMD.x

f = open("Times.txt", "w")
f.write(str(Times))
f.close()

f = open("Crits.txt", "w")
f.write(str(Crits))
f.close()

f = open("ratio.txt", "w")
f.write(str(ratio))
f.close()

f = open("SNR.txt", "w")
f.write(str(SNR))
f.close()


savemat(os.getcwd() + "/outputs/y_restored_long_BP3MG.mat", {"y_BP_long": x})

