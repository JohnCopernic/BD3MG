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

from threeMG.MajorizeMinimizeMemoryGradient3D import MajorizeMinimizeMemoryGradient3D
from PSF_tools.gaussian_kernel_3D import gaussian_kernel_3D
from BD3MG.blur_alt_z import blur_alt_z
from BD3MG.adjblur_alt_z import adjblur_alt_z



print("imports ok")
# Loading image
y = loadmat(os.getcwd() + "/Images/aneurysm1.mat")["I"]
Nx, Ny, Nz = y.shape
###ATTENTION : Il faut Nz=nb cpu
# Slicing the layers along the z-axis
y = y[:, :, 80:100]
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
Timemax = 600
NbIt = 10  # Max iterations number
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


MMMG = MajorizeMinimizeMemoryGradient3D(y, h, Hty, eta, lambda_, delta,tau,phi,I,x0,x0, xmin, xmax,NbIt, Timemax)

start=time.time()
MMMG.optimize()
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
Crits[cores] = MMMG.Crit
Times[cores] = np.cumsum(MMMG.Time)
SNR[cores] = MMMG.SNR

layer = 20
y = MMMG.y
x = MMMG.x

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


savemat(os.getcwd() + "/outputs/y_MMG.mat", {"y_MMG": x})

