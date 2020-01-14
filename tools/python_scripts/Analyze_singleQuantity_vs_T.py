#!/usr/local/bin/python
import numpy as np
from numpy import *
import h5py
import sys
import os
import pandas as p

sys.path.append("/Projects/ETH_Reps/")
import solveBSE_fromG4

# fileBase = 'dca_PARTICLE_PARTICLE_UP_DOWN.hdf5'
fileBase = './dca_PARTICLE_HOLE_MAGNETIC_[0, 0].hdf5'
ALLQ = True; iq = 0
# fileBase = 'dca_PARTICLE_HOLE_MAGNETIC_[0, 0].hdf5'

# temps = [0.2,0.15,0.125,0.1,0.08,0.06,0.04]
# temps = [0.4,0.3,0.25,0.225,0.2,0.175,0.15,0.125,0.1]
temps = [0.4 ,  0.3, 0.25 , 0.225, 0.2, 0.175, 0.15, 0.125,0.1,0.08,0.06]
# temps = [0.4 ,  0.3, 0.25 , 0.225, 0.2, 0.175, 0.15, 0.125,0.1,0.08]
# temps = [0.4,0.3,0.25,0.225,0.2,0.175,0.15,0.125,0.08,0.06]
# temps = [0.04]
nTemps = temps.__len__()

# Parameters:
U     = zeros((nTemps))
tp    = zeros((nTemps))
cs    = zeros((nTemps))
dens  = zeros((nTemps))
sign  = zeros((nTemps))
# Observable
chiC  = zeros((nTemps))



for T_ind, T in enumerate(temps):
    dir_str = "./T=" + str(T)
    print("Analyzing directory ",dir_str)

    b=solveBSE_fromG4.BSE(fileG4=dir_str+"/" + fileBase,symmetrize_G4=True,newMaster=True,oldFormat=False,allq=ALLQ,evenFreqOnly=False,iq=iq)
    # b=solveBSE_fromG4.BSE(fileG4=dir_str+"/" + fileBase,symmetrize_G4=True,newMaster=False,oldFormat=False,allq=ALLQ,evenFreqOnly=False,iq=0)

    U[T_ind]     = b.U
    tp[T_ind]    = b.tp
    cs[T_ind]    = int(b.Nc)
    dens[T_ind]  = round(b.dens[0],3)
    sign[T_ind]  = b.qmcSign[0]

    chiC[T_ind]    = real(b.chiC)

# data = np.column_stack((temps,Pd,Pd0,Pd02,Pd03,ld,Pxs))
# np.savetxt("Pairing.txt", data, header = "T, Pd, Pd0, Pd02, Pd03, lambdad, Pxs")

# df = p.DataFrame({'T'         : temps, 
#                                 'chiC'      : chiC,
#                                 'chiL'      : chiL,
#                                 'U'         : U,
#                                 'tp'        : tp,
#                                 'Nc'        : cs,
#                                 'dens'          : dens
#                                 })

df = p.DataFrame({'T'         : temps, 
                  'chiC(q=0)'        : chiC,
                  'U'         : U,
                  'tp'        : tp,
                  'Nc'        : cs,
                  'dens'      : dens
                                  })

df.to_csv("Analysis_chiCQ0_vs_T.txt",index=False)

# Plot some data and fit

# import seaborn as sns
# sns.set_style('darkgrid',{'xtick.bottom': True, 'ytick.left': True})

# sns.scatterplot(x=df['T'],y=1-df['lambda_d'])

# from scipy.optimize import curve_fit
# def fitfunc(x,a,b):
# 	return a*log(x/b)

# popt,pcov = curve_fit(fitfunc,df["T"][-5:],1-df["lambda_d"][-5:])
# Tnew = linspace(0.00001,0.2,100)
# plot(Tnew,fitfunc(Tnew,*popt),'k--')
# ylim(-0.05,1.05)
# ylabel(r"$1-\lambda_d(T)$")
