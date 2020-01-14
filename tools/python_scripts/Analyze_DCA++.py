#!/usr/local/bin/python3.7
import numpy as np
from numpy import *
import h5py
import sys
import os
import pandas as p

sys.path.append("/Users/t7m/Projects/DCAPP/DCAPP_maierta_fork/DCA/tools/python_scripts/")
import solveBSE_fromG4


# fileBase = 'dca_PARTICLE_PARTICLE_UP_DOWN.hdf5'
# fileBase = 'dca_PARTICLE_PARTICLE_UP_DOWN_[0, 0].hdf5'
# newMaster = False
ALLQ = False
# fileBase = 'dca_PARTICLE_HOLE_MAGNETIC_[0, 0].hdf5'

# temps = [0.2,0.15,0.125,0.1,0.08,0.06,0.04]
# temps = [0.4,0.3,0.25,0.225,0.2,0.175,0.15,0.125,0.1]
# temps = [0.4 ,  0.3, 0.25 , 0.225, 0.2, 0.175, 0.15, 0.125,0.1,0.08,0.06]
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
# Observables
Pd    = zeros((nTemps))
Pd01  = zeros((nTemps))
Pd02  = zeros((nTemps))
Pd03  = zeros((nTemps))
Pxs   = zeros((nTemps))
Ppp   = zeros((nTemps))
ld    = zeros((nTemps))
Vd1   = zeros((nTemps))
Vd2   = zeros((nTemps))
Vd3   = zeros((nTemps))

chiC = zeros((nTemps))
chiL = zeros((nTemps))



for T_ind, T in enumerate(temps):
    dir_str = "./T=" + str(T)
    print("Analyzing directory ",dir_str)

    b=solveBSE_fromG4.BSE(fileG4=dir_str+"/" + fileBase,symmetrize_G4=True,newMaster=newMaster,oldFormat=False,allq=ALLQ,evenFreqOnly=False,iq=0)

    U[T_ind]     = b.U
    tp[T_ind]    = b.tp
    cs[T_ind]    = int(b.Nc)
    dens[T_ind]  = round(b.dens[0],3)
    sign[T_ind]  = b.qmcSign[0]

    Pd[T_ind]    = b.Pd
    Pd01[T_ind]  = b.Pd01
    Pd02[T_ind]  = b.Pd02
    Pd03[T_ind]  = b.Pd03
    Vd1[T_ind]   = b.Vd01
    Vd2[T_ind]   = b.Vd02
    Vd3[T_ind]   = b.Vd03
    # Pxs[T_ind]   = b.Pxs
    # Ppp[T_ind]   = b.Ppx
    ld[T_ind]    = b.lambdad

    # chiC[T_ind] = b.chiC
    # chiL[T_ind] = b.chiL

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
                  'Pd'        : Pd,
                  'Pd01'      : Pd01,
                  'Pd02'      : Pd02,
                  'Pd03'      : Pd03,
                  'Vd1'       : Vd1,
                  'Vd2'       : Vd2,
                  'Vd3'       : Vd3,
                  'lambda_d'  : ld,
                  'U'         : U,
                  'tp'        : tp,
                  'Nc'        : cs,
                  'dens'      : dens
                                  })

df.to_csv("Analysis.txt",index=False)

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
