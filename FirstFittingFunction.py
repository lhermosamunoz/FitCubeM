import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy.constants import c
from pyspeckit import speclines as pylines
import Ofuncts
import lmfit
import scipy.stats as stats
import os

from plot_broad import broad_plot
from DefRegionLines import DefRegionLines
from plot_refer_lines import refer_plot
from FitThirdComp import FitThirdComp


##################################################################################################################
# DEFINE THE FITTING FUNCTION FOR 1, 1+B, 2, 2+B COMPONENTS USING LMFIT AND THE EPSILON AND AIC CRITERIUM
#
def FirstFittingFunction(galaxy,galaxy2,l,data_cor,mu0,sig0,amp0,sig20,amp20,sig30,amp30,mu1,sig1,amp1,sig21,amp21,sig31,amp31,l2,l3,l5,l6,l9,l10):
    # 
    # Start the parameters for the LINEAR fit
    in_slope = 0.
    in_intc  = data_cor[0]

    # Redefine the lambda zone with the first and last zones
    newl = l[np.where(l<l[1])[0][-1]+5:np.where(l>6250.)[0][0]-10]
    zone_O_N = l[np.where(l<6450.)[0][-1]+10:np.where(l>l9)[0][0]-10]
    zone_N_S = l[np.where(l<l6)[0][-1]+30:np.where(l>l3)[0][0]-30]
    newl = np.append(newl,zone_O_N)
    newl = np.append(newl,zone_N_S)
    newl = np.append(newl,l[-100:-1])
    # now we do the same but with the flux data (y vector)
    newflux = data_cor[np.where(l<l[1])[0][-1]+5:np.where(l>6250.)[0][0]-10]
    zon_O_N = data_cor[np.where(l<6450.)[0][-1]+10:np.where(l>l9)[0][0]-10]
    zon_N_S = data_cor[np.where(l<l6)[0][-1]+30:np.where(l>l3)[0][0]-30]
    newflux = np.append(newflux,zon_O_N)
    newflux = np.append(newflux,zon_N_S)
    newflux = np.append(newflux,data_cor[-100:-1])

    ############## Standard deviation of the continuum ###############
    #
    # In order to determine if the lines need one more gaussian to be fit correctly, we apply the 
    # condition that the std dev of the continuum should be higher than 3 times the std dev of 
    # the residuals of the fit of the line. 
    # We have to calculate the stddev of the continuum in a place where there are no lines 
    # Calculate the std dev of a part of the continuum without lines nor contribution of them
    std0 = np.where(l>6450.)[0][0]
    std1 = np.where(l<6500.)[0][-1]
    stadev = np.std(data_cor[std0:std1])

    ############################# Start the fit and the MODEL ####################################
    #
    # First we have to initialise the model in the Halpha+NII lines by doing
    lin_mod = lmfit.Model(Ofuncts.linear)
    one_mod = lmfit.Model(Ofuncts.twogaussian)
    two_mod = lmfit.Model(Ofuncts.funcSII2comp)
    three_mod = lmfit.Model(Ofuncts.funcSII3comp)

    # and initialise the model in the whole spectra for several different models
    comp_mod = lmfit.Model(Ofuncts.funcgauslin)
    broad_mod = lmfit.Model(Ofuncts.funcbroad)
    twocomp_mod = lmfit.Model(Ofuncts.func2com)
    threecomp_mod = lmfit.Model(Ofuncts.func3com)
    twobroadcomp_mod = lmfit.Model(Ofuncts.func2bcom)
