import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.constants import c
from pyspeckit import speclines as pylines
import Ofuncts
import lmfit
import os

from DefRegionLines import DefRegionLines
from plot_refer_lines_Mix import refer_plot
from FitThirdComp import FitThirdComp
import configMix


##############################################################################################################
# Define the parameters to enter in lmfit
def DefComp(lineSO,mu0,sig0,amp0,mu1,sig1,amp1,sig20,amp20,sig21,amp21,sig30,amp30,sig31,amp31):
    if lineSO == 'S':
        cd1 = lmfit.Parameter('mu_0', value=mu0)
        de = lmfit.Parameter('sig_0', value=sig0,min=0.)#,max=6.)
        ef = lmfit.Parameter('amp_0', value=amp0,min=0.)
        fg = lmfit.Parameter('mu_1', value=mu1,expr='mu_0*(6716.44/6730.82)')
        gh = lmfit.Parameter('sig_1', value=sig1,expr='sig_0')
        hi = lmfit.Parameter('amp_1', value=amp1,min=0.)
        # second components
        aaa = lmfit.Parameter('mu_20', value=mu0)
        aab = lmfit.Parameter('sig_20', value=sig20,min=0)#6.,max=12.)
        aac = lmfit.Parameter('amp_20', value=amp20,min=0.)
        aad = lmfit.Parameter('mu_21', value=mu1,expr='mu_20*(6716.44/6730.82)')
        aae = lmfit.Parameter('sig_21', value=sig21,expr='sig_20')
        aaf = lmfit.Parameter('amp_21', value=amp21,min=0.)
        # third components 
        aba = lmfit.Parameter('mu_30',value=mu0)
        abb = lmfit.Parameter('sig_30', value=sig30,min=0.)
        abc = lmfit.Parameter('amp_30',value=amp30,min=0.)
        abd = lmfit.Parameter('mu_31',value=mu1,expr='mu_30*(6716.44/6730.82)')
        abe = lmfit.Parameter('sig_31', value=sig31,expr='sig_30')
        abf = lmfit.Parameter('amp_31',value=amp31,min=0.)
    elif lineSO == 'O':
        cd1 = lmfit.Parameter('mu_5', value=mu0)
        de = lmfit.Parameter('sig_5', value=sig0,min=0.)#,max=6.)
        ef = lmfit.Parameter('amp_5', value=amp0,min=0.)
        fg = lmfit.Parameter('mu_6', value=mu1,expr='mu_5*(6363.77/6300.30)')
        gh = lmfit.Parameter('sig_6', value=sig1,expr='sig_5')
        hi = lmfit.Parameter('amp_6', value=amp1,min=0.)
        # second components
        aaa = lmfit.Parameter('mu_25', value=mu0)
        aab = lmfit.Parameter('sig_25', value=sig20,min=0)#6.,max=12.)
        aac = lmfit.Parameter('amp_25', value=amp20,min=0.)
        aad = lmfit.Parameter('mu_26', value=mu1,expr='mu_25*(6363.77/6300.30)')
        aae = lmfit.Parameter('sig_26', value=sig21,expr='sig_25')
        aaf = lmfit.Parameter('amp_26', value=amp21,min=0.)
        # third components 
        aba = lmfit.Parameter('mu_35',value=mu0)
        abb = lmfit.Parameter('sig_35', value=sig30,min=0.)
        abc = lmfit.Parameter('amp_35',value=amp30,min=0.)
        abd = lmfit.Parameter('mu_36',value=mu1,expr='mu_35*(6363.77/6300.30)')
        abe = lmfit.Parameter('sig_36', value=sig31,expr='sig_35')
        abf = lmfit.Parameter('amp_36',value=amp31,min=0.)

    return cd1,de,ef,fg,gh,hi,aaa,aab,aac,aad,aae,aaf,aba,abb,abc,abd,abe,abf

##################################################################################################################
# DEFINE THE FITTING FUNCTION FOR 1, 1+B, 2, 2+B COMPONENTS USING LMFIT AND THE EPSILON AND AIC CRITERIUM
#
def FirstFittingFunction(galaxy,galaxy2,l,data_cor,mu0,sig0,amp0,sig20,amp20,sig30,amp30,mu1,sig1,amp1,sig21,amp21,sig31,amp31,mu5,sig5,amp5,sig25,amp25,sig35,amp35,mu6,sig6,amp6,sig26,amp26,sig36,amp36,l2,l3,l5,l6,l9,l10,l11,l12,l13,l14):
    path = configMix.path
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
    #one_mod = lmfit.Model(Ofuncts.twogaussian)
    #two_mod = lmfit.Model(Ofuncts.funcSII2comp)
    one_mod = lmfit.Model(Ofuncts.fourgaussian)
    two_mod = lmfit.Model(Ofuncts.eightgaussian)
    three_mod = lmfit.Model(Ofuncts.twelvegaussian)

    # and initialise the model in the whole spectra for several different models
    comp_mod = lmfit.Model(Ofuncts.funcgauslin)
    broad_mod = lmfit.Model(Ofuncts.funcbroad)
    twocomp_mod = lmfit.Model(Ofuncts.func2com)
    threecomp_mod = lmfit.Model(Ofuncts.func3com)
    twobroadcomp_mod = lmfit.Model(Ofuncts.func2bcom)
    threebroadcomp_mod = lmfit.Model(Ofuncts.func3bcom)

    # We make the linear fit only with some windows of the spectra and calculate the 
    # straight line to introduce it in the formula
    linresu  = lin_mod.fit(newflux,slope=in_slope,intc=in_intc,x=newl)
    new_slop = linresu.values['slope']
    new_intc = linresu.values['intc']
    lin_data_fin = (linresu.values['slope']*l+linresu.values['intc'])

    # Now we define the initial guesses and the constraints
    params1 = lmfit.Parameters()
    params2 = lmfit.Parameters()
    params3 = lmfit.Parameters()

    sl = lmfit.Parameter('slop', value=new_slop,vary=False)
    it = lmfit.Parameter('intc', value=new_intc,vary=False)

    meth = 'M1'#input('Which method to be applied? ("S"/"O", not "M1"/"M2"): ')  # Method to fit

    # Define the parent Folder to save the data depending on the method
    if not os.path.exists(path+str(galaxy)+'_'+str(galaxy2)+'_results/'+meth+'_method'):
        os.mkdir(path+str(galaxy)+'_'+str(galaxy2)+'_results/'+meth+'_method')

    parentFold = path+str(galaxy)+'_'+str(galaxy2)+'_results/'+meth+'_method/'

    # -------------
    # Define the parameters using the function DefComp.
    # If neccesary, other parameters can be applied
    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    cds,des,efs,fgs,ghs,his,aaas,aabs,aacs,aads,aaes,aafs,abas,abbs,abcs,abds,abes,abfs = DefComp('S',mu0,sig0,amp0,mu1,sig1,amp1,sig20,amp20,sig21,amp21,sig30,amp30,sig31,amp31)
    cdo,deo,efo,fgo,gho,hio,aaao,aabo,aaco,aado,aaeo,aafo,abao,abbo,abco,abdo,abeo,abfo = DefComp('O',mu5,sig5,amp5,mu6,sig6,amp6,sig25,amp25,sig26,amp26,sig35,amp35,sig36,amp36)
    
    params1.add_many(sl,it,cds,des,efs,fgs,ghs,his,cdo,deo,efo,fgo,gho,hio)
    params2.add_many(sl,it,cds,des,efs,fgs,ghs,his,aaas,aabs,aacs,aads,aaes,aafs,cdo,deo,efo,fgo,gho,hio,aaao,aabo,aaco,aado,aaeo,aafo)
    params3.add_many(sl,it,cds,des,efs,fgs,ghs,his,aaas,aabs,aacs,aads,aaes,aafs,abas,abbs,abcs,abds,abes,abfs,cdo,deo,efo,fgo,gho,hio,aaao,aabo,aaco,aado,aaeo,aafo,abao,abbo,abco,abdo,abeo,abfo)

    ###############################################################################################
    # Make the fit using lmfit
    # simultanously to SII and OI lines (but independently)
    oneresu = one_mod.fit(data_cor,params1,x=l)
    tworesu = two_mod.fit(data_cor,params2,x=l)
    threeresu = three_mod.fit(data_cor,params3,x=l)

    with open(parentFold+'fit_oneC_result.txt', 'w') as fh:
        fh.write(oneresu.fit_report())
    with open(parentFold+'fit_twoC_result.txt', 'w') as fh:
        fh.write(tworesu.fit_report())
    with open(parentFold+'fit_threeC_result.txt', 'w') as fh:
        fh.write(threeresu.fit_report())

    ###############################################################################################
    # Select if one or two components in the SII lines and then apply to the rest
    #ep_1,ep_2,ep2_1,ep2_2,ep3_1,ep3_2 = refer_plot(parentFold,l,data_cor,meth,linresu,oneresu,tworesu,threeresu,l5,l6,l9,l10,l11,l12,l13,l14,std0,std1)
    # APPLY THIS LINE and obtain epsilons PER LINE ([SII] AND [OI])
    ep_1s,ep_1o,ep2_1s,ep2_1o,ep3_1s,ep3_1o = refer_plot(parentFold,l,data_cor,meth,linresu,oneresu,tworesu,threeresu,l5,l6,l9,l10,l11,l12,l13,l14,std0,std1)
    #if oneresu.chisqr < tworesu.chisqr:
    #    trigger = 'Y'
    #else:
    if ep_1s > 3 and ep_1o > 3:
        trigger = 'N'
    else:
        trigger = 'Y'
    #trigger='Y'   #input('Is the fit good enough with one component? ("Y"/"N"): ')

    return std0,std1,stadev,linresu,lin_data_fin,sl,it,meth,oneresu,tworesu,threeresu,trigger,comp_mod,broad_mod,twocomp_mod,twobroadcomp_mod,threecomp_mod,threebroadcomp_mod,ep_1s,ep2_1s,ep3_1s,ep_1o,ep2_1o,ep3_1o


