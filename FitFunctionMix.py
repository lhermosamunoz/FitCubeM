'''
This script makes a gaussian fit to the emission lines of AGN spectra
It is needed a parentFold, the spectrum in which the fit is going to be made and the initial estimation of the fit
'''

from configMix import * 
import configMix

from plot_broad_Mix import broad_plot
from DefRegionLines import DefRegionLines
from plot_refer_lines_Mix import refer_plot
from FitThirdComp import FitThirdComp
from FirstFittingFunctionMix import FirstFittingFunction
 
def FitPoint(listxy):
    path = configMix.path
    meth = configMix.meth
    x,y,nana=listxy
    galaxy2 = y
    galaxy  = x
    fitoutput = FitOutput(x,y)
    print('')
    print('Starting the fit for spec_'+str(galaxy)+'_'+str(galaxy2))
    #print('Done {:.2f} %'.format(100*((galaxy+1)+(galaxy2*33))/990.))
    print('')
    # Define the path to save the data
    l = np.copy(l_init)#*10000)-red_lambda_cor_SII2
    data = datos[:,galaxy2,galaxy]
    beg_point = np.where(l>6200)[0][0]
    end_point = np.where(l<6850)[0][-1]
    l = l[beg_point:end_point]
    data_cor = data[beg_point:end_point]
    ########################## Transform data and plot the spectra ##########################
    # Plot the spectra and check that everything is ok
    plt.close('all')
    plt.subplot(211)
    plt.plot(l,data_cor,'k',linewidth=1.)
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel(r'Flux (erg $\rm s^{-1} cm^{-2} \AA^{-1}$)')
    plt.xlim(l[0],l[-1])
    plt.subplot(212)
    plt.plot(l_init[400:-400],datos[400:-400,galaxy2,galaxy],'k',linewidth=1.)
    plt.xlim(l_init[400],l_init[-400])
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux ($erg/s/cm^{2} / \AA$)')
    plt.savefig(path+str(galaxy)+'_'+str(galaxy2)+'_results/spec.png')

    # some points may be nan, we have to create a mask to avoid them
    # when applying the mask, all nan values will turn to 0
    # if all points are 0, then break the loop and take the next spec
    mask = np.isnan(data_cor)
    data_cor[mask] = 0
    if np.any(data_cor) == 0:
        print('For the spectra in '+str(galaxy)+'_'+str(galaxy2)+' all the datapoints where nan')
        print('Continue...')
        fitoutput.prov_VN = (np.nan);fitoutput.prov_VN2 = (np.nan);fitoutput.prov_VN_2 = (np.nan);fitoutput.prov_VS = (np.nan);fitoutput.prov_VS_2 = (np.nan);fitoutput.prov_VB = (np.nan);fitoutput.prov_eVN = (np.nan);fitoutput.prov_eVN2 = (np.nan);fitoutput.prov_eVN_2 = (np.nan);fitoutput.prov_eVS = (np.nan);fitoutput.prov_eVS_2 = (np.nan);fitoutput.prov_eVB = (np.nan);fitoutput.prov_SN = (np.nan);fitoutput.prov_SN2 = (np.nan);fitoutput.prov_SN_2 = (np.nan);fitoutput.prov_SS = (np.nan);fitoutput.prov_SS_2 = (np.nan);fitoutput.prov_SB = (np.nan);fitoutput.prov_eSN = (np.nan);fitoutput.prov_eSN2 = (np.nan);fitoutput.prov_eSN_2 = (np.nan);fitoutput.prov_eSS = (np.nan);fitoutput.prov_eSS_2 = (np.nan);fitoutput.prov_eSB = (np.nan);fitoutput.prov_notSN = (np.nan);fitoutput.prov_noteSN = (np.nan);fitoutput.prov_notVN = (np.nan);fitoutput.prov_noteVN = (np.nan);fitoutput.prov_fHa = (np.nan);fitoutput.prov_fN1 = (np.nan);fitoutput.prov_fN2 = (np.nan);fitoutput.prov_fS1 = (np.nan);fitoutput.prov_fS2 = (np.nan);fitoutput.prov_fO1 = (np.nan);fitoutput.prov_fO2 = (np.nan);fitoutput.prov_f2Ha = (np.nan);fitoutput.prov_f2N1 = (np.nan);fitoutput.prov_f2N2 = (np.nan);fitoutput.prov_f2S1 = (np.nan);fitoutput.prov_f2S2 = (np.nan);fitoutput.prov_f2O1 = (np.nan);fitoutput.prov_f2O2 = (np.nan);fitoutput.prov_fSHa = (np.nan);fitoutput.prov_fSN1 = (np.nan);fitoutput.prov_fSN2 = (np.nan);fitoutput.prov_fSS1 = (np.nan);fitoutput.prov_fSS2 = (np.nan);fitoutput.prov_fSO1 = (np.nan);fitoutput.prov_fSO2 = (np.nan);fitoutput.prov_e1 = (np.nan);fitoutput.prov_e2 = (np.nan);fitoutput.prov_e3 = (np.nan);fitoutput.prov_a1 = (np.nan);fitoutput.prov_a2 = (np.nan);fitoutput.prov_a3 = (np.nan);fitoutput.prov_b1 = (np.nan);fitoutput.prov_b2 = (np.nan);fitoutput.prov_b3 = (np.nan);fitoutput.prov_e1_2 = (np.nan);fitoutput.prov_e2_2 = (np.nan);fitoutput.prov_e3_2 = (np.nan)
        return  

    parentFold = GetParentFold(path,galaxy,galaxy2)

    ############################ MAIN ################################
    ##################################################################
    l1 = 6765. #input('lambda inf for SII 2 (angs)?: ') 6725.
    l2 = 6790. #input('lambda sup for SII 2 (angs)?: ') 6740.
    l3 = 6730. #input('lambda inf for SII 1 (angs)?: ') 6710.
    l4 = 6765. #input('lambda sup for SII 1 (angs)?: ') 6725.
    l5 = 6610. #input('lambda inf for NII 2 (angs)?: ') 6575.
    l6 = 6635. #input('lambda sup for NII 2 (angs)?: ') 6595.
    l7 = 6590. #input('lambda inf for Halpha (angs)?: ') 6555.
    l8 = 6610. #input('lambda sup for Halpha (angs)?: ') 6575.
    l9 = 6560. #input('lambda inf for NII 1 (angs)?: ') 6535.
    l10 = 6590. #input('lambda sup for NII 1 (angs)?: ') 6555.
    l11 = 6300.	# lambda inf for OI 1 (angs)
    l12 = 6350.	# lambda sup for OI 1 (angs)
    l13 = 6380.	# lambda inf for OI 2 (angs)
    l14 = 6410.	# lambda sup for OI 2 (angs)
    ##################################################################
    # Retrieve the inital conditions and values for the lines
    # the function selects the regions in which the lines are expected
    # given the redshift of the source
    sig0,sig20,sig30,mu0,amp0,amp20,amp30,sig1,sig21,sig31,mu1,amp1,amp21,amp31,sig2,sig22,sig32,mu2,amp2,amp22,amp32,sig3,sig23,sig33,mu3,amp3,amp23,amp33,sig4,sig24,sig34,mu4,amp4,amp24,amp34,sig5,sig25,sig35,mu5,amp5,amp25,amp35,sig6,sig26,sig36,mu6,amp6,amp26,amp36 = DefRegionLines(l,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,sig_inicial,data_cor)
    # Fit the SII and OI lines and initialise all the models to be fitted during the process
    std0,std1,stadev,linresu,lin_data_fin,sl,it,meth,oneresu,tworesu,threeresu,trigger,comp_mod,broad_mod,twocomp_mod,twobroadcomp_mod,threecomp_mod,threebroadcomp_mod,ep_1s,ep2_1s,ep3_1s,ep_1o,ep2_1o,ep3_1o = FirstFittingFunction(galaxy,galaxy2,l,data_cor,mu0,sig0,amp0,sig20,amp20,sig30,amp30,mu1,sig1,amp1,sig21,amp21,sig31,amp31,mu5,sig5,amp5,sig25,amp25,sig35,amp35,mu6,sig6,amp6,sig26,amp26,sig36,amp36,l2,l3,l5,l6,l9,l10,l11,l12,l13,l14)
    # Save the epsilon and aic information from the models 
    fitoutput.prov_e1 = (ep_1s);fitoutput.prov_e2 = (ep2_1s);fitoutput.prov_e3 = (ep3_1s)
    fitoutput.prov_e1_2 = (ep_1o);fitoutput.prov_e2_2 = (ep2_1o);fitoutput.prov_e3_2 = (ep3_1o)
    fitoutput.prov_a1 = (oneresu.aic);fitoutput.prov_a2 = (tworesu.aic);fitoutput.prov_a3 = (threeresu.aic)
    fitoutput.prov_b1 = (oneresu.bic);fitoutput.prov_b2 = (tworesu.bic);fitoutput.prov_b3 = (threeresu.bic)
    #trigger = 'Y' # EDITAR ESTO! BUSCAR ESTO TODo 
    #############################################################################################################
    # Now we start the fit to the complete spectra depending on the selected method
    #############################################################################################################
    if trigger == 'Y':
        # There is no secondary or third component
        # thus one can add a 0 value directly here to the second comp list
        fitoutput.prov_SS = (np.nan);fitoutput.prov_SN2 = (np.nan)
        fitoutput.prov_VS = (np.nan);fitoutput.prov_VN2 = (np.nan)
        fitoutput.prov_eSS = (np.nan);fitoutput.prov_eSN2 = (np.nan)
        fitoutput.prov_eVS = (np.nan);fitoutput.prov_eVN2 = (np.nan)
        fitoutput.prov_fSHa = (np.nan);fitoutput.prov_fSN1 = (np.nan);fitoutput.prov_fSN2 = (np.nan);fitoutput.prov_fSS1 = (np.nan);fitoutput.prov_fSS2 = (np.nan);fitoutput.prov_fSO1 = (np.nan);fitoutput.prov_fSO2 = (np.nan);fitoutput.prov_f2Ha = (np.nan);fitoutput.prov_f2N1 = (np.nan);fitoutput.prov_f2N2 = (np.nan);fitoutput.prov_f2S1 = (np.nan);fitoutput.prov_f2S2 = (np.nan);fitoutput.prov_f2O1 = (np.nan);fitoutput.prov_f2O2 = (np.nan)
        ########################################################################################################
        params = lmfit.Parameters() # initialise the parameters
        # Now we define the initial guesses and the constraints
        # as both have already S and O lines fitted, we can take them out of the if
        cd1 = lmfit.Parameter('mu_0', value=oneresu.values['mu_0'],vary=False)
        de = lmfit.Parameter('sig_0', value=oneresu.values['sig_0'],vary=False)
        ef = lmfit.Parameter('amp_0', value=oneresu.values['amp_0'],vary=False)
        fg = lmfit.Parameter('mu_1', value=oneresu.values['mu_1'],vary=False)
        gh = lmfit.Parameter('sig_1', value=oneresu.values['sig_1'],vary=False)
        hi = lmfit.Parameter('amp_1', value=oneresu.values['amp_1'],vary=False)
        rs = lmfit.Parameter('mu_5', value=oneresu.values['mu_5'],vary=False)# mu5,expr='mu_0*(6300.304/6730.82)')
        st = lmfit.Parameter('sig_5', value=oneresu.values['sig_5'],vary=False)
        tu = lmfit.Parameter('amp_5', value=oneresu.values['amp_5'],vary=False)
        uv = lmfit.Parameter('mu_6', value=oneresu.values['mu_6'],vary=False)
        vw = lmfit.Parameter('sig_6', value=oneresu.values['sig_6'],vary=False)
        wy = lmfit.Parameter('amp_6', value=oneresu.values['amp_6'],vary=False)
        if meth == 'M1':
            ij = lmfit.Parameter('mu_2', value=mu2,expr='mu_0*(6583.45/6730.82)')
            jk = lmfit.Parameter('sig_2', value=sig2,expr='sig_0')
            kl = lmfit.Parameter('amp_2', value=amp2,min=0)
            lm = lmfit.Parameter('mu_3', value=mu3,expr='mu_5*(6562.801/6300.304)')
            mn = lmfit.Parameter('sig_3', value=sig3,expr='sig_5')
            no = lmfit.Parameter('amp_3', value=amp3,min=0.)
            op = lmfit.Parameter('mu_4', value=mu4,expr='mu_0*(6548.05/6730.82)')
            pq = lmfit.Parameter('sig_4', value=sig4,expr='sig_0')
            qr = lmfit.Parameter('amp_4', value=amp4,min=0.,expr='amp_2*(1./3.)')

        elif meth == 'M2':
            ij = lmfit.Parameter('mu_2', value=mu2,expr='mu_5*(6583.45/6300.304)')
            jk = lmfit.Parameter('sig_2', value=sig2,expr='sig_5')
            kl = lmfit.Parameter('amp_2', value=amp2,min=0)
            lm = lmfit.Parameter('mu_3', value=mu3,expr='mu_5*(6562.801/6300.304)')
            mn = lmfit.Parameter('sig_3', value=sig3,expr='sig_5')
            no = lmfit.Parameter('amp_3', value=amp3,min=0.)
            op = lmfit.Parameter('mu_4', value=mu4,expr='mu_5*(6548.05/6300.304)')
            pq = lmfit.Parameter('sig_4', value=sig4,expr='sig_5')
            qr = lmfit.Parameter('amp_4', value=amp4,min=0.,expr='amp_2*(1./3.)')

        params.add_many(sl,it,cd1,de,ef,fg,gh,hi,rs,st,tu,uv,vw,wy,ij,jk,kl,lm,mn,no,op,pq,qr)        
 
        # Initial guesses and fit
        resu1 = comp_mod.fit(data_cor,params,x=l)
        #lmfit.model.save_modelresult(resu1, parentFold+'one_modelresult.sav')
        with open(parentFold+'fitone_result.txt', 'w') as fh:
            fh.write(resu1.fit_report())
        ########################### Calculate gaussians and final fit #############################
        # Now we create and plot the individual gaussians of the fit
        gaus1 = Ofuncts.gaussian(l,resu1.values['mu_0'],resu1.values['sig_0'],resu1.values['amp_0'])
        gaus2 = Ofuncts.gaussian(l,resu1.values['mu_1'],resu1.values['sig_1'],resu1.values['amp_1'])
        gaus3 = Ofuncts.gaussian(l,resu1.values['mu_2'],resu1.values['sig_2'],resu1.values['amp_2'])
        gaus4 = Ofuncts.gaussian(l,resu1.values['mu_3'],resu1.values['sig_3'],resu1.values['amp_3'])
        gaus5 = Ofuncts.gaussian(l,resu1.values['mu_4'],resu1.values['sig_4'],resu1.values['amp_4'])
        gaus6 = Ofuncts.gaussian(l,resu1.values['mu_5'],resu1.values['sig_5'],resu1.values['amp_5'])
        gaus7 = Ofuncts.gaussian(l,resu1.values['mu_6'],resu1.values['sig_6'],resu1.values['amp_6'])
        fin_fit = resu1.best_fit

        ########################################################################################################
        # TODo LISTO HASTA AQUI A 17 DE NOVIEMBRE: hay que arreglar esta parte del epsilon como en referplot
        ########################################################################################################

        # one component
        stdf_s2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10]-fin_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10])
        stdf_s1 = np.std(data_cor[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]]-fin_fit[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]])
        stdf_n2 = np.std(data_cor[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10]-fin_fit[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10])
        stdf_ha = np.std(data_cor[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]]-fin_fit[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]])
        stdf_n1 = np.std(data_cor[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]]-fin_fit[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]])
        stdf_o1 = np.std(data_cor[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]]-fin_fit[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]])
        stdf_o2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-fin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 components is... ')
        print('         For SII2: '+str(stdf_s2/stadev)+' < 3')
        print('         For SII1: '+str(stdf_s1/stadev)+' < 3')
        print('         For NII2: '+str(stdf_n2/stadev)+' < 3')
        print('         For Halpha: '+str(stdf_ha/stadev)+' < 3')
        print('         For NII1: '+str(stdf_n1/stadev)+' < 3')

        if os.path.exists(parentFold+'eps_adj'+str(meth)+'_1.txt'): os.remove(parentFold+'eps_adj'+str(meth)+'_1.txt')
        np.savetxt(parentFold+'eps_adj'+str(meth)+'_1.txt',np.c_[stdf_s2/stadev,stdf_s1/stadev,stdf_n2/stadev,stdf_ha/stadev,stdf_n1/stadev,stdf_o2/stadev,stdf_o1/stadev,resu1.chisqr], ('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('SII2\tSII1\tNII2\tHa\tNII1\tOI1\tOI2\tChi2'))

        if resu1.values['mu_0'] > 6800:
            maxfS1,maxfS2,maxfN1,maxfHa,maxfN2,maxfO1,maxfO2 = 0.,0.,0.,0.,0.,0.,0.
        else:
            maxfS1 = fin_fit[np.where(abs(resu1.values['mu_0']-l)<0.5)[0][0]]
            maxfS2 = fin_fit[np.where(abs(resu1.values['mu_1']-l)<0.5)[0][0]]
            maxfN1 = fin_fit[np.where(abs(resu1.values['mu_2']-l)<0.5)[0][0]]
            maxfHa = fin_fit[np.where(abs(resu1.values['mu_3']-l)<0.5)[0][0]]
            maxfN2 = fin_fit[np.where(abs(resu1.values['mu_4']-l)<0.5)[0][0]]
            maxfO1 = fin_fit[np.where(abs(resu1.values['mu_5']-l)<0.5)[0][0]]
            maxfO2 = fin_fit[np.where(abs(resu1.values['mu_6']-l)<0.5)[0][0]]

        # Estimate the sigma for both SII and OI lines
        sigS2 = pix_to_v*np.sqrt(oneresu.values['sig_0']**2-sig_inst**2)
        sigO2 = pix_to_v*np.sqrt(oneresu.values['sig_5']**2-sig_inst**2)

        if oneresu.params['sig_0'].stderr == None:
            print('Problem determining the errors! First component sigma ')
            esigS2 = 0.
        elif oneresu.params['sig_0'].stderr != None:
            esigS2 = pix_to_v*(2*oneresu.values['sig_0']*oneresu.params['sig_0'].stderr)/(np.sqrt(oneresu.values['sig_0']**2-sig_inst**2))
        if oneresu.params['sig_5'].stderr == None:
            print('Problem determining the errors! First component sigma ')
            esigO2 = 0.
        elif oneresu.params['sig_5'].stderr != None:
            esigO2 = pix_to_v*(2*oneresu.values['sig_5']*oneresu.params['sig_5'].stderr)/(np.sqrt(oneresu.values['sig_5']**2-sig_inst**2))

        # Estimate the velocities for both SII and OI lines
        vS2 = v_luz*((resu1.values['mu_0']-l_SII_2)/l_SII_2)
        vO2 = v_luz*((resu1.values['mu_5']-l_OI_1)/l_OI_1)
        if oneresu.params['mu_0'].stderr == None:
            print('Problem determining the errors! First component [SII] lines')
            evS2 = 0.
        elif oneresu.params['mu_0'].stderr != None:
            evS2 = ((v_luz/l_SII_2)*oneresu.params['mu_0'].stderr)
        if oneresu.params['mu_5'].stderr == None:
            print('Problem determining the errors! First component [OI] line')
            evO2 = 0.
        elif oneresu.params['mu_5'].stderr != None:
            evO2 = ((v_luz/l_OI_1)*oneresu.params['mu_5'].stderr)

        # Save all the v and sigma in txt files in the parent folder
        if os.path.exists(parentFold+'v_sig_adj_1.txt'): os.remove(parentFold+'v_sig_adj_1.txt')
        np.savetxt(parentFold+'v_sig_adj'+str(meth)+'_1.txt',np.c_[vS2,evS2,sigS2,esigS2,vO2,evO2,sigO2,esigO2],('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('v_refS\tev_refS\tsig_refS\tesig_refS\tv_refO\tev_refO\tsig_refO\tesig_refO'))
        if os.path.exists(parentFold+'fluxes_1C.txt'): os.remove(parentFold+'fluxes_1C.txt')
        np.savetxt(parentFold+'fluxes_1C.txt',np.c_[sum(gaus1),sum(gaus2),sum(gaus3),sum(gaus4),sum(gaus5),sum(gaus6),sum(gaus7)],('%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))

        # Saving the max flux of each line to do the final flux ma;
        #fitoutput.prov_fS2 = (sum(gaus1));fitoutput.prov_fS1 = (sum(gaus2));fitoutput.prov_fN2 = (sum(gaus3));fitoutput.prov_fHa = (sum(gaus4));fitoutput.prov_fN1 = (sum(gaus5));fitoutput.prov_fO1 = (sum(gaus6));fitoutput.prov_fO2 = (sum(gaus7))

        ############################ PLOT ############################
        plt.close('all')
        # MAIN plot
        fig1   = plt.figure(1,figsize=(10, 9))
        frame1 = fig1.add_axes((.1,.25,.85,.65))             # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
        plt.plot(l,data_cor,'k',linewidth=2)                 # Initial data
        plt.plot(l[std0:std1],data_cor[std0:std1],'y',linewidth=4.)  # Zone where the stddev is calculated
        plt.plot(l[std0:std1],data_cor[std0:std1],'k',linewidth=1)                   # Initial data
        #plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='y',linestyle=(0, (5, 8)))#,label='Linear fit')
        plt.plot(l,gaus1+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,gaus2+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,gaus6+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,gaus7+lin_data_fin,c='goldenrod',linestyle='-')
        if meth == 'M1': 
            plt.plot(l,gaus3+lin_data_fin,c='darkgreen',linestyle='-')
            plt.plot(l,gaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,gaus5+lin_data_fin,c='darkgreen',linestyle='-')
        elif meth == 'M2': 
            plt.plot(l,gaus3+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,gaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,gaus5+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,fin_fit,'r-')
        textstr = '\n'.join((r'$V_{SII_{2}}$ = '+ '{:.2f} +- {:.2f}'.format(vS2,evS2),
                         r'$V_{OI_{1}}$ = '+ '{:.2f} +- {:.2f}'.format(vO2,evO2),
                         r'$\sigma_{SII_{2}}$ = '+ '{:.2f} +- {:.2f}'.format(sigS2,esigS2),
                         r'$\sigma_{OI_{1}}$ = '+ '{:.2f} +- {:.2f}'.format(sigO2,esigO2),
                         r'$\frac{F_{SII_{2}}}{F_{SII_{1}}}$ = '+ '{:.3f}'.format(maxfS2/maxfS1)))

        frame1.set_xticklabels([])                      # Remove x-tic labels for the first frame
        plt.ylabel(r'Flux (erg $\rm cm^{-2} s^{-1} \AA^{-1}$)',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(prune='lower'))
        #frame1.yaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
        plt.text(0.87,0.94,'N-SII',color='darkgreen',transform=frame1.transAxes,fontsize=21)
        plt.text(0.87,0.9,'N-OI',color='goldenrod',transform=frame1.transAxes,fontsize=21)

        # RESIDUAL plot
        frame2 = fig1.add_axes((.1,.1,.85,.15))
        plt.plot(l,np.zeros(len(l)),c='orange',linestyle='-')           # Line around zero
        plt.plot(l,data_cor-fin_fit,c='k')              # Main
        plt.xlabel(r'Wavelength ($\rm \AA$)',fontsize=19)
        plt.ylabel('Residuals',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.plot(l,np.zeros(len(l))+1.5*stadev,c='orange',linestyle=(0,(5,8)))  # 3 sigma upper limit
        plt.plot(l,np.zeros(len(l))-1.5*stadev,c='orange',linestyle=(0,(5,8)))  # 3 sigma down limit
        plt.ylim(-(3*stadev)*3,(3*stadev)*3)

        plt.savefig(parentFold+'adj_full_1comp.pdf',format='pdf',bbox_inches='tight',pad_inches=0.2)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        frame1.text(6250,max(data_cor), textstr, fontsize=12,verticalalignment='top', bbox=props)
        plt.savefig(parentFold+'adj_full_1comp.png',bbox_inches='tight',pad_inches=0.2)

        fitoutput.prov_VN = (vS2);fitoutput.prov_eVN = (evS2);fitoutput.prov_SN = (sigS2);fitoutput.prov_eSN = (esigS2)
        fitoutput.prov_VN_2 = (vO2);fitoutput.prov_eVN_2 = (evO2);fitoutput.prov_SN_2 = (sigO2);fitoutput.prov_eSN_2 = (esigO2)

        ##############################################################
	##############################################################
        # Ask if a broad component is needed 
        trigger2 = 'N'
        if galaxy >= 15 and galaxy <= 20 and galaxy2 >= 12 and galaxy2 <= 16: 
            if stdf_n2/stadev > 3 or stdf_ha/stadev > 3 or stdf_n1/stadev > 3: trigger2 = 'Y'
            else: trigger2 = 'N'
        trigger2 = 'N'
        #if galaxy >= 15 and galaxy <= 20 and galaxy2 >= 12 and galaxy2 <= 16: trigger2 = 'Y' #33  30
        # Ask if a broad component is needed 
        if trigger2 == 'N':
            print('The final plots are already printed and have been already saved!')
            print('No broad component from the BLR added to the fit of this spaxel')
            # Append all the values to the corresponding vector
            fitoutput.prov_fS2 = (sum(gaus1));fitoutput.prov_fS1 = (sum(gaus2));fitoutput.prov_fN2 = (sum(gaus3));fitoutput.prov_fHa = (sum(gaus4));fitoutput.prov_fN1 = (sum(gaus5));fitoutput.prov_fO1 = (sum(gaus6));fitoutput.prov_fO2 = (sum(gaus7))
            fitoutput.prov_SB = (np.nan);fitoutput.prov_eSB = (np.nan);fitoutput.prov_VB = (np.nan);fitoutput.prov_eVB = (np.nan)
            np.savetxt(parentFold+'fit1comp_best_values.txt',np.c_[l,resu1.data,resu1.best_fit,lin_data_fin,gaus1,gaus2,gaus3,gaus4,gaus5,gaus6,gaus7],fmt=('%5.6f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Wavelength\tReal_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tNarrow_SII6716\tNarrow_NII6584\tNarrow_Halpha\tNarrow_NII6548\tNarrow_OI6300\tNarrow_OI6363'))

        elif trigger2 == 'Y':
            # Now we define the initial guesses and the constraints
            newxb = l[np.where(l>l9)[0][0]:np.where(l<l6)[0][-1]]
            newyb = data_cor[np.where(l>l9)[0][0]:np.where(l<l6)[0][-1]]
            sigb = 16.
            mub  = mu3
            ampb = amp3/3.
            paramsbH = lmfit.Parameters()
            # broad components
            ab = lmfit.Parameter('mu_b',value=mub)#6605.,vary=False)
            bc = lmfit.Parameter('sig_b',value=sigb,min=12.)#,vary=False) sig_inst)minbroad
            yz = lmfit.Parameter('amp_b',value=ampb,min=0.)
            paramsbH.add_many(sl,it,ab,bc,yz,cd1,de,ef,fg,gh,hi,rs,st,tu,uv,vw,wy,ij,jk,kl,lm,mn,no,op,pq,qr)
  
            broadresu = broad_mod.fit(data_cor,paramsbH,x=l)
            lmfit.model.save_modelresult(broadresu, parentFold+'broadone_modelresult.sav')
            with open(parentFold+'fitbroad_result.txt', 'w') as fh:
                fh.write(broadresu.fit_report())
 
            # PLOT AND PRINT THE RESULTS 
            refer2 = broad_plot(parentFold,l,data_cor,meth,trigger,linresu,oneresu,fin_fit,broadresu,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,std0,std1)
            broad_fit = broadresu.best_fit
            # Save the results in the corresponding array
            fitoutput.prov_VB = refer2[0]
            fitoutput.prov_eVB = refer2[1]
            fitoutput.prov_SB = refer2[2]
            fitoutput.prov_eSB = refer2[3]
            fitoutput.prov_fS2 = sum(refer2[4]);fitoutput.prov_fS1 = (sum(refer2[5]));fitoutput.prov_fN2 = (sum(refer2[6]));fitoutput.prov_fHa = (sum(refer2[7]));fitoutput.prov_fN1 = (sum(refer2[8]));fitoutput.prov_fO1 = (sum(refer2[9]));fitoutput.prov_fO2 = (sum(refer2[10]))
            print('Info saved for broad component with 1 Narrow fit of this spaxel')

    ##############################################################
    ##############################################################

    elif trigger == 'N':
        params2c = lmfit.Parameters()
        cd1 = lmfit.Parameter('mu_0', value=tworesu.values["mu_0"],vary=False)
        de = lmfit.Parameter('sig_0', value=tworesu.values["sig_0"],vary=False)
        ef = lmfit.Parameter('amp_0', value=tworesu.values["amp_0"],vary=False)
        fg = lmfit.Parameter('mu_1', value=tworesu.values["mu_1"],vary=False)
        gh = lmfit.Parameter('sig_1', value=tworesu.values["sig_1"],vary=False)
        hi = lmfit.Parameter('amp_1', value=tworesu.values["amp_1"],vary=False)
        rs = lmfit.Parameter('mu_5', value=tworesu.values["mu_5"],vary=False)
        st = lmfit.Parameter('sig_5', value=tworesu.values["sig_5"],vary=False)
        tu = lmfit.Parameter('amp_5', value=tworesu.values["amp_5"],vary=False)
        uv = lmfit.Parameter('mu_6', value=tworesu.values["mu_6"],vary=False)
        vw = lmfit.Parameter('sig_6', value=tworesu.values["sig_6"],vary=False)
        wy = lmfit.Parameter('amp_6', value=tworesu.values["amp_6"],vary=False)
        aaa = lmfit.Parameter('mu_20', value=tworesu.values["mu_20"],vary=False)
        aab = lmfit.Parameter('sig_20', value=tworesu.values["sig_20"],vary=False)
        aac = lmfit.Parameter('amp_20', value=tworesu.values["amp_20"],vary=False)
        aad = lmfit.Parameter('mu_21', value=tworesu.values["mu_21"],vary=False)
        aae = lmfit.Parameter('sig_21', value=tworesu.values["sig_21"],vary=False)
        aaf = lmfit.Parameter('amp_21', value=tworesu.values["amp_21"],vary=False)
        aap = lmfit.Parameter('mu_25', value=tworesu.values["mu_25"],vary=False)
        aaq = lmfit.Parameter('sig_25', value=tworesu.values["sig_25"],vary=False)
        aar = lmfit.Parameter('amp_25', value=tworesu.values["amp_25"],vary=False)
        aas = lmfit.Parameter('mu_26', value=tworesu.values["mu_26"],vary=False)
        aat = lmfit.Parameter('sig_26', value=tworesu.values["sig_26"],vary=False)
        aau = lmfit.Parameter('amp_26', value=tworesu.values["amp_26"],vary=False)
        if meth == 'M1':
            # Now we define the initial guesses and the constraints
            ij = lmfit.Parameter('mu_2', value=mu2,expr='mu_0*(6584./6731.)')
            jk = lmfit.Parameter('sig_2', value=sig2,expr='sig_0')
            kl = lmfit.Parameter('amp_2', value=amp2,min=0.)
            lm = lmfit.Parameter('mu_3', value=mu3,expr='mu_5*(6563./6300.3)')
            mn = lmfit.Parameter('sig_3', value=sig3,expr='sig_5')
            no = lmfit.Parameter('amp_3', value=amp3,min=0.)
            op = lmfit.Parameter('mu_4', value=mu4,expr='mu_0*(6548./6731.)')
            pq = lmfit.Parameter('sig_4', value=sig4,expr='sig_0')
            qr = lmfit.Parameter('amp_4', value=amp4,min=0.,expr='amp_2*(1./3.)')
            aag = lmfit.Parameter('mu_22', value=mu2,expr='mu_20*(6584./6731.)')
            aah = lmfit.Parameter('sig_22', value=sig22,expr='sig_20')
            aai = lmfit.Parameter('amp_22', value=amp22,min=0.)
            aaj = lmfit.Parameter('mu_23', value=mu3,expr='mu_25*(6563./6300.3)')
            aak = lmfit.Parameter('sig_23', value=sig23,expr='sig_25')
            aal = lmfit.Parameter('amp_23', value=amp23,min=0.)
            aam = lmfit.Parameter('mu_24', value=mu4,expr='mu_20*(6548./6731.)')
            aan = lmfit.Parameter('sig_24', value=sig24,expr='sig_20')
            aao = lmfit.Parameter('amp_24', value=amp24,min=0.,expr='amp_22*(1./3.)')

        elif meth == 'M2':
            # Now we define the initial guesses and the constraints
            ij = lmfit.Parameter('mu_2', value=mu2,expr='mu_5*(6584./6300.30)')
            jk = lmfit.Parameter('sig_2', value=sig2,expr='sig_5')
            kl = lmfit.Parameter('amp_2', value=amp2,min=0.)
            lm = lmfit.Parameter('mu_3', value=mu3,expr='mu_5*(6563./6300.30)')
            mn = lmfit.Parameter('sig_3', value=sig3,expr='sig_5')
            no = lmfit.Parameter('amp_3', value=amp3,min=0.)
            op = lmfit.Parameter('mu_4', value=mu4,expr='mu_5*(6548./6300.30)')
            pq = lmfit.Parameter('sig_4', value=sig4,expr='sig_5')
            qr = lmfit.Parameter('amp_4', value=amp4,min=0.,expr='amp_2*(1./3.)')
            aag = lmfit.Parameter('mu_22', value=mu2,expr='mu_25*(6584./6300.30)')
            aah = lmfit.Parameter('sig_22', value=sig22,expr='sig_25')
            aai = lmfit.Parameter('amp_22', value=amp22,min=0.)
            aaj = lmfit.Parameter('mu_23', value=mu3,expr='mu_25*(6563./6300.30)')
            aak = lmfit.Parameter('sig_23', value=sig23,expr='sig_25')
            aal = lmfit.Parameter('amp_23', value=amp23,min=0.)
            aam = lmfit.Parameter('mu_24', value=mu4,expr='mu_25*(6548./6300.30)')
            aan = lmfit.Parameter('sig_24', value=sig24,expr='sig_25')
            aao = lmfit.Parameter('amp_24', value=amp24,min=0.,expr='amp_22*(1./3.)')

        params2c.add_many(sl,it,cd1,de,ef,fg,gh,hi,rs,st,tu,uv,vw,wy,aaa,aab,aac,aad,aae,aaf,aap,aaq,aar,aas,aat,aau,ij,jk,kl,lm,mn,no,op,pq,qr,aag,aah,aai,aaj,aak,aal,aam,aan,aao)

	# FIT
        twocompresu = twocomp_mod.fit(data_cor,params2c,x=l)
	
        lmfit.model.save_modelresult(twocompresu, parentFold+'two_modelresult.sav')
        with open(parentFold+'fit_'+str(meth)+'_two_result.txt', 'w') as fh:
            fh.write(twocompresu.fit_report())

	####################### Calculate gaussians and final fit #######################
	# Now we create and plot the individual gaussians of the fit
        tgaus1 = Ofuncts.gaussian(l,twocompresu.values['mu_0'],twocompresu.values['sig_0'],twocompresu.values['amp_0']) 
        tgaus2 = Ofuncts.gaussian(l,twocompresu.values['mu_1'],twocompresu.values['sig_1'],twocompresu.values['amp_1'])
        tgaus3 = Ofuncts.gaussian(l,twocompresu.values['mu_2'],twocompresu.values['sig_2'],twocompresu.values['amp_2']) 
        tgaus4 = Ofuncts.gaussian(l,twocompresu.values['mu_3'],twocompresu.values['sig_3'],twocompresu.values['amp_3'])
        tgaus5 = Ofuncts.gaussian(l,twocompresu.values['mu_4'],twocompresu.values['sig_4'],twocompresu.values['amp_4'])
        tgaus6 = Ofuncts.gaussian(l,twocompresu.values['mu_5'],twocompresu.values['sig_5'],twocompresu.values['amp_5'])
        tgaus7 = Ofuncts.gaussian(l,twocompresu.values['mu_6'],twocompresu.values['sig_6'],twocompresu.values['amp_6'])
        tgaus8 = Ofuncts.gaussian(l,twocompresu.values['mu_20'],twocompresu.values['sig_20'],twocompresu.values['amp_20']) 
        tgaus9 = Ofuncts.gaussian(l,twocompresu.values['mu_21'],twocompresu.values['sig_21'],twocompresu.values['amp_21'])
        tgaus10 = Ofuncts.gaussian(l,twocompresu.values['mu_22'],twocompresu.values['sig_22'],twocompresu.values['amp_22']) 
        tgaus11 = Ofuncts.gaussian(l,twocompresu.values['mu_23'],twocompresu.values['sig_23'],twocompresu.values['amp_23'])
        tgaus12 = Ofuncts.gaussian(l,twocompresu.values['mu_24'],twocompresu.values['sig_24'],twocompresu.values['amp_24'])
        tgaus13 = Ofuncts.gaussian(l,twocompresu.values['mu_25'],twocompresu.values['sig_25'],twocompresu.values['amp_25'])
        tgaus14 = Ofuncts.gaussian(l,twocompresu.values['mu_26'],twocompresu.values['sig_26'],twocompresu.values['amp_26'])
        fin2_fit = twocompresu.best_fit

    ########################################################################################################
    # TODo LISTO HASTA AQUI A 24 DE NOVIEMBRE EN TODOS LOS CODIGOS A LOS QUE SE LLAMAN.
    ########################################################################################################
	# two components
        stdf2_s2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10]-fin2_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10])
        stdf2_s1 = np.std(data_cor[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]]-fin2_fit[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]])
        stdf2_n2 = np.std(data_cor[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10]-fin2_fit[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10])
        stdf2_ha = np.std(data_cor[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]]-fin2_fit[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]])
        stdf2_n1 = np.std(data_cor[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]]-fin2_fit[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]])
        stdf2_o1 = np.std(data_cor[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]]-fin2_fit[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]])
        stdf2_o2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-fin2_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 2 components is... ')
        print('		For SII2: '+str(stdf2_s2/stadev)+' < 3')
        print('		For SII1: '+str(stdf2_s1/stadev)+' < 3')
        print('		For NII2: '+str(stdf2_n2/stadev)+' < 3')
        print('		For Halpha: '+str(stdf2_ha/stadev)+' < 3')
        print('		For NII1: '+str(stdf2_n1/stadev)+' < 3')
        print('         For OI1: '+str(stdf2_o1/stadev)+' < 3')
        print('         For OI2: '+str(stdf2_o2/stadev)+' < 3')
	
        if os.path.exists(parentFold+'eps_adj_2C.txt'): os.remove(parentFold+'eps_adj_2C.txt')
	np.savetxt(parentFold+'eps_adj'+str(meth)+'_2.txt',np.c_[stdf2_s2/stadev,stdf2_s1/stadev,stdf2_n2/stadev,stdf2_ha/stadev,stdf2_n1/stadev,stdf2_o1/stadev,stdf2_o2/stadev,twocompresu.chisqr], ('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('SII2\tSII1\tNII2\tHa\tNII1\tOI1\tOI2\tChi2'))

        try:
	    # We determine the maximum flux of the fit for all the lines, and the velocity and sigma components
            max2S1 = fin2_fit[np.where(abs(twocompresu.values['mu_0']-l)<0.5)[0][0]] 
            max2S2 = fin2_fit[np.where(abs(twocompresu.values['mu_1']-l)<0.5)[0][0]] 
            max2N1 = fin2_fit[np.where(abs(twocompresu.values['mu_2']-l)<0.5)[0][0]] 
            max2Ha = fin2_fit[np.where(abs(twocompresu.values['mu_3']-l)<0.5)[0][0]] 
            max2N2 = fin2_fit[np.where(abs(twocompresu.values['mu_4']-l)<0.5)[0][0]] 
            max2O1 = fin2_fit[np.where(abs(twocompresu.values['mu_5']-l)<0.5)[0][0]]
            max2O2 = fin2_fit[np.where(abs(twocompresu.values['mu_6']-l)<0.5)[0][0]]
        except IndexError: 
            print('ERROR: index out of range. Setting the flux values of the OI 1 line to 0.')
	
        # two comps
	sig2S2 = pix_to_v*np.sqrt(twocompresu.values['sig_0']**2-sig_inst**2)
	sig2O2 = pix_to_v*np.sqrt(twocompresu.values['sig_5']**2-sig_inst**2)
	sig20S2 = pix_to_v*np.sqrt(twocompresu.values['sig_20']**2-sig_inst**2)
	sig20O2 = pix_to_v*np.sqrt(twocompresu.values['sig_25']**2-sig_inst**2)
	
        if tworesu.params['sig_0'].stderr == None: 
             print('Problem determining the errors! First component sigma SII')
             esig2S2 = 0.
        elif tworesu.params['sig_0'].stderr != None: 
             esig2S2 = pix_to_v*(2*twocompresu.values['sig_0']*tworesu.params['sig_0'].stderr)/(np.sqrt(twocompresu.values['sig_0']**2-sig_inst**2))
        if tworesu.params['sig_5'].stderr == None: 
             print('Problem determining the errors! First component sigma OI')
             esig2O2 = 0.
        elif tworesu.params['sig_5'].stderr != None: 
             esig2O2 = pix_to_v*(2*twocompresu.values['sig_5']*tworesu.params['sig_5'].stderr)/(np.sqrt(twocompresu.values['sig_5']**2-sig_inst**2))
	     
        if tworesu.params['sig_20'].stderr == None:
            print('Problem determining the errors! Second component sigma ')
            esig20S2 = 0.
        elif tworesu.params['sig_20'].stderr != None:
            esig20S2 = pix_to_v*(2*twocompresu.values['sig_20']*tworesu.params['sig_20'].stderr)/(np.sqrt(twocompresu.values['sig_20']**2-sig_inst**2))
        if tworesu.params['sig_25'].stderr == None:
            print('Problem determining the errors! Second component sigma ')
            esig20O2 = 0.
        elif tworesu.params['sig_25'].stderr != None:
            esig20O2 = pix_to_v*(2*twocompresu.values['sig_25']*tworesu.params['sig_25'].stderr)/(np.sqrt(twocompresu.values['sig_25']**2-sig_inst**2))
	    
        v2S2 = v_luz*((twocompresu.values['mu_0']-l_SII_2)/l_SII_2)
        v2O2 = v_luz*((twocompresu.values['mu_5']-l_OI_1)/l_OI_1)
        v20S2 = v_luz*((twocompresu.values['mu_20']-l_SII_2)/l_SII_2)
        v20O2 = v_luz*((twocompresu.values['mu_25']-l_OI_1)/l_OI_1)
        if tworesu.params['mu_0'].stderr == None: 
            print('Problem determining the errors! First component ')
            ev2S2= 0.
        elif tworesu.params['mu_0'].stderr != None:
            print('Problem determining the errors! Second component ')
            ev2S2 = ((v_luz/l_SII_2)*tworesu.params['mu_0'].stderr)
        if tworesu.params['mu_5'].stderr == None: 
            print('Problem determining the errors! First component ')
            ev2O2= 0.
        elif tworesu.params['mu_5'].stderr != None:
            print('Problem determining the errors! Second component ')
            ev2O2 = ((v_luz/l_OI_1)*tworesu.params['mu_5'].stderr)
        if tworesu.params['mu_20'].stderr == None: 
            ev20S2 = 0.
        elif tworesu.params['mu_20'].stderr != None: 
            ev20S2 = ((v_luz/l_SII_2)*tworesu.params['mu_20'].stderr)
        if tworesu.params['mu_25'].stderr == None: 
            ev20O2 = 0.
        elif tworesu.params['mu_25'].stderr != None: 
            ev20O2 = ((v_luz/l_OI_q)*tworesu.params['mu_25'].stderr)

        # Save txt with the information from velocity, sigma and fluxes
        #	
        # Save the velocity and sigma for all components
        if os.path.exists(parentFold+'v_sig_adj_2C.txt'): os.remove(parentFold+'v_sig_adj_2C.txt')
	np.savetxt(parentFold+'v_sig_adj_2C.txt',np.c_[v2S2,ev2S2,v20S2,ev20S2,sig2S2,esig2S2,sig20S2,esig20S2,v2O2,ev2O2,v20O2,ev20O2,sig2O2,esig2O2,sig20O2,esig20O2],('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('v_refS2\tev_refS2\tv_2refS2\tev_2refS2\tsig_refS2\tesig_refS2\tsig_2refS2\tesig_2refS2\tv_refO2\tev_refO2\tv_2refO2\tev_2refO2\tsig_refO2\tesig_refO2\tsig_2refO2\tesig_2refO2'))
        
        # Save the fluxes for all components
        if os.path.exists(parentFold+'fluxes_'+str(meth)+'_2_Ncomp.txt'): os.remove(parentFold+'fluxes_2C_Ncomp.txt')
        np.savetxt(parentFold+'fluxes_2C_Ncomp.txt',np.c_[sum(tgaus1),sum(tgaus2),sum(tgaus3),sum(tgaus4),sum(tgaus5),sum(tgaus6),sum(tgaus7)],fmt=('%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))
        if os.path.exists(parentFold+'fluxes_2C_Scomp.txt'): os.remove(parentFold+'fluxes_2C_Scomp.txt')
        np.savetxt(parentFold+'fluxes_2C_Scomp.txt',np.c_[sum(tgaus8),sum(tgaus9),sum(tgaus10),sum(tgaus11),sum(tgaus12),sum(tgaus13),sum(tgaus14)],('%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))

	########################### PLOT #############################
        plt.close('all')
	# MAIN plot
        fig1   = plt.figure(1,figsize=(10, 9))
        frame1 = fig1.add_axes((.1,.25,.85,.65)) 	     # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
        plt.plot(l,data_cor,'k',linewidth=2)		     # Initial data
        plt.plot(l[std0:std1],data_cor[std0:std1],c='y',linewidth=4)  # Zone where the stddev is calculated
        plt.plot(l[std0:std1],data_cor[std0:std1],'k',linewidth=1)		     # Initial data
        #plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='y',linestyle=(0, (5, 8)),label='Linear fit')
        plt.plot(l,tgaus1+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,tgaus2+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,tgaus6+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,tgaus7+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,tgaus8+lin_data_fin,c='dodgerblue',linestyle='-')
        plt.plot(l,tgaus9+lin_data_fin,c='dodgerblue',linestyle='-')
        plt.plot(l,tgaus13+lin_data_fin,c='darkred',linestyle='-')
        plt.plot(l,tgaus14+lin_data_fin,c='darkred',linestyle='-')
        if meth == 'M1':
            plt.plot(l,tgaus3+lin_data_fin,c='darkgreen',linestyle='-')
            plt.plot(l,tgaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,tgaus5+lin_data_fin,c='darkgreen',linestyle='-')
            plt.plot(l,tgaus10+lin_data_fin,c='dodgerblue',linestyle='-')
            plt.plot(l,tgaus11+lin_data_fin,c='darkred',linestyle='-')
            plt.plot(l,tgaus12+lin_data_fin,c='dodgerblue',linestyle='-')
        if meth == 'M2':
            plt.plot(l,tgaus3+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,tgaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,tgaus5+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,tgaus10+lin_data_fin,c='darkred',linestyle='-')
            plt.plot(l,tgaus11+lin_data_fin,c='darkred',linestyle='-')
            plt.plot(l,tgaus12+lin_data_fin,c='darkred',linestyle='-')
        plt.plot(l,fin2_fit,'r-')
        textstr = '\n'.join((r'$V_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(v2S2,ev2S2),
                            r'$V_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(v2O2,ev2O2),
			    r'$V_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v20S2,ev20S2),
			    r'$V_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v20O2,ev20O2),
			    r'$\sigma_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig2S2,esig2S2),
			    r'$\sigma_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig2O2,esig2O2),
			    r'$\sigma_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig20S2,esig20S2),
			    r'$\sigma_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig20O2,esig20O2),
			    r'$F_{H_{\alpha}}$ = '+ '{:.3f}'.format(max2Ha)+' $10^{-14}$'))

        frame1.set_xticklabels([]) 			# Remove x-tic labels for the first frame
        plt.ylabel(r'Flux (A.U.)',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(prune='lower'))
        plt.text(0.84,0.9,'N$_{SII}$',color='darkgreen',transform=frame1.transAxes,fontsize=21)
        plt.text(0.84,0.94,'N$_{OI}$',color='goldenrod',transform=frame1.transAxes,fontsize=21)
        plt.text(0.87,0.92,'+',color='k',transform=frame1.transAxes,fontsize=21)
        plt.text(0.9,0.9,'S$_{SII}$',color='dodgerblue',transform=frame1.transAxes,fontsize=21)
        plt.text(0.9,0.94,'S$_{OI}$',color='darkred',transform=frame1.transAxes,fontsize=21)

	# RESIDUAL plot
        frame2 = fig1.add_axes((.1,.1,.85,.15))
        plt.plot(l,np.zeros(len(l)),c='orange',linestyle='-')         	# Line around zero
        plt.plot(l,data_cor-fin2_fit,color='k')		# Main
        plt.xlabel('Wavelength ($\AA$)',fontsize=19)
        plt.ylabel('Residuals',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.plot(l,np.zeros(len(l))+1.5*stadev,c='orange',linestyle=(0,(5,8)))	# 3 sigma upper limit
        plt.plot(l,np.zeros(len(l))-1.5*stadev,c='orange',linestyle=(0,(5,8))) 	# 3 sigma down limit
        plt.ylim(-(3*stadev)*3,(3*stadev)*3)

        plt.savefig(parentFold+'adj_full_2comp.pdf',format='pdf',bbox_inches='tight',pad_inches=0.2)
        props = dict(boxstyle='round',facecolor='white', alpha=0.5)
        frame1.text(6250.,max(data_cor),textstr,fontsize=12,verticalalignment='top', bbox=props)
        plt.savefig(parentFold+'adj_full_2comp.png',bbox_inches='tight',pad_inches=0.2)

        ########################################################################################################
        # TODo LISTO HASTA AQUI A 24 DE NOVIEMBRE EN TODOS LOS CODIGOS A LOS QUE SE LLAMAN.
        ########################################################################################################

        ##############################################################
        ##############################################################
        # Now we consider if a third narrow component is necessary:
        if stdf2_s2/stadev > 3 and stdf2_s1/stadev > 3:
            trigger3 = 'Y'
        else:
            trigger3 = 'N'

        trigger3 = 'N'
        # trigger3 = input('Do the fit needs a second narrow component? ("Y"/"N"): ')
        if trigger3 == 'N':
            print('No second narrow component is needed')
            fitoutput.prov_eSN = (esig2S2);fitoutput.prov_SN = (sig2S2);fitoutput.prov_VN = (v2S2);fitoutput.prov_eVN = (ev2S2);fitoutput.prov_notSN = (twocompresu.best_values['sig_3']);fitoutput.prov_noteSN = (twocompresu.params['sig_3'].stderr);fitoutput.prov_notVN = (twocompresu.params['mu_3']);fitoutput.prov_noteVN = (twocompresu.params['mu_3'].stderr);fitoutput.prov_eSS = (esig20S2);fitoutput.prov_SS = (sig20S2);fitoutput.prov_VS = (v20S2);fitoutput.prov_eVS = (ev20S2)
            fitoutput.prov_eSN_2 = (esig2O2);fitoutput.prov_SN_2 = (sig2O2);fitoutput.prov_VN_2 = (v2O2);fitoutput.prov_eVN_2 = (ev2O2);fitoutput.prov_eSS_2 = (esig20O2);fitoutput.prov_SS_2 = (sig20O2);fitoutput.prov_VS_2 = (v20O2);fitoutput.prov_eVS_2 = (ev20O2)
            fitoutput.prov_VN2 = (np.nan);fitoutput.prov_eVN2 = (np.nan);fitoutput.prov_SN2 = (np.nan);fitoutput.prov_eSN2 = (np.nan)
            fitoutput.prov_f2Ha = (np.nan);fitoutput.prov_f2N1 = (np.nan);fitoutput.prov_f2N2 = (np.nan);fitoutput.prov_f2S1 = (np.nan);fitoutput.prov_f2S2 = (np.nan);fitoutput.prov_f2O1 = (np.nan);fitoutput.prov_f2O2 = (np.nan)
            np.savetxt(parentFold+'fittwo_best_values.txt',np.c_[l,twocompresu.data,twocompresu.best_fit,lin_data_fin,tgaus1,tgaus8,tgaus2,tgaus9,tgaus3,tgaus10,tgaus4,tgaus11,tgaus5,tgaus12,tgaus6,tgaus13,tgaus7,tgaus14],fmt=('%5.6f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Wavelength\tReal_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tSecond_SII1\tNarrow_SII6716\tSecond_SII2\tNarrow_NII6584\tSecond_NII2\tNarrow_Halpha\tSecond_Halpha\tNarrow_NII6548\tSecond_NII1\tNarrow_OI6300\tSecond_OI6300\tNarrow_OI6363\tSecond_OI6363'))

        elif trigger3 == 'Y':
            print('A second narrow component is needed')
            Sgaus1,Sgaus2,Sgaus3,Sgaus4,Sgaus5,Sgaus6,Sgaus7,Sgaus11,Sgaus12,Sgaus13,Sgaus14,Sgaus15,Sgaus16,Sgaus17,Sgaus21,Sgaus22,Sgaus23,Sgaus24,Sgaus25,Sgaus26,Sgaus27,v30S2,ev30S2,v31S2,ev31S2,v32S2,ev32S2,sig30S2,esig30S2,sig31S2,esig31S2,sig30S2,esig30S2,T3CResu = FitThirdComp(parentFold,l,data_cor,lin_data_fin,threeresu,meth,mu0,sig0,amp0,mu1,sig1,amp1,mu2,sig2,amp2,mu3,sig3,amp3,mu4,sig4,amp4,mu5,sig5,amp5,mu6,sig6,amp6,sig20,amp20,sig21,amp21,sig22,amp22,sig23,amp23,sig24,amp24,sig25,amp25,sig26,amp26,sl,it,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,stadev)
            np.savetxt(parentFold+'fitthree_best_values.txt',np.c_[l,twocompresu.data,twocompresu.best_fit,lin_data_fin,Sgaus1,Sgaus11,Sgaus21,Sgaus2,Sgaus12,Sgaus22,Sgaus3,Sgaus13,Sgaus23,Sgaus4,Sgaus14,Sgaus24,Sgaus5,Sgaus15,Sgaus25,Sgaus6,Sgaus16,Sgaus26,Sgaus7,Sgaus17,Sgaus27],fmt=('%5.6f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Wavelength\tReal_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tSecond_SII1\tThird_SII1\tNarrow_SII6716\tSecond_SII2\tThird_SII2\tNarrow_NII6584\tSecond_NII2\tThird_NII2\tNarrow_Halpha\tSecond_Halpha\tThird_Halpha\tNarrow_NII6548\tSecond_NII1\tThird_NII1\tNarrow_OI6300\tSecond_OI6300\tThird_OI6300\tNarrow_OI6363\tSecond_OI6363\tThird_OI6363'))
            fitoutput.prov_fS2 = sum(Sgaus1);fitoutput.prov_fS1 = sum(Sgaus2);fitoutput.prov_fN2 = (sum(Sgaus3));fitoutput.prov_fHa = (sum(Sgaus4));fitoutput.prov_fN1 = (sum(Sgaus5));fitoutput.prov_fO1 = (sum(Sgaus6));fitoutput.prov_fO2 = (sum(Sgaus7));fitoutput.prov_fSHa = (sum(Sgaus14));fitoutput.prov_fSN1 = (sum(Sgaus15));fitoutput.prov_fSN2 = (sum(Sgaus13));fitoutput.prov_fSS1 = (sum(Sgaus12));fitoutput.prov_fSS2 = (sum(Sgaus11));fitoutput.prov_fSO1 = (sum(Sgaus16));fitoutput.prov_fSO2 = (sum(Sgaus17));fitoutput.prov_f2S2 = (sum(Sgaus21));fitoutput.prov_f2S1 = (sum(Sgaus22));fitoutput.prov_f2N2 = (sum(Sgaus23));fitoutput.prov_f2Ha = (sum(Sgaus24));fitoutput.prov_f2N1 = (sum(Sgaus25));fitoutput.prov_f2O1 = (sum(Sgaus26));fitoutput.prov_f2O2 = (sum(Sgaus27))
            fitoutput.prov_VN2 = (v32S2);fitoutput.prov_eVN2 = (ev32S2);fitoutput.prov_SN2 = (sig32S2);fitoutput.prov_eSN2 = (esig32S2);fitoutput.prov_VN = (v30S2);fitoutput.prov_eVN = (ev30S2);fitoutput.prov_SN = (sig30S2);fitoutput.prov_eSN = (esig30S2);fitoutput.prov_VS = (v31S2);fitoutput.prov_eVS = (ev31S2);fitoutput.prov_SS = (sig31S2);fitoutput.prov_eSS = (esig31S2);fitoutput.prov_notSN = (T3CResu.params['sig_3']);fitoutput.prov_noteSN = (T3CResu.params['sig_3'].stderr);fitoutput.prov_notVN = (T3CResu.params['mu_3']);fitoutput.prov_noteVN = (T3CResu.params['mu_3'].stderr)

	##############################################################
        trigger2 = 'N'
        if galaxy >= 15 and galaxy <= 20 and galaxy2 >= 12 and galaxy2 <= 16:
            if stdf2_n2/stadev > 3 or stdf2_ha/stadev > 3 or stdf2_n1/stadev > 3: trigger2 = 'Y'
            else: trigger2 = 'N'
        #if galaxy >= 15 and galaxy <= 20 and galaxy2 >= 12 and galaxy2 <= 16: trigger2 = 'Y' #33  30
	#trigger2 = input('Do the fit needs a broad Halpha component? ("Y"/"N"): ')
        #fitoutput.prov_eSN = (esig2S2),fitoutput.prov_SN = (sig2S2),fitoutput.prov_VN = (v2S2),fitoutput.prov_eVN = (ev2S2),fitoutput.prov_eSS = (esig20S2),fitoutput.prov_SS = (sig20S2),fitoutput.prov_VS = (v20S2),fitoutput.prov_eVS = (ev20S2),fitoutput.prov_notSN = (twocompresu.best_values['sig_3']),fitoutput.prov_noteSN = (twocompresu.params['sig_3'].stderr),prov_notVN.append(twocompresu.params['mu_3']),fitoutput.prov_noteVN = (twocompresu.params['mu_3'].stderr)

        if trigger2 == 'N': 
            print('The final plots are already printed and have been already saved!')
            fitoutput.prov_eSB = (np.nan)
            fitoutput.prov_SB = (np.nan)
            fitoutput.prov_VB = (np.nan)
            fitoutput.prov_eVB = (np.nan)
            if trigger3 == 'N':
                np.savetxt(parentFold+'fittwo_best_values.txt',np.c_[l,twocompresu.data,twocompresu.best_fit,lin_data_fin,tgaus1,tgaus8,tgaus2,tgaus9,tgaus3,tgaus10,tgaus4,tgaus11,tgaus5,tgaus12,tgaus6,tgaus13,tgaus7,tgaus14],fmt=('%5.6f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Wavelength\tReal_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tSecond_SII1\tNarrow_SII6716\tSecond_SII2\tNarrow_NII6584\tSecond_NII2\tNarrow_Halpha\tSecond_Halpha\tNarrow_NII6548\tSecond_NII1\tNarrow_OI6300\tSecond_OI6300\tNarrow_OI6363\tSecond_OI6363'))
                # Saving the max flux of each line to do the final flux map
                fitoutput.prov_fS2 = (sum(tgaus1));fitoutput.prov_fS1 = (sum(tgaus2));fitoutput.prov_fN2 = (sum(tgaus3));fitoutput.prov_fHa = (sum(tgaus4));fitoutput.prov_fN1 = (sum(tgaus5));fitoutput.prov_fO1 = (sum(tgaus6));fitoutput.prov_fO2 = (sum(tgaus7));fitoutput.prov_fSHa = (sum(tgaus11));fitoutput.prov_fSN1 = (sum(tgaus12));fitoutput.prov_fSN2 = (sum(tgaus10));fitoutput.prov_fSS1 = (sum(tgaus9));fitoutput.prov_fSS2 = (sum(tgaus8));fitoutput.prov_fSO1 = (sum(tgaus13));fitoutput.prov_fSO2 = (sum(tgaus14))


        elif trigger2 == 'Y':
            # Now we define the initial guesses and the constraints
            newxb = l[np.where(l>l9)[0][0]:np.where(l<l6)[0][-1]]
            newyb = data_cor[np.where(l>l9)[0][0]:np.where(l<l6)[0][-1]]
            sigb = 16.
            mub  = mu3
            ampb = amp3/3.
            paramsbH = lmfit.Parameters()
	        # broad components
            ab = lmfit.Parameter('mu_b',value=mub)#6605.,vary=False)
            bc = lmfit.Parameter('sig_b',value=sigb,min=12.)#twocompresu.values['sig_23'])#29.51,vary=False) minbroad
            yz = lmfit.Parameter('amp_b',value=ampb,min=0.)
            #aramsbH.add_many(sl,it,ab,bc,yz,cd1,de,ef,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr,aaa,aab,aac,aad,aae,aaf,aag,aah,aai,aaj,aak,aal,aam,aan,aao)
            paramsbH.add_many(sl,it,ab,bc,yz,cd1,de,ef,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr,rs,st,tu,uv,vw,wy,aaa,aab,aac,aad,aae,aaf,aag,aah,aai,aaj,aak,aal,aam,aan,aao,aap,aaq,aar,aas,aat,aau)

            twobroadresu = twobroadcomp_mod.fit(data_cor,paramsbH,x=l)
            lmfit.model.save_modelresult(twobroadresu, parentFold+'broadtwo_modelresult.sav')
            with open(parentFold+'fit_twobroad_result.txt', 'w') as fh:
                fh.write(twobroadresu.fit_report())

	    # PLOT AND PRINT THE RESULTS 
            refer2 = broad_plot(parentFold,l,data_cor,meth,trigger,linresu,tworesu,fin2_fit,twobroadresu,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,std0,std1)
            twobroad_fit = twobroadresu.best_fit
            # Save the results in the corresponding array
            fitoutput.prov_VB = (refer2[0])
            fitoutput.prov_eVB = (refer2[1])
            fitoutput.prov_SB = (refer2[2])
            fitoutput.prov_eSB = (refer2[3])
            fitoutput.prov_fS2 = (sum(refer2[4]));fitoutput.prov_fS1 = (sum(refer2[5]));fitoutput.prov_fN2 = (sum(refer2[6]));fitoutput.prov_fHa = (sum(refer2[7]));fitoutput.prov_fN1 = (sum(refer2[8]));fitoutput.prov_fO1 = (sum(refer2[9]));fitoutput.prov_fO2 = (sum(refer2[10]));fitoutput.prov_fSHa = (sum(refer2[14]));fitoutput.prov_fSN1 = (sum(refer2[15]));fitoutput.prov_fSN2 = (sum(refer2[13]));fitoutput.prov_fSS1 = (sum(refer2[12]));fitoutput.prov_fSS2 = (sum(refer2[11]));fitoutput.prov_fSO1 = (sum(refer2[16]));fitoutput.prov_fSO2 = (sum(refer2[17]))

    else: 
        print('Please use "Y" or "N"')

    
    #### Save
    outname = 'fit_%i_%i.p'%(x,y)
    fout = open(parentFold+outname,'w')
    pickle.dump(fitoutput, fout)
