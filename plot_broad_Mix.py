import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import Ofuncts
import os

def broad_plot(parentFold,l,data_cor,meth,trigger,linresu,refresu,fullresu,broadresu,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,std0,std1):
    '''
	It gives the plots for one and two components + a broad component in the whole spectra

	The parameters needed are:
	parentFold:          Path to the data
	l:             Wavelength range
	data_cor:      Flux for each wavelength
	meth:          Method to be applied (S/O)
	trigger:       This had to be said to the program to decide whether 1 or 2 components
	linresu:       Result of the linear fit of the spectra
	refresu:       Result of the linear+gaussian fit for the reference lines with one component or two components
	fullresu:      Result of the linear+gaussian fit for the spectra with one component or two components
	broadresu:     Result of the linear+gaussian+broad Ha fit for the spectra with one or two components
	l1-l14:        Parts of the spectra where the lines are located
	std0/std1:     Where the standard deviation of the continuum is calculated
    '''
    # Rest values of the line wavelengths 
    l_OI_1  = 6300.30
    l_OI_2  = 6363.30
    l_Halpha = 6562.801
    l_NII_1  = 6548.05
    l_NII_2  = 6583.45
    l_SII_1  = 6716.44
    l_SII_2  = 6730.82
	
    # Constants and STIS parameters
    v_luz = 299792.458 # km/s
    plate_scale = 0.58
    fwhm = 2*np.sqrt(2*np.log(2)) # times sigma
    pix_to_v = 50.
    sig_inst = 2.2

    new_slop = linresu.values['slope']
    new_intc = linresu.values['intc']
    lin_data_fin = (linresu.values['slope']*l+linresu.values['intc'])
    stadev = np.std(data_cor[std0:std1])

    if trigger=='Y':
        ################################# Calculate gaussians and final fit ###################################
        # Now we create and plot the individual gaussians of the fit
        bgaus1 = Ofuncts.gaussian(l,broadresu.values['mu_0'],broadresu.values['sig_0'],broadresu.values['amp_0']) 
        bgaus2 = Ofuncts.gaussian(l,broadresu.values['mu_1'],broadresu.values['sig_1'],broadresu.values['amp_1'])
        bgaus3 = Ofuncts.gaussian(l,broadresu.values['mu_2'],broadresu.values['sig_2'],broadresu.values['amp_2']) 
        bgaus4 = Ofuncts.gaussian(l,broadresu.values['mu_3'],broadresu.values['sig_3'],broadresu.values['amp_3'])
        bgaus5 = Ofuncts.gaussian(l,broadresu.values['mu_4'],broadresu.values['sig_4'],broadresu.values['amp_4'])
        bgaus6 = Ofuncts.gaussian(l,broadresu.values['mu_5'],broadresu.values['sig_5'],broadresu.values['amp_5'])
        bgaus7 = Ofuncts.gaussian(l,broadresu.values['mu_6'],broadresu.values['sig_6'],broadresu.values['amp_6'])
        bgaus8 = Ofuncts.gaussian(l,broadresu.values['mu_b'],broadresu.values['sig_b'],broadresu.values['amp_b'])
        broad_fit = broadresu.best_fit

        # We have to calculate the contribution of each component to the global fit
        # Lets define the linear fit data to add to each individual gaussian
        bgaus_total = broad_fit - lin_data_fin
        np.savetxt(parentFold+'fitbroad_best_values.txt',np.c_[broadresu.data,broadresu.best_fit,lin_data_fin,bgaus1,bgaus2,bgaus3,bgaus4,bgaus5,bgaus6,bgaus7,bgaus8],fmt=('%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Real_data\tBest_fit\tLineal_fit\tNarrow_SII6741\tNarrow_SII6716\tNarrow_NII6584\tNarrow_Halpha\tNarrow_NII6548\tNarrow_OI6300\tNarrow_OI6363\tBroad_Halpha'))
        # Now lets determine the contribution of the individual components as follows:
        contr_HaN = sum(bgaus4)
        contr_HaB = sum(bgaus8)
        ix_Br_sup = np.where(bgaus8 > 10**-5)[0][-1]
        ix_Br_inf = np.where(bgaus8 > 10**-5)[0][0]
        contr_NII2N = sum(bgaus3)
        contr_NII1N = sum(bgaus5)
        total_flux_NII_Halp = sum(bgaus_total[ix_Br_inf:ix_Br_sup])

        contr_HaBtoNHa = 100*(contr_HaB/total_flux_NII_Halp)
        contr_HaNtoNHa = 100*(contr_HaN/total_flux_NII_Halp)
        contr_NII2NtoNHa = 100*(contr_NII2N/total_flux_NII_Halp)
        contr_NII1NtoNHa = 100*(contr_NII1N/total_flux_NII_Halp)

        print('The contribution of the broad component to the total Halpha+N flux is: '+'{:.2f}'.format(contr_HaBtoNHa)+'%')
        np.savetxt(parentFold+'1cBroad_N+Ha_indivcontr.txt',np.c_[contr_NII1NtoNHa,contr_HaNtoNHa,contr_HaBtoNHa,contr_NII2NtoNHa],fmt=('%10.7f','%10.7f','%10.7f','%10.7f'),header=('Narrow_NII1(%)\tNarrow_Halpha(%)\tBroad_Halpha(%)\tNarrow_NII2(%)'))

        # Now we calculate the epsilon values derived from the fit
        stdb_s2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10]-broad_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10])
        stdb_s1 = np.std(data_cor[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]]-broad_fit[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]])
        stdb_n2 = np.std(data_cor[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10]-broad_fit[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10])
        stdb_ha = np.std(data_cor[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]]-broad_fit[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]])
        stdb_n1 = np.std(data_cor[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]]-broad_fit[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]])
        stdb_o1 = np.std(data_cor[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]]-broad_fit[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]])
        stdb_o2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-broad_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 component + Ha is... ')
        print('		For SII2: '+str(stdb_s2/stadev)+' < 3')
        print('		For SII1: '+str(stdb_s1/stadev)+' < 3')
        print('		For NII2: '+str(stdb_n2/stadev)+' < 3')
        print('		For Halp: '+str(stdb_ha/stadev)+' < 3')
        print('		For NII1: '+str(stdb_n1/stadev)+' < 3')
        print('		For OI1 : '+str(stdb_o1/stadev)+' < 3')
        print('		For OI2 : '+str(stdb_o2/stadev)+' < 3')
    	    
        if os.path.exists(parentFold+'eps_adj'+str(meth)+'_1Cb.txt'): os.remove(parentFold+'eps_adj'+str(meth)+'_1Cb.txt')
        np.savetxt(parentFold+'eps_adj_1Cb.txt',np.c_[stdb_s2/stadev,stdb_s1/stadev,stdb_n2/stadev,stdb_ha/stadev,stdb_n1/stadev,stdb_o1/stadev,stdb_o2/stadev,broadresu.chisqr],('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('SII2\tSII1\tNII2\tHa\tNII1\tOI1\tOI2\tChi2'))

   	    # We determine the maximum flux of the fit for all the lines, and the velocity and sigma components
        try:
            maxbS1 = broad_fit[np.where(abs(broadresu.values['mu_0']-l)<0.5)[0][0]] #max(broad_fit[np.where(l>l3)[0][0]:np.where(l<l4)[0][-1]])
            maxbS2 = broad_fit[np.where(abs(broadresu.values['mu_1']-l)<0.5)[0][0]] #max(broad_fit[np.where(l>l1)[0][0]:np.where(l<l2)[0][-1]])
            maxbN1 = broad_fit[np.where(abs(broadresu.values['mu_2']-l)<0.5)[0][0]] #max(broad_fit[np.where(l>l9)[0][0]:np.where(l<l10)[0][-1]])
            maxbHa = broad_fit[np.where(abs(broadresu.values['mu_3']-l)<0.5)[0][0]] #max(broad_fit[np.where(l>l7)[0][0]:np.where(l<l8)[0][-1]])
            maxbN2 = broad_fit[np.where(abs(broadresu.values['mu_4']-l)<0.5)[0][0]] #max(broad_fit[np.where(l>l5)[0][0]:np.where(l<l6)[0][-1]])
            maxbO1 = broad_fit[np.where(abs(broadresu.values['mu_5']-l)<0.5)[0][0]] #max(broad_fit[np.where(l>l5)[0][0]:np.where(l<l6)[0][-1]])
            maxbO2 = broad_fit[np.where(abs(broadresu.values['mu_6']-l)<0.5)[0][0]] #max(broad_fit[np.where(l>l5)[0][0]:np.where(l<l6)[0][-1]])
        except IndexError:
            print('ERROR: index out of range. Setting the flux values of the OI 1 line to 0.')

	# one component + Halpha
        sigS2 = pix_to_v*np.sqrt(broadresu.values['sig_0']**2-sig_inst**2)
        sigO2 = pix_to_v*np.sqrt(broadresu.values['sig_5']**2-sig_inst**2)
        sigbS2 = pix_to_v*np.sqrt(broadresu.values['sig_b']**2-sig_inst**2)
        if refresu.params['sig_0'].stderr == None:
            esigS2 = 0.
        else: 
            esigS2 = pix_to_v*(2*broadresu.values['sig_0']*refresu.params['sig_0'].stderr)/(np.sqrt(broadresu.values['sig_0']**2-sig_inst**2))
        if refresu.params['sig_5'].stderr == None:
            esigO2 = 0.
        else: 
            esigO2 = pix_to_v*(2*broadresu.values['sig_5']*refresu.params['sig_5'].stderr)/(np.sqrt(broadresu.values['sig_5']**2-sig_inst**2))
        if broadresu.params['sig_b'].stderr == None:
            esigbS2 = 0.
        else: 
            esigbS2 = pix_to_v*(2*broadresu.values['sig_b']*broadresu.params['sig_b'].stderr)/(np.sqrt(broadresu.values['sig_b']**2-sig_inst**2))

        vS2 = v_luz*((broadresu.values['mu_0']-l_SII_2)/l_SII_2)
        vO2 = v_luz*((broadresu.values['mu_5']-l_OI_1)/l_OI_1)
        vbS2 = v_luz*((broadresu.values['mu_b']-l_Halpha)/l_Halpha)
        if refresu.params['mu_0'].stderr == None: 
            print('Problem determining the errors! First component ')
            evS2 = 0.
        elif refresu.params['mu_0'].stderr != None: 
            evS2 = ((v_luz/l_SII_2)*refresu.params['mu_0'].stderr)
        if refresu.params['mu_5'].stderr == None: 
            print('Problem determining the errors! First component ')
            evO2 = 0.
        elif refresu.params['mu_5'].stderr != None: 
            evO2 = ((v_luz/l_OI_1)*refresu.params['mu_5'].stderr)
        if broadresu.params['mu_b'].stderr == None:
            evbS2 = 0.
        else:
            evbS2 = ((v_luz/l_Halpha)*broadresu.params['mu_b'].stderr)
        textstr = '\n'.join((r'$V_{SII_{2}}$ = '+ '{:.2f} +- {:.2f}'.format(vS2,evS2),
		r'$V_{OI_{1}}$ = '+ '{:.2f} +- {:.2f}'.format(vO2,evO2),
		r'$V_{H_{\alpha-broad}}$ = '+ '{:.2f} +- {:.2f}'.format(vbS2,evbS2),
	    	r'$\sigma_{SII_{2}}$ = '+ '{:.2f} +- {:.2f}'.format(sigS2,esigS2),
	    	r'$\sigma_{OI_{1}}$ = '+ '{:.2f} +- {:.2f}'.format(sigO2,esigO2),
	    	r'$\sigma_{H_{\alpha-broad}}$ = '+ '{:.2f} +- {:.2f}'.format(sigbS2,esigbS2),
		r'$F_{SII_{2}}/F_{SII_{1}}$ = '+ '{:.3f}'.format(maxbS2/maxbS1),
	    	r'$F_{H_{\alpha}}$ = '+ '{:.3f}'.format(maxbHa)+' $10^{-14}$'))


        # Save the velocity and sigma for all components
        if os.path.exists(parentFold+'v_sig_adj_1Cb.txt'): os.remove(parentFold+'v_sig_adj_1Cb.txt')
        np.savetxt(parentFold+'v_sig_adj_1Cb.txt',
                       np.c_[vS2,evS2,vO2,evO2,vbS2,evbS2,sigS2,esigS2,sigO2,esigO2,sigbS2,esigbS2], 
                       ('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),
                       header=('v_refS\tev_refS\tv_refO\tev_refO\tv_broadH\tev_broadH\tsig_refS\tesig_refS\tsig_refO\tesig_refO\tsig_broadH\tesig_broadH'))
        if os.path.exists(parentFold+'fluxes_'+str(meth)+'_1Cb.txt'): os.remove(parentFold+'fluxes_'+str(meth)+'_1Cb.txt')
        np.savetxt(parentFold+'fluxes_'+str(meth)+'_1Cb.txt', 
                       np.c_[(sum(bgaus1)+sum(bgaus2))/sum(bgaus4),sum(bgaus3)/sum(bgaus4),sum(bgaus4),sum(bgaus5)/sum(bgaus4)],
                       ('%8.6f','%8.6f','%8.15f','%8.6f'),
                       header=('SII_6731+6716\tNII_6584\tHalpha\tNII_6548'))

   	########################### PLOT #############################
        plt.close('all')
        # MAIN plot
        fig1   = plt.figure(1,figsize=(10, 9))
        frame1 = fig1.add_axes((.1,.25,.85,.65)) 	     # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
        plt.plot(l,data_cor,'k',linewidth=2)	     # Initial data
        plt.plot(l[std0:std1],data_cor[std0:std1],'y',linewidth=4.)	# Zone where the stddev is calculated
        plt.plot(l[std0:std1],data_cor[std0:std1],'k',linewidth=1.)	     # Initial data
#    	    plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='gold',linestyle=(0, (5, 8)),label='Linear fit')
        plt.plot(l,bgaus1+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,bgaus2+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,bgaus6+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,bgaus7+lin_data_fin,c='goldenrod',linestyle='-')
        if meth == 'M1':
            plt.plot(l,bgaus3+lin_data_fin,c='darkgreen',linestyle='-')
            plt.plot(l,bgaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,bgaus5+lin_data_fin,c='darkgreen',linestyle='-')
        if meth == 'M2':
            plt.plot(l,bgaus3+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,bgaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,bgaus5+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,bgaus8+lin_data_fin,c='darkviolet',linestyle='-')#,label='Broad component')
        plt.plot(l,broad_fit,'r-')

        frame1.set_xticklabels([]) 			# Remove x-tic labels for the first frame
        plt.ylabel('Flux (erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(prune='lower'))
        plt.text(0.84,0.94,'N$_{S}$',color='darkgreen',transform=frame1.transAxes,fontsize=21)
        plt.text(0.84,0.9,'N$_{O}',color='goldenrod',transform=frame1.transAxes,fontsize=21)
        plt.text(0.87,0.92,'+',color='k',transform=frame1.transAxes,fontsize=21)
        plt.text(0.9,0.92,'B',color='darkviolet',transform=frame1.transAxes,fontsize=21)

        # RESIDUAL plot
        frame2 = fig1.add_axes((.1,.1,.85,.15))
        plt.plot(l,np.zeros(len(l)),c='orange',linestyle='-')         	# Line around zero
        plt.plot(l,data_cor-broad_fit,c='k')		# Main
        plt.xlabel('Wavelength ($\AA$)',fontsize=19)
        plt.ylabel('Residuals',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.plot(l,np.zeros(len(l))+1.5*stadev,c='orange',linestyle=(0,(5,8)))	# 3 sigma upper limit
        plt.plot(l,np.zeros(len(l))-1.5*stadev,c='orange',linestyle=(0,(5,8))) 	# 3 sigma down limit
        plt.ylim(-(3*stadev)*3,(3*stadev)*3)

        plt.savefig(parentFold+'adj_met'+str(meth)+'_full_1comp_broadH.pdf',format='pdf',bbox_inches='tight',pad_inches=0.2)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        frame1.text(2.2,max(data_cor)+0.05, textstr, fontsize=12,verticalalignment='top', bbox=props)
        plt.savefig(parentFold+'adj_met'+str(meth)+'_full_1comp_broadH.png',bbox_inches='tight',pad_inches=0.2)
        
        return [vbS2,evbS2,sigbS2,esigbS2,bgaus1,bgaus2,bgaus3,bgaus4,bgaus5,bgaus6,bgaus7,bgaus8]
        
############################################################################################################################################
############################################################################################################################################

    elif trigger=='N':
        ############# Calculate gaussians and final fit #################
        # Now we create and plot the individual gaussians of the fit
        b2gaus1 = Ofuncts.gaussian(l,broadresu.values['mu_0'],broadresu.values['sig_0'],broadresu.values['amp_0']) 
        b2gaus2 = Ofuncts.gaussian(l,broadresu.values['mu_1'],broadresu.values['sig_1'],broadresu.values['amp_1'])
        b2gaus3 = Ofuncts.gaussian(l,broadresu.values['mu_2'],broadresu.values['sig_2'],broadresu.values['amp_2']) 
        b2gaus4 = Ofuncts.gaussian(l,broadresu.values['mu_3'],broadresu.values['sig_3'],broadresu.values['amp_3'])
        b2gaus5 = Ofuncts.gaussian(l,broadresu.values['mu_4'],broadresu.values['sig_4'],broadresu.values['amp_4'])
        b2gaus6 = Ofuncts.gaussian(l,broadresu.values['mu_5'],broadresu.values['sig_5'],broadresu.values['amp_5'])
        b2gaus7 = Ofuncts.gaussian(l,broadresu.values['mu_6'],broadresu.values['sig_6'],broadresu.values['amp_6'])
        b2gaus8 = Ofuncts.gaussian(l,broadresu.values['mu_20'],broadresu.values['sig_20'],broadresu.values['amp_20']) 
        b2gaus9 = Ofuncts.gaussian(l,broadresu.values['mu_21'],broadresu.values['sig_21'],broadresu.values['amp_21'])
        b2gaus10 = Ofuncts.gaussian(l,broadresu.values['mu_22'],broadresu.values['sig_22'],broadresu.values['amp_22']) 
        b2gaus11 = Ofuncts.gaussian(l,broadresu.values['mu_23'],broadresu.values['sig_23'],broadresu.values['amp_23'])
        b2gaus12 = Ofuncts.gaussian(l,broadresu.values['mu_24'],broadresu.values['sig_24'],broadresu.values['amp_24'])
        b2gaus13 = Ofuncts.gaussian(l,broadresu.values['mu_25'],broadresu.values['sig_25'],broadresu.values['amp_25'])
        b2gaus14 = Ofuncts.gaussian(l,broadresu.values['mu_26'],broadresu.values['sig_26'],broadresu.values['amp_26'])
        b2gausb = Ofuncts.gaussian(l,broadresu.values['mu_b'],broadresu.values['sig_b'],broadresu.values['amp_b'])
        twobroad_fit = broadresu.best_fit

        # We have to calculate the contribution of each component to the global fit
        # Lets define the linear fit data to add to each individual gaussian
        b2gaus_total = twobroad_fit - lin_data_fin
        np.savetxt(parentFold+'fitbroadtwo_best_values.txt',np.c_[broadresu.data,broadresu.best_fit,lin_data_fin,b2gaus1,b2gaus2,b2gaus3,b2gaus4,b2gaus5,b2gaus6,b2gaus7,b2gaus8,b2gaus9,b2gaus10,b2gaus11,b2gaus12,b2gaus13,b2gaus14,b2gausb],fmt=('%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Real_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tNarrow_SII6716\tNarrow_NII6584\tNarrow_Halpha\tNarrow_NII6548\tNarrow_OI6300\tNarrow_OI6363\tSecond_SII1\tSecond_SII2\tSecond_NII2\tSecond_Halpha\tSecond_NII1\tSecond_OI1\tSecond_OI2\tBroad_Halpha'))
        # Now lets determine the contribution of the individual components as follows:
        contr_HaN = sum(b2gaus4)
        contr_HaB = sum(b2gausb)
        contr_HaS = sum(b2gaus11)
        ix_Br_sup = np.where(b2gausb > 10**-5)[0][-1]
        ix_Br_inf = np.where(b2gausb > 10**-5)[0][0]
        contr_NII2N = sum(b2gaus3)
        contr_NII1N = sum(b2gaus5)
        contr_NII2S = sum(b2gaus10)
        contr_NII1S = sum(b2gaus12)
        total_flux_NII_Halp = sum(b2gaus_total[ix_Br_inf:ix_Br_sup])

        contr_HaBtoNHa = 100*(contr_HaB/total_flux_NII_Halp)
        contr_HaNtoNHa = 100*(contr_HaN/total_flux_NII_Halp)
        contr_HaStoNHa = 100*(contr_HaS/total_flux_NII_Halp)
        contr_NII2NtoNHa = 100*(contr_NII2N/total_flux_NII_Halp)
        contr_NII2StoNHa = 100*(contr_NII2S/total_flux_NII_Halp)
        contr_NII1NtoNHa = 100*(contr_NII1N/total_flux_NII_Halp)
        contr_NII1StoNHa = 100*(contr_NII1S/total_flux_NII_Halp)

        print('The contribution of the broad component to the total Halpha+N flux is: '+'{:.2f}'.format(contr_HaBtoNHa)+'%')
        np.savetxt(parentFold+'2cBroad_N+Ha_indivcontr.txt',np.c_[contr_NII1NtoNHa,contr_NII1StoNHa,contr_HaNtoNHa,contr_HaStoNHa,contr_HaBtoNHa,contr_NII2NtoNHa,contr_NII2StoNHa],fmt=('%10.7f','%10.7f','%10.7f','%10.7f','%10.7f','%10.7f','%10.7f'),header=('Narrow_NII1(%)\tSecond_NII1(%)\tNarrow_Halpha(%)\tSecond_Halpha(%)\tBroad_Halpha(%)\tNarrow_NII2(%)\tSecond_NII2(%)'))

        # Now we calculate the epsilon under the lines
        stdb2_s2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10]-twobroad_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10])
        stdb2_s1 = np.std(data_cor[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]]-twobroad_fit[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]])
        stdb2_n2 = np.std(data_cor[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10]-twobroad_fit[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10])
        stdb2_ha = np.std(data_cor[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]]-twobroad_fit[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]])
        stdb2_n1 = np.std(data_cor[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]]-twobroad_fit[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]])
        stdb2_o1 = np.std(data_cor[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]]-twobroad_fit[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]])
        stdb2_o2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-twobroad_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])

        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 component + Ha is... ')
        print('		For SII2: '+str(stdb2_s2/stadev)+' < 3')
        print('		For SII1: '+str(stdb2_s1/stadev)+' < 3')
        print('		For NII2: '+str(stdb2_n2/stadev)+' < 3')
        print('		For Halp: '+str(stdb2_ha/stadev)+' < 3')
        print('		For NII1: '+str(stdb2_n1/stadev)+' < 3')
        print('		For OI1 : '+str(stdb2_o1/stadev)+' < 3')
        print('		For OI2 : '+str(stdb2_o2/stadev)+' < 3')
	    
        if os.path.exists(parentFold+'eps_adj_2Cb.txt'): os.remove(parentFold+'eps_adj_2Cb.txt')
        np.savetxt(parentFold+'eps_adj_2Cb.txt',np.c_[stdb2_s2/stadev,stdb2_s1/stadev,stdb2_n2/stadev,stdb2_ha/stadev,stdb2_n1/stadev,stdb2_o1/stadev,stdb2_o2/stadev,broadresu.chisqr],('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('SII2\tSII1\tNII2\tHa\tNII1\tOI1\tOI2\tChi2'))

        # We determine the maximum flux of the fit for all the lines, and the velocity and sigma components
        try:
            maxfbS1 = twobroad_fit[np.where(abs(broadresu.values['mu_0']-l)<0.5)[0][0]]
            maxfbS2 = twobroad_fit[np.where(abs(broadresu.values['mu_1']-l)<0.5)[0][0]] 
            maxfbN1 = twobroad_fit[np.where(abs(broadresu.values['mu_2']-l)<0.5)[0][0]] 
            maxfbHa = twobroad_fit[np.where(abs(broadresu.values['mu_3']-l)<0.5)[0][0]] 
            maxfbN2 = twobroad_fit[np.where(abs(broadresu.values['mu_4']-l)<0.5)[0][0]] 
            maxfbO1 = twobroad_fit[np.where(abs(broadresu.values['mu_5']-l)<0.5)[0][0]] 
            maxfbO2 = twobroad_fit[np.where(abs(broadresu.values['mu_6']-l)<0.5)[0][0]] 
        except IndexError:
            if broadresu.values['mu_0']>l[-1]:
                print('ERROR: index out of range. Setting the flux values of the SI 1 line to 0.')
                maxfbS1 = 0.

	# two comps + Halpha
        sigS2 = pix_to_v*np.sqrt(broadresu.values['sig_0']**2-sig_inst**2) 
        sigO2 = pix_to_v*np.sqrt(broadresu.values['sig_5']**2-sig_inst**2) 
        sig2S2 = pix_to_v*np.sqrt(broadresu.values['sig_20']**2-sig_inst**2)
        sig2O2 = pix_to_v*np.sqrt(broadresu.values['sig_25']**2-sig_inst**2)
        sigbS2 = pix_to_v*np.sqrt(broadresu.values['sig_b']**2-sig_inst**2)
        if refresu.params['sig_0'].stderr == None and refresu.params['sig_20'].stderr == None:
            esigS2,esig2S2 = 0.,0.
        else: 
            esigS2 = pix_to_v*(2*broadresu.values['sig_0']*refresu.params['sig_0'].stderr)/(np.sqrt(broadresu.values['sig_0']**2-sig_inst**2))
            esig2S2 = pix_to_v*(2*broadresu.values['sig_20']*refresu.params['sig_20'].stderr)/(np.sqrt(broadresu.values['sig_20']**2-sig_inst**2))
        if refresu.params['sig_5'].stderr == None and refresu.params['sig_25'].stderr == None:
            esigO2,esig2O2 = 0.,0.
        else: 
            esigO2 = pix_to_v*(2*broadresu.values['sig_5']*refresu.params['sig_5'].stderr)/(np.sqrt(broadresu.values['sig_5']**2-sig_inst**2))
            esig2O2 = pix_to_v*(2*broadresu.values['sig_25']*refresu.params['sig_25'].stderr)/(np.sqrt(broadresu.values['sig_25']**2-sig_inst**2))
        if broadresu.params['sig_b'].stderr == None:
            esigbS2 = 0.
        else: 
            esigbS2 = pix_to_v*(2*broadresu.values['sig_b']*broadresu.params['sig_b'].stderr)/(np.sqrt(broadresu.values['sig_b']**2-sig_inst**2))

        # Calculate velocity and sigma
        vS2 = v_luz*((broadresu.values['mu_0']-l_SII_2)/l_SII_2)
        vO2 = v_luz*((broadresu.values['mu_5']-l_OI_1)/l_OI_1)
        v2S2 = v_luz*((broadresu.values['mu_20']-l_SII_2)/l_SII_2)
        v2O2 = v_luz*((broadresu.values['mu_25']-l_OI_1)/l_OI_1)
        vbS2 = v_luz*((broadresu.values['mu_b']-l_Halpha)/l_Halpha)
        if refresu.params['mu_0'].stderr == None: 
            print('Problem determining the errors! First component ')
            evS2,ev2S2 = 0.,0.
        elif refresu.params['mu_0'].stderr != None: 
            evS2 = ((v_luz/l_SII_2)*refresu.params['mu_0'].stderr)
            ev2S2 = ((v_luz/l_SII_2)*refresu.params['mu_20'].stderr)
        if refresu.params['mu_5'].stderr == None: 
            print('Problem determining the errors! First component ')
            evO2,ev2O2 = 0.,0.
        elif refresu.params['mu_5'].stderr != None: 
            evO2 = ((v_luz/l_OI_1)*refresu.params['mu_5'].stderr)
            ev2O2 = ((v_luz/l_OI_1)*refresu.params['mu_25'].stderr)
        if broadresu.params['mu_b'].stderr == None:
            evbS2 = 0.
        else:
            evbS2 = ((v_luz/l_Halpha)*broadresu.params['mu_b'].stderr)
		
        textstr = '\n'.join((r'$V_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(vS2,evS2),
			r'$V_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v2S2,ev2S2),
			r'$V_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(vO2,evO2),
			r'$V_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v2O2,ev2O2),
			r'$V_{H_{\alpha-broad}}$ = '+ '{:.2f} +- {:.2f}'.format(vbS2,evbS2),
		    	r'$\sigma_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sigS2,esigS2),
		    	r'$\sigma_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig2S2,esig2S2),
		    	r'$\sigma_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sigO2,esigO2),
		    	r'$\sigma_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig2O2,esig2O2),
		    	r'$\sigma_{H_{\alpha-broad}}$ = '+ '{:.2f} +- {:.2f}'.format(sigbS2,esigbS2),
			r'$F_{SII_{2}}/F_{SII_{1}}$ = '+ '{:.3f}'.format(maxfbS2/maxfbS1),
		    	r'$F_{H_{\alpha}}$ = '+ '{:.3f}'.format(maxfbHa)+' $10^{-14}$'))
	    
        # Save the velocity and sigma for all components
        if os.path.exists(parentFold+'v_sig_adj_2Cb.txt'): os.remove(parentFold+'v_sig_adj_2Cb.txt')
        np.savetxt(parentFold+'v_sig_adj_2Cb.txt',
                       np.c_[vS2,evS2,v2S2,ev2S2,vO2,evO2,v2O2,ev2O2,vbS2,evbS2,sigS2,esigS2,sig2S2,esig2S2,sigO2,esigO2,sig2O2,esig2O2,sigbS2,esigbS2],
                       ('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),
                       header=('v_ref2S\tev_ref2S\tv_2ref2S\tev_2ref2S\tv_ref2O\tev_ref2O\tv_2ref2O\tev_2ref2O\tv_broadH\tev_broadH\tsig_ref2S\tesig_ref2S\tsig_2ref2S\tesig_2ref2S\tsig_ref2O\tesig_ref2O\tsig_2ref2O\tesig_2ref2O\tsig_broadH\tesig_broadH'))

        if os.path.exists(parentFold+'fluxes_2Cb_Ncomp.txt'): os.remove(parentFold+'fluxes_2Cb_Ncomp.txt')
        np.savetxt(parentFold+'fluxes_2Cb_Ncomp.txt',np.c_[(b2gaus1,b2gaus2,b2gaus3,b2gaus4,b2gaus5,b2gaus6,b2gaus7)],('%8.6f','%8.6f','%8.6f','%8.6f','%8.6f','%8.15f','%8.6f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))
        if os.path.exists(parentFold+'fluxes_2Cb_Scomp.txt'): os.remove(parentFold+'fluxes_2Cb_Scomp.txt')
        np.savetxt(parentFold+'fluxes_2Cb_Scomp.txt',np.c_[(b2gaus8,b2gaus9,b2gaus10,b2gaus11,b2gaus12,b2gaus13,b2gaus14)],('%8.6f','%8.6f','%8.6f','%8.6f','%8.6f','%8.15f','%8.6f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))


        ########################### PLOT #############################
        plt.close('all')
        # MAIN plot
        fig1   = plt.figure(1,figsize=(10, 9))
        frame1 = fig1.add_axes((.1,.25,.85,.65)) 	     # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
        plt.plot(l,data_cor,'k',linewidth=2)	     # Initial data
        plt.plot(l[std0:std1],data_cor[std0:std1],'y',linewidth=4.)	# Zone where the stddev is calculated
        plt.plot(l[std0:std1],data_cor[std0:std1],'k',linewidth=1.)	     # Initial data
#	    plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='gold',linestyle=(0, (5, 8)),label='Linear fit')
        plt.plot(l,b2gaus1+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,b2gaus2+lin_data_fin,c='darkgreen',linestyle='-')
        plt.plot(l,b2gaus6+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,b2gaus7+lin_data_fin,c='goldenrod',linestyle='-')
        plt.plot(l,b2gaus8+lin_data_fin,c='dodgerblue',linestyle='-')
        plt.plot(l,b2gaus9+lin_data_fin,c='dodgerblue',linestyle='-')
        plt.plot(l,b2gaus13+lin_data_fin,c='darkred',linestyle='-')
        plt.plot(l,b2gaus14+lin_data_fin,c='darkred',linestyle='-')
        if meth == 'M1':
            plt.plot(l,b2gaus3+lin_data_fin,c='darkgreen',linestyle='-')
            plt.plot(l,b2gaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,b2gaus5+lin_data_fin,c='darkgreen',linestyle='-')
            plt.plot(l,b2gaus10+lin_data_fin,c='dodgerblue',linestyle='-')
            plt.plot(l,b2gaus11+lin_data_fin,c='darkred',linestyle='-')
            plt.plot(l,b2gaus12+lin_data_fin,c='dodgerblue',linestyle='-')
        if meth == 'M2':
            plt.plot(l,b2gaus3+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,b2gaus4+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,b2gaus5+lin_data_fin,c='goldenrod',linestyle='-')
            plt.plot(l,b2gaus10+lin_data_fin,c='darkred',linestyle='-')
            plt.plot(l,b2gaus11+lin_data_fin,c='darkred',linestyle='-')
            plt.plot(l,b2gaus12+lin_data_fin,c='darkred',linestyle='-')
        plt.plot(l,b2gausb+lin_data_fin,c='darkviolet',linestyle='-')#,label='Broad component')
        plt.plot(l,twobroad_fit,'r-')
	    
        frame1.set_xticklabels([]) 			# Remove x-tic labels for the first frame
        plt.ylabel('Flux (erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(prune='lower'))
        plt.text(0.81,0.94,'N$_{S}$',color='darkgreen',transform=frame1.transAxes,fontsize=21)
        plt.text(0.81,0.9,'N$_{O}$',color='goldenrod',transform=frame1.transAxes,fontsize=21)
        plt.text(0.84,0.92,'+',color='k',transform=frame1.transAxes,fontsize=21)
        plt.text(0.87,0.94,'S$_{S}$',color='dodgerblue',transform=frame1.transAxes,fontsize=21)
        plt.text(0.87,0.9,'S$_{O}$',color='darkred',transform=frame1.transAxes,fontsize=21)
        plt.text(0.90,0.92,'+',color='k',transform=frame1.transAxes,fontsize=21)
        plt.text(0.93,0.92,'B',color='darkviolet',transform=frame1.transAxes,fontsize=21)

        # RESIDUAL plot
        frame2 = fig1.add_axes((.1,.1,.85,.15))
        plt.plot(l,np.zeros(len(l)),c='orange',linestyle='-')         	# Line around zero
        plt.plot(l,data_cor-twobroad_fit,c='k')		# Main
        plt.xlabel('Wavelength ($\AA$)',fontsize=19)
        plt.ylabel('Residuals',fontsize=19)
        plt.tick_params(axis='both', labelsize=17)
        plt.xlim(l[0],l[-1])
        plt.plot(l,np.zeros(len(l))+1.5*stadev,c='orange',linestyle=(0,(5,8)))	# 3 sigma upper limit
        plt.plot(l,np.zeros(len(l))-1.5*stadev,c='orange',linestyle=(0,(5,8))) 	# 3 sigma down limit
        plt.ylim(-(3*stadev)*3,(3*stadev)*3)
	    
        plt.savefig(parentFold+'adj_met'+str(meth)+'_full_2comp_broadH.pdf',format='pdf',bbox_inches='tight',pad_inches=0.2)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        frame1.text(2.2,max(data_cor)+0.05, textstr, fontsize=12,verticalalignment='top', bbox=props)
        plt.savefig(parentFold+'adj_met'+str(meth)+'_full_2comp_broadH.png',bbox_inches='tight',pad_inches=0.2)
   
        return [vbS2,evbS2,sigbS2,esigbS2,b2gaus1,b2gaus2,b2gaus3,b2gaus4,b2gaus5,b2gaus6,b2gaus7,b2gaus8,b2gaus9,b2gaus10,b2gaus11,b2gaus12,b2gaus13,b2gaus14,b2gausb]

    else: 
        print('ERROR in the plot! ')
        return 0
