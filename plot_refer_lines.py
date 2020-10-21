import numpy as np
import matplotlib.pyplot as plt
from pyspeckit import speclines as pylines
from astropy.constants import c
import Ofuncts

def refer_plot(parentFold,l,data_cor,meth,linresu,oneresu,tworesu,threeresu,l1,l2,l3,l4,l11,l12,l13,l14,std0,std1):
	'''
	It gives the plots for one and two components in the reference lines NII+Halpha

	The parameters needed are:
	parentFold:      Path to the data
	l:         Wavelength range
	data_cor:  Flux for each wavelength
	meth:      Method to be applied (S/O)
	linresu:   Result of the linear fit of the spectra
	oneresu:   Result of the linear+gaussian fit for the reference lines with one component
	tworesu:   Result of the linear+gaussian fit for the reference lines with two components
	threeresu: Result of the linear+gaussian fit for the reference lines with three components
	l1-l14:    Parts of the spectra where the lines are located
	std0/std1: Where the standard deviation of the continuum is calculated
	'''
	# Rest values of the line wavelengths 
	l_Halpha = 6562.801
	l_SII_1  = pylines.optical.lines['SIIa'][0]     # 6716.44
	l_SII_2  = pylines.optical.lines['SIIb'][0]     # 6730.82
	l_OI_1  = pylines.optical.lines['OI'][0]     # 6300.30
	
	# Constants and STIS parameters
	v_luz = c.value/10**3 # km/s
	plate_scale = 0.2
	fwhm = 2*np.sqrt(2*np.log(2)) # times sigma
	pix_to_v = 49.	# km/s
	sig_inst = 1.0
	
	# Parameters of the linear fit and the std of the continuum	
	new_slop = linresu.values['slope']
	new_intc = linresu.values['intc']
	stadev = np.std(data_cor[std0:std1])
	
	############################### PLOT and PRINT for the SII lines ##########################################
	#
	# Now we create the individual gaussians in order to plot and print the results for only 1 component
	#print('				RESULTS OF THE FIT: ')
	#print('Linear fit equation: {:.5f}*x + {:.5f}'.format(linresu.values['slope'], linresu.values['intc']))
	#print('')
	#print('The rest of the results can be displayed all together with two/oneresu.params; the data can be accesed with two/oneresu.values['']')
	#print('')
	#print('The chi-square of the fit for 1 gaussian for the reference line is: {:.5f}'.format(oneresu.chisqr))
	#print('The chi-square of the fit for 2 gaussian for the reference line is: {:.5f}'.format(tworesu.chisqr))
	#print('')
	
	# Now we create and plot the individual gaussians of the fit
	gaus1 = Ofuncts.gaussian(l,oneresu.values['mu_0'],oneresu.values['sig_0'],oneresu.values['amp_0']) 
	gaus2 = Ofuncts.gaussian(l,oneresu.values['mu_1'],oneresu.values['sig_1'],oneresu.values['amp_1'])
	gaus21 = Ofuncts.gaussian(l,tworesu.values['mu_0'],tworesu.values['sig_0'],tworesu.values['amp_0']) 
	gaus22 = Ofuncts.gaussian(l,tworesu.values['mu_1'],tworesu.values['sig_1'],tworesu.values['amp_1'])
	gaus24 = Ofuncts.gaussian(l,tworesu.values['mu_20'],tworesu.values['sig_20'],tworesu.values['amp_20'])
	gaus25 = Ofuncts.gaussian(l,tworesu.values['mu_21'],tworesu.values['sig_21'],tworesu.values['amp_21'])
	gaus31 = Ofuncts.gaussian(l,threeresu.values['mu_0'],threeresu.values['sig_0'],tworesu.values['amp_0']) 
	gaus32 = Ofuncts.gaussian(l,threeresu.values['mu_1'],threeresu.values['sig_1'],threeresu.values['amp_1'])
	gaus33 = Ofuncts.gaussian(l,threeresu.values['mu_20'],threeresu.values['sig_20'],threeresu.values['amp_20'])
	gaus34 = Ofuncts.gaussian(l,threeresu.values['mu_21'],threeresu.values['sig_21'],threeresu.values['amp_21'])
	gaus35 = Ofuncts.gaussian(l,threeresu.values['mu_30'],threeresu.values['sig_30'],threeresu.values['amp_30'])
	gaus36 = Ofuncts.gaussian(l,threeresu.values['mu_31'],threeresu.values['sig_31'],threeresu.values['amp_31'])
	onefin_fit = Ofuncts.twogaussian(l,new_slop,new_intc,
					 oneresu.values['mu_0'],oneresu.values['sig_0'],oneresu.values['amp_0'],
					 oneresu.values['mu_1'],oneresu.values['sig_1'],oneresu.values['amp_1'])
	twofin_fit = Ofuncts.funcSII2comp(l,new_slop,new_intc,
					 tworesu.values['mu_0'],tworesu.values['sig_0'],tworesu.values['amp_0'],
					 tworesu.values['mu_1'],tworesu.values['sig_1'],tworesu.values['amp_1'],
					 tworesu.values['mu_20'],tworesu.values['sig_20'],tworesu.values['amp_20'],
					 tworesu.values['mu_21'],tworesu.values['sig_21'],tworesu.values['amp_21'])
	threefin_fit = Ofuncts.funcSII3comp(l,new_slop,new_intc,
					 threeresu.values['mu_0'],threeresu.values['sig_0'],threeresu.values['amp_0'],
					 threeresu.values['mu_1'],threeresu.values['sig_1'],threeresu.values['amp_1'],
					 threeresu.values['mu_20'],threeresu.values['sig_20'],threeresu.values['amp_20'],
					 threeresu.values['mu_21'],threeresu.values['sig_21'],threeresu.values['amp_21'],
					 threeresu.values['mu_30'],threeresu.values['sig_30'],threeresu.values['amp_30'],
					 threeresu.values['mu_31'],threeresu.values['sig_31'],threeresu.values['amp_31'])
	if meth == 'S':
    	# one component
	    ixSTD1 = l[np.where(l > oneresu.values['mu_0'])[0][0]] - 1.5*oneresu.values['sig_0']
	    ixSTD2 = l[np.where(l > oneresu.values['mu_0'])[0][0]] + 1.5*oneresu.values['sig_0']
	    ixSTD1_2 = l[np.where(l > oneresu.values['mu_1'])[0][0]] - 1.5*oneresu.values['sig_1']
	    ixSTD2_2 = l[np.where(l > oneresu.values['mu_1'])[0][0]] + 1.5*oneresu.values['sig_1']
	    std_2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]]-onefin_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]])
	    std_1 = np.std(data_cor[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]]-onefin_fit[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]])
	    #std_2 = np.std(data_cor[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]]-onefin_fit[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]])
	    #std_1 = np.std(data_cor[np.where(l<ixSTD1_2)[0][-1]:np.where(l>ixSTD2_2)[0][0]]-onefin_fit[np.where(l<ixSTD1_2)[0][-1]:np.where(l>ixSTD2_2)[0][0]])
	    ep_1 = std_1/stadev
	    ep_2 = std_2/stadev
	    print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 component is... ')
	    print('	For the SII2 line: '+str(ep_2)+' < 3')
	    print('	For the SII1 line: '+str(ep_1)+' < 3')
	    # two components
	    ixSTD1 = l[np.where(l > tworesu.values['mu_0'])[0][0]] - 1.5*tworesu.values['sig_0']
	    ixSTD2 = l[np.where(l > tworesu.values['mu_0'])[0][0]] + 1.5*tworesu.values['sig_0']
	    ixSTD1_2 = l[np.where(l > tworesu.values['mu_1'])[0][0]] - 1.5*tworesu.values['sig_1']
	    ixSTD2_2 = l[np.where(l > tworesu.values['mu_1'])[0][0]] + 1.5*tworesu.values['sig_1']
	    #std2_2 = np.std(data_cor[np.where(l>ixSTD1)[0][0]:np.where(l<ixSTD2)[0][-1]]-twofin_fit[np.where(l>ixSTD1)[0][0]:np.where(l<ixSTD2)[0][-1]])
	    #std2_1 = np.std(data_cor[np.where(l>ixSTD1_2)[0][0]:np.where(l<ixSTD2_2)[0][-1]]-twofin_fit[np.where(l>ixSTD1_2)[0][0]:np.where(l<ixSTD2_2)[0][-1]])
	    std2_2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]]-twofin_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]])
	    std2_1 = np.std(data_cor[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]]-twofin_fit[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]])
	    ep2_1 = std2_1/stadev
	    ep2_2 = std2_2/stadev
	    print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 2 components is... ')
	    print('	For the SII2 line: '+str(ep2_2)+' < 3')
	    print('	For the SII1 line: '+str(ep2_1)+' < 3')

            # three components
	    ixSTD1 = l[np.where(l > threeresu.values['mu_0'])[0][0]] - 1.5*threeresu.values['sig_0']
	    ixSTD2 = l[np.where(l > threeresu.values['mu_0'])[0][0]] + 1.5*threeresu.values['sig_0']
	    ixSTD1_2 = l[np.where(l > threeresu.values['mu_1'])[0][0]] - 1.5*threeresu.values['sig_1']
	    ixSTD2_2 = l[np.where(l > threeresu.values['mu_1'])[0][0]] + 1.5*threeresu.values['sig_1']
	    #std3_2 = np.std(data_cor[np.where(l>ixSTD1)[0][0]:np.where(l<ixSTD2)[0][-1]]-threefin_fit[np.where(l>ixSTD1)[0][0]:np.where(l<ixSTD2)[0][-1]])
	    #std3_1 = np.std(data_cor[np.where(l>ixSTD1_2)[0][0]:np.where(l<ixSTD2_2)[0][-1]]-threefin_fit[np.where(l>ixSTD1_2)[0][0]:np.where(l<ixSTD2_2)[0][-1]])
            std3_2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]]-threefin_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]])
            std3_1 = np.std(data_cor[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]]-threefin_fit[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]])
            ep3_1 = std3_1/stadev
            ep3_2 = std3_2/stadev
            print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 3 components is... ')
            print('     For the SII2 line: '+str(ep3_2)+' < 3')
            print('     For the SII1 line: '+str(ep3_1)+' < 3')

        elif meth == 'O':
	    # one component
            std_1 = np.std(data_cor[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10]-onefin_fit[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10])
            std_2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-onefin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
            ep_1 = std_1/stadev
            ep_2 = std_2/stadev
            print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 component is... ')
            print('     For the OI2 line: '+str(ep_2)+' < 3')
            print('     For the OI1 line: '+str(ep_1)+' < 3')
            # two components
            std2_1 = np.std(data_cor[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10]-twofin_fit[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10])
            std2_2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-twofin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
            ep2_1 = std2_1/stadev
            ep2_2 = std2_2/stadev
            print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 2 components is... ')
            print('     For the OI2 line: '+str(ep2_2)+' < 3')
            print('     For the OI1 line: '+str(ep2_1)+' < 3')

            # three components
            std3_1 = np.std(data_cor[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10]-threefin_fit[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10])
            std3_2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-threefin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
            ep3_1 = std3_1/stadev
            ep3_2 = std3_2/stadev
            print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 3 components is... ')
            print('     For the OI2 line: '+str(ep3_2)+' < 3')
            print('     For the OI1 line: '+str(ep3_1)+' < 3')
            l_SII_2 = l_OI_1 


	# We determine the maximum flux of the fit for all the lines, and the velocity and sigma components
	maxN1 = onefin_fit[np.where(abs(oneresu.values['mu_1']-l)<1.0)[0][0]]
	maxN2 = onefin_fit[np.where(abs(oneresu.values['mu_0']-l)<1.0)[0][0]]
	max2N1 = twofin_fit[np.where(abs(tworesu.values['mu_1']-l)<1.0)[0][0]]
	max2N2 = twofin_fit[np.where(abs(tworesu.values['mu_0']-l)<1.0)[0][0]]
	max3N1 = threefin_fit[np.where(abs(threeresu.values['mu_1']-l)<1.0)[0][0]]
	max3N2 = threefin_fit[np.where(abs(threeresu.values['mu_0']-l)<1.0)[0][0]]
	# one component
	vS2 = v_luz*((oneresu.values['mu_0']-l_SII_2)/l_SII_2)
	sigS2 = pix_to_v*np.sqrt(oneresu.values['sig_0']**2-sig_inst**2)
	# two comps
	v2S2 = v_luz*((tworesu.values['mu_0']-l_SII_2)/l_SII_2)
	v20S2 = v_luz*((tworesu.values['mu_20']-l_SII_2)/l_SII_2)
	sig2S2 = pix_to_v*np.sqrt(tworesu.values['sig_0']**2-sig_inst**2)
	sig20S2 = pix_to_v*np.sqrt(tworesu.values['sig_20']**2-sig_inst**2)
	# three comps
	v3S2 = v_luz*((threeresu.values['mu_0']-l_SII_2)/l_SII_2)
	v30S2 = v_luz*((threeresu.values['mu_20']-l_SII_2)/l_SII_2)
	v300S2 = v_luz*((threeresu.values['mu_30']-l_SII_2)/l_SII_2)
	sig3S2 = pix_to_v*np.sqrt(threeresu.values['sig_0']**2-sig_inst**2)
	sig30S2 = pix_to_v*np.sqrt(threeresu.values['sig_20']**2-sig_inst**2)
	sig300S2 = pix_to_v*np.sqrt(threeresu.values['sig_30']**2-sig_inst**2)

	if oneresu.params['mu_0'].stderr == None: 
	    print('Problem determining the errors!')
	    evS2,esigS2 = 0.,0.
	elif oneresu.params['mu_0'].stderr != None: 
	    evS2 = ((v_luz/l_SII_2)*oneresu.params['mu_0'].stderr)
	    esigS2 = pix_to_v*(oneresu.values['sig_0']*oneresu.params['sig_0'].stderr)/(np.sqrt(oneresu.values['sig_0']**2-sig_inst**2))

	if tworesu.params['mu_20'].stderr == None:
	    print('Problem determining the errors!')
	    ev20S2,ev2S2,esig2S2,esig20S2 = 0.,0.,0.,0.
	elif tworesu.params['mu_20'].stderr != None:
	    ev2S2 = ((v_luz/l_Halpha)*tworesu.params['mu_0'].stderr)
	    ev20S2 = ((v_luz/l_Halpha)*tworesu.params['mu_20'].stderr)
	    esig2S2 = pix_to_v*(tworesu.values['sig_0']*tworesu.params['sig_0'].stderr)/(np.sqrt(tworesu.values['sig_0']**2-sig_inst**2))
	    esig20S2 = pix_to_v*(tworesu.values['sig_20']*tworesu.params['sig_20'].stderr)/(np.sqrt(tworesu.values['sig_20']**2-sig_inst**2))

        if threeresu.params['mu_30'].stderr == None:
            print('Problem determining the errors!')
            ev30S2,ev300S2,ev3S2,esig3S2,esig30S2,esig300S2 = 0.,0.,0.,0.,0.,0.
        elif threeresu.params['mu_30'].stderr != None:
            ev3S2 = ((v_luz/l_Halpha)*threeresu.params['mu_0'].stderr)
            ev30S2 = ((v_luz/l_Halpha)*threeresu.params['mu_20'].stderr)
            ev300S2 = ((v_luz/l_Halpha)*threeresu.params['mu_30'].stderr)
            esig3S2 = pix_to_v*(threeresu.values['sig_0']*threeresu.params['sig_0'].stderr)/(np.sqrt(threeresu.values['sig_0']**2-sig_inst**2))
            esig30S2 = pix_to_v*(threeresu.values['sig_20']*threeresu.params['sig_20'].stderr)/(np.sqrt(threeresu.values['sig_20']**2-sig_inst**2))
            esig300S2 = pix_to_v*(threeresu.values['sig_30']*threeresu.params['sig_30'].stderr)/(np.sqrt(threeresu.values['sig_30']**2-sig_inst**2))


	textstr = '\n'.join((r'$V_{Ha_{2}}$ = '+ '{:.2f} +- {:.2f}'.format(vS2,evS2),
			    r'$\sigma_{Ha_{2}}$ = '+ '{:.2f} +- {:.2f}'.format(sigS2,esigS2),
			    r'$\frac{F_{NII_{2}}}{F_{NII_{1}}}$ = '+ '{:.3f}'.format(maxN2/maxN1)))
	textstr2 = '\n'.join((r'$V_{Ha_{2-1comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v2S2,ev2S2),
			    r'$V_{Ha_{2-2comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v20S2,ev20S2),
			    r'$\sigma_{Ha_{2-1comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig2S2,esig2S2),
			    r'$\sigma_{Ha_{2-2comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig20S2,esig20S2),
			    r'$\frac{F_{NII_{2}}}{F_{NII_{1}}}$ = '+ '{:.3f}'.format(max2N2/max2N1)))
#			    r'$F_{SII_{1}}$ = '+ '{:.3f}'.format(max2N1)+' $10^{-14}$'))
	textstr3 = '\n'.join((r'$V_{Ha_{3-1comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v3S2,ev3S2),
			    r'$V_{Ha_{3-2comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v30S2,ev30S2),
			    r'$V_{Ha_{3-3comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v300S2,ev300S2),
			    r'$\sigma_{Ha_{3-1comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig3S2,esig3S2),
			    r'$\sigma_{Ha_{3-2comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig30S2,esig30S2),
			    r'$\sigma_{Ha_{3-3comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig300S2,esig300S2),
			    r'$\frac{F_{NII_{2}}}{F_{NII_{1}}}$ = '+ '{:.3f}'.format(max3N2/max3N1)))
#			    r'$F_{SII_{1}}$ = '+ '{:.3f}'.format(max2N1)+' $10^{-14}$'))
	
	
	################################################ PLOT ######################################################
	plt.close()
	# MAIN plot
	fig1   = plt.figure(1,figsize=(10, 9))
	frame1 = fig1.add_axes((.1,.25,.85,.65)) 	     # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
	plt.plot(l,data_cor,'k')			     # Initial data
	plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='y',linestyle='-.',label='Linear fit')
	plt.plot(l,gaus1,c='darkgreen',linestyle='-')
	plt.plot(l,gaus2,c='darkgreen',linestyle='-')
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	frame1.text(6350.,max(data_cor), textstr, fontsize=12,verticalalignment='top', bbox=props)
	plt.plot(l[std0:std1],data_cor[std0:std1],c='darkorange')	# Zone where the stddev is calculated
	plt.plot(l,onefin_fit,'r-')

	frame1.set_xticklabels([]) 			# Remove x-tic labels for the first frame
	plt.ylabel('Flux ($\mathrm{erg/s/cm^{2} / \AA}$)',fontsize=14)
	plt.tick_params(axis='both', labelsize=12)
	plt.xlim(l[0],l[-1])
	plt.legend(loc='best')
	
	# RESIDUAL plot
	frame2 = fig1.add_axes((.1,.1,.85,.15))
	plt.plot(l,data_cor-onefin_fit,c='k')		# Main
	plt.xlabel('Wavelength ($\AA$)',fontsize=14)
	plt.ylabel('Residuals',fontsize=14)
	plt.tick_params(axis='both', labelsize=12)
	plt.xlim(l[0],l[-1])
	plt.plot(l,np.zeros(len(l)),c='grey',linestyle='--')         	# Line around zero
	plt.plot(l,np.zeros(len(l))+1.5*stadev,c='grey',linestyle='--')	# 3 sigma upper limit
	plt.plot(l,np.zeros(len(l))-1.5*stadev,c='grey',linestyle='--') 	# 3 sigma down limit
	plt.ylim(-(3*stadev)*2,(3*stadev)*2)
	
	plt.savefig(parentFold+'adj_met'+str(meth)+'_ref_1comp.png')
	
	#######################################################################################
	# Two components in reference line
	# MAIN plot
	fig2   = plt.figure(2,figsize=(10, 9))
	frame3 = fig2.add_axes((.1,.25,.85,.65)) 	     # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
	plt.plot(l,data_cor,'k')		     # Initial data
	plt.plot(l,twofin_fit,'r-')
	plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='y',linestyle='-.',label='Linear fit')
	plt.plot(l,gaus21,c='darkgreen',linestyle='-')
	plt.plot(l,gaus22,c='darkgreen',linestyle='-')
	plt.plot(l,gaus24,c='dodgerblue',linestyle='-')
	plt.plot(l,gaus25,c='dodgerblue',linestyle='-')
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	frame3.text(6350.,max(data_cor), textstr2, fontsize=12,verticalalignment='top', bbox=props)
	plt.plot(l[std0:std1],data_cor[std0:std1],c='darkorange')	# Zone where the stddev is calculated
	
	frame3.set_xticklabels([]) 			# Remove x-tic labels for the first frame
	plt.ylabel('Flux ($\mathrm{erg/s/cm^{2} / \AA}$)',fontsize=14)
	plt.tick_params(axis='both', labelsize=12)
	plt.xlim(l[0],l[-1])
	plt.legend(loc='best')
	
	# RESIDUAL plot
	frame4 = fig2.add_axes((.1,.1,.85,.15))
	plt.plot(l,data_cor-twofin_fit,c='k')		# Main
	plt.xlabel('Wavelength ($\AA$)',fontsize=14)
	plt.ylabel('Residuals',fontsize=14)
	plt.tick_params(axis='both', labelsize=12)
	plt.xlim(l[0],l[-1])
	plt.plot(l,np.zeros(len(l)),c='grey',linestyle='--')         	# Line around zero
	plt.plot(l,np.zeros(len(l))+1.5*stadev,c='grey',linestyle='--')	# 3 sigma upper limit
	plt.plot(l,np.zeros(len(l))-1.5*stadev,c='grey',linestyle='--') 	# 3 sigma down limit
	plt.ylim(-(3*stadev)*2,(3*stadev)*2)
	
	plt.savefig(parentFold+'adj_met'+str(meth)+'_ref_2comp.png')


        #######################################################################################
        # Three components in reference line
        # MAIN plot
        fig3   = plt.figure(3,figsize=(10, 9))
        frame5 = fig3.add_axes((.1,.25,.85,.65))             # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
        plt.plot(l,data_cor,'k')                     # Initial data
        plt.plot(l,threefin_fit,'r-')
        plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='y',linestyle='-.',label='Linear fit')
        plt.plot(l,gaus31,c='darkgreen',linestyle='-')
        plt.plot(l,gaus32,c='darkgreen',linestyle='-')
        plt.plot(l,gaus33,c='dodgerblue',linestyle='-')
        plt.plot(l,gaus34,c='dodgerblue',linestyle='-')
        plt.plot(l,gaus35,c='magenta',linestyle='-')
        plt.plot(l,gaus36,c='magenta',linestyle='-')
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        frame5.text(6350.,max(data_cor), textstr3, fontsize=12,verticalalignment='top', bbox=props)
        plt.plot(l[std0:std1],data_cor[std0:std1],c='darkorange')       # Zone where the stddev is calculated

        frame5.set_xticklabels([])                      # Remove x-tic labels for the first frame
        plt.ylabel('Flux ($\mathrm{erg/s/cm^{2} / \AA}$)',fontsize=14)
        plt.tick_params(axis='both', labelsize=12)
        plt.xlim(l[0],l[-1])
        plt.legend(loc='best')

        # RESIDUAL plot
        frame6 = fig3.add_axes((.1,.1,.85,.15))
        plt.plot(l,data_cor-threefin_fit,c='k')           # Main
        plt.xlabel('Wavelength ($\AA$)',fontsize=14)
        plt.ylabel('Residuals',fontsize=14)
        plt.tick_params(axis='both', labelsize=12)
        plt.xlim(l[0],l[-1])
        plt.plot(l,np.zeros(len(l)),c='grey',linestyle='--')            # Line around zero
        plt.plot(l,np.zeros(len(l))+1.5*stadev,c='grey',linestyle='--') # 3 sigma upper limit
        plt.plot(l,np.zeros(len(l))-1.5*stadev,c='grey',linestyle='--')         # 3 sigma down limit
        plt.ylim(-(3*stadev)*2,(3*stadev)*2)

        plt.savefig(parentFold+'adj_met'+str(meth)+'_ref_3comp.png')


	return ep_1,ep_2,ep2_1,ep2_2,ep3_1,ep3_2
