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
	gaus3 = Ofuncts.gaussian(l,oneresu.values['mu_5'],oneresu.values['sig_5'],oneresu.values['amp_5'])
	gaus4 = Ofuncts.gaussian(l,oneresu.values['mu_6'],oneresu.values['sig_6'],oneresu.values['amp_6'])
	gaus21 = Ofuncts.gaussian(l,tworesu.values['mu_0'],tworesu.values['sig_0'],tworesu.values['amp_0']) 
	gaus22 = Ofuncts.gaussian(l,tworesu.values['mu_1'],tworesu.values['sig_1'],tworesu.values['amp_1'])
	gaus23 = Ofuncts.gaussian(l,tworesu.values['mu_5'],tworesu.values['sig_5'],tworesu.values['amp_5'])
	gaus24 = Ofuncts.gaussian(l,tworesu.values['mu_6'],tworesu.values['sig_6'],tworesu.values['amp_6'])
	gaus25 = Ofuncts.gaussian(l,tworesu.values['mu_20'],tworesu.values['sig_20'],tworesu.values['amp_20'])
	gaus26 = Ofuncts.gaussian(l,tworesu.values['mu_21'],tworesu.values['sig_21'],tworesu.values['amp_21'])
	gaus27 = Ofuncts.gaussian(l,tworesu.values['mu_25'],tworesu.values['sig_25'],tworesu.values['amp_25'])
	gaus28 = Ofuncts.gaussian(l,tworesu.values['mu_26'],tworesu.values['sig_26'],tworesu.values['amp_26'])
	gaus31 = Ofuncts.gaussian(l,threeresu.values['mu_0'],threeresu.values['sig_0'],tworesu.values['amp_0']) 
	gaus32 = Ofuncts.gaussian(l,threeresu.values['mu_1'],threeresu.values['sig_1'],threeresu.values['amp_1'])
	gaus33 = Ofuncts.gaussian(l,threeresu.values['mu_5'],threeresu.values['sig_5'],threeresu.values['amp_5'])
	gaus34 = Ofuncts.gaussian(l,threeresu.values['mu_6'],threeresu.values['sig_6'],threeresu.values['amp_6'])
	gaus35 = Ofuncts.gaussian(l,threeresu.values['mu_20'],threeresu.values['sig_20'],threeresu.values['amp_20'])
	gaus36 = Ofuncts.gaussian(l,threeresu.values['mu_21'],threeresu.values['sig_21'],threeresu.values['amp_21'])
	gaus37 = Ofuncts.gaussian(l,threeresu.values['mu_25'],threeresu.values['sig_25'],threeresu.values['amp_25'])
	gaus38 = Ofuncts.gaussian(l,threeresu.values['mu_26'],threeresu.values['sig_26'],threeresu.values['amp_26'])
	gaus39 = Ofuncts.gaussian(l,threeresu.values['mu_30'],threeresu.values['sig_30'],threeresu.values['amp_30'])
	gaus310 = Ofuncts.gaussian(l,threeresu.values['mu_31'],threeresu.values['sig_31'],threeresu.values['amp_31'])
	gaus311 = Ofuncts.gaussian(l,threeresu.values['mu_35'],threeresu.values['sig_35'],threeresu.values['amp_35'])
	gaus312 = Ofuncts.gaussian(l,threeresu.values['mu_36'],threeresu.values['sig_36'],threeresu.values['amp_36'])
	onefin_fit = Ofuncts.funcSII2comp(l,new_slop,new_intc,
					 oneresu.values['mu_0'],oneresu.values['sig_0'],oneresu.values['amp_0'],
					 oneresu.values['mu_1'],oneresu.values['sig_1'],oneresu.values['amp_1'],
					 oneresu.values['mu_5'],oneresu.values['sig_5'],oneresu.values['amp_5'],
					 oneresu.values['mu_6'],oneresu.values['sig_6'],oneresu.values['amp_6'])
	twofin_fit = Ofuncts.eightgaussian(l,new_slop,new_intc,
					 tworesu.values['mu_0'],tworesu.values['sig_0'],tworesu.values['amp_0'],
					 tworesu.values['mu_1'],tworesu.values['sig_1'],tworesu.values['amp_1'],
					 tworesu.values['mu_5'],tworesu.values['sig_5'],tworesu.values['amp_5'],
					 tworesu.values['mu_6'],tworesu.values['sig_6'],tworesu.values['amp_6'],
					 tworesu.values['mu_20'],tworesu.values['sig_20'],tworesu.values['amp_20'],
					 tworesu.values['mu_25'],tworesu.values['sig_25'],tworesu.values['amp_25'],
					 tworesu.values['mu_26'],tworesu.values['sig_26'],tworesu.values['amp_26'],
					 tworesu.values['mu_21'],tworesu.values['sig_21'],tworesu.values['amp_21'])
	threefin_fit = Ofuncts.twelvegaussian(l,new_slop,new_intc,
					 threeresu.values['mu_0'],threeresu.values['sig_0'],threeresu.values['amp_0'],
					 threeresu.values['mu_1'],threeresu.values['sig_1'],threeresu.values['amp_1'],
					 threeresu.values['mu_5'],threeresu.values['sig_5'],threeresu.values['amp_5'],
					 threeresu.values['mu_6'],threeresu.values['sig_6'],threeresu.values['amp_6'],
					 threeresu.values['mu_20'],threeresu.values['sig_20'],threeresu.values['amp_20'],
					 threeresu.values['mu_21'],threeresu.values['sig_21'],threeresu.values['amp_21'],
					 threeresu.values['mu_25'],threeresu.values['sig_25'],threeresu.values['amp_25'],
					 threeresu.values['mu_26'],threeresu.values['sig_26'],threeresu.values['amp_26'],
					 threeresu.values['mu_30'],threeresu.values['sig_30'],threeresu.values['amp_30'],
					 threeresu.values['mu_31'],threeresu.values['sig_31'],threeresu.values['amp_31'],
					 threeresu.values['mu_35'],threeresu.values['sig_35'],threeresu.values['amp_35'],
					 threeresu.values['mu_36'],threeresu.values['sig_36'],threeresu.values['amp_36'])

        residuals1 = data_cor - oneresu.best_fit
        residuals2 = data_cor - tworesu.best_fit
        residuals3 = data_cor - threeresu.best_fit

        ##########################################################################################################################
        # Now we have to calculate the epsilon value for both SII and OI (mu_5) lines
        ep_1,ep_2,ep2_1,ep2_2,ep3_1,ep3_2 = 0.,0.,0.,0.,0.,0.

	# [SII] lines

    	# one component
	ixSTD1 = l[np.where(l > np.mean([oneresu.values['mu_0'],oneresu.values['mu_1']]))[0][0]] - 3.0*oneresu.values['sig_0']
	ixSTD2 = l[np.where(l > np.mean([oneresu.values['mu_0'],oneresu.values['mu_1']]))[0][0]] + 3.0*oneresu.values['sig_0']
        DA = residuals1[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]]
        std_1 = np.std(DA)
        # Fixed range 
        ep_1 = std_1/stadev
        #ep_1 = np.std(residuals1[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]])/stadev
        #ep_2 = np.std(residuals1[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]])/stadev
	print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 component is... ')
	print('	For the SII2 lines: '+str(ep_1)+' < 3 -> '+str(ep_1<3))
	#print('	For the SII1 line: '+str(ep_2)+' < 3 -> '+str(ep_2<3))
	print('	For the SII1 line: '+str(np.std(residuals1[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]]))+' < 3') 

	# two components
	ixSTD1 = l[np.where(l > np.mean([tworesu.values['mu_0'],tworesu.values['mu_1']]))[0][0]] - 3.0*tworesu.values['sig_20']
	ixSTD2 = l[np.where(l > np.mean([tworesu.values['mu_0'],tworesu.values['mu_1']]))[0][0]] + 3.0*tworesu.values['sig_20']
        DA = residuals2[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]]
        std2_1 = np.std(DA)
	ep2_1 = std2_1/stadev
	#ep2_2 = std2_2/stadev
	#ep2_1 = np.std(residuals2[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]])/stadev
	#ep2_2 = np.std(residuals2[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]])/stadev
	print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 2 components is... ')
	print('	For the SII2 lines: '+str(ep2_1)+' < 3')
	#print('	For the SII1 line: '+str(ep2_2)+' < 3')

        # three components
	ixSTD1 = l[np.where(l > np.mean([threeresu.values['mu_0'],threeresu.values['mu_1']]))[0][0]] - 3.0*threeresu.values['sig_20']
	ixSTD2 = l[np.where(l > np.mean([threeresu.values['mu_0'],threeresu.values['mu_1']]))[0][0]] + 3.0*threeresu.values['sig_20']
        DA = residuals3[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]]
        #DA = np.append(DA1,DA2)
        std3_1 = np.std(DA)
        #std3_2 = np.std(residuals3[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]])
        #std3_1 = np.std(residuals3[np.where(l<l3)[0][-1]:np.where(l>l4)[0][0]])
        ep3_1 = std3_1/stadev
        #ep3_2 = std3_2/stadev
        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 3 components is... ')
        print('     For the SII2 lines: '+str(ep3_1)+' < 3')
        #print('     For the SII1 line: '+str(ep3_2)+' < 3')


        # [OI] 6300 AA line: 

	# one component
	ixSTD1 = l[np.where(l > oneresu.values['mu_5'])[0][0]] - 2.0*oneresu.values['sig_5']
	ixSTD2 = l[np.where(l > oneresu.values['mu_5'])[0][0]] + 2.0*oneresu.values['sig_5']
        DA = residuals1[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]]
        std_1 = np.std(DA)
        # Fixed range 
        #std_1 = np.std(data_cor[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10]-onefin_fit[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10])
        #std_2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-onefin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
        ep_2 = std_1/stadev
        #ep_2 = std_2/stadev
        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 component is... ')
        print('     For the OI1 line: '+str(ep_2)+' < 3')
        
        # two components
	ixSTD1 = l[np.where(l > tworesu.values['mu_5'])[0][0]] - 3.0*tworesu.values['sig_25']
	ixSTD2 = l[np.where(l > tworesu.values['mu_5'])[0][0]] + 3.0*tworesu.values['sig_25']
        DA = residuals2[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]]
        std2_1 = np.std(DA)
        #std2_1 = np.std(data_cor[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10]-twofin_fit[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10])
        #std2_2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-twofin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
        ep2_2 = std2_1/stadev
        #ep2_2 = std2_2/stadev
        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 2 components is... ')
        print('     For the OI1 line: '+str(ep2_2)+' < 3')

        # three components
	ixSTD1 = l[np.where(l > threeresu.values['mu_5'])[0][0]] - 3.0*threeresu.values['sig_25']
	ixSTD2 = l[np.where(l > threeresu.values['mu_5'])[0][0]] + 3.0*threeresu.values['sig_25']
        DA = residuals3[np.where(l<ixSTD1)[0][-1]:np.where(l>ixSTD2)[0][0]]
        std3_1 = np.std(DA)
        #std3_1 = np.std(data_cor[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10]-threefin_fit[np.where(l<l11)[0][-1]:np.where(l>l12)[0][0]+10])
        #std3_2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-threefin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
        #ep3_1 = std3_1/stadev
        ep3_2 = std3_1/stadev
        print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 3 components is... ')
        print('     For the OI1 line: '+str(ep3_2)+' < 3')

        ####################################################################################################################################
	# We determine the velocity and sigma components for each line
	# one component
	vS2 = v_luz*((oneresu.values['mu_0']-l_SII_2)/l_SII_2)
	vO2 = v_luz*((oneresu.values['mu_5']-l_OI_1)/l_OI_1)
	sigS2 = pix_to_v*np.sqrt(oneresu.values['sig_0']**2-sig_inst**2)
	sigO2 = pix_to_v*np.sqrt(oneresu.values['sig_5']**2-sig_inst**2)
	# two comps
	v2S2 = v_luz*((tworesu.values['mu_0']-l_SII_2)/l_SII_2)
	v20S2 = v_luz*((tworesu.values['mu_20']-l_SII_2)/l_SII_2)
	v2O2 = v_luz*((tworesu.values['mu_5']-l_OI_1)/l_OI_1)
	v20O2 = v_luz*((tworesu.values['mu_25']-l_OI_1)/l_OI_1)
	sig2S2 = pix_to_v*np.sqrt(tworesu.values['sig_0']**2-sig_inst**2)
	sig20S2 = pix_to_v*np.sqrt(tworesu.values['sig_20']**2-sig_inst**2)
	sig2O2 = pix_to_v*np.sqrt(tworesu.values['sig_5']**2-sig_inst**2)
	sig20O2 = pix_to_v*np.sqrt(tworesu.values['sig_25']**2-sig_inst**2)
	# three comps
	v3S2 = v_luz*((threeresu.values['mu_0']-l_SII_2)/l_SII_2)
	v30S2 = v_luz*((threeresu.values['mu_20']-l_SII_2)/l_SII_2)
	v300S2 = v_luz*((threeresu.values['mu_30']-l_SII_2)/l_SII_2)
	v3O2 = v_luz*((threeresu.values['mu_5']-l_OI_1)/l_OI_1)
	v30O2 = v_luz*((threeresu.values['mu_25']-l_OI_1)/l_OI_1)
	v300O2 = v_luz*((threeresu.values['mu_35']-l_OI_1)/l_OI_1)
	sig3S2 = pix_to_v*np.sqrt(threeresu.values['sig_0']**2-sig_inst**2)
	sig30S2 = pix_to_v*np.sqrt(threeresu.values['sig_20']**2-sig_inst**2)
	sig300S2 = pix_to_v*np.sqrt(threeresu.values['sig_30']**2-sig_inst**2)
	sig3O2 = pix_to_v*np.sqrt(threeresu.values['sig_5']**2-sig_inst**2)
	sig30O2 = pix_to_v*np.sqrt(threeresu.values['sig_25']**2-sig_inst**2)
	sig300O2 = pix_to_v*np.sqrt(threeresu.values['sig_35']**2-sig_inst**2)

	if oneresu.params['mu_0'].stderr == None: 
	    print('Problem determining the errors!')
	    evS2,esigS2 = 0.,0.
	elif oneresu.params['mu_0'].stderr != None: 
	    evS2 = ((v_luz/l_SII_2)*oneresu.params['mu_0'].stderr)
	    esigS2 = pix_to_v*(oneresu.values['sig_0']*oneresu.params['sig_0'].stderr)/(np.sqrt(oneresu.values['sig_0']**2-sig_inst**2))
	if oneresu.params['mu_5'].stderr == None: 
	    print('Problem determining the errors!')
	    evO2,esigO2 = 0.,0.
	elif oneresu.params['mu_5'].stderr != None: 
	    evO2 = ((v_luz/l_OI_1)*oneresu.params['mu_5'].stderr)
	    esigO2 = pix_to_v*(oneresu.values['sig_0']*oneresu.params['sig_0'].stderr)/(np.sqrt(oneresu.values['sig_0']**2-sig_inst**2))

	if tworesu.params['mu_20'].stderr == None:
	    print('Problem determining the errors!')
	    ev20S2,ev2S2,esig2S2,esig20S2 = 0.,0.,0.,0.
	elif tworesu.params['mu_20'].stderr != None:
	    ev2S2 = ((v_luz/l_SII_2)*tworesu.params['mu_0'].stderr)
	    ev20S2 = ((v_luz/l_SII_2)*tworesu.params['mu_20'].stderr)
	    esig2S2 = pix_to_v*(tworesu.values['sig_0']*tworesu.params['sig_0'].stderr)/(np.sqrt(tworesu.values['sig_0']**2-sig_inst**2))
	    esig20S2 = pix_to_v*(tworesu.values['sig_20']*tworesu.params['sig_20'].stderr)/(np.sqrt(tworesu.values['sig_20']**2-sig_inst**2))
	if tworesu.params['mu_25'].stderr == None:
	    print('Problem determining the errors!')
	    ev20O2,ev2O2,esig2O2,esig20O2 = 0.,0.,0.,0.
	elif tworesu.params['mu_25'].stderr != None:
	    ev2O2 = ((v_luz/l_OI_1)*tworesu.params['mu_5'].stderr)
	    ev20O2 = ((v_luz/l_OI_1)*tworesu.params['mu_25'].stderr)
	    esig2O2 = pix_to_v*(tworesu.values['sig_5']*tworesu.params['sig_5'].stderr)/(np.sqrt(tworesu.values['sig_5']**2-sig_inst**2))
	    esig20O2 = pix_to_v*(tworesu.values['sig_25']*tworesu.params['sig_25'].stderr)/(np.sqrt(tworesu.values['sig_25']**2-sig_inst**2))

        if threeresu.params['mu_30'].stderr == None:
            print('Problem determining the errors!')
            ev30S2,ev300S2,ev3S2,esig3S2,esig30S2,esig300S2 = 0.,0.,0.,0.,0.,0.
        elif threeresu.params['mu_30'].stderr != None:
            ev3S2 = ((v_luz/l_SII_2)*threeresu.params['mu_0'].stderr)
            ev30S2 = ((v_luz/l_SII_2)*threeresu.params['mu_20'].stderr)
            ev300S2 = ((v_luz/l_SII_2)*threeresu.params['mu_30'].stderr)
            esig3S2 = pix_to_v*(threeresu.values['sig_0']*threeresu.params['sig_0'].stderr)/(np.sqrt(threeresu.values['sig_0']**2-sig_inst**2))
            esig30S2 = pix_to_v*(threeresu.values['sig_20']*threeresu.params['sig_20'].stderr)/(np.sqrt(threeresu.values['sig_20']**2-sig_inst**2))
            esig300S2 = pix_to_v*(threeresu.values['sig_30']*threeresu.params['sig_30'].stderr)/(np.sqrt(threeresu.values['sig_30']**2-sig_inst**2))
        if threeresu.params['mu_35'].stderr == None:
            print('Problem determining the errors!')
            ev30O2,ev300O2,ev3O2,esig3O2,esig30O2,esig300O2 = 0.,0.,0.,0.,0.,0.
        elif threeresu.params['mu_35'].stderr != None:
            ev3O2 = ((v_luz/l_OI_1)*threeresu.params['mu_5'].stderr)
            ev30O2 = ((v_luz/l_OI_1)*threeresu.params['mu_25'].stderr)
            ev300O2 = ((v_luz/l_OI_1)*threeresu.params['mu_35'].stderr)
            esig3O2 = pix_to_v*(threeresu.values['sig_5']*threeresu.params['sig_5'].stderr)/(np.sqrt(threeresu.values['sig_5']**2-sig_inst**2))
            esig30O2 = pix_to_v*(threeresu.values['sig_25']*threeresu.params['sig_25'].stderr)/(np.sqrt(threeresu.values['sig_25']**2-sig_inst**2))
            esig300O2 = pix_to_v*(threeresu.values['sig_35']*threeresu.params['sig_35'].stderr)/(np.sqrt(threeresu.values['sig_35']**2-sig_inst**2))


	textstr = '\n'.join((r'$V_{SII}$ = '+ '{:.2f} +- {:.2f}'.format(vS2,evS2),
			    r'$V_{OI}$ = '+ '{:.2f} +- {:.2f}'.format(vO2,evO2),
			    r'$\sigma_{SII}$ = '+ '{:.2f} +- {:.2f}'.format(sigS2,esigS2),
			    r'$\sigma_{OI}$ = '+ '{:.2f} +- {:.2f}'.format(sigO2,esigO2)))
	textstr2 = '\n'.join((r'$V_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(v2S2,ev2S2),
			    r'$V_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v20S2,ev20S2),
			    r'$V_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(v2O2,ev2O2),
			    r'$V_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v20O2,ev20O2),
			    r'$\sigma_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig2S2,esig2S2),
			    r'$\sigma_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig20S2,esig20S2),
			    r'$\sigma_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig2O2,esigO2),
			    r'$\sigma_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig20O2,esig20O2)))
	textstr3 = '\n'.join((r'$V_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(v3S2,ev3S2),
			    r'$V_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v30S2,ev30S2),
			    r'$V_{SII_{3c}}$ = '+ '{:.2f} +- {:.2f}'.format(v300S2,ev300S2),
			    r'$V_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(v3O2,ev3O2),
			    r'$V_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(v30O2,ev30O2),
			    r'$V_{OI_{3c}}$ = '+ '{:.2f} +- {:.2f}'.format(v300O2,ev300O2),
			    r'$\sigma_{SII_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig3S2,esig3S2),
			    r'$\sigma_{SII_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig30S2,esig30S2),
			    r'$\sigma_{SII_{3c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig300S2,esig300S2),
			    r'$\sigma_{OI_{1c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig3O2,esig3O2),
			    r'$\sigma_{OI_{2c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig30O2,esig30O2),
			    r'$\sigma_{OI_{3c}}$ = '+ '{:.2f} +- {:.2f}'.format(sig300O2,esig300O2)))

        lin_data_fin = (linresu.values['slope']*l+linresu.values['intc'])
        # save the data of the individual fits in txt files
        np.savetxt(parentFold+'fit1comp_values.txt',np.c_[l,oneresu.data,oneresu.best_fit,lin_data_fin,gaus1,gaus2,gaus3,gaus4],fmt=('%5.6f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Wavelength\tReal_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tNarrow_SII6716\tNarrow_OI6300\tNarrow_OI6363'))
        np.savetxt(parentFold+'fit2comp_values.txt',np.c_[l,tworesu.data,tworesu.best_fit,lin_data_fin,gaus21,gaus22,gaus23,gaus24,gaus25,gaus26,gaus27,gaus28],fmt=('%5.6f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Wavelength\tReal_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tNarrow_SII6716\tNarrow_OI6300\tNarrow_OI6363\tSecond_SII6731\tSecond_SII6716\tSecond_OI6300\tSecond_OI6363'))
        np.savetxt(parentFold+'fit3comp_values.txt',np.c_[l,threeresu.data,threeresu.best_fit,lin_data_fin,gaus31,gaus32,gaus33,gaus34,gaus35,gaus36,gaus37,gaus38,gaus39,gaus310,gaus311,gaus312],fmt=('%5.6f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f','%5.10f'),header=('Wavelength\tReal_data\tBest_fit\tLineal_fit\tNarrow_SII6731\tNarrow_SII6716\tNarrow_OI6300\tNarrow_OI6363\tSecond_SII6731\tSecond_SII6716\tSecond_OI6300\tSecond_OI6363\tThird_SII6731\tThird_SII6716\tThird_OI6300\tThird_OI6363'))
	
	################################################ PLOT ######################################################
	plt.close()
	# MAIN plot
	fig1   = plt.figure(1,figsize=(10, 9))
	frame1 = fig1.add_axes((.1,.25,.85,.65)) 	     # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
	plt.plot(l,data_cor,'k')			     # Initial data
	plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='y',linestyle='-.',label='Linear fit')
	plt.plot(l,gaus1,c='darkgreen',linestyle='-')
	plt.plot(l,gaus2,c='darkgreen',linestyle='-',label='Narrow S')
	plt.plot(l,gaus3,c='goldenrod',linestyle='-')
	plt.plot(l,gaus4,c='goldenrod',linestyle='-',label='Narrow O')
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	frame1.text(6350.,max(data_cor), textstr, fontsize=12,verticalalignment='top', bbox=props)
	plt.plot(l[std0:std1],data_cor[std0:std1],c='darkorange')	# Zone where the stddev is calculated
	plt.plot(l,onefin_fit,'r-')

	frame1.set_xticklabels([]) 			# Remove x-tic labels for the first frame
	plt.ylabel('Flux (A.U.)',fontsize=14)#($\mathrm{erg/s/cm^{2} / \AA}$)'
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
	plt.plot(l,gaus22,c='darkgreen',linestyle='-',label='Narrow S')
	plt.plot(l,gaus23,c='goldenrod',linestyle='-')
	plt.plot(l,gaus24,c='goldenrod',linestyle='-',label='Narrow O')
	plt.plot(l,gaus25,c='dodgerblue',linestyle='-')
	plt.plot(l,gaus26,c='dodgerblue',linestyle='-',label='Second S')
	plt.plot(l,gaus27,c='darkred',linestyle='-')
	plt.plot(l,gaus28,c='darkred',linestyle='-',label='Second O')
	props = dict(boxstyle='round', facecolor='white', alpha=0.5)
	frame3.text(6350.,max(data_cor), textstr2, fontsize=12,verticalalignment='top', bbox=props)
	plt.plot(l[std0:std1],data_cor[std0:std1],c='darkorange')	# Zone where the stddev is calculated
	
	frame3.set_xticklabels([]) 			# Remove x-tic labels for the first frame
	plt.ylabel('Flux (A.U.)',fontsize=14) #($\mathrm{erg/s/cm^{2} / \AA}$)'
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
        plt.plot(l,gaus32,c='darkgreen',linestyle='-',label='Narrow S')
	plt.plot(l,gaus33,c='goldenrod',linestyle='-')
	plt.plot(l,gaus34,c='goldenrod',linestyle='-',label='Narrow O')
        plt.plot(l,gaus35,c='dodgerblue',linestyle='-')
        plt.plot(l,gaus36,c='dodgerblue',linestyle='-',label='Second S')
	plt.plot(l,gaus37,c='darkred',linestyle='-')
	plt.plot(l,gaus38,c='darkred',linestyle='-',label='Second O')
        plt.plot(l,gaus39,c='magenta',linestyle='-')
        plt.plot(l,gaus310,c='magenta',linestyle='-',label='Third S')
        plt.plot(l,gaus311,c='springgreen',linestyle='-')
        plt.plot(l,gaus312,c='springgreen',linestyle='-',label='Third O')
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        frame5.text(6350.,max(data_cor), textstr3, fontsize=12,verticalalignment='top', bbox=props)
        plt.plot(l[std0:std1],data_cor[std0:std1],c='darkorange')       # Zone where the stddev is calculated

        frame5.set_xticklabels([])                      # Remove x-tic labels for the first frame
        plt.ylabel('Flux (A.U.)',fontsize=14) #($\mathrm{erg/s/cm^{2} / \AA}$)'
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
