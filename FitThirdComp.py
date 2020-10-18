import numpy as np
import lmfit 
from pyspeckit import speclines as pylines
from astropy.constants import c
import Ofuncts
import os

def FitThirdComp(parentFold,l,data_cor,lin_data_fin,threeresu,meth,mu0,sig0,amp0,mu1,sig1,amp1,mu2,sig2,amp2,mu3,sig3,amp3,mu4,sig4,amp4,mu5,sig5,amp5,mu6,sig6,amp6,sig20,amp20,sig21,amp21,sig22,amp22,sig23,amp23,sig24,amp24,sig25,amp25,sig26,amp26,sl,it,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,stadev):
    # Constants to be used
    v_luz = c.value/10**3 # km/s
    pix_to_v = 50.      # km/s per arcsec BASADO EN R = c/deltaV

    l_Halpha = pylines.optical.lines['H_alpha'][0]
    l_OI_1 = pylines.optical.lines['OI'][0] # 6300.304
    l_SII_2  = pylines.optical.lines['SIIb'][0]     # 6730.82

    threecomp_mod = lmfit.Model(Ofuncts.func3com)
    params3c = lmfit.Parameters()

    if meth == 'S':
        cd1 = lmfit.Parameter('mu_0', value=threeresu.values["mu_0"],vary=False)
        de = lmfit.Parameter('sig_0', value=threeresu.values["sig_0"],vary=False)
        ef = lmfit.Parameter('amp_0', value=threeresu.values["amp_0"],vary=False)
        fg = lmfit.Parameter('mu_1', value=threeresu.values["mu_1"],vary=False)
        gh = lmfit.Parameter('sig_1', value=threeresu.values["sig_1"],vary=False)
        hi = lmfit.Parameter('amp_1', value=threeresu.values["amp_1"],vary=False)
        ij = lmfit.Parameter('mu_2', value=mu2,expr='mu_0*(6584./6731.)')
        jk = lmfit.Parameter('sig_2', value=sig2,expr='sig_0')
        kl = lmfit.Parameter('amp_2', value=amp2,min=0.)
        lm = lmfit.Parameter('mu_3', value=mu3,expr='mu_0*(6563./6731.)')
        mn = lmfit.Parameter('sig_3', value=sig3,expr='sig_0')
        no = lmfit.Parameter('amp_3', value=amp3,min=0.)
        op = lmfit.Parameter('mu_4', value=mu4,expr='mu_0*(6548./6731.)')
        pq = lmfit.Parameter('sig_4', value=sig4,expr='sig_0')
        qr = lmfit.Parameter('amp_4', value=amp4,min=0.,expr='amp_2*(1./3.)')
        rs = lmfit.Parameter('mu_5', value=mu5,expr='mu_0*(6300.304/6730.82)')
        st = lmfit.Parameter('sig_5', value=sig5,expr='sig_0')
        tu = lmfit.Parameter('amp_5', value=amp5,min=0.)
        uv = lmfit.Parameter('mu_6', value=mu6,expr='mu_0*(6363.77/6730.82)')
        vw = lmfit.Parameter('sig_6', value=sig6,expr='sig_0')
        wy = lmfit.Parameter('amp_6', value=amp6,min=0.,expr='amp_5*(1./3.)')
        aaa = lmfit.Parameter('mu_20', value=threeresu.values["mu_20"],vary=False)
        aab = lmfit.Parameter('sig_20', value=threeresu.values["sig_20"],vary=False)
        aac = lmfit.Parameter('amp_20', value=threeresu.values["amp_20"],vary=False)
        aad = lmfit.Parameter('mu_21', value=threeresu.values["mu_21"],vary=False)
        aae = lmfit.Parameter('sig_21', value=threeresu.values["sig_21"],vary=False)
        aaf = lmfit.Parameter('amp_21', value=threeresu.values["amp_21"],vary=False)
        aag = lmfit.Parameter('mu_22', value=mu2,expr='mu_20*(6584./6731.)')
        aah = lmfit.Parameter('sig_22', value=sig22,expr='sig_20')
        aai = lmfit.Parameter('amp_22', value=amp22,min=0.)
        aaj = lmfit.Parameter('mu_23', value=mu3,expr='mu_20*(6563./6731.)')
        aak = lmfit.Parameter('sig_23', value=sig23,expr='sig_20')
        aal = lmfit.Parameter('amp_23', value=amp23,min=0.)
        aam = lmfit.Parameter('mu_24', value=mu4,expr='mu_20*(6548./6731.)')
        aan = lmfit.Parameter('sig_24', value=sig24,expr='sig_20')
        aao = lmfit.Parameter('amp_24', value=amp24,min=0.,expr='amp_22*(1./3.)')
        aap = lmfit.Parameter('mu_25', value=mu5,expr='mu_20*(6300.304/6730.82)')
        aaq = lmfit.Parameter('sig_25', value=sig25,expr='sig_20')
        aar = lmfit.Parameter('amp_25', value=amp25,min=0.)
        aas = lmfit.Parameter('mu_26', value=mu6,expr='mu_20*(6363.77/6730.82)')
        aat = lmfit.Parameter('sig_26', value=sig26,expr='sig_20')
        aau = lmfit.Parameter('amp_26', value=amp26,min=0.,expr='amp_25*(1./3.)')
        aba = lmfit.Parameter('mu_30', value=threeresu.values["mu_20"],vary=False)
        abb = lmfit.Parameter('sig_30', value=threeresu.values["sig_20"],vary=False)
        abc = lmfit.Parameter('amp_30', value=threeresu.values["amp_20"],vary=False)
        abd = lmfit.Parameter('mu_31', value=threeresu.values["mu_21"],vary=False)
        abe = lmfit.Parameter('sig_31', value=threeresu.values["sig_21"],vary=False)
        abf = lmfit.Parameter('amp_31', value=threeresu.values["amp_21"],vary=False)
        abg = lmfit.Parameter('mu_32', value=mu2,expr='mu_30*(6584./6731.)')
        abh = lmfit.Parameter('sig_32', value=sig22,expr='sig_30')
        abi = lmfit.Parameter('amp_32', value=amp22,min=0.)
        abj = lmfit.Parameter('mu_33', value=mu3,expr='mu_30*(6563./6731.)')
        abk = lmfit.Parameter('sig_33', value=sig23,expr='sig_30')
        abl = lmfit.Parameter('amp_33', value=amp23,min=0.)
        abm = lmfit.Parameter('mu_34', value=mu4,expr='mu_30*(6548./6731.)')
        abn = lmfit.Parameter('sig_34', value=sig24,expr='sig_30')
        abo = lmfit.Parameter('amp_34', value=amp24,min=0.,expr='amp_32*(1./3.)')
        abp = lmfit.Parameter('mu_35', value=mu5,expr='mu_30*(6300.304/6730.82)')
        abq = lmfit.Parameter('sig_35', value=sig25,expr='sig_30')
        abr = lmfit.Parameter('amp_35', value=amp25,min=0.)
        abt = lmfit.Parameter('mu_36', value=mu6,expr='mu_30*(6363.77/6730.82)')
        abu = lmfit.Parameter('sig_36', value=sig26,expr='sig_30')
        abv = lmfit.Parameter('amp_36', value=amp26,min=0.,expr='amp_35*(1./3.)')
        params3c.add_many(sl,it,cd1,de,ef,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr,rs,st,tu,uv,vw,wy,aaa,aab,aac,aad,aae,aaf,aag,aah,aai,aaj,aak,aal,aam,aan,aao,aap,aaq,aar,aas,aat,aau,aba,abb,abc,abd,abe,abf,abg,abh,abi,abj,abk,abl,abm,abn,abo,abp,abq,abr,abt,abu,abv)


    elif meth == 'O':
        cd1 = lmfit.Parameter('mu_0', value=mu0,expr='mu_5*(6730.82/6300.30)')
        de = lmfit.Parameter('sig_0', value=sig0,expr='sig_5')
        ef = lmfit.Parameter('amp_0', value=amp0,min=0.)
        fg = lmfit.Parameter('mu_1', value=mu1,expr='mu_5*(6716.44/6300.30)')
        gh = lmfit.Parameter('sig_1', value=sig1,expr='sig_5')
        hi = lmfit.Parameter('amp_1', value=amp1,min=0.)
        ij = lmfit.Parameter('mu_2', value=mu2,expr='mu_5*(6584./6300.30)')
        jk = lmfit.Parameter('sig_2', value=sig2,expr='sig_5')
        kl = lmfit.Parameter('amp_2', value=amp2,min=0.)
        lm = lmfit.Parameter('mu_3', value=mu3,expr='mu_5*(6563./6300.30)')
        mn = lmfit.Parameter('sig_3', value=sig3,expr='sig_5')
        no = lmfit.Parameter('amp_3', value=amp3,min=0.)
        op = lmfit.Parameter('mu_4', value=mu4,expr='mu_5*(6548./6300.30)')
        pq = lmfit.Parameter('sig_4', value=sig4,expr='sig_5')
        qr = lmfit.Parameter('amp_4', value=amp4,min=0.,expr='amp_2*(1./3.)')
        rs = lmfit.Parameter('mu_5', value=threeresu.values["mu_0"],vary=False)#mu5,expr='mu_0*(6300.304/6730.82)')
        st = lmfit.Parameter('sig_5', value=threeresu.values["sig_0"],vary=False)
        tu = lmfit.Parameter('amp_5', value=threeresu.values["amp_0"],vary=False)
        uv = lmfit.Parameter('mu_6', value=threeresu.values["mu_1"],vary=False)#mu6,expr='mu_0*(6363.77/6730.82)')
        vw = lmfit.Parameter('sig_6', value=threeresu.values["sig_1"],vary=False)
        wy = lmfit.Parameter('amp_6', value=threeresu.values["amp_1"],vary=False)

        aaa = lmfit.Parameter('mu_20', value=mu0,expr='mu_25*(6730.82/6300.30)')
        aab = lmfit.Parameter('sig_20', value=sig20,expr='sig_25')
        aac = lmfit.Parameter('amp_20', value=amp20,min=0.)
        aad = lmfit.Parameter('mu_21', value=mu1,expr='mu_25*(6716.44/6300.30)')
        aae = lmfit.Parameter('sig_21', value=sig21,expr='sig_25')
        aaf = lmfit.Parameter('amp_21', value=amp21,min=0.)
        aag = lmfit.Parameter('mu_22', value=mu2,expr='mu_25*(6584./6731.)')
        aah = lmfit.Parameter('sig_22', value=sig22,expr='sig_25')
        aai = lmfit.Parameter('amp_22', value=amp22,min=0.)
        aaj = lmfit.Parameter('mu_23', value=mu3,expr='mu_25*(6563./6731.)')
        aak = lmfit.Parameter('sig_23', value=sig23,expr='sig_25')
        aal = lmfit.Parameter('amp_23', value=amp23,min=0.)
        aam = lmfit.Parameter('mu_24', value=mu4,expr='mu_20*(6548./6731.)')
        aan = lmfit.Parameter('sig_24', value=sig24,expr='sig_25')
        aao = lmfit.Parameter('amp_24', value=amp24,min=0.,expr='amp_22*(1./3.)')
        aap = lmfit.Parameter('mu_25', value=threeresu.values["mu_20"],vary=False)#mu5,expr='mu_20*(6300.304/6730.82)')
        aaq = lmfit.Parameter('sig_25', value=threeresu.values["sig_20"],vary=False)
        aar = lmfit.Parameter('amp_25', value=threeresu.values["amp_20"],vary=False)
        aas = lmfit.Parameter('mu_26', value=threeresu.values["mu_21"],vary=False)#mu6,expr='mu_20*(6363.77/6730.82)')
        aat = lmfit.Parameter('sig_26', value=threeresu.values["sig_21"],vary=False)
        aau = lmfit.Parameter('amp_26', value=threeresu.values["mu_21"],vary=False)

        aba = lmfit.Parameter('mu_30', value=mu0,expr='mu_35*(6730.82/6300.30)')
        abb = lmfit.Parameter('sig_30', value=sig30,expr='sig_35')
        abc = lmfit.Parameter('amp_30', value=amp30,min=0.)
        abd = lmfit.Parameter('mu_31', value=mu1,expr='mu_35*(6716.44/6300.30)')
        abe = lmfit.Parameter('sig_31', value=sig31,expr='sig_35')
        abf = lmfit.Parameter('amp_31', value=amp31,min=0.)
        abg = lmfit.Parameter('mu_32', value=mu2,expr='mu_30*(6584./6731.)')
        abh = lmfit.Parameter('sig_32', value=sig22,expr='sig_30')
        abi = lmfit.Parameter('amp_32', value=amp22,min=0.)
        abj = lmfit.Parameter('mu_33', value=mu3,expr='mu_30*(6563./6731.)')
        abk = lmfit.Parameter('sig_33', value=sig23,expr='sig_30')
        abl = lmfit.Parameter('amp_33', value=amp23,min=0.)
        abm = lmfit.Parameter('mu_34', value=mu4,expr='mu_30*(6548./6731.)')
        abn = lmfit.Parameter('sig_34', value=sig24,expr='sig_30')
        abo = lmfit.Parameter('amp_34', value=amp24,min=0.,expr='amp_32*(1./3.)')
        abp = lmfit.Parameter('mu_35', value=threeresu.values["mu_30"],vary=False)#mu5,expr='mu_30*(6300.304/6730.82)')
        abq = lmfit.Parameter('sig_35', value=threeresu.values["sig_30"],vary=False)
        abr = lmfit.Parameter('amp_35', value=threeresu.values["amp_30"],vary=False)
        abt = lmfit.Parameter('mu_36', value=threeresu.values["mu_31"],vary=False)#mu6,expr='mu_30*(6363.77/6730.82)')
        abu = lmfit.Parameter('sig_36', value=threeresu.values["sig_31"],vary=False)
        abv = lmfit.Parameter('amp_36', value=threeresu.values["amp_31"],vary=False)
        params3c.add_many(sl,it,rs,st,tu,uv,vw,wy,cd1,de,ef,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr,aap,aaq,aar,aas,aat,aau,aaa,aab,aac,aad,aae,aaf,aag,aah,aai,aaj,aak,aal,aam,aan,aao,abp,abq,abr,abt,abu,abv,aba,abb,abc,abd,abe,abf,abg,abh,abi,abj,abk,abl,abm,abn,abo)

    T3CResu = threecomp_mod.fit(data_cor,params=params3c,x=l)
    
    with open(parentFold+'fit3CAll_result.txt', 'w') as fh: fh.write(T3CResu.fit_report())

    ########################## Calculate gaussians and final fit ##########################
    # Now we create and plot the individual gaussians of the fit
    gaus1 = Ofuncts.gaussian(l,T3CResu.values['mu_0'],T3CResu.values['sig_0'],T3CResu.values['amp_0'])
    gaus2 = Ofuncts.gaussian(l,T3CResu.values['mu_1'],T3CResu.values['sig_1'],T3CResu.values['amp_1'])
    gaus3 = Ofuncts.gaussian(l,T3CResu.values['mu_2'],T3CResu.values['sig_2'],T3CResu.values['amp_2'])
    gaus4 = Ofuncts.gaussian(l,T3CResu.values['mu_3'],T3CResu.values['sig_3'],T3CResu.values['amp_3'])
    gaus5 = Ofuncts.gaussian(l,T3CResu.values['mu_4'],T3CResu.values['sig_4'],T3CResu.values['amp_4'])
    gaus6 = Ofuncts.gaussian(l,T3CResu.values['mu_5'],T3CResu.values['sig_5'],T3CResu.values['amp_5'])
    gaus7 = Ofuncts.gaussian(l,T3CResu.values['mu_6'],T3CResu.values['sig_6'],T3CResu.values['amp_6'])
    gaus21 = Ofuncts.gaussian(l,T3CResu.values['mu_20'],T3CResu.values['sig_20'],T3CResu.values['amp_20'])
    gaus22 = Ofuncts.gaussian(l,T3CResu.values['mu_21'],T3CResu.values['sig_21'],T3CResu.values['amp_21'])
    gaus23 = Ofuncts.gaussian(l,T3CResu.values['mu_22'],T3CResu.values['sig_22'],T3CResu.values['amp_22'])
    gaus24 = Ofuncts.gaussian(l,T3CResu.values['mu_23'],T3CResu.values['sig_23'],T3CResu.values['amp_23'])
    gaus25 = Ofuncts.gaussian(l,T3CResu.values['mu_24'],T3CResu.values['sig_24'],T3CResu.values['amp_24'])
    gaus26 = Ofuncts.gaussian(l,T3CResu.values['mu_25'],T3CResu.values['sig_25'],T3CResu.values['amp_25'])
    gaus27 = Ofuncts.gaussian(l,T3CResu.values['mu_26'],T3CResu.values['sig_26'],T3CResu.values['amp_26'])
    gaus31 = Ofuncts.gaussian(l,T3CResu.values['mu_30'],T3CResu.values['sig_30'],T3CResu.values['amp_30'])
    gaus32 = Ofuncts.gaussian(l,T3CResu.values['mu_31'],T3CResu.values['sig_31'],T3CResu.values['amp_31'])
    gaus33 = Ofuncts.gaussian(l,T3CResu.values['mu_32'],T3CResu.values['sig_32'],T3CResu.values['amp_32'])
    gaus34 = Ofuncts.gaussian(l,T3CResu.values['mu_33'],T3CResu.values['sig_33'],T3CResu.values['amp_33'])
    gaus35 = Ofuncts.gaussian(l,T3CResu.values['mu_34'],T3CResu.values['sig_34'],T3CResu.values['amp_34'])
    gaus36 = Ofuncts.gaussian(l,T3CResu.values['mu_35'],T3CResu.values['sig_35'],T3CResu.values['amp_35'])
    gaus37 = Ofuncts.gaussian(l,T3CResu.values['mu_36'],T3CResu.values['sig_36'],T3CResu.values['amp_36'])
    T3Cfin_fit = T3CResu.best_fit

    # one component
    stdf_s2 = np.std(data_cor[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10]-T3Cfin_fit[np.where(l<l1)[0][-1]:np.where(l>l2)[0][0]+10])
    stdf_s1 = np.std(data_cor[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]]-T3Cfin_fit[np.where(l<l3)[0][-1]-10:np.where(l>l4)[0][0]])
    stdf_n2 = np.std(data_cor[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10]-T3Cfin_fit[np.where(l<l5)[0][-1]:np.where(l>l6)[0][0]+10])
    stdf_ha = np.std(data_cor[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]]-T3Cfin_fit[np.where(l<l7)[0][-1]:np.where(l>l8)[0][0]])
    stdf_n1 = np.std(data_cor[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]]-T3Cfin_fit[np.where(l<l9)[0][-1]-10:np.where(l>l10)[0][0]])
    stdf_o1 = np.std(data_cor[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]]-T3Cfin_fit[np.where(l<l11)[0][-1]-10:np.where(l>l12)[0][0]])
    stdf_o2 = np.std(data_cor[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]]-T3Cfin_fit[np.where(l<l13)[0][-1]-10:np.where(l>l14)[0][0]])
    print('The condition for each line (in the same order as before) needs to be std_line < 3*std_cont --> for 1 components is... ')
    print('         For SII2: '+str(stdf_s2/stadev)+' < 3')
    print('         For SII1: '+str(stdf_s1/stadev)+' < 3')
    print('         For NII2: '+str(stdf_n2/stadev)+' < 3')
    print('         For Halpha: '+str(stdf_ha/stadev)+' < 3')
    print('         For NII1: '+str(stdf_n1/stadev)+' < 3')
    print('         For OI1: '+str(stdf_o1/stadev)+' < 3')
    print('         For OI2: '+str(stdf_o2/stadev)+' < 3')

    if os.path.exists(parentFold+'eps_adj'+str(meth)+'_3C.txt'): os.remove(parentFold+'eps_adj'+str(meth)+'_3C.txt')
    np.savetxt(parentFold+'eps_adj'+str(meth)+'_3C.txt',np.c_[stdf_s2/stadev,stdf_s1/stadev,stdf_n2/stadev,stdf_ha/stadev,stdf_n1/stadev,stdf_o2/stadev,stdf_o1/stadev,T3CResu.chisqr], ('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('SII2\tSII1\tNII2\tHa\tNII1\tOI1\tOI2\tChi2'))


    try:
        # We determine the maximum flux of the fit for all the lines, and the velocity and sigma components
        max2S1 = T3Cfin_fit[np.where(abs(T3CResu.values['mu_0']-l)<0.5)[0][0]]
        max2S2 = T3Cfin_fit[np.where(abs(T3CResu.values['mu_1']-l)<0.5)[0][0]]
        max2N1 = T3Cfin_fit[np.where(abs(T3CResu.values['mu_2']-l)<0.5)[0][0]]
        max2Ha = T3Cfin_fit[np.where(abs(T3CResu.values['mu_3']-l)<0.5)[0][0]]
        max2N2 = T3Cfin_fit[np.where(abs(T3CResu.values['mu_4']-l)<0.5)[0][0]]
        max2O1 = T3Cfin_fit[np.where(abs(T3CResu.values['mu_5']-l)<0.5)[0][0]]
        max2O2 = T3Cfin_fit[np.where(abs(T3CResu.values['mu_6']-l)<0.5)[0][0]]
    except IndexError:
        print('ERROR: index out of range. Setting the flux values of the OI 1 line to 0.')

    # Calculus of the velocity and sigma for the three components
    sig30S2 = pix_to_v*np.sqrt(T3CResu.values['sig_3']**2-sig_inst**2)
    sig31S2 = pix_to_v*np.sqrt(T3CResu.values['sig_23']**2-sig_inst**2)
    sig32S2 = pix_to_v*np.sqrt(T3CResu.values['sig_33']**2-sig_inst**2)
    
    v30S2 = v_luz*((T3CResu.values['mu_3']-l_Halpha)/l_Halpha)
    v31S2 = v_luz*((T3CResu.values['mu_23']-l_Halpha)/l_Halpha)
    v32S2 = v_luz*((T3CResu.values['mu_33']-l_Halpha)/l_Halpha)

    # Using SII lines as reference
    if meth == 'S':
        if threeresu.params['sig_0'].stderr == None:
            print('Problem determining the errors! First component sigma ')
            esig30S2 = 0.
        elif threeresu.params['sig_0'].stderr != None:
            esig30S2 = pix_to_v*(2*T3CResu.values['sig_0']*threeresu.params['sig_0'].stderr)/(np.sqrt(T3CResu.values['sig_0']**2-sig_inst**2))

        if threeresu.params['sig_20'].stderr == None:
            print('Problem determining the errors! Second component sigma ')
            esig31S2 = 0.
        elif threeresu.params['sig_20'].stderr != None:
            esig31S2 = pix_to_v*(2*T3CResu.values['sig_20']*threeresu.params['sig_20'].stderr)/(np.sqrt(T3CResu.values['sig_20']**2-sig_inst**2))

        if threeresu.params['sig_30'].stderr == None:
            print('Problem determining the errors! Second component sigma ')
            esig32S2 = 0.
        elif threeresu.params['sig_30'].stderr != None:
            esig32S2 = pix_to_v*(2*T3CResu.values['sig_30']*threeresu.params['sig_30'].stderr)/(np.sqrt(T3CResu.values['sig_30']**2-sig_inst**2))

        if threeresu.params['mu_0'].stderr == None:
            print('Problem determining the errors! First component ')
            ev30S2= 0.
        elif tworesu.params['mu_0'].stderr != None:
            print('Problem determining the errors! Second component ')
            ev30S2 = ((v_luz/l_SII_2)*T3CResu.params['mu_0'].stderr)

        if threeresu.params['mu_20'].stderr == None:
            ev31S2 = 0.
        elif threeresu.params['mu_20'].stderr != None:
            ev31S2 = ((v_luz/l_SII_2)*threeresu.params['mu_20'].stderr)

        if threeresu.params['mu_30'].stderr == None:
            ev32S2 = 0.
        elif threeresu.params['mu_30'].stderr != None:
            ev32S2 = ((v_luz/l_SII_2)*threeresu.params['mu_30'].stderr)

    # Using OI lines as reference
    elif meth == 'O':
        if threeresu.params['sig_5'].stderr == None:
            print('Problem determining the errors! First component sigma ')
            esig30S2 = 0.
        elif threeresu.params['sig_5'].stderr != None:
            esig30S2 = pix_to_v*(2*T3CResu.values['sig_5']*threeresu.params['sig_5'].stderr)/(np.sqrt(T3CResu.values['sig_5']**2-sig_inst**2))

        if threeresu.params['sig_25'].stderr == None:
            print('Problem determining the errors! Second component sigma ')
            esig31S2 = 0.
        elif threeresu.params['sig_25'].stderr != None:
            esig31S2 = pix_to_v*(2*T3CResu.values['sig_25']*threeresu.params['sig_25'].stderr)/(np.sqrt(T3CResu.values['sig_25']**2-sig_inst**2))

        if threeresu.params['sig_35'].stderr == None:
            print('Problem determining the errors! Third component sigma ')
            esig32S2 = 0.
        elif threeresu.params['sig_35'].stderr != None:
            esig32S2 = pix_to_v*(2*T3CResu.values['sig_35']*threeresu.params['sig_35'].stderr)/(np.sqrt(T3CResu.values['sig_35']**2-sig_inst**2))

        if threeresu.params['mu_5'].stderr == None:
            print('Problem determining the errors! First component ')
            ev30S2 = 0.
        elif tworesu.params['mu_5'].stderr != None:
            print('Problem determining the errors! Second component ')
            ev30S2 = ((v_luz/l_OI_1)*T3CResu.params['mu_5'].stderr)

        if threeresu.params['mu_25'].stderr == None:
            ev31S2 = 0.
        elif threeresu.params['mu_25'].stderr != None:
            ev31S2 = ((v_luz/l_OI_1)*threeresu.params['mu_25'].stderr)

        if threeresu.params['mu_35'].stderr == None:
            ev32S2 = 0.
        elif threeresu.params['mu_35'].stderr != None:
            ev32S2 = ((v_luz/l_OI_1)*threeresu.params['mu_35'].stderr)

    # Save the velocity and sigma for all components
    if os.path.exists(parentFold+'v_sig_adj_3C.txt'): os.remove(parentFold+'v_sig_adj_3C.txt')
    np.savetxt(parentFold+'v_sig_adj_3C.txt',np.c_[v30S2,ev30S2,v31S2,ev31S2,v32S2,ev32S2,sig30S2,esig30S2,sig31S2,esig31S2,sig32S2,esig32S2],('%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f','%8.5f'),header=('v_ref2\tev_ref2\tv_2ref2\tev_2ref2\tv_3ref2\tev_3ref2\tsig_ref2\tesig_ref2\tsig_2ref2\tesig_2ref2\tsig_3ref2\tesig_3ref2'))

    # Save the fluxes for all components
    if os.path.exists(parentFold+'fluxes_'+str(meth)+'_3C_Ncomp.txt'): os.remove(parentFold+'fluxes_'+str(meth)+'_3C_Ncomp.txt')
    np.savetxt(parentFold+'fluxes_3C_Ncomp.txt',np.c_[sum(gaus1),sum(gaus2),sum(gaus3),sum(gaus4),sum(gaus5),sum(gaus6),sum(gaus7)],fmt=('%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))
    if os.path.exists(parentFold+'fluxes_'+str(meth)+'_3C_Scomp.txt'): os.remove(parentFold+'fluxes_'+str(meth)+'_3C_Scomp.txt')
    np.savetxt(parentFold+'fluxes_3C_Scomp.txt',np.c_[sum(gaus11),sum(gaus12),sum(gaus13),sum(gaus14),sum(gaus15),sum(gaus16),sum(gaus17)],('%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))
    if os.path.exists(parentFold+'fluxes_'+str(meth)+'_3C_N2comp.txt'): os.remove(parentFold+'fluxes_'+str(meth)+'_3C_N2comp.txt')
    np.savetxt(parentFold+'fluxes_3C_N2comp.txt',np.c_[sum(gaus21),sum(gaus22),sum(gaus23),sum(gaus24),sum(gaus25),sum(gaus26),sum(gaus27)],('%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f','%8.16f'),header=('SII_6731\tSII_6716\tNII_6584\tHalpha\tNII_6548\tOI_6300\tOI_6363'))


    ########################### PLOT #############################
    plt.close('all')
    # MAIN plot
    fig1   = plt.figure(1,figsize=(10, 9))
    frame1 = fig1.add_axes((.1,.25,.85,.65))             # xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
    plt.plot(l,data_cor,'k',linewidth=2)                 # Initial data
    plt.plot(l[std0:std1],data_cor[std0:std1],c='y',linewidth=4)  # Zone where the stddev is calculated
    plt.plot(l[std0:std1],data_cor[std0:std1],'k',linewidth=1)                   # Initial data
    #plt.plot(l,(linresu.values['slope']*l+linresu.values['intc']),c='y',linestyle=(0, (5, 8)),label='Linear fit')
    plt.plot(l,gaus1+lin_data_fin,c='g',linestyle='-')
    plt.plot(l,gaus2+lin_data_fin,c='g',linestyle='-')
    plt.plot(l,gaus3+lin_data_fin,c='g',linestyle='-')
    plt.plot(l,gaus4+lin_data_fin,c='g',linestyle='-')
    plt.plot(l,gaus5+lin_data_fin,c='g',linestyle='-')
    plt.plot(l,gaus6+lin_data_fin,c='g',linestyle='-')
    plt.plot(l,gaus7+lin_data_fin,c='g',linestyle='-')
    plt.plot(l,gaus11+lin_data_fin,c='dodgerblue',linestyle='-')
    plt.plot(l,gaus12+lin_data_fin,c='dodgerblue',linestyle='-')
    plt.plot(l,gaus13+lin_data_fin,c='dodgerblue',linestyle='-')
    plt.plot(l,gaus14+lin_data_fin,c='dodgerblue',linestyle='-')
    plt.plot(l,gaus15+lin_data_fin,c='dodgerblue',linestyle='-')
    plt.plot(l,gaus16+lin_data_fin,c='dodgerblue',linestyle='-')
    plt.plot(l,gaus17+lin_data_fin,c='dodgerblue',linestyle='-')
    plt.plot(l,gaus21+lin_data_fin,c='magenta',linestyle='-')
    plt.plot(l,gaus22+lin_data_fin,c='magenta',linestyle='-')
    plt.plot(l,gaus23+lin_data_fin,c='magenta',linestyle='-')
    plt.plot(l,gaus24+lin_data_fin,c='magenta',linestyle='-')
    plt.plot(l,gaus25+lin_data_fin,c='magenta',linestyle='-')
    plt.plot(l,gaus26+lin_data_fin,c='magenta',linestyle='-')
    plt.plot(l,gaus27+lin_data_fin,c='magenta',linestyle='-')
    plt.plot(l,T3Cfin_fit,'r-')
    textstr = '\n'.join((r'$V_{SII_{3-1comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v30S2,ev30S2),
                         r'$V_{SII_{3-2comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v31S2,ev31S2),
                         r'$V_{SII_{3-3comp}}$ = '+ '{:.2f} +- {:.2f}'.format(v32S2,ev32S2),
                         r'$\sigma_{SII_{3-1comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig30S2,esig30S2),
                         r'$\sigma_{SII_{3-2comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig31S2,esig31S2),
                         r'$\sigma_{SII_{3-3comp}}$ = '+ '{:.2f} +- {:.2f}'.format(sig32S2,esig32S2),
                         #r'$\frac{F_{SII_{2}}}{F_{SII_{1}}}$ = '+ '{:.3f}'.format(max2S2/max2S1),
                         r'$F_{H_{\alpha}}$ = '+ '{:.3f}'.format(max2Ha)+' $10^{-14}$'))

    frame1.set_xticklabels([])                      # Remove x-tic labels for the first frame
    plt.ylabel(r'Flux (erg $\rm cm^{-2} s^{-1} \AA^{-1}$)',fontsize=19)
    plt.tick_params(axis='both', labelsize=17)
    plt.xlim(l[0],l[-1])
    plt.gca().yaxis.set_major_locator(plt.MaxNLocator(prune='lower'))
    plt.text(0.81,0.9,'N1',color='g',transform=frame1.transAxes,fontsize=21)
    plt.text(0.84,0.9,'+',color='k',transform=frame1.transAxes,fontsize=21)
    plt.text(0.87,0.9,'S',color='dodgerblue',transform=frame1.transAxes,fontsize=21)
    plt.text(0.9,0.9,'+',color='k',transform=frame1.transAxes,fontsize=21)
    plt.text(0.93,0.9,'N2',color='magenta',transform=frame1.transAxes,fontsize=21)

    # RESIDUAL plot
    frame2 = fig1.add_axes((.1,.1,.85,.15))
    plt.plot(l,np.zeros(len(l)),c='orange',linestyle='-')           # Line around zero
    plt.plot(l,data_cor-T3Cfin_fit,color='k')         # Main
    plt.xlabel('Wavelength ($\AA$)',fontsize=19)
    plt.ylabel('Residuals',fontsize=19)
    plt.tick_params(axis='both', labelsize=17)
    plt.xlim(l[0],l[-1])
    plt.plot(l,np.zeros(len(l))+1.5*stadev,c='orange',linestyle=(0,(5,8)))  # 3 sigma upper limit
    plt.plot(l,np.zeros(len(l))-1.5*stadev,c='orange',linestyle=(0,(5,8)))  # 3 sigma down limit
    plt.ylim(-(3*stadev)*3,(3*stadev)*3)

    plt.savefig(parentFold+'adj_full_3comp.pdf',format='pdf',bbox_inches='tight',pad_inches=0.2)
    props = dict(boxstyle='round',facecolor='white', alpha=0.5)
    frame1.text(6250.,max(data_cor),textstr,fontsize=12,verticalalignment='top', bbox=props)
    plt.savefig(parentFold+'adj_full_3comp.png',bbox_inches='tight',pad_inches=0.2)


    return gaus1,gaus2,gaus3,gaus4,gaus5,gaus6,gaus7,gaus11,gaus12,gaus13,gaus14,gaus15,gaus16,gaus17,gaus21,gaus22,gaus23,gaus24,gaus25,gaus26,gaus27,v30S2,ev30S2,v31S2,ev31S2,v32S2,ev32S2,sig30S2,esig30S2,sig31S2,esig31S2,sig30S2,esig30S2,T3CResu
