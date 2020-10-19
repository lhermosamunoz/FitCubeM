from config import *
import sys

checkAll = True
n = 0
for y in np.arange(0,30,1):
  for x in np.arange(0,33,1):
    if not CheckFit(x,y): 
        print('Missing fit [x, y] = [%i, %i]'%(x,y))
        checkAll = False
        n += 1

if not checkAll and n > 8: sys.exit()


sig_narrow,sig_narrow2,sig_second,sig_broad = [],[],[],[]
esig_narrow,esig_narrow2,esig_second,esig_broad = [],[],[],[]
vel_narrow,vel_narrow2,vel_second,vel_broad = [],[],[],[]
evel_narrow,evel_narrow2,evel_second,evel_broad = [],[],[],[]
vel_narrow_NoCor,evel_narrow_NoCor,sig_narrow_NoCor,esig_narrow_NoCor = [],[],[],[]
flux_Ha,flux_NII1,flux_NII2,flux_SII1,flux_SII2 = [],[],[],[],[]
flux2_Ha,flux2_NII1,flux2_NII2,flux2_SII1,flux2_SII2 = [],[],[],[],[]
flux_second_Ha,flux_second_NII1,flux_second_NII2,flux_second_SII1,flux_second_SII2 = [],[],[],[],[]
flux_OI1,flux_OI2,flux2_OI1,flux2_OI2,flux_second_OI1,flux_second_OI2 = [],[],[],[],[],[]
epsilon1,epsilon2,epsilon3,AIC1,AIC2,AIC3,BIC1,BIC2,BIC3 = [],[],[],[],[],[],[],[],[]

for y in np.arange(0,30,1):
  prov_SN,prov_eSN,prov_SN2,prov_eSN2,prov_SS,prov_eSS,prov_SB,prov_eSB = [],[],[],[],[],[],[],[]
  prov_notSN,prov_noteSN,prov_notVN,prov_noteVN = [],[],[],[]
  prov_VS,prov_eVS,prov_VN,prov_eVN,prov_VN2,prov_eVN2,prov_VB,prov_eVB = [],[],[],[],[],[],[],[]
  prov_fHa,prov_fN1,prov_fN2,prov_fS1,prov_fS2 = [],[],[],[],[]
  prov_f2Ha,prov_f2N1,prov_f2N2,prov_f2S1,prov_f2S2 = [],[],[],[],[]
  prov_fSHa,prov_fSN1,prov_fSN2,prov_fSS1,prov_fSS2 = [],[],[],[],[]
  prov_fO1,prov_fO2,prov_f2O1,prov_f2O2,prov_fSO1,prov_fSO2 = [],[],[],[],[],[]
  prov_e1,prov_e2,prov_e3,prov_a1,prov_a2,prov_a3,prov_b1,prov_b2,prov_b3 = [],[],[],[],[],[],[],[],[]
  for x in np.arange(0,33,1):
    name = 'fit_%i_%i.p'%(x,y)
    parentFold = GetParentFold(path, x, y)
    try: 
        output = pickle.load(open(parentFold+name,'rb'))
        prov_SN.append(output.prov_SN),prov_eSN.append(output.prov_eSN),prov_SN2.append(output.prov_SN2),prov_eSN2.append(output.prov_eSN2),prov_SS.append(output.prov_SS),prov_eSS.append(output.prov_eSS),prov_SB.append(output.prov_SB),prov_eSB.append(output.prov_eSB),prov_notSN.append(output.prov_notSN),prov_noteSN.append(output.prov_noteSN),prov_notVN.append(output.prov_notVN),prov_noteVN.append(output.prov_noteVN),prov_VS.append(output.prov_VS),prov_eVS.append(output.prov_eVS),prov_VN.append(output.prov_VN),prov_eVN.append(output.prov_eVN),prov_VN2.append(output.prov_VN2),prov_eVN2.append(output.prov_eVN2),prov_VB.append(output.prov_VB),prov_eVB.append(output.prov_eVB),prov_fHa.append(output.prov_fHa),prov_fN1.append(output.prov_fN1),prov_fN2.append(output.prov_fN2),prov_fS1.append(output.prov_fS1),prov_fS2.append(output.prov_fS2),prov_f2Ha.append(output.prov_f2Ha),prov_f2N1.append(output.prov_f2N1),prov_f2N2.append(output.prov_f2N2),prov_f2S1.append(output.prov_f2S1),prov_f2S2.append(output.prov_f2S2),prov_fSHa.append(output.prov_fSHa),prov_fSN1.append(output.prov_fSN1),prov_fSN2.append(output.prov_fSN2),prov_fSS1.append(output.prov_fSS1),prov_fSS2.append(output.prov_fSS2),prov_fO1.append(output.prov_fO1),prov_fO2.append(output.prov_fO2),prov_f2O1.append(output.prov_f2O1),prov_f2O2.append(output.prov_f2O2),prov_fSO1.append(output.prov_fSO1),prov_fSO2.append(output.prov_fSO2),prov_e1.append(output.prov_e1),prov_e2.append(output.prov_e2),prov_e3.append(output.prov_e3),prov_a1.append(output.prov_a1),prov_a2.append(output.prov_a2),prov_a3.append(output.prov_a3),prov_b1.append(output.prov_b1),prov_b2.append(output.prov_b2),prov_b3.append(output.prov_b3)
    except IOError:
        prov_SN.append(np.nan),prov_eSN.append(np.nan),prov_SN2.append(np.nan),prov_eSN2.append(np.nan),prov_SS.append(np.nan),prov_eSS.append(np.nan),prov_SB.append(np.nan),prov_eSB.append(np.nan),prov_notSN.append(np.nan),prov_noteSN.append(np.nan),prov_notVN.append(np.nan),prov_noteVN.append(np.nan),prov_VS.append(np.nan),prov_eVS.append(np.nan),prov_VN.append(np.nan),prov_eVN.append(np.nan),prov_VN2.append(np.nan),prov_eVN2.append(np.nan),prov_VB.append(np.nan),prov_eVB.append(np.nan),prov_fHa.append(np.nan),prov_fN1.append(np.nan),prov_fN2.append(np.nan),prov_fS1.append(np.nan),prov_fS2.append(np.nan),prov_f2Ha.append(np.nan),prov_f2N1.append(np.nan),prov_f2N2.append(np.nan),prov_f2S1.append(np.nan),prov_f2S2.append(np.nan),prov_fSHa.append(np.nan),prov_fSN1.append(np.nan),prov_fSN2.append(np.nan),prov_fSS1.append(np.nan),prov_fSS2.append(np.nan),prov_fO1.append(np.nan),prov_fO2.append(np.nan),prov_f2O1.append(np.nan),prov_f2O2.append(np.nan),prov_fSO1.append(np.nan),prov_fSO2.append(np.nan),prov_e1.append(np.nan),prov_e2.append(np.nan),prov_e3.append(np.nan),prov_a1.append(np.nan),prov_a2.append(np.nan),prov_a3.append(np.nan),prov_b1.append(np.nan),prov_b2.append(np.nan),prov_b3.append(np.nan)
        

  vel_narrow.append(prov_VN),evel_narrow.append(prov_eVN),vel_second.append(prov_VS),evel_second.append(prov_eVS),vel_broad.append(prov_VB),evel_broad.append(prov_eVB),sig_narrow.append(prov_SN),esig_narrow.append(prov_eSN),sig_second.append(prov_SS),esig_second.append(prov_eSS),sig_broad.append(prov_SB),esig_broad.append(prov_eSB),vel_narrow_NoCor.append(prov_notVN),evel_narrow_NoCor.append(prov_noteVN),sig_narrow_NoCor.append(prov_notSN),esig_narrow_NoCor.append(prov_noteSN),flux_Ha.append(prov_fHa),flux_SII1.append(prov_fS1),flux_SII2.append(prov_fS2),flux_NII1.append(prov_fN1),flux_NII2.append(prov_fN2),flux_OI1.append(prov_fO1),flux_OI2.append(prov_fO2),flux_second_Ha.append(prov_fSHa),flux_second_SII1.append(prov_fSS1),flux_second_SII2.append(prov_fSN2),flux_second_NII1.append(prov_fSN1),flux_second_NII2.append(prov_fSN2),flux_second_OI1.append(prov_fSO1),flux_second_OI2.append(prov_fSO2),epsilon1.append(prov_e1),epsilon2.append(prov_e2),epsilon3.append(prov_e3),AIC1.append(prov_a1),AIC2.append(prov_a2),AIC3.append(prov_a3),BIC1.append(prov_b1),BIC2.append(prov_b2),BIC3.append(prov_b3)

# Now we save a txt for each of the components (v and sigma)
# It would be ordered as: each line [i,:], each column [:,i]
# Watch out which are the real positions of each vector in order 
# to create the final maps!!!
np.savetxt(path+'epsilon1c.txt',np.matrix(epsilon1),fmt='%.6f')
np.savetxt(path+'AIC1c.txt',np.matrix(AIC1),fmt='%.6f')
np.savetxt(path+'BIC1c.txt',np.matrix(BIC1),fmt='%.6f')
np.savetxt(path+'epsilon2c.txt',np.matrix(epsilon2),fmt='%.6f')
np.savetxt(path+'AIC2c.txt',np.matrix(AIC2),fmt='%.6f')
np.savetxt(path+'BIC2c.txt',np.matrix(BIC2),fmt='%.6f')
np.savetxt(path+'epsilon3c.txt',np.matrix(epsilon3),fmt='%.6f')
np.savetxt(path+'AIC3c.txt',np.matrix(AIC3),fmt='%.6f')
np.savetxt(path+'BIC3c.txt',np.matrix(BIC3),fmt='%.6f')

np.savetxt(path+'velnarrow.txt',np.matrix(vel_narrow),fmt='%.4f')
np.savetxt(path+'evelnarrow.txt',np.matrix(evel_narrow),fmt='%.4f')
#np.savetxt(path+'velnarrow_notCor.txt',np.matrix(vel_narrow_NoCor),fmt='%.8f')
np.savetxt(path+'signarrow.txt',np.matrix(sig_narrow),fmt='%.4f')
np.savetxt(path+'esignarrow.txt',np.matrix(esig_narrow),fmt='%.4f')
#np.savetxt(path+'signarrow_notCor.txt',np.matrix(sig_narrow_NoCor),fmt='%.8f')

np.savetxt(path+'velsecond.txt',np.matrix(vel_second),fmt='%.4f')
np.savetxt(path+'evelsecond.txt',np.matrix(evel_second),fmt='%.4f')
np.savetxt(path+'sigsecond.txt',np.matrix(sig_second),fmt='%.4f')
np.savetxt(path+'esigsecond.txt',np.matrix(esig_second),fmt='%.4f')
np.savetxt(path+'velbroad.txt',np.matrix(vel_broad),fmt='%.4f')
np.savetxt(path+'evelbroad.txt',np.matrix(evel_broad),fmt='%.4f')
np.savetxt(path+'sigbroad.txt',np.matrix(sig_broad),fmt='%.4f')
np.savetxt(path+'esigbroad.txt',np.matrix(esig_broad),fmt='%.4f')
np.savetxt(path+'flux_Halpha.txt',np.matrix(flux_Ha),fmt='%.8f')
np.savetxt(path+'flux_NII6584.txt',np.matrix(flux_NII2),fmt='%.8f')
np.savetxt(path+'flux_NII6548.txt',np.matrix(flux_NII1),fmt='%.8f')
np.savetxt(path+'flux_SII6730.txt',np.matrix(flux_SII2),fmt='%.8f')
np.savetxt(path+'flux_SII6716.txt',np.matrix(flux_SII1),fmt='%.8f')
np.savetxt(path+'flux_OI6300.txt',np.matrix(flux_OI1),fmt='%.8f')
np.savetxt(path+'flux_OI6363.txt',np.matrix(flux_OI2),fmt='%.8f')
np.savetxt(path+'flux_second_Halpha.txt',np.matrix(flux_second_Ha),fmt='%.8f')
np.savetxt(path+'flux_second_NII6584.txt',np.matrix(flux_second_NII2),fmt='%.8f')
np.savetxt(path+'flux_second_NII6548.txt',np.matrix(flux_second_NII1),fmt='%.8f')
np.savetxt(path+'flux_second_SII6730.txt',np.matrix(flux_second_SII2),fmt='%.8f')
np.savetxt(path+'flux_second_SII6716.txt',np.matrix(flux_second_SII1),fmt='%.8f')
np.savetxt(path+'flux_second_OI6300.txt',np.matrix(flux_second_OI1),fmt='%.8f')
np.savetxt(path+'flux_second_OI6363.txt',np.matrix(flux_second_OI2),fmt='%.8f')

# Plot the maps
plt.figure(figsize=(18,8))

ax = plt.subplot(131)
masked_array = np.ma.array(vel_narrow,mask=np.isnan(vel_narrow))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax.invert_yaxis()
plt.title('V (km/s)')

ax1 = plt.subplot(132)
masked_array = np.ma.array(sig_narrow,mask=np.isnan(sig_narrow))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax1)
ax1.invert_yaxis()
plt.title('$\sigma$ (km/s)')

ax2 = plt.subplot(133)
masked_array = np.ma.array(flux_Ha,mask=np.isnan(flux_Ha))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax2.invert_yaxis()
plt.title(r'Flux$_{H\alpha}$')

plt.savefig(path+'maps_'+str(cubo)+'.png',bbox_inches='tight',pad_inches=0.2)

#######################################################################
# Plot the flux maps 
#######################################################################
plt.figure(figsize=(18,8))
ax = plt.subplot(231)
masked_array = np.ma.array(flux_Ha,mask=np.isnan(flux_Ha))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax.invert_yaxis()
plt.title(r'Flux$_{H\alpha}$')

ax1 = plt.subplot(232)
masked_array = np.ma.array(flux_NII2,mask=np.isnan(flux_NII2))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax1)
ax1.invert_yaxis()
plt.title(r'Flux$_{[NII]6584}$')

ax2 = plt.subplot(233)
masked_array = np.ma.array(flux_NII1,mask=np.isnan(flux_NII1))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax2.invert_yaxis()
plt.title(r'Flux$_{[NII]6548}$')

ax3 = plt.subplot(234)
masked_array = np.ma.array(flux_SII2,mask=np.isnan(flux_SII2))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax3.invert_yaxis()
plt.title(r'Flux$_{[SII]6731}$')

ax4 = plt.subplot(235)
masked_array = np.ma.array(flux_SII1,mask=np.isnan(flux_SII1))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax4)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax4.invert_yaxis()
plt.title(r'Flux$_{[SII]6716}$')

ax6 = plt.subplot(236)
flux_NII_Ha = np.log10(np.array(flux_NII2)/np.array(flux_Ha))
masked_array = np.ma.array(flux_NII_Ha,mask=np.isnan(flux_NII_Ha))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax6)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax6.invert_yaxis()
plt.title(r'log($[NII]/H_{\alpha}$)')

plt.savefig(path+'maps_flux_allLines_'+str(cubo)+'.png',bbox_inches='tight',pad_inches=0.2)
plt.savefig(path+'maps_flux_allLines_'+str(cubo)+'.pdf',bbox_inches='tight',pad_inches=0.2)

##########################################
# Secondary component
plt.figure(figsize=(18,8))

ax = plt.subplot(131)
masked_array = np.ma.array(vel_second,mask=np.isnan(vel_second))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax.ticklabel_format(fontsize='large')
ax.invert_yaxis()
plt.title('V (km/s)')

ax1 = plt.subplot(132)
masked_array = np.ma.array(sig_second,mask=np.isnan(sig_second))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax1)
ax1.ticklabel_format(fontsize='large')
ax1.invert_yaxis()
plt.title('$\sigma$ (km/s)')

ax2 = plt.subplot(133)
masked_array = np.ma.array(flux_second_Ha,mask=np.isnan(flux_second_Ha))
im = plt.imshow(masked_array,cmap='jet')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax=cax)
ax2.invert_yaxis()
plt.title(r'Flux$_{H\alpha}$')

plt.savefig(path+'maps_notCor_second_'+str(cubo)+'.png',bbox_inches='tight',pad_inches=0.2)
