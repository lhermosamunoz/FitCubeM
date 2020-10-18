'''
This script makes a gaussian fit to the emission lines of AGN spectra
It is needed a parentFold, the spectrum in which the fit is going to be made and the initial estimation of the fit
'''

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
import pickle

from plot_broad import broad_plot
from DefRegionLines import DefRegionLines
from plot_refer_lines import refer_plot
from FitThirdComp import FitThirdComp
from FirstFittingFunction import FirstFittingFunction

fordername = 'results_nolim_SII_1610/'


# Conf
CubePath  = "/mnt/data/lhermosa/MEGARA/DATA/NGC1052/MixCubeRV/PPXF/FinalSubtract/"
cubePath  = "/mnt/data/lhermosa/MEGARA/DATA/NGC1052/MixCubeRV/EmissionLines/"
if not os.path.exists(cubePath+foldername): os.mkdir(cubePath+foldername)
path = cubePath+foldername #_2comp1broad_lim300/'
cubo = 'cubo_NGC1052'
z = 0.005

# Rest values of the line wavelengths 
l_Halpha = pylines.optical.lines['H_alpha'][0]
l_NII_1  = pylines.optical.lines['NIIa'][0]     # 6548.05
l_NII_2  = pylines.optical.lines['NIIb'][0]     # 6583.45
l_SII_1  = pylines.optical.lines['SIIa'][0]     # 6716.44
l_SII_2  = pylines.optical.lines['SIIb'][0]     # 6730.82
l_OI_1 = pylines.optical.lines['OI'][0] # 6300.304
l_OI_2 = 6363.77

# Constants and parameters
v_luz = c.value/10**3 # km/s
plate_scale = 0.2  #arcsec per pix
ang_to_pix = 0.28   # Angstrom per pixel - LR-R = 0.28 AA/pix
fwhm = 2*np.sqrt(2*np.log(2)) # por sigma
sigma_per_arc = 82  # km/s per arcsec derived from the FWHM per arcsec
pix_to_v = 50.      # km/s per arcsec BASADO EN R = c/deltaV

sig_inst = 0.654   # A  delta lambda (FWHM) = 3.6 A
sig_inicial = 1.0

meth = 'S'#input('Which method to be applied? ("S"/"O", not "M1"/"M2"): ')  # Method to fit

GetParentFold = lambda path,x,y : path+str(galaxy)+'_'+str(galaxy2)+'_results/'

'''
# Select all the files (sustitute of the previous lines)
cube = fits.open(cubePath+cubo)
datos = cube[0].data
header = cube[0].header
init_lambda = header['CRVAL3']
fin_lambda = header['CRVAL3']+(len(datos)*header['CDELT3'])
step = header['CDELT3']

cube.close()
# Define the wavelength range
l_init = np.linspace(init_lambda,fin_lambda,len(datos))
gal_i,x_i = np.shape(datos[0])
'''
listDir = os.listdir(CubePath)
headName = listDir[0].split('_')[0]

x = np.arange(0,33,1)
gal = np.arange(0,30,1)
l_init = np.genfromtxt(CubePath+listDir[0])[:,0]

def CreateFolders():
  datos = np.zeros((len(l_init),30,33))
  for i in x:
    for ia in gal:
      try: 
          DatProv = np.genfromtxt(CubePath+headName+str(i)+'_'+str(ia)+'.txt')[:,1]
      except IOError: 
          DatProv = np.zeros(len(l_init))
      datos[:,ia,i] = DatProv
      if not os.path.exists(path+str(i)+'_'+str(ia)+'_results'):
          os.mkdir(path+str(i)+'_'+str(ia)+'_results')
          print('Folder '+path+str(i)+'_'+str(ia)+'_results created in path')
  return datos

datos = CreateFolders

class FitOutput:
  def __init__(self, x, y):
    self.x = x
    self.y = y
    self.prov_SN,self.prov_eSN,self.prov_SN2,self.prov_eSN2,self.prov_SS,self.prov_eSS,self.prov_SB,self.prov_eSB = 0,0,0,0,0,0,0,0
    self.prov_notSN,self.prov_noteSN,self.prov_notVN,self.prov_noteVN = 0,0,0,0
    self.prov_VS,self.prov_eVS,self.prov_VN,self.prov_eVN,self.prov_VN2,self.prov_eVN2,self.prov_VB,self.prov_eVB = 0,0,0,0,0,0,0,0
    self.prov_fHa,self.prov_fN1,self.prov_fN2,self.prov_fS1,self.prov_fS2 = 0,0,0,0,0
    self.prov_f2Ha,self.prov_f2N1,self.prov_f2N2,self.prov_f2S1,self.prov_f2S2 = 0,0,0,0,0
    self.prov_fSHa,self.prov_fSN1,self.prov_fSN2,self.prov_fSS1,self.prov_fSS2 = 0,0,0,0,0
    self.prov_fO1,self.prov_fO2,self.prov_f2O1,self.prov_f2O2,self.prov_fSO1,self.prov_fSO2 = 0,0,0,0,0,0
    self.prov_e1,self.prov_e2,self.prov_e3,self.prov_a1,self.prov_a2,self.prov_a3,self.prov_b1,self.prov_b2,self.prov_b3 = 0,0,0,0,0,0,0,0,0

def CheckFit(x, y):
  name = 'fit_%i_%i.p'%(x,y)
  parentFold = GetParentFold(path, x, y)
  return os.path.isfile(parentFold+name)
