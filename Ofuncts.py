import numpy as np

# Create a function to fit the data to a Gaussian given some initial values
def gaussian(x,mu,sigm,amp):
    '''
    Gaussian distribution
    
    x - values for the fit
    p[0]: mu - mean of the distribution --> wavelength for fitting the lines
    p[1]: sigma - stddev		--> standard deviation
    p[2]: amplitude			--> defined as the inverse of the square root of 2pi times the variance (sigma**2)
    '''
    return amp*np.exp(-(x-mu)**2/(2*sigm**2))

def linear(x,slope,intc):
    '''
    Linear equation
    '''
    y = slope*x + intc
    return y

# Function to create the gaussian and the linear one component fit for the reference lines only
def twogaussian(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1):
    '''
    Function to fit the reference lines to a gaussian + linear with only one component.
    The parameters to introduce have to be the initial guesses. 
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1)
    return y

# Function to create the gaussian and the linear one component fit for the reference lines only
def threegaussian(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_2,sig_2,amp_2):
    '''
    Function to fit the reference lines to a gaussian + linear with only one component.
    The parameters to introduce have to be the initial guesses. 
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_2,sig_2,amp_2)
    return y


# Function to create the gaussian and the linear two component fit for the reference lines only
def funcSII2comp(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_20,sig_20,amp_20,mu_21,sig_21,amp_21):
    '''
    Function to fit the reference lines to a gaussian + linear with two components.
    The parameters to introduce have to be the initial guesses. 
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_20,sig_20,amp_20) + gaussian(x,mu_21,sig_21,amp_21)
    return y

# Function to create the gaussian and the linear three component fit for the reference lines only
def funcSII3comp(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_20,sig_20,amp_20,mu_21,sig_21,amp_21,mu_30,sig_30,amp_30,mu_31,sig_31,amp_31):
    '''
    Function to fit the reference lines to a gaussian + linear with three components.
    The parameters to introduce have to be the initial guesses. 
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_20,sig_20,amp_20) + gaussian(x,mu_21,sig_21,amp_21) + gaussian(x,mu_30,sig_30,amp_30) + gaussian(x,mu_31,sig_31,amp_31)
    return y

# Function to create the gaussian and linear fit for all the lines with one component
def funcgauslin(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_2,sig_2,amp_2,mu_3,sig_3,amp_3,mu_4,sig_4,amp_4,mu_5,sig_5,amp_5,mu_6,sig_6,amp_6):
    '''
    Function to fit the spectra to a gaussian + linear.

    The parameters to introduce have to be the initial guesses for both components. 
    The first two values need to be the slope and the intercept, and then the rest 
    will be the parameters for fitting the gaussians.
    x - values for the fit
    params: The first two have to be the slope and the intercept of the linear fit
	    1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + slop*x + intc
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_2,sig_2,amp_2) + gaussian(x,mu_3,sig_3,amp_3) + gaussian(x,mu_4,sig_4,amp_4) + gaussian(x,mu_5,sig_5,amp_5) + gaussian(x,mu_6,sig_6,amp_6)
    return y

# Broad component of Halpha + full one component fit
def funcbroad(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_2,sig_2,amp_2,mu_3,sig_3,amp_3,mu_4,sig_4,amp_4,mu_5,sig_5,amp_5,mu_6,sig_6,amp_6,mu_b,sig_b,amp_b):
    '''
    Function to fit the spectra to a broad Halpha component.
    The parameters to introduce have to be the initial guesses. 
    It is necesary to have made the linear fit first
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_2,sig_2,amp_2) + gaussian(x,mu_3,sig_3,amp_3) + gaussian(x,mu_4,sig_4,amp_4) + gaussian(x,mu_5,sig_5,amp_5) + gaussian(x,mu_6,sig_6,amp_6) + gaussian(x,mu_b,sig_b,amp_b)
    return y

# Gaussian and linear fit for all the lines including a second component
def func2com(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_2,sig_2,amp_2,mu_3,sig_3,amp_3,mu_4,sig_4,amp_4,mu_5,sig_5,amp_5,mu_6,sig_6,amp_6,mu_20,sig_20,amp_20,mu_21,sig_21,amp_21,mu_22,sig_22,amp_22,mu_23,sig_23,amp_23,mu_24,sig_24,amp_24,mu_25,sig_25,amp_25,mu_26,sig_26,amp_26):
    '''
    Function to fit the lines to a second component.
    The parameters to introduce have to be the initial guesses. 
    It is necesary to have made the linear fit first
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_2,sig_2,amp_2) + gaussian(x,mu_3,sig_3,amp_3) + gaussian(x,mu_4,sig_4,amp_4) + gaussian(x,mu_5,sig_5,amp_5) + gaussian(x,mu_6,sig_6,amp_6) + gaussian(x,mu_20,sig_20,amp_20) + gaussian(x,mu_21,sig_21,amp_21) + gaussian(x,mu_22,sig_22,amp_22) + gaussian(x,mu_23,sig_23,amp_23) + gaussian(x,mu_24,sig_24,amp_24) + gaussian(x,mu_25,sig_25,amp_25) + gaussian(x,mu_26,sig_26,amp_26)
    return y

# Second component + broad Halpha of the lines
def func2bcom(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_2,sig_2,amp_2,mu_3,sig_3,amp_3,mu_4,sig_4,amp_4,mu_5,sig_5,amp_5,mu_6,sig_6,amp_6,mu_20,sig_20,amp_20,mu_21,sig_21,amp_21,mu_22,sig_22,amp_22,mu_23,sig_23,amp_23,mu_24,sig_24,amp_24,mu_25,sig_25,amp_25,mu_26,sig_26,amp_26,mu_b,sig_b,amp_b):
    '''
    Function to fit the lines to a second component + a broad Halpha component.
    The parameters to introduce have to be the initial guesses. 
    It is necesary to have made the linear fit first
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_2,sig_2,amp_2) + gaussian(x,mu_3,sig_3,amp_3) + gaussian(x,mu_4,sig_4,amp_4) + gaussian(x,mu_5,sig_5,amp_5) + gaussian(x,mu_6,sig_6,amp_6) + gaussian(x,mu_20,sig_20,amp_20) + gaussian(x,mu_21,sig_21,amp_21) + gaussian(x,mu_22,sig_22,amp_22) + gaussian(x,mu_23,sig_23,amp_23) + gaussian(x,mu_24,sig_24,amp_24) + gaussian(x,mu_25,sig_25,amp_25) + gaussian(x,mu_26,sig_26,amp_26) + gaussian(x,mu_b,sig_b,amp_b)
    return y

# Gaussian and linear fit for all the lines including a narrow, a second narrow and an intermediate component
def func3com(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_2,sig_2,amp_2,mu_3,sig_3,amp_3,mu_4,sig_4,amp_4,mu_5,sig_5,amp_5,mu_6,sig_6,amp_6,mu_20,sig_20,amp_20,mu_21,sig_21,amp_21,mu_22,sig_22,amp_22,mu_23,sig_23,amp_23,mu_24,sig_24,amp_24,mu_25,sig_25,amp_25,mu_26,sig_26,amp_26,mu_30,sig_30,amp_30,mu_31,sig_31,amp_31,mu_32,sig_32,amp_32,mu_33,sig_33,amp_33,mu_34,sig_34,amp_34,mu_35,sig_35,amp_35,mu_36,sig_36,amp_36):
    '''
    Function to fit the lines to a second component.
    The parameters to introduce have to be the initial guesses. 
    It is necesary to have made the linear fit first
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_2,sig_2,amp_2) + gaussian(x,mu_3,sig_3,amp_3) + gaussian(x,mu_4,sig_4,amp_4) + gaussian(x,mu_5,sig_5,amp_5) + gaussian(x,mu_6,sig_6,amp_6) + gaussian(x,mu_20,sig_20,amp_20) + gaussian(x,mu_21,sig_21,amp_21) + gaussian(x,mu_22,sig_22,amp_22) + gaussian(x,mu_23,sig_23,amp_23) + gaussian(x,mu_24,sig_24,amp_24) + gaussian(x,mu_25,sig_25,amp_25) + gaussian(x,mu_26,sig_26,amp_26) + gaussian(x,mu_30,sig_30,amp_30) + gaussian(x,mu_31,sig_31,amp_31) + gaussian(x,mu_32,sig_32,amp_32) + gaussian(x,mu_33,sig_33,amp_33) + gaussian(x,mu_34,sig_34,amp_34) + gaussian(x,mu_35,sig_35,amp_35) + gaussian(x,mu_36,sig_36,amp_36)
    return y

# Gaussian and linear fit for all the lines including a narrow, a second narrow, an intermediate component and a broad component
def func3bcom(x,slop,intc,mu_0,sig_0,amp_0,mu_1,sig_1,amp_1,mu_2,sig_2,amp_2,mu_3,sig_3,amp_3,mu_4,sig_4,amp_4,mu_5,sig_5,amp_5,mu_6,sig_6,amp_6,mu_20,sig_20,amp_20,mu_21,sig_21,amp_21,mu_22,sig_22,amp_22,mu_23,sig_23,amp_23,mu_24,sig_24,amp_24,mu_25,sig_25,amp_25,mu_26,sig_26,amp_26,mu_30,sig_30,amp_30,mu_31,sig_31,amp_31,mu_32,sig_32,amp_32,mu_33,sig_33,amp_33,mu_34,sig_34,amp_34,mu_35,sig_35,amp_35,mu_36,sig_36,amp_36,mu_b,sig_b,amp_b):
    '''
    Function to fit the lines to a second component.
    The parameters to introduce have to be the initial guesses. 
    It is necesary to have made the linear fit first
    x - values for the fit
    params: 1. mu - mean of the distribution
	    2. sigma - stddev
	    3. amplitude
    '''
    y = np.zeros_like(x)
    y = y + (slop*x+intc)
    y = y + gaussian(x,mu_0,sig_0,amp_0) + gaussian(x,mu_1,sig_1,amp_1) + gaussian(x,mu_2,sig_2,amp_2) + gaussian(x,mu_3,sig_3,amp_3) + gaussian(x,mu_4,sig_4,amp_4) + gaussian(x,mu_5,sig_5,amp_5) + gaussian(x,mu_6,sig_6,amp_6) + gaussian(x,mu_20,sig_20,amp_20) + gaussian(x,mu_21,sig_21,amp_21) + gaussian(x,mu_22,sig_22,amp_22) + gaussian(x,mu_23,sig_23,amp_23) + gaussian(x,mu_24,sig_24,amp_24) + gaussian(x,mu_25,sig_25,amp_25) + gaussian(x,mu_26,sig_26,amp_26) + gaussian(x,mu_30,sig_30,amp_30) + gaussian(x,mu_31,sig_31,amp_31) + gaussian(x,mu_32,sig_32,amp_32) + gaussian(x,mu_33,sig_33,amp_33) + gaussian(x,mu_34,sig_34,amp_34) + gaussian(x,mu_35,sig_35,amp_35) + gaussian(x,mu_36,sig_36,amp_36) + gaussian(x,mu_b,sig_b,amp_b)
    return y

