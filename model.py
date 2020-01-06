import numpy as np


def beta_model(So, r, rc, beta):
    ''' Beta model with 3 parameters. 
    So -- normalization
    rc -- core radius
    beta -- powerlaw slope

    returns the flux (Cnts/pixel)

    '''
    
    return So * ( 1.0 + (r / rc)**2)**(-3.0 * beta + 0.5)

def integ_beta_model(r, rc, beta, r_start=0):
    ''' 2pi*r integral of the above Beta model with 3 parameters. 
    r -- r to integrate to
    rc -- core radius
    beta -- powerlaw slope
    r_start -- starting point for integration
    
    This equation comes from Moretti et al. 2005. Their eq. 2.
    
    It's been updated with our beta.
    
    '''

    rc2 = rc * rc
    
    beta = -3.0 * beta + 0.5
    
    area1 = np.pi * rc2 / (1 - beta) * ((1 + r_start**2 / rc2)**(1 - beta) - 1)
    area2 = np.pi * rc2 / (1 - beta) * ((1 + r**2 / rc2)**(1 - beta) - 1)
    
    return area2 - area1

def inv_beta_model(So, rc, beta, bkg):
    ''' Solves for where (radius) the beta model is equal to the background

    Make sure you keep track of the units going into all this. You wanna make 
    sure the radius comes out as something that makes sense.

    '''

    a = -3 * beta + 0.5

    b = (bkg / So)**(1 / a) - 1

    c = b * rc**2

    return np.sqrt(c)

def chi2(model, y, y_err):
    '''Chi error. We are going to square this to make it the chi2 error.'''
    return np.sum(((model - y) / y_err)**2)

def like(theta, x, y, yerr):
    So, rc, beta, bg = theta
    model = beta_model(So, x, rc, beta) + bg
    #return -0.5 * np.sum(np.log(2 * np.pi * yerr) + (y - model)**2 / yerr)
    return -chi2(model, y, yerr)

def prior(theta):
    So, rc, beta, bg = theta
    if So < 0:
        return -np.inf
    elif rc < 0 or rc > 10:
        return -np.inf
    elif beta < 0 or beta > 3: 
        return -np.inf
    elif bg < 0 or bg > 1:
        return -np.inf
    else:
        return 0.0

def prob(theta, x, y, yerr):
    lp = prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + like(theta, x, y, yerr)