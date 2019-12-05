import numpy as np


def beta_model(So, r, rc, beta):
    ''' Beta model with 3 parameters. 
    So -- normalization
    rc -- core radius
    beta -- powerlaw slope

    returns the flux (Cnts/pixel)

    '''
    
    return So * ( 1.0 + (r / rc)**2)**(-3.0 * beta + 0.5)

def integ_beta_model(r, rc, beta):
    ''' 2pi*r integral of the above Beta model with 3 parameters. 
    r -- r
    rc -- core radius
    beta -- powerlaw slope

    '''

    rc2 = rc * rc
    
    return np.pi * rc2 / (1 - beta) * ((1 + r**2 / rc2)**(1 - beta) - 1)

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