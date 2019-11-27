#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#%matplotlib notebook
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column
from scipy.stats import ks_2samp
from scipy.stats import ksone

import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

# parallel processor
from utilities import parallel_process

# In[ ]:


def ks_critical_value(alpha, n_trials):
    ''' returns the critial value for a given alpha and number of trials'''
    return ksone.ppf(1 - alpha / 2, n_trials)


def ks_alpha(cv, n_trials):
    ''' returns the alpha for a given critical value and number of trials'''
    return ksone.sf(cv, n_trials) * 2


def test_ks_cv(alphas=None, n_trials=None):

    if not n_trials:
        n_trials = range(5, 15)
    if not alphas:
        alphas = [0.5, 0.1, 0.05, 0.02, 0.01]

    if isinstance(n_trials, int) or isinstance(n_trials, float):
        n_trials = range(n_trials - 2, n_trials + 2)

    # Print table headers
    print('{:<6}|{:<6} Level of significance, alpha'.format(' ', ' '))
    print('{:<6}|{:>8} {:>8} {:>8} {:>8} {:>8}'.format(*['Trials'] + alphas))
    print('-' * 42)
    # Print critical values for each n_trials x alpha combination
    for t in n_trials:
        print('{:6d}|{:>8.5f} {:>8.5f} {:>8.5f} {:>8.5f} {:>8.5f}'.format(
            *[t] + [ks_critical_value(a, t) for a in alphas]))


def ks_test_sources(name, outpath, plotting=False):

    # detections
    if not os.path.isfile(f'{outpath}/{name}/{name}_vtp.detect'):
        return
    else:
        srcs = f'{outpath}/{name}/{name}_vtp.detect'

    # mcmc fit data
    if not os.path.isfile(f'{outpath}/{name}/{name}_mcmcfits.txt'):
        return
    else:
        mcmcfits = f'{outpath}/{name}/{name}_mcmcfits.txt'
        fit = Table.read(mcmcfits, format='ascii', header_start=0)

    # now we need to read the individual detections
    detects = Table.read(srcs, hdu=1)

    # add a column for classifications!
    try:
        # add the columns to the detection catalog
        detects.add_column(
            Column(data=np.ones(len(detects)) * -1,
                   name='Extended',
                   dtype='>i8'))
    except ValueError:
        pass

    # loop over the sources -- only the first 3
    for indx, i in enumerate(detects['INDEX']):
        if os.path.isfile(f'{outpath}/{name}/{name}_vtp_{i}.radprof'):
            data = Table.read(f'{outpath}/{name}/{name}_vtp_{i}.radprof',
                              format='ascii',
                              header_start=2)
        else:
            continue

        # where is it in the mcmc file?
        loc = fit['ID'] == i
        if not loc.any():  # we didn't fit that one with mcmc
            continue
        bkg = fit['bg_50'][loc]  # cnts/s/arcmin2

        # convert to cnts/s
        pixscale = 2.36  # arcs/pix
        bkg *= pixscale**2 / 60**2  # cnts/s/arcmin2 * arcmin2/pix2
        source_exp = data['x'] / data['y'] - bkg * data[
            'w']  # cnts/s - cnts/s/pix2 * pix2
        psf_exp = data['psf'] / data['y']  # cnts/s

        # cummulative distributions
        source_cum = np.cumsum(source_exp)
        psf_cum = np.cumsum(psf_exp)

        # how far out do we want to go?
        bins = 12

        # normalized
        source_cum_n = source_cum / source_cum[bins - 1]
        psf_cum_n = psf_cum / psf_cum[bins - 1]

        cv = ks_critical_value(0.05, bins)

        #         print((psf_cum_n - source_cum_n)[:bins])

        #         print(i, (psf_cum_n[:bins] - source_cum_n[:bins]).max(), cv)

        if (psf_cum_n[:bins] - source_cum_n[:bins]).max() > cv:
            detects['Extended'][indx] = 1
#             print(i, 'NOT POINT SOURCE')
        else:
            detects['Extended'][indx] = 0
#             print(i, 'POINT SOURCE')

# plotting
        if plotting:
            # x-axis, in arcminutes
            x = (data['r1'] + data['r2']) / 2. / 60. * pixscale
            plt.plot(x[:bins], source_cum_n[:bins], label=i)

# ## if we want to try to use scipy's kstest
#         result = ks_2samp(source_cum_n[:bins], psf_cum_n[:bins])
#         print(i, (psf_cum_n - source_cum_n).max())
#         print(i, result)

    detects.write(srcs, format='fits', overwrite=True)

    # write out the regions
    with open(f'{outpath}/{name}/{name}_vtp.reg', 'w') as reg:
        reg.write("# Region file format: DS9 version 4.1\n")
        reg.write(
            'global color=cyan dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 '
            'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        )
        reg.write('fk5\n')
        for j, xc, yc, rc, rotc, extd in detects[[
                'INDEX', 'RA', 'DEC', 'R', 'ROTANG', 'Extended'
        ]]:
            reg.write(
                f'ellipse({xc},{yc},{(rc[0] * 2.36):.3f}",{(rc[1] * 2.36):.3f}",{rotc:.3f}) '
            )
            if extd > 0:
                reg.write(f'# color=magenta ')
            elif extd < 0:
                reg.write('# color=cyan ')
            else:
                reg.write(f'# color=yellow ')
            reg.write(f'text={{{j}}}\n')

    if plotting:
        plt.plot(x[:bins], psf_cum_n[:bins], label='PSF')
        plt.legend()


# In[ ]:

# get file data
data = load_PSZcatalog()
data = data.sort_index(axis=1)

outpath = './data_full'

arr = [{'name': n.replace(' ', '_'), 'outpath': outpath} for n in data['NAME']]

###
# This is the order you should call the functions.
# There are other functions in this notebook, but these are the only ones
# that should be called directly.
###

parallel_process(arr, ks_test_sources, use_kwargs=True, n_jobs=1)

# In[ ]:

### For testing ###

outpath = './data_full'
name = 'PSZ2_G094.00+27.41'

# In[ ]:

ks_test_sources(name, outpath, True)

# In[ ]:

# detections
if not os.path.isfile(f'{outpath}/{name}/{name}_vtp.detect'):
    pass
else:
    srcs = f'{outpath}/{name}/{name}_vtp.detect'

# mcmc fit data
if not os.path.isfile(f'{outpath}/{name}/{name}_mcmcfits.txt'):
    pass
else:
    mcmcfits = f'{outpath}/{name}/{name}_mcmcfits.txt'
    fit = Table.read(mcmcfits, format='ascii', header_start=0)

# now we need to read the individual detections
detects = Table.read(srcs, hdu=1)

# In[ ]:

detects

# In[ ]:

test_ks_cv(n_trials=30)

# In[ ]:

##############
### Legacy code I didn't want to throw away
##############

# def integ_beta_model(r, rc, beta):
#     # from notebook 4a

#     ''' 2pi*r integral of the above Beta model with 3 parameters.
#     r -- r
#     rc -- core radius
#     beta -- powerlaw slope

#     '''

#     rc2 = rc * rc

#     return np.pi * rc2 / (1 - beta) * ((1 + r**2 / rc2)**(1 - beta) - 1)

#         targetr1 = integ_beta_model(data['r1'], rc=(fit['rc_50'][loc] * 60/2.36), beta=fit['beta_50'][loc])
#         targetr2 = integ_beta_model(data['r2'], rc=(fit['rc_50'][loc] * 60/2.36), beta=fit['beta_50'][loc])

#         target = targetr2 - targetr1
#         source_exp = target / data['y'] # y is the exposure time in each bin

#         # same for the psf -- I've already computed the profiles and stored them in the file
#         # see notebook 4a
#         psf_exp = data['psf'] / data['y']

#         # cummulative distributions
#         source_cum = np.cumsum(source_exp)
#         psf_cum = np.cumsum(psf_exp)

#         # how far out do we want to go?
#         bins = 30

#         # normalized
#         source_cum_n = source_cum / source_cum[bins - 1]
#         psf_cum_n = psf_cum / psf_cum[bins - 1]

# #         # make CDFs
#         source_cum_n = source_cum / source_cum[bins - 1]
#         psf_cum_n = psf_cum / psf_cum[bins - 1]

#         pixscale = 2.36 # arcs/pix

# In[ ]:

np.arange(10)[:4]

# In[ ]:
