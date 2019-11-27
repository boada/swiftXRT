#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#%matplotlib notebook
import matplotlib.pyplot as plt

from astropy import table
from astropy.table import Table, Column
import tempfile
import numpy.ma as ma
import numpy as np
import emcee
import corner
import subprocess
import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

from astropy.io import fits
from astropy import wcs
import astropy.units as u

from regions import PixCoord, EllipsePixelRegion, CircleAnnulusPixelRegion, CirclePixelRegion

from math import sqrt, sin, cos, pow

# parallel processor
from utilities import parallel_process

# In[ ]:

# Utility functions


def compound_regions(region_list, img_size=1000):

    if not isinstance(region_list, list):
        m = region_list.to_mask()
        compound_region = m.to_image((img_size, img_size))
    else:
        for r in region_list:
            # handle the first pass when we only have one part
            try:
                m = r.to_mask()
                compound_region = compound_region + m.to_image(
                    (img_size, img_size))
            except NameError:
                m = r.to_mask()
                compound_region = m.to_image((img_size, img_size))

    return compound_region


def check_exe(exe, verb=False):
    ''' Checks to make sure we have the appropriate system command available.
    If we don't it raises an exception.

    '''

    path = os.environ['PATH'].split(':')
    for p in path:
        f = os.path.join(p, exe)
        if os.path.isfile(f):
            if verb:
                print("# Found %s in %s" % (exe, f), file=sys.stderr)
            return True
    # it wasn't found
    print("# ERROR: Couldn't find %s" % exe, file=sys.stderr)
    raise FileNotFoundError(exe)
    return False


def clean_files(outpath):
    print('Cleaning files...')

    os.system(f'find {outpath}/PSZ* -name "*corner_*.png" -type f -delete')
    os.system(f'find {outpath}/PSZ* -name "*_mcmcfits.txt" -type f -delete')
    os.system(
        f'find {outpath}/PSZ* -name "*_radproffit_*.png" -type f -delete')
    os.system(f'find {outpath}/PSZ* -name "*.radprof" -type f -delete')


def cp_results(outpath):
    print('Copying files...')

    os.system(
        f'find {outpath}/PSZ* -name "*corner_*.png" -exec cp -t {outpath}/pngs/fits/ '
        '{} \;')
    os.system(
        f'find {outpath}/PSZ* -name "*radproffit_*.png" -exec cp -t {outpath}/pngs/fits/ '
        '{} \;')


# In[ ]:


def beta_model(So, r, rc, beta):
    ''' Beta model with 3 parameters.
    So -- normalization
    rc -- core radius
    beta -- powerlaw slope

    returns the flux (Cnts/pixel)

    '''

    return So * (1.0 + (r / rc)**2)**(-3.0 * beta + 0.5)


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


def fit_radprof(name, outpath):

    # swift pixelscale in degrees
    pixsc = 6.548089E-04 * 3600

    if os.path.isfile(f'{outpath}/{name}/{name}_vtp.detect'):
        srcs = f'{outpath}/{name}/{name}_vtp.detect'
        # now we need to read the individual detections
        detects = Table.read(srcs, hdu=1)
    else:
        return

    # we are going to fit profiles to the top 3 biggest
    detects.sort(['NET_COUNTS'], reverse=True)

    # open a file to write out the results
    outfile = open(f'{outpath}/{name}/{name}_mcmcfits.txt', 'w')
    outfile.write(
        '# Field ID So_50 So_84 So_16 '
        'rc_50 rc_84 rc_16 beta_50 beta_84 beta_16 bg_50 bg_84 bg_16\n')

    # loop over the sources -- only the first 20
    for i in detects['INDEX'][:20]:

        if os.path.isfile(f'{outpath}/{name}/{name}_vtp_{i}.radprof'):
            data = Table.read(f'{outpath}/{name}/{name}_vtp_{i}.radprof',
                              format='ascii',
                              header_start=2)
        else:
            continue

        # don't try to fit things if there isn't any data
        if len(data) <= 60:
            continue

        # mask out any points where the sb is zero
        mask = data['sb'] != 0
        data = data[mask]

        # x-axis, in arcminutes
        x = (data['r1'] + data['r2']) / 2. / 60. * pixsc

        # this is the parameter fitting
        ndim = 4  # number of parameters in the model
        nwalkers = 100  # number of MCMC walkers
        nburn = 100  # "burn-in" period to let chains stabilize
        nsteps = 500  # number of MCMC steps to take

        pos = [
            np.array([0.001, 1., 1., 0.001]) + 1e-4 * np.random.randn(ndim)
            for i in range(nwalkers)
        ]
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        prob,
                                        args=(x, data['sb'], data['sb_err']))
        sampler.run_mcmc(pos, nsteps)
        emcee_trace = sampler.get_chain(discard=nburn, thin=15, flat=True)

        # plot the result
        fignum = np.random.randint(1, 1000)
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6.5, 5), num=fignum)
        # data
        ax.errorbar(x, data['sb'], yerr=data['sb_err'], fmt='o', label='data')
        for So, rc, beta, bg in emcee_trace[np.random.randint(len(emcee_trace),
                                                              size=100)]:
            ax.plot(x, beta_model(So, x, rc, beta) + bg, color='k', alpha=0.1)
        ylims = ax.get_ylim()

        # fit -- median (and error) values
        So, rc, beta, bg = map(
            lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
            zip(*np.percentile(emcee_trace, [16, 50, 84], axis=0)))
        ax.plot(x,
                beta_model(So[0], x, rc[0], beta[0]) + bg[0],
                label=r'$\beta$ + bg')
        ax.plot(x, beta_model(So[0], x, rc[0], beta[0]), label=r'$\beta$')

        # add the psf
        ax.plot(x, data['psfsb'] + bg[0], label='psf', color='orange')

        # write out the fit parameters
        outfile.write(f'{name} ')
        outfile.write(f'{i} ')
        outfile.write(f'{So[0]:.5f} {So[1]:.5f} {So[2]:.5f} ')
        outfile.write(f'{rc[0]:.5f} {rc[1]:.5f} {rc[2]:.5f} ')
        outfile.write(f'{beta[0]:.5f} {beta[1]:.5f} {beta[2]:.5f} ')
        outfile.write(f'{bg[0]:.5f} {bg[1]:.5f} {bg[2]:.5f}\n')

        # rc and bg lines
        ax.axhline(bg[0], label='bg', zorder=0, lw=1)
        ax.axvline(rc[0], zorder=0, lw=1)

        ax.semilogy()

        if ylims[0] < 1e-5:
            ax.text(rc[0] + 0.1,
                    2e-5,
                    f"rc = {rc[0]:.2f}'",
                    rotation='vertical',
                    ha='left')
            ax.set_ylim(1e-5, ylims[1])
        else:
            ax.text(rc[0] + 0.1,
                    ylims[0],
                    f"rc = {rc[0]:.2f}'",
                    rotation='vertical',
                    ha='left')
            ax.set_ylim(ylims)
        ax.set_xlabel('Radius [arcmin]')
        ax.set_ylabel('Flux $[cnts/s/arcmin^2]$')
        ax.legend(loc='upper right')

        fig.savefig(f'{outpath}/{name}/{name}_radproffit_{i}.png',
                    bbox='tight',
                    dpi=180)
        plt.close(fig)

        # add a corner plot of the fits.
        # plot the bg and So in log to make them fit on the plot better
        emcee_trace[:, 0] = np.log10(emcee_trace[:, 0])
        emcee_trace[:, 3] = np.log10(emcee_trace[:, 3])
        f = corner.corner(
            emcee_trace,
            labels=['log So', 'rc', r'$\beta$', 'log bg'],
            bins=[50, 50, 50, 50],
            smooth=True,
            fill_contours=True,
            truths=[np.log10(So[0]), rc[0], beta[0],
                    np.log10(bg[0])])

        f.savefig(f'{outpath}/{name}/{name}_corner_{i}.png', bbox='tight')
        plt.close(f)
    outfile.close()


# In[ ]:


def get_radprof(name, outpath):
    pixscale = 6.548089E-04 * 3600

    # model parameters
    params = {}
    params['n'] = 75
    params['r0'] = 0
    params['dr'] = 4

    # check for files
    # events image
    if not os.path.isfile(f'{outpath}/{name}/{name}_img_50-200.fits'):
        return
    else:
        evnts = f'{outpath}/{name}/{name}_img_50-200.fits'
    # expo map
    if not os.path.isfile(f'{outpath}/{name}/{name}_exp.fits'):
        return
    else:
        expmap = f'{outpath}/{name}/{name}_exp.fits'
    # detections
    if not os.path.isfile(f'{outpath}/{name}/{name}_vtp.detect'):
        return
    else:
        srcs = f'{outpath}/{name}/{name}_vtp.detect'

    # load the data
    evnt_data = fits.getdata(evnts)
    expmap_data = fits.getdata(expmap)

    # now we need to read the individual detections
    detects = Table.read(srcs, hdu=1)

    # we are going to compute profiles for the 10 biggest -- this is for speed.
    detects.sort(['NET_COUNTS'], reverse=True)

    for i in range(len(detects[:20])):
        # mask out all sources except the source we are working with
        if len(detects) < 2:
            # here there is only one source in the field. Don't mask anything.
            mask_sources = False
        else:
            mask_sources = True
            src_mask = [j for j in range(len(detects)) if j != i]

        if mask_sources:
            regs = []
            ### have to subtract 1 from the region positions, because FITS are 1-index'd and python
            ### is 0-index'd.
            for idx, xc1, yc1, rc, rotc in detects[src_mask][[
                    'INDEX', 'X', 'Y', 'R', 'ROTANG'
            ]]:
                ellipse_region = EllipsePixelRegion(center=PixCoord(x=xc1 - 1,
                                                                    y=yc1 - 1),
                                                    width=2 * rc[0],
                                                    height=2 * rc[1],
                                                    angle=rotc * u.deg,
                                                    meta={'ID': idx})
                regs.append(ellipse_region)

            # create the compound regions
            cmp_reg = compound_regions(regs, evnt_data.shape[0])
            # invert the comp region -- we want the places where there aren't regions
            cmp_regmask = np.where(cmp_reg == 0, 1, 0)

        else:
            cmp_regmask = np.ones_like(evnt_data)

        with open(f'{outpath}/{name}/{name}_vtp_{detects["INDEX"][i]}.radprof',
                  'w') as radprof:
            radprof.write(
                f"# source number and position: {name}_{detects['INDEX'][i]} "
                f"{detects['X'][i]:.3f} {detects['Y'][i]:.3f}\n")
            radprof.write(
                f"# profile parameters: {params['n']} {params['r0']} {params['dr']}\n"
            )
            radprof.write(f"# bin r1 r2 x y w sb sb_err psf area psfsb\n")

            for irad in range(params['n']):
                r1 = params['r0'] + irad * params['dr']
                r2 = params['r0'] + (irad + 1) * params['dr']

                center = PixCoord(x=detects['X'][i] - 1, y=detects['Y'][i] - 1)
                if r1 == 0:
                    ann_reg = CirclePixelRegion(center=center, radius=r2)
                else:
                    ann_reg = CircleAnnulusPixelRegion(center=center,
                                                       inner_radius=r1,
                                                       outer_radius=r2)
                # make it into an image
                ann_regmask = compound_regions([ann_reg], evnt_data.shape[0])

                # combine the regions with the annulus.
                final_regmask = (ann_regmask * cmp_regmask)
                final_regmask_inv = np.where(
                    final_regmask == 0, 1,
                    0)  # this is the final mask for the data

                evnt_data_masked = ma.masked_array(evnt_data,
                                                   mask=final_regmask_inv)
                expmap_data_masked = ma.masked_array(expmap_data,
                                                     mask=final_regmask_inv)

                x = evnt_data_masked.sum()
                y = expmap_data_masked.mean()
                w = expmap_data_masked.count()

                sb = x / (y * w * pixscale**2 / 3600)  #  cts/s/arcmin^2
                err_gehrels = sqrt(x + 0.75)
                sbe = (1. + err_gehrels) / (y * w * pixscale**2 / 3600
                                            )  #  cts/s/arcmin^2

                # build the model of the psf
                psfr1 = integ_beta_model(r1, rc=(5.6 / 2.36), beta=1.5)  # cnts
                psfr2 = integ_beta_model(r2, rc=(5.6 / 2.36), beta=1.5)  # cnts
                psf = psfr2 - psfr1  # total flux in the annulus

                area = np.pi * r2**2 - np.pi * r1**2  # number of square pixels

                psfsb = psf / (y * area * pixscale**2 / 3600)

                radprof.write(
                    f'{irad:3d} {r1:7d} {r2:7d} {x:9.1f} {y:9.1f} {w:9.1f}'
                    f'{sb:12.4e} {sbe:12.4e} '
                    f'{psf:9.1f} {area:9.1f} {psfsb:12.4e} \n')


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

#clean_files(outpath)

#parallel_process(arr, get_radprof, use_kwargs=True, n_jobs=6)
parallel_process(arr, fit_radprof, use_kwargs=True, n_jobs=6)

cp_results(outpath)

# In[ ]:

### For testing ###

outpath = './data_full'
name = 'PSZ1_G221.64-58.20'

# In[ ]:

get_radprof(name, outpath)

# In[ ]:

fit_radprof(name, outpath)

# In[ ]:

pixscale = 6.548089E-04 * 3600

# model parameters
params = {}
params['n'] = 75
params['r0'] = 0
params['dr'] = 4

# check for files
# events image
if not os.path.isfile(f'{outpath}/{name}/{name}_img_50-200.fits'):
    pass
else:
    evnts = f'{outpath}/{name}/{name}_img_50-200.fits'
# expo map
if not os.path.isfile(f'{outpath}/{name}/{name}_exp.fits'):
    pass
else:
    expmap = f'{outpath}/{name}/{name}_exp.fits'
# detections
if not os.path.isfile(f'{outpath}/{name}/{name}_vtp.detect'):
    pass
else:
    srcs = f'{outpath}/{name}/{name}_vtp.detect'

# load the data
evnt_data = fits.getdata(evnts)
expmap_data = fits.getdata(expmap)

# now we need to read the individual detections
detects = Table.read(srcs, hdu=1)

# we are going to compute profiles for the 10 biggest -- this is for speed.
detects.sort(['NET_COUNTS'], reverse=True)

i = 1
print(i)
# mask out all sources except the source we are working with
if len(detects) < 2:
    mask_sources = False
    src_mask = [0]
else:
    mask_sources = True
    src_mask = [j for j in range(len(detects)) if j != i]

if mask_sources:
    regs = []
    ### have to subtract 1 from the region positions, because FITS are 1-index'd and python
    ### is 0-index'd.
    for idx, xc1, yc1, rc, rotc in detects[src_mask][[
            'INDEX', 'X', 'Y', 'R', 'ROTANG'
    ]]:
        ellipse_region = EllipsePixelRegion(center=PixCoord(x=xc1 - 1,
                                                            y=yc1 - 1),
                                            width=2 * rc[0],
                                            height=2 * rc[1],
                                            angle=rotc * u.deg,
                                            meta={'ID': idx})
        regs.append(ellipse_region)

    # create the compound regions
    cmp_reg = compound_regions(regs, evnt_data.shape[0])
    # invert the comp region -- we want the places where there aren't regions
    cmp_regmask = np.where(cmp_reg == 0, 1, 0)

else:
    cmp_regmask = np.ones_like(evnt_data)

for irad in range(params['n']):
    r1 = params['r0'] + irad * params['dr']
    r2 = params['r0'] + (irad + 1) * params['dr']

    center = PixCoord(x=detects['X'][i] - 1, y=detects['Y'][i] - 1)
    if r1 == 0:
        ann_reg = CirclePixelRegion(center=center, radius=r2)
        ann_reg = CircleAnnulusPixelRegion(center=center,
                                           inner_radius=r1,
                                           outer_radius=r2)
    else:
        ann_reg = CircleAnnulusPixelRegion(center=center,
                                           inner_radius=r1,
                                           outer_radius=r2)
    # make it into an image
    ann_regmask = compound_regions([ann_reg], evnt_data.shape[0])

    # combine the regions with the annulus.
    final_regmask = (ann_regmask * cmp_regmask)
    final_regmask_inv = np.where(final_regmask == 0, 1,
                                 0)  # this is the final mask for the data

    evnt_data_masked = ma.masked_array(evnt_data, mask=final_regmask_inv)
    expmap_data_masked = ma.masked_array(expmap_data, mask=final_regmask_inv)

    x = evnt_data_masked.sum()
    y = expmap_data_masked.mean()
    w = expmap_data_masked.count()

    sb = x / (y * w * pixscale**2 / 3600)  #  cts/s/arcmin^2
    err_gehrels = sqrt(x + 0.75)
    sbe = (1. + err_gehrels) / (y * w * pixscale**2 / 3600)  #  cts/s/arcmin^2

    # build the model of the psf
    psfr1 = integ_beta_model(r1, rc=(5.6 / 2.36), beta=1.5)  # cnts
    psfr2 = integ_beta_model(r2, rc=(5.6 / 2.36), beta=1.5)  # cnts
    psf = psfr2 - psfr1  # total flux in the annulus

    area = np.pi * r2**2 - np.pi * r1**2  # number of pixels

    psfsb = psf / (y * area * pixscale**2 / 3600)

    print(f'{irad:3d} {r1:7d} {r2:7d} {x:9.1f} {y:9.1f} {w:9.1f}'
          f'{sb:12.4e} {sbe:12.4e} '
          f'{psf:9.1f} {area:9.1f} {psfsb:12.4e} \n')
    if irad > 4:
        break

# In[ ]:

# swift pixelscale in degrees
pixsc = 6.548089E-04 * 3600

if os.path.isfile(f'{outpath}/{name}/{name}_vtp.detect'):
    srcs = f'{outpath}/{name}/{name}_vtp.detect'
    # now we need to read the individual detections
    detects = Table.read(srcs, hdu=1)

# we are going to fit profiles to the top 3 biggest
detects.sort(['NET_COUNTS'], reverse=True)

i = 0
print(name)
if os.path.isfile(f'{outpath}/{name}/{name}_vtp_{i}.radprof'):
    data = Table.read(f'{outpath}/{name}/{name}_vtp_{i}.radprof',
                      format='ascii',
                      header_start=2)

# x-axis, in arcminutes
x = (data['r1'] + data['r2']) / 2. / 60. * pixsc

# this is the parameter fitting
ndim = 4  # number of parameters in the model
nwalkers = 100  # number of MCMC walkers
nburn = 100  # "burn-in" period to let chains stabilize
nsteps = 500  # number of MCMC steps to take

pos = [
    np.array([0.001, 1., 1., 0.001]) + 1e-4 * np.random.randn(ndim)
    for i in range(nwalkers)
]
sampler = emcee.EnsembleSampler(nwalkers,
                                ndim,
                                prob,
                                args=(x, data['sb'], data['sb_err']))
sampler.run_mcmc(pos, nsteps)
emcee_trace = sampler.get_chain(discard=nburn, thin=15, flat=True)

# plot the result
fignum = np.random.randint(1, 1000)
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6.5, 5), num=fignum)
# data
ax.errorbar(x, data['sb'], yerr=data['sb_err'], fmt='o', label='data')
#ax.scatter(x, data['psfsb'], label='psf', color='orange')
for So, rc, beta, bg in emcee_trace[np.random.randint(len(emcee_trace),
                                                      size=100)]:
    ax.plot(x, beta_model(So, x, rc, beta) + bg, color='k', alpha=0.1)
ylims = ax.get_ylim()

#fit -- median (and error) values
So, rc, beta, bg = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                       zip(*np.percentile(emcee_trace, [16, 50, 84], axis=0)))
ax.plot(x, beta_model(So[0], x, rc[0], beta[0]) + bg[0], label=r'$\beta$ + bg')
ax.plot(x, beta_model(So[0], x, rc[0], beta[0]), label=r'$\beta$')

ax.plot(x, data['psfsb'] + bg[0], label='psf', color='orange')

# rc and bg lines
ax.axhline(bg[0], label='bg', zorder=0, lw=1)
ax.axvline(rc[0], zorder=0, lw=1)

ax.semilogy()

if ylims[0] < 1e-5:
    ax.text(rc[0] + 0.1,
            2e-5,
            f"rc = {rc[0]:.2f}'",
            rotation='vertical',
            ha='left')
    ax.set_ylim(1e-5, ylims[1])
else:
    ax.text(rc[0] + 0.1,
            ylims[0],
            f"rc = {rc[0]:.2f}'",
            rotation='vertical',
            ha='left')
    ax.set_ylim(ylims)
ax.set_xlabel('Radius [arcmin]')
ax.set_ylabel('Flux [cnts/s/$arcmin^2$]')
ax.legend(loc='upper right')

# plt.close()

# # add a corner plot of the fits.
# # plot the bg and So in log to make them fit on the plot better
# emcee_trace[:,0] = np.log10(emcee_trace[:,0])
# emcee_trace[:,3] = np.log10(emcee_trace[:,3])
# f = corner.corner(emcee_trace, labels=['log So', 'rc', r'$\beta$', 'log bg'],
#             bins=[50, 50, 50, 50],
#             smooth=True,
#             fill_contours=True,
#             truths=[np.log10(So[0]), rc[0], beta[0], np.log10(bg[0])])

# In[ ]:

fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = labels = ['log So', 'rc', r'$\beta$', 'log bg']
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")

# In[ ]:

emcee_trace_lp = sampler.get_log_prob(discard=nburn, thin=15, flat=True)

# In[ ]:

emcee_trace = sampler.get_chain(discard=nburn, thin=15, flat=True)

# In[ ]:

emcee_trace.shape

# In[ ]:

emcee_trace_lp.shape

# In[ ]:

emcee_trace[:, 0] = np.log10(emcee_trace[:, 0])
emcee_trace[:, 3] = np.log10(emcee_trace[:, 3])
f = corner.corner(emcee_trace,
                  labels=['log So', 'rc', r'$\beta$', 'log bg'],
                  bins=[50, 50, 50, 50],
                  smooth=True,
                  fill_contours=True,
                  truths=[np.log10(So[0]), rc[0], beta[0],
                          np.log10(bg[0])],
                  weights=-1 / emcee_trace_lp,
                  quantiles=[.16, .84])

# In[ ]:

get_ipython().run_line_magic('matplotlib', 'notebook')
n, xe, ye, _ = plt.hist2d(emcee_trace[:, 1],
                          emcee_trace[:, 2],
                          weights=-1 / emcee_trace_lp,
                          bins=25)

# In[ ]:

plt.hist2d

# In[ ]:

np.unravel_index(np.argmax(n, axis=None), n.shape)

# In[ ]:

xe[19]

# In[ ]:

ye[24]

# In[ ]:

fignum = np.random.randint(1, 1000)
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6.5, 5), num=fignum)
# data
ax.errorbar(x, data['sb'], yerr=data['sb_err'], fmt='o', label='data')
for So, rc, beta, bg in emcee_trace[np.random.randint(len(emcee_trace),
                                                      size=100)]:
    ax.plot(x, beta_model(So, x, rc, beta) + bg, color='k', alpha=0.1)
ylims = ax.get_ylim()

#fit -- median (and error) values
So, rc, beta, bg = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                       zip(*np.percentile(emcee_trace, [16, 50, 84], axis=0)))
ax.plot(x, beta_model(So[0], x, rc[0], beta[0]) + bg[0], label=r'$\beta$ + bg')
ax.plot(x, beta_model(So[0], x, rc[0], beta[0]), label=r'$\beta$')

# rc and bg lines
ax.axhline(bg[0], label='bg', zorder=0, lw=1)
ax.axvline(rc[0], zorder=0, lw=1)

ax.semilogy()

# In[ ]:
