#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#%matplotlib notebook
import matplotlib.pyplot as plt

from astropy import table
from astropy.table import Table, Column
from tqdm import tqdm_notebook
import tempfile
import emcee
import numpy.ma as ma
import numpy as np
import subprocess
import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy.coordinates import SkyCoord

from regions import read_ds9, write_ds9
from regions import PixCoord, EllipsePixelRegion, CircleAnnulusPixelRegion

from math import sqrt, sin, cos, pow

# parallel processor
from utilities import parallel_process

# In[ ]:

# Utility functions


def pointInEllipse(xo, yo, xp, yp, d, D, angle):
    #tests if a point[xp,yp] is within
    #boundaries defined by the ellipse
    #of center[x,y], diameter d D, and tilted at angle

    #### the angle must be in radians!!!! ####

    ## This is for the detect_vtp function #

    cosa = cos(angle)
    sina = sin(angle)

    dd = d**2
    DD = D**2

    a = pow(cosa * (xp - xo) + sina * (yp - yo), 2)
    b = pow(sina * (xp - xo) - cosa * (yp - yo), 2)
    ellipse = (a / dd) + (b / DD)

    if ellipse <= 1:
        return True
    else:
        return False


def compound_regions(region_list, img_size=1000):
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


# In[ ]:


def detect_vtp(name, outpath):
    # set up all the non-file-specific parameters
    params = {}
    params['limit'] = 1E-6
    params['coarse'] = 25
    params['maxiter'] = 10
    params['clobber'] = 'yes'
    params['verbose'] = 1

    evts = f'{outpath}/{name}/{name}_events.fits'
    expmap = f'{outpath}/{name}/{name}_exp.fits'

    # check to make sure the files exist.
    if not os.path.isfile(evts) and not os.path.isfile(expmap):
        return 0

    # We are going to run vtpdetect twice... once with a low scale and once with a high(er) scale
    params['regfile'] = f'{outpath}/{name}/{name}_vtp_low.reg'
    params['log'] = f'{outpath}/{name}/{name}_vtp_low.log'
    params['scale'] = 1
    outfits = f'{outpath}/{name}/{name}_vtp_low.detect'
    with open(f'{outpath}/{name}/{name}_vtp_low.in', 'w') as f:
        # build the cmd
        cmd = f'vtpdetect {evts}[pi=50:600] {expmap} {outfits} '
        for param, value in list(params.items()):
            cmd += f'{param}={value} '
        f.writelines(f'{cmd}\n')

    os.system(cmd)

    # now for the higher level
    params['regfile'] = f'{outpath}/{name}/{name}_vtp_high.reg'
    params['log'] = f'{outpath}/{name}/{name}_vtp_high.log'
    params['scale'] = 1.8
    outfits = f'{outpath}/{name}/{name}_vtp_high.detect'
    with open(f'{outpath}/{name}/{name}_vtp_high.in', 'w') as f:
        # build the cmd
        cmd = f'vtpdetect {evts}[pi=50:600] {expmap} {outfits} '
        for param, value in list(params.items()):
            cmd += f'{param}={value} '
        f.writelines(f'{cmd}\n')

    os.system(cmd)

    ###
    # Now we are building the final catalog by comparing the low and high scale catalogs.
    ###
    try:
        low = Table.read(f'{outpath}/{name}/{name}_vtp_low.detect', hdu=1)
        high = Table.read(f'{outpath}/{name}/{name}_vtp_high.detect', hdu=1)
    except FileNotFoundError:
        return

    # create a new table to store our results -- just makes an empty copy of the original
    final = Table(dtype=low.dtype)
    final.add_column(Column(name='INDEX', dtype='>i4'), index=0)
    final.add_column(Column(name='HIGH', dtype='>i4'))

    # have to add columns to the original tables to make them match
    low.add_column(Column(data=np.zeros(len(low)), name='INDEX', dtype='>i4'),
                   index=0)
    low.add_column(Column(data=np.zeros(len(low)), name='HIGH', dtype='>i4'))
    high.add_column(Column(data=np.zeros(len(high)), name='INDEX',
                           dtype='>i4'),
                    index=0)
    high.add_column(Column(data=np.zeros(len(high)), name='HIGH', dtype='>i4'))

    index = 1  # vtp is normally 1 indexed
    for x_l, y_l, rad_l, rot_l, comp_l in low[[
            'X', 'Y', 'R', 'ROTANG', 'COMPONENT'
    ]]:

        added_high = False  # keep track if we found a high source

        for x_h, y_h, comp_h in high[['X', 'Y', 'COMPONENT']]:
            if pointInEllipse(x_l, y_l, x_h, y_h, rad_l[0], rad_l[1],
                              np.radians(rot_l)):
                final.add_row(high[comp_h - 1])  # add the row
                final['INDEX'][index - 1] = index  # add the index to the row
                final['HIGH'][
                    index -
                    1] = 1  # say that the source came from high catalog
                index += 1
                added_high = True

        if not added_high:
            final.add_row(low[comp_l - 1])  # add the row
            final['INDEX'][index - 1] = index  # add the index to the row
            final['HIGH'][index -
                          1] = 0  # say that the source came from low catalog
            index += 1

    # no sources detected
    if not len(final):
        return final

    # for some reason there are duplicate sources
    final = table.unique(final, keys=['X', 'Y'])
    final['INDEX'] = [i + 1 for i in range(len(final))]

    # remove any source where the center of the region lies inside another region
    # gonna do two loops
    remove_idx = []
    for x1, y1, rad1, rot1, idx1, area1 in final[[
            'X', 'Y', 'R', 'ROTANG', 'INDEX', 'SRC_AREA'
    ]]:
        for x2, y2, rad2, rot2, idx2, area2 in final[[
                'X', 'Y', 'R', 'ROTANG', 'INDEX', 'SRC_AREA'
        ]]:
            if idx1 == idx2:
                continue
            elif any(rad2 < 1e-3):  # some sources have radii = 0
                remove_idx.append(
                    idx2 - 1)  # remember that our source index is 1 indexed
            elif pointInEllipse(x1, y1, x2, y2, rad1[0], rad1[1],
                                np.radians(rot1)):
                remove_idx.append(
                    idx2 - 1)  # remember that our source index is 1 indexed

    # remove sources!
    final.remove_rows(remove_idx)

    #reindex the table since we've removed sources
    final['INDEX'] = [i + 1 for i in range(len(final))]

    # write detections!
    final.write(f'{outpath}/{name}/{name}_vtp.detect',
                format='fits',
                overwrite=True)

    # write out the regions
    with open(f'{outpath}/{name}/{name}_vtp.reg', 'w') as reg:
        reg.write("# Region file format: DS9 version 4.1\n")
        reg.write(
            'global color=cyan dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 '
            'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        )
        reg.write('fk5\n')
        for j, xc, yc, rc, rotc in final[['INDEX', 'RA', 'DEC', 'R',
                                          'ROTANG']]:
            reg.write(
                f'ellipse({xc},{yc},{(rc[0] * 2.36):.3f}",{(rc[1] * 2.36):.3f}",{rotc:.3f}) '
            )
            reg.write(f'# text={{{j}}}\n')

    return final


# In[ ]:


def centroid_ciao(name, outpath):

    # check to make sure we've started the ciaotools
    if not check_exe('dmstat'):
        print('dmstat not available. Have you started ciaotools?')
        return

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

    # WCS info to convert pixels to world coords
    hdr = fits.getheader(evnts)
    w = wcs.WCS(hdr)

    # now we need to read the individual detections
    detects = Table.read(srcs, hdu=1)

    # rename columns
    try:
        detects['X'].name = 'X_VTP'
        detects['Y'].name = 'Y_VTP'
        # add the columns to the detection catalog
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='X', dtype='>f8'))
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='Y', dtype='>f8'))
    except KeyError:
        pass

    try:
        detects['X_ERR'].name = 'X_ERR_VTP'
        detects['Y_ERR'].name = 'Y_ERR_VTP'
        # add the columns to the detection catalog
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='X_ERR', dtype='>f8'))
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='Y_ERR', dtype='>f8'))
    except KeyError:
        pass

    try:
        detects['RA'].name = 'RA_VTP'
        detects['DEC'].name = 'DEC_VTP'
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='RA', dtype='>f8'))
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='DEC', dtype='>f8'))
    except KeyError:
        pass

    try:
        detects['RA_ERR'].name = 'RA_ERR_VTP'
        detects['DEC_ERR'].name = 'DEC_ERR_VTP'
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='RA_ERR', dtype='>f8'))
        detects.add_column(
            Column(data=np.zeros(len(detects)), name='DEC_ERR', dtype='>f8'))
    except KeyError:
        pass

    for i, (xc, yc, rc,
            rotc) in enumerate(detects[['X_VTP', 'Y_VTP', 'R', 'ROTANG']]):
        fd, path = tempfile.mkstemp()
        with os.fdopen(fd, 'w') as reg:
            reg.write("# Region file format: CIAO version 1.0\n")
            reg.write(f"+ellipse({xc},{yc},{rc[0]},{rc[1]},{rotc})\n")

        cmd = f"dmstat '{evnts}[(x,y)=region({path})]' centroid=yes clip=yes"

        result = subprocess.Popen(cmd, shell=True,
                                  stdout=subprocess.PIPE).communicate()[0]

        #print(result.decode())

        center = result.decode().split('cntrd[phys]')[1].split('\n')[0].split(
            ' ')[1:3]
        errors = result.decode().split('sigma_cntrd')[1].split('\n')[0].split(
            ' ')[1:3]

        # write the new coordinates to the table
        detects['X'][i], detects['Y'][i] = float(center[0]), float(center[1])
        detects['X_ERR'][i], detects['Y_ERR'][i] = float(errors[0]), float(
            errors[1])

        ### have to subtract 1 from the region positions, because FITS are 1-index'd and python
        ### is 0-index'd.
        world = w.wcs_pix2world(float(center[0]) - 1, float(center[1]) - 1, 0)
        detects['RA'][i], detects['DEC'][i] = world[0], world[1]
        detects['RA_ERR'][i], detects['DEC_ERR'][i] = (detects['X_ERR'][i] *
                                                       2.36 / 3600,
                                                       detects['Y_ERR'][i] *
                                                       2.36 / 3600)

        os.remove(path)

    detects.write(f'{outpath}/{name}/{name}_vtp.detect',
                  format='fits',
                  overwrite=True)

    with open(f'{outpath}/{name}/{name}_vtp.reg', 'w') as reg:
        reg.write("# Region file format: DS9 version 4.1\n")
        reg.write(
            'global color=cyan width=1 font="helvetica 10 normal roman" \n')
        reg.write('fk5\n')
        for j, xc, yc, rc, rotc in detects[[
                'INDEX', 'RA', 'DEC', 'R', 'ROTANG'
        ]]:
            reg.write(
                f'ellipse({xc},{yc},{(rc[0] * 2.36):.3f}",{(rc[1] * 2.36):.3f}",{rotc:.3f}) '
            )
            reg.write(f'# text={{{j}}}\n')

    return detects


# In[ ]:


def beta_model(So, r, rc, beta):
    ''' Beta model with 3 parameters.
    So -- normalization
    rc -- core radius
    beta -- powerlaw slope

    '''

    return So * (1.0 + (r / rc)**2)**(-3.0 * beta + 0.5)


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


def fit_profile(name, outpath, src_num=0, srcs=None):

    # swift pixelscale in degrees
    pixsc = 6.548089E-04 * 3600

    # detections
    if not srcs:
        srcs = f'{outpath}/{name}/{name}_vtp.detect'
        # now we need to read the individual detections
        detects = Table.read(srcs, hdu=1)
    else:
        detects = srcs

    # get the source index
    i = detects['INDEX'][src_num]

    if os.path.isfile(f'{outpath}/{name}/{name}_vtp_{i}.tmpprof'):
        data = Table.read(f'{outpath}/{name}/{name}_vtp_{i}.tmpprof',
                          format='ascii',
                          header_start=2)
    else:
        return

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

    #fit -- median (and error) values
    So, rc, beta, bg = map(
        lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
        zip(*np.percentile(emcee_trace, [16, 50, 84], axis=0)))

    # open a file to write out the results
    with open(f'{outpath}/{name}/{name}_{i}_mcmcfits.tmp', 'w') as outfile:
        outfile.write('# Field ID So_50 So_84 So_16 rc_50 rc_84 rc_16 '
                      'beta_50 beta_84 beta_16 bg_50 bg_84 bg_16\n')

        # write out the fit parameters
        outfile.write(f'{name} ')
        outfile.write(f'{i} ')
        outfile.write(f'{So[0]:.5f} {So[1]:.5f} {So[2]:.5f} ')
        outfile.write(f'{rc[0]:.5f} {rc[1]:.5f} {rc[2]:.5f} ')
        outfile.write(f'{beta[0]:.5f} {beta[1]:.5f} {beta[2]:.5f} ')
        outfile.write(f'{bg[0]:.5f} {bg[1]:.5f} {bg[2]:.5f}\n')


# In[ ]:


def get_radprof(name, outpath, src_num=0, srcs=None):
    pixscale = 6.548089E-04 * 3600  # arcsec/pixel

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
    if not srcs:
        srcs = f'{outpath}/{name}/{name}_vtp.detect'
        # now we need to read the individual detections
        detects = Table.read(srcs, hdu=1)
    else:
        detects = srcs

    i = detects['INDEX'][src_num]
    xc = detects['X'][src_num]
    yc = detects['Y'][src_num]

    # load the data
    evnt_data = fits.getdata(evnts)
    expmap_data = fits.getdata(expmap)

    # mask out all sources except the source we are working with
    if len(detects) < 2:
        # here there is only one source in the field. Don't mask anything.
        mask_sources = False
    else:
        mask_sources = True
        src_mask = [j for j in range(len(detects)) if j != src_num]

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

    # output file
    with open(f'{outpath}/{name}/{name}_vtp_{i}.tmpprof', 'w') as radprof:
        radprof.write(
            f"# source number and position: {name}_{i} {xc:.3f} {yc:.3f}\n")
        radprof.write(
            f"# profile parameters: {params['n']} {params['r0']} {params['dr']}\n"
        )
        radprof.write(f"# bin r1 r2 x y w sb sb_err\n")

        for irad in range(params['n']):
            r1 = params['r0'] + irad * params['dr']
            r2 = params['r0'] + (irad + 1) * params['dr']

            center = PixCoord(x=detects['X'][src_num] - 1,
                              y=detects['Y'][src_num] - 1)
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

            radprof.write(
                f'{irad:3d} {r1:7d} {r2:7d} {x:9.1f} {y:9.1f} {w:9.1f} {sb:12.4e} {sbe:12.4e}\n'
            )


# In[ ]:


def clean_sources(name, outpath):

    itter = 1
    srcs = Table.read(f'{outpath}/{name}/{name}_vtp.detect', hdu=1)
    while True:
        # print(itter, len(srcs))
        src = np.argmax(srcs['NET_COUNTS'])
        get_radprof(name, outpath, src, srcs)
        fit_profile(name, outpath, src, srcs)
        fit = Table.read(
            f'{outpath}/{name}/{name}_{srcs["INDEX"][src]}_mcmcfits.tmp',
            format='ascii',
            header_start=0)

        weak_srcs = []
        for i in range(len(srcs)):
            if i == src:
                continue

            clust = SkyCoord(srcs['RA'][src],
                             srcs['DEC'][src],
                             unit='deg',
                             frame='icrs')
            src1 = SkyCoord(srcs['RA'][i],
                            srcs['DEC'][i],
                            unit='deg',
                            frame='icrs')

            sep = clust.separation(src1)
            r = sep.arcmin

            flux = beta_model(fit['So_50'], r, fit['rc_50'],
                              fit['beta_50']) + fit['bg_50']
            flux = flux.data[0]

            cnts = srcs['NET_COUNTS'][i] + srcs['BKG_COUNTS'][i]
            exptime = srcs['EXPTIME'][i]
            src_area = srcs['SRC_AREA'][i] * 2.36**2 / 60**2

            src_flux = cnts / exptime / src_area
            src_flux_err = sqrt(cnts) / exptime / src_area

            sigma = (src_flux - flux) / src_flux_err
            if sigma < 5:
                weak_srcs.append(i)
            #print(srcs['INDEX'][i], sigma)

        # print('bad sources', len(weak_srcs))
        if not len(weak_srcs):
            # remove tmp files.
            # print('remove files', srcs["INDEX"][src])
            os.remove(
                f'{outpath}/{name}/{name}_vtp_{srcs["INDEX"][src]}.tmpprof')
            os.remove(
                f'{outpath}/{name}/{name}_{srcs["INDEX"][src]}_mcmcfits.tmp')
            break
        else:
            srcs.remove_rows(weak_srcs)

        itter += 1

    srcs['INDEX'] = [i + 1 for i in range(len(srcs))]
    srcs.write(f'{outpath}/{name}/{name}_vtp.detect',
               format='fits',
               overwrite=True)

    with open(f'{outpath}/{name}/{name}_vtp.reg', 'w') as reg:
        #with open(f'{outpath}/{name}/test.reg', 'w') as reg:
        reg.write("# Region file format: DS9 version 4.1\n")
        reg.write(
            'global color=cyan width=1 font="helvetica 10 normal roman" \n')
        reg.write('fk5\n')
        for j, xc, yc, rc, rotc in srcs[['INDEX', 'RA', 'DEC', 'R', 'ROTANG']]:
            reg.write(
                f'ellipse({xc},{yc},{(rc[0] * 2.36):.3f}",{(rc[1] * 2.36):.3f}",{rotc:.3f}) '
            )
            reg.write(f'# text={{{j}}}\n')

    return srcs


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

parallel_process(arr, detect_vtp, use_kwargs=True, n_jobs=1)
parallel_process(arr, centroid_ciao, use_kwargs=True, n_jobs=6)
parallel_process(arr, clean_sources, use_kwargs=True, n_jobs=6)

# In[ ]:

### For testing ###

outpath = './data_full'
name = 'PSZ1_G023.53-36.53'

# In[ ]:

detect_vtp(name, outpath)

# In[ ]:

centroid_ciao(name, outpath)

# In[ ]:

clean_sources(name, outpath)

# In[ ]:
