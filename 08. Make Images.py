#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#%matplotlib notebook
import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog
import numpy as np
import subprocess
import aplpy
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from astropy.table import Table
from astropy.io.fits import getheader
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from tqdm import tqdm_notebook

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy import log
log.setLevel('WARN')

# # Add extra code to run things in parallel.
# The combining exposure maps takes a REALLY long time. So I added this bit so I can run it all in parallel.
#
# You don't need to run any of this if you don't want to run things in paralle. I'll try to comment out code at the bottom so you can run things in parallel or not depending on what you want to do.

# In[ ]:

from concurrent.futures import ThreadPoolExecutor, as_completed


def parallel_process(array,
                     function,
                     n_jobs=None,
                     use_kwargs=False,
                     front_num=0):
    """
        A parallel version of the map function with a progress bar.

        Args:
            array (array-like): An array to iterate over.
            function (function): A python function to apply to the elements of array
            n_jobs (int, default=16): The number of cores to use
            use_kwargs (boolean, default=False): Whether to consider the elements of array as dictionaries of
                keyword arguments to function
            front_num (int, default=3): The number of iterations to run serially before kicking off the
                parallel job. This can be useful for catching bugs
        Returns:
            [function(array[0]), function(array[1]), ...]
    """
    #We run the first few iterations serially to catch bugs
    if front_num > 0:
        front = [
            function(**a) if use_kwargs else function(a)
            for a in array[:front_num]
        ]
    #If we set n_jobs to 1, just run a list comprehension. This is useful for benchmarking and debugging.
    if n_jobs == 1:
        return [
            function(**a) if use_kwargs else function(a)
            for a in tqdm_notebook(array[front_num:])
        ]
    #Assemble the workers
    with ThreadPoolExecutor(max_workers=n_jobs) as pool:
        #Pass the elements of array into function
        if use_kwargs:
            futures = [pool.submit(function, **a) for a in array[front_num:]]
        else:
            futures = [pool.submit(function, a) for a in array[front_num:]]
        kwargs = {
            'total': len(futures),
            'unit': 'it',
            'unit_scale': False,
            'leave': True
        }
        #Print out the progress as tasks complete
        for f in tqdm_notebook(as_completed(futures), **kwargs):
            pass


# In[ ]:


def cp_results(outpath):
    print('Copying files...')

    os.system(
        f'find {outpath}/PSZ* -name "*vtp.png" -exec cp -t {outpath}/pngs/images/ '
        '{} \;')


# In[ ]:


def show_sources(name, outpath):
    # check for files
    # check for detections
    if os.path.isfile(f'{outpath}/{name}/{name}_vtp.detect'):
        srcs = f'{outpath}/{name}/{name}_vtp.detect'
    else:
        return
    # events image
    if os.path.isfile(f'{outpath}/{name}/{name}_img_50-600.fits'):
        evnts = f'{outpath}/{name}/{name}_img_50-600.fits'
    else:
        return

    # check for imaging -- sdss first
    if os.path.isfile(f'{outpath}/{name}/{name}_DESstack_r.fits'):
        survey = 'DES'
        optimage = True
        img = f'{outpath}/{name}/{name}_DESstack_r.fits'
    elif os.path.isfile(f'{outpath}/{name}/{name}_DECaLSstack_r.fits'):
        survey = 'DECaLS'
        optimage = True
        img = f'{outpath}/{name}/{name}_DECaLSstack_r.fits'
    elif os.path.isfile(f'{outpath}/{name}/{name}_SDSSstack_i.fits'):
        survey = 'SDSS'
        optimage = True
        img = f'{outpath}/{name}/{name}_SDSSstack_i.fits'
    elif os.path.isfile(f'{outpath}/{name}/{name}_PS1stack_i.fits'):
        survey = 'PS1'
        optimage = True
        img = f'{outpath}/{name}/{name}_PS1stack_i.fits'
    else:
        optimage = False

    # now we need to read the individual detections
    detects = Table.read(srcs, hdu=1)
    info = Table.read(f'{outpath}/{name}/{name}.info', format='ascii.fast_csv')

    # show the figure
    gc = aplpy.FITSFigure(evnts, figsize=(10, 10))
    gc.show_grayscale(vmin=0, pmax=100, stretch='linear', interpolation='none')

    # add all the sources -- the float is the pixscale in deg
    #     r1 = detects['R'].data[:,0] * 6.548089E-04
    #     r2 = detects['R'].data[:,1] * 6.548089E-04
    #     gc.show_ellipses(detects['RA'], detects['DEC'], r1, r2, detects['ROTANG'],
    #                      coords_frame='world', edgecolor='cyan')

    gc.show_regions(f'{outpath}/{name}/{name}_vtp.reg')

    # add PSZ info and circles
    gc.show_circles(info['RA'],
                    info['DEC'],
                    2 / 60,
                    linestyle='--',
                    edgecolor='#e24a33',
                    facecolor='none',
                    path_effects=[
                        pe.Stroke(linewidth=1.2, foreground='white'),
                        pe.Normal()
                    ])
    gc.show_circles(info['RA'],
                    info['DEC'],
                    5 / 60,
                    linestyle='-',
                    edgecolor='#e24a33',
                    facecolor='none',
                    path_effects=[
                        pe.Stroke(linewidth=1.2, foreground='white'),
                        pe.Normal()
                    ])
    gc.show_markers(info['RA'],
                    info['DEC'],
                    marker='*',
                    s=150,
                    layer='psz',
                    edgecolor='#e24a33',
                    path_effects=[
                        pe.Stroke(linewidth=1.2, foreground='white'),
                        pe.Normal()
                    ])

    # write the exposure time
    exp_time = getheader(evnts)['EXPOSURE']
    text = f'exp time: {exp_time:.2f}s'
    xo, yo = (0.05, 0.05)
    gc.add_label(xo,
                 yo,
                 text,
                 relative=True,
                 fontsize=18,
                 color='white',
                 horizontalalignment='left')

    # write redshift
    ztext = f'z: {info["REDSHIFT"][0]:.3f}'
    gc.add_label(xo,
                 yo + 0.03,
                 ztext,
                 relative=True,
                 fontsize=18,
                 color='white',
                 horizontalalignment='left')

    # write legend
    xo, yo = (0.95, 0.05)
    gc.add_label(xo,
                 yo,
                 'Extended',
                 relative=True,
                 fontsize=18,
                 color='magenta',
                 horizontalalignment='right')
    gc.add_label(xo,
                 yo + 0.03,
                 'P-Source',
                 relative=True,
                 fontsize=18,
                 color='yellow',
                 horizontalalignment='right')

    gc.save(f'{outpath}/{name}/{name}_XRT_vtp.png', dpi=90)

    gc.close()

    ### optical imaging ###
    if optimage:
        # make sure the links aren't broken
        if os.path.exists(f'{outpath}/{name}/{name}_{survey}stack.jpg'):
            ending = 'stack.jpg'
        elif os.path.exists(f'{outpath}/{name}/{name}_{survey}stack_irg.tiff'):
            ending = 'stack_irg.tiff'
        else:
            return

        # show the figure
        gc = aplpy.FITSFigure(img, figsize=(10, 10))
        try:
            gc.show_rgb(f'{outpath}/{name}/{name}_{survey}{ending}')
        except FileNotFoundError:
            gc.show_grayscale(stretch='arcsinh', pmin=1, pmax=98)
            gc.set_theme('publication')

        #gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm')
        #gc.set_tick_labels_size('small')

        # add all the sources -- the float is the pixscale in deg
#         r1 = detects['R'].data[:,0] * 6.548089E-04
#         r2 = detects['R'].data[:,1] * 6.548089E-04
#         gc.show_ellipses(detects['RA'], detects['DEC'], r1, r2, detects['ROTANG'],
#                          coords_frame='world', edgecolor='cyan')

        gc.show_regions(f'{outpath}/{name}/{name}_vtp.reg')

        # add PSZ info and circles
        gc.show_circles(info['RA'],
                        info['DEC'],
                        2 / 60,
                        linestyle='--',
                        edgecolor='#e24a33',
                        facecolor='none',
                        path_effects=[
                            pe.Stroke(linewidth=1.2, foreground='white'),
                            pe.Normal()
                        ])
        gc.show_circles(info['RA'],
                        info['DEC'],
                        5 / 60,
                        linestyle='-',
                        edgecolor='#e24a33',
                        facecolor='none',
                        path_effects=[
                            pe.Stroke(linewidth=1.2, foreground='white'),
                            pe.Normal()
                        ])
        gc.show_markers(info['RA'],
                        info['DEC'],
                        marker='*',
                        s=150,
                        layer='psz',
                        edgecolor='#e24a33',
                        path_effects=[
                            pe.Stroke(linewidth=1.2, foreground='white'),
                            pe.Normal()
                        ])

        xo, yo = (0.05, 0.05)
        # write the exposure time
        gc.add_label(xo,
                     yo,
                     text,
                     relative=True,
                     fontsize=18,
                     color='white',
                     horizontalalignment='left')

        # write redshift
        gc.add_label(xo,
                     yo + 0.03,
                     ztext,
                     relative=True,
                     fontsize=18,
                     color='white',
                     horizontalalignment='left')

        # write legend
        xo, yo = (0.95, 0.05)
        gc.add_label(xo,
                     yo,
                     'Extended',
                     relative=True,
                     fontsize=18,
                     color='magenta',
                     horizontalalignment='right')
        gc.add_label(xo,
                     yo + 0.03,
                     'P-Source',
                     relative=True,
                     fontsize=18,
                     color='yellow',
                     horizontalalignment='right')

        gc.save(f'{outpath}/{name}/{name}_OP_vtp.png', dpi=90)

        gc.close()

    return


# In[ ]:

# get file data
data = load_PSZcatalog()
data = data.sort_index(axis=1)

outpath = './data_full'

arr = [{'name': n.replace(' ', '_'), 'outpath': outpath} for n in data['NAME']]
parallel_process(arr, show_sources, use_kwargs=True, n_jobs=10)

cp_results(outpath)

# In[ ]:

outpath = './data_full'
name = 'PSZ2_G287.96-32.99'

# In[ ]:

show_sources(name, outpath)

# In[ ]:

data

# In[ ]:
