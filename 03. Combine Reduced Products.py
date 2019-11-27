#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
from glob import glob
from tqdm import tqdm_notebook

import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

# parallel processor
from utilities import parallel_process

# In[ ]:


def combineXRT(name, outpath):
    ''' This combines the individual observations

    IMPORTANT! The scripts that this (and other functions) creat are designed
    to be run from the same directory as this notebook. They WILL NOT work
    if you try to run them from the individual data directories.

    '''

    if not os.path.isdir(f'{outpath}/{name}'):
        return

    # find the x-ray files
    files = glob(f'{outpath}/{name}/reduced/**/*xpcw[1-4]po_cl.evt',
                 recursive=True)

    # sort the files so we can control their input order
    # this is used to make sure all of the data products are the same
    # when we are combining observations

    files = np.sort(files)

    if len(files) < 1:
        return

    # write xsel.in
    with open(f'{outpath}/{name}/{name}_xsel.in', 'w') as f:
        for i, f_in in enumerate(files):
            f_parts = f_in.split('/')
            if i == 0:
                # this is the session name... specify random letters to run many at once.
                f.writelines(f'{name}\n')
                f.writelines('read events\n')
                # set the data directory
                f.writelines('/'.join(f_parts[:3]) + '\n')
                # first entry
                f.writelines('/'.join(f_parts[3:]) + '\n')
                f.writelines('yes\n')
                continue

            f.writelines('read events\n')
            f.writelines('/'.join(f_parts[3:]) + '\n')
            # if you try to read more than 20 exposures, it says "more?"
            if i >= 19:
                f.writelines('\n')
            if i >= 42:
                f.writelines('\n')
            if i >= 65:
                f.writelines('\n')
            if i >= 88:
                f.writelines('\n')
            if i >= 111:
                f.writelines('\n')
            if i >= 134:
                f.writelines('\n')
            if i >= 157:
                f.writelines('\n')
            if i >= 180:
                f.writelines('\n')
            if i >= 203:
                f.writelines('\n')
            if i >= 226:
                f.writelines('\n')
            if i >= 249:
                f.writelines('\n')
        f.writelines('extract events\n')
        f.writelines(f'save events {outpath}/{name}/{name}_events.fits\n')
        if os.path.isfile(f'{outpath}/{name}/{name}_events.fits'):
            f.writelines('yes\n')
        f.writelines('yes\n')

        f.writelines('set phaname PI\n')
        # here we are going to make a few binned images for a few different energy ranges
        # energies in loop
        for eng in [200, 300, 400, 500, 600]:
            f.writelines(f'filter pha_cutoff 50 {eng}\n')

            # save non-binned image -- the yes's are to overwrite if file is already there
            f.writelines('set xybinsize 1\n')
            f.writelines('extract image\n')
            f.writelines(
                f'save image {"/".join(f_parts[:3])}/{name}_img_50-{eng}.fits\n'
            )
            if os.path.isfile(f'{outpath}/{name}/{name}_img_50-{eng}.fits'):
                f.writelines('yes\n')

            # save binned image -- see above
            f.writelines('set xybinsize 8\n')
            f.writelines('extract image\n')
            f.writelines(
                f'save image {"/".join(f_parts[:3])}/{name}_img_50-{eng}_bl8.fits\n'
            )
            if os.path.isfile(
                    f'{outpath}/{name}/{name}_img_50-{eng}_bl8.fits'):
                f.writelines('yes\n')

            f.writelines('set xybinsize 4\n')
            f.writelines('extract image\n')
            f.writelines(
                f'save image {"/".join(f_parts[:3])}/{name}_img_50-{eng}_bl4.fits\n'
            )
            if os.path.isfile(
                    f'{outpath}/{name}/{name}_img_50-{eng}_bl4.fits'):
                f.writelines('yes\n')

        f.writelines('exit\n')
        f.writelines('no\n')

    # log the output
    log_file = f'{outpath}/{name}/{name}_xsel.log'

    # call xselect
    os.system(f'xselect < {outpath}/{name}/{name}_xsel.in > {log_file}')

    return


# In[ ]:


def combineXRT_exp(name, outpath):
    ''' This combines the exposure maps'''

    if not os.path.isdir(f'{outpath}/{name}'):
        return

    # find the x-ray files
    files = glob(f'{outpath}/{name}/reduced/**/*xpcw[1-4]po_ex.img',
                 recursive=True)

    # sort the observations -- see above
    files = np.sort(files)

    if len(files) < 1:
        return name

    # remove the old file if it is there
    if os.path.isfile(f'{outpath}/{name}/{name}_exp.fits'):
        os.remove(f'{outpath}/{name}/{name}_exp.fits')

    # write xsel.in
    with open(f'{outpath}/{name}/{name}_ximg_exp.in', 'w') as f:
        for i, f_in in enumerate(files):
            f_parts = f_in.split('/')
            f.writelines(f'read {f_in}\n')
            if i == 0:
                continue
            f.writelines('sum\n')
            f.writelines('save\n')

        f.writelines(f'write/fits {"/".join(f_parts[:3])}/{name}_exp.fits\n')

        f.writelines('exit\n')

    # log the output
    log_file = f'{outpath}/{name}/{name}_ximg_exp.log'
    # call ximage
    os.system(f'ximage < {outpath}/{name}/{name}_ximg_exp.in > {log_file}')

    return name


# In[ ]:

# get file data
data = load_PSZcatalog()
data = data.sort_index(axis=1)

outpath = './data_full'

arr = [{'name': n.replace(' ', '_'), 'outpath': outpath} for n in data['NAME']]
parallel_process(arr, combineXRT, use_kwargs=True, n_jobs=6)
parallel_process(arr, combineXRT_exp, use_kwargs=True, n_jobs=6)

# In[ ]:

outpath = './data_full'
name = 'PSZ2_G027.99-69.85'

# In[ ]:

combineXRT_exp(name, outpath)

# In[ ]:
