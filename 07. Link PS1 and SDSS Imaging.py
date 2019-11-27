#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog
from tqdm import tqdm_notebook

# In[5]:

data = load_PSZcatalog()

PS1_dir = f'{os.environ["HOME"]}/Projects/planckClusters/data/extern/PS1'
SDSS_dir = f'{os.environ["HOME"]}/Projects/planckClusters/data/extern/SDSS'
DECaLS_dir = f'{os.environ["HOME"]}/Projects/planckClusters/data/extern/DECaLS'
DES_dir = f'{os.environ["HOME"]}/Projects/planckClusters/data/extern/DES'
outpath = './data_full'

for name in tqdm_notebook(data['NAME'], total=len(data['NAME'])):
    name = name.replace(' ', '_')

    if not os.path.exists(f'{outpath}/{name}'):
        continue

    if os.path.isdir(f'{PS1_dir}/{name}'):
        relpath = os.path.relpath(f'{PS1_dir}/{name}', f'{outpath}/{name}')

        target_files = [
            '_PS1stack_g.fits', '_PS1stack_r.fits', '_PS1stack_i.fits',
            '_PS1stack_z.fits', '_PS1stack_y.fits', '_PS1stack_irg.tiff'
        ]

        for file in target_files:
            try:
                os.symlink(f'{PS1_dir}/{name}/{name}{file}',
                           f'{outpath}/{name}/{name}{file}')
            except FileExistsError:
                pass

    if os.path.isdir(f'{SDSS_dir}/{name}'):
        relpath = os.path.relpath(f'{SDSS_dir}/{name}', f'{outpath}/{name}')

        target_files = [
            '_SDSSstack_g.fits', '_SDSSstack_r.fits', '_SDSSstack_i.fits',
            '_SDSSstack_z.fits', '_SDSSstack_irg.tiff'
        ]

        for file in target_files:
            try:
                os.symlink(f'{SDSS_dir}/{name}/{name}{file}',
                           f'{outpath}/{name}/{name}{file}')
            except FileExistsError:
                pass

    if os.path.isdir(f'{DECaLS_dir}/{name}'):
        relpath = os.path.relpath(f'{DECaLS_dir}/{name}', f'{outpath}/{name}')

        target_files = ['_DECaLSstack_r.fits', '_DECaLSstack.jpg']

        for file in target_files:
            try:
                os.symlink(f'{DECaLS_dir}/{name}/{name}{file}',
                           f'{outpath}/{name}/{name}{file}')
            except FileExistsError:
                pass

    if os.path.isdir(f'{DES_dir}/{name}'):
        relpath = os.path.relpath(f'{DES_dir}/{name}', f'{outpath}/{name}')

        target_files = ['_DESstack_r.fits', '_DESstack.jpg']

        for file in target_files:
            try:
                os.symlink(f'{DES_dir}/{name}/{name}{file}',
                           f'{outpath}/{name}/{name}{file}')
            except FileExistsError:
                pass

# In[ ]:
