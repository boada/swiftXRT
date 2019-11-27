#!/usr/bin/env python
# coding: utf-8

import numpy as np
from time import sleep
from astropy.io import fits
from astropy.time import Time

import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

# parallel processor
from utilities import parallel_process


def get_catalog(ra, dec, name, outpath, table='swiftmastr', oformat='fits'):

    # build the command
    cmd = (
        f'./browse_extract_wget.pl table={table} position={ra},{dec} '
        f'outfile={outpath}/{name}/{name}xrt.fits format={oformat} radius=10')

    if not os.path.isdir(f'{outpath}/{name}'):
        os.makedirs(f'{outpath}/{name}')

    #print(cmd)
    os.system(cmd)

    sleep(0.5)

def get_files(name, outpath):

    # at least one file was 0 bytes...
    try:
        cat = fits.getdata(f'{outpath}/{name}/{name}xrt.fits')
    except OSError:
        return

    if cat.shape[0] < 1:
        os.remove(f'{outpath}/{name}/{name}xrt.fits')
        os.rmdir(f'{outpath}/{name}')
        return

    for i in range(cat.shape[0]):
        obsid = cat['obsid'][i]
        t = Time(cat['start_time'][i], format='mjd')
        y, m = t.iso.split('-')[:2]

        wget_base = (
            f'wget -P {outpath}/{name} -q -nH --no-check-certificate --cut-dirs=5 '
            '-r -l0 -c -nc -np -R "index*" -erobots=off --retr-symlinks '
            f'https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{y}_{m}/')

        for prod in ['xrt', 'auxil', 'log']:
            wget = f'{wget_base}{obsid}/{prod}/'
            os.system(wget)

def write_info(outpath, idx, data):
    target = data.iloc[idx]['NAME'].replace(' ', '_')

    if not os.path.isdir(f'{outpath}/{target}'):
        return

    cols = data.columns.tolist()
    vals = data.iloc[idx]
    # convert the values into all strings
    vals_str = [
        f'{i:.8}' if isinstance(i, float) else i for i in vals.values.tolist()
    ]
    with open(f'{outpath}/{target}/{target}.info', 'w') as f:
        f.write(f'{",".join(cols)}\n')
        f.write(f'{",".join(vals_str)}\n')


# get file data
data = load_PSZcatalog()
data = data.sort_index(axis=1)

outpath = './data_full'

arr = [{
    'name': n.replace(' ', '_'),
    'outpath': outpath,
    'ra': ra,
    'dec': dec
} for n, ra, dec in zip(data['NAME'], data['RA'], data['DEC'])]

parallel_process(arr, get_catalog, use_kwargs=True, n_jobs=6)

arr = [{'name': n.replace(' ', '_'), 'outpath': outpath} for n in data['NAME']]
parallel_process(arr, get_files, use_kwargs=True, n_jobs=6)

for i, n in enumerate(data['NAME']):
    write_info(outpath, i, data)