import os
import numpy as np
from time import sleep
from astropy.io import fits
from astropy.time import Time
import pandas as pd

def get_catalog(ra, dec, target, table='swiftmastr', oformat='fits'):

    # build the command
    cmd = "./browse_extract_wget.pl table={} position={},{} ".format(
                table, ra, dec)

    cmd += 'outfile=./data/{}/{}xrt.fits format={} radius=10'.format(target,
                                                                     target,
                                                                     oformat)

    print(cmd)
    os.system(cmd)

def get_files(target):

    cat = fits.getdata('./data/{}/{}xrt.fits'.format(target, target))

    if cat.shape[0] < 1:
        return

    for i in range(cat.shape[0]):
        obsid = cat['obsid'][i]
        t = Time(cat['start_time'][i], format='mjd')
        y, m = t.iso.split('-')[:2]

        wget = ('wget -P ./data/{} -q -nH --no-check-certificate --cut-dirs=5 '
                '-r -l0 -c -N -np -R "index*" -erobots=off --retr-symlinks '
                'https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{}_{}//'
                '{}/xrt/products/sw{}xpc_ex.img.gz'.format(target, y, m, obsid,
                                                           obsid))

        print(wget)
        os.system(wget)

        wget = ('wget -P ./data/{} -q -nH --no-check-certificate --cut-dirs=5 '
                '-r -l0 -c -N -np -R "index*" -erobots=off --retr-symlinks '
                'https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{}_{}//'
                '{}/xrt/event/sw{}xpcw3po_cl.evt.gz'.format(target, y, m,
                                                            obsid, obsid))

        print(wget)
        os.system(wget)


# get file data
data = np.genfromtxt('./PSZcatalogs/PSZ2_unconfirmed_catalog - current.csv',
           delimiter=',', names=True, dtype=None)

data = pd.read_csv('./PSZ_unconfirmed.csv')

for i, (ra, dec,
        name) in enumerate(zip(data['RA'], data['DEC'], data['NAME'])):
    #print(data['Name'][i])

    print(name)

    name = name.replace(' ', '_')

    #name = ''.join(e for e in name.decode() if e.isalnum())
    if not os.path.isdir('./data/{}'.format(name)):
        os.makedirs('./data/{}'.format(name))

    get_catalog(ra, dec, '{}'.format(name))
    get_files('{}'.format(name))

    sleep(1)
