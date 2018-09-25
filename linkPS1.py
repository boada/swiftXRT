import os
import numpy as np
import pandas as pd

# get file data
#data = np.genfromtxt('./PSZcatalogs/PSZ2_unconfirmed_catalog - Master.csv',
#           delimiter=',', names=True, dtype=None)

data = pd.read_csv('./PSZ_unconfirmed.csv')

PS1_dir = '/home/boada/Projects/planckClusters/data/extern/PS1'

for name in data['NAME']:
    print(name)
    name = name.replace(' ', '_')

    #name = ''.join(e for e in name if e.isalnum())
    if not os.path.isdir('./data/{}'.format(name)):
        os.makedirs('./data/{}'.format(name))

    relpath = os.path.relpath('{}/{}'.format(PS1_dir, name),
                              './data/{}'.format(name))

    target_files = ['_PS1stack_g.fits', '_PS1stack_r.fits', '_PS1stack_i.fits',
                    '_PS1stack_z.fits', '_PS1stack_y.fits', 'irg.tiff']

    for file in target_files:
        try:
            os.symlink('{}/{}/{}{}'.format(PS1_dir, name,
                                           name, file),
                       './data/{}/{}{}'.format(name, name,
                       file))
        except FileExistsError:
            pass


