import os
import numpy as np

# get file data
data = np.genfromtxt('./PSZcatalogs/PSZ2_unconfirmed_catalog - Master.csv',
           delimiter=',', names=True, dtype=None)

PS1_dir = '/home/boada/Projects/planckClusters/data/extern/PS1'


for name in data['Name']:
    print(name.decode())

    #name = ''.join(e for e in name.decode() if e.isalnum())
    if not os.path.isdir('./data/{}'.format(name.decode())):
        os.makedirs('./data/{}'.format(name.decode()))

    relpath = os.path.relpath('{}/{}'.format(PS1_dir, name.decode()),
                              './data/{}'.format(name.decode()))

    target_files = ['_PS1stack_g.fits', '_PS1stack_r.fits', '_PS1stack_i.fits',
                    '_PS1stack_z.fits', '_PS1stack_y.fits', 'irg.tiff']

    for file in target_files:
        try:
            os.symlink('{}/{}/{}{}'.format(PS1_dir, name.decode(),
                                           name.decode(), file),
                       './data/{}/{}{}'.format(name.decode(), name.decode(),
                       file))
        except FileExistsError:
            pass


