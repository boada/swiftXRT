import os
import numpy as np
import subprocess
from glob import glob
import aplpy
from time import sleep
import matplotlib.pyplot as plt

def combineXRT(name):
    # find the x-ray files
    files = glob('./data/{}/**/*xpcw3po_cl.evt.gz'.format(name), recursive=True)

    if len(files) < 1:
        return

    # remove the old file if it is there
    if os.path.isfile('./data/{}/{}_xsel.in'.format(name, name)):
        os.remove('./data/{}/{}_xsel.in'.format(name, name))

    # write xsel.in
    with open('./data/{}/{}_xsel.in'.format(name, name), 'w') as f:
        for i, f_in in enumerate(files):
            f_parts = f_in.split('/')
            if i == 0:
                f.writelines('xsel\n')
                f.writelines('read events\n')
                # set the data directory
                f.writelines('/'.join(f_parts[:3]) + '\n')
                # first entry
                f.writelines('/'.join(f_parts[3:]) + '\n')
                f.writelines('yes\n')
                continue

            f.writelines('read events\n')
            f.writelines('/'.join(f_parts[3:]) + '\n')

        f.writelines('extract events\n')
        f.writelines('save events {}/{}_events.fits\n'.format('/'.join(
                                                        f_parts[:3]), name))
        f.writelines('yes\n')
        f.writelines('exit\n')
        f.writelines('no\n')

    # call xselect
    os.system('xselect < ./data/{}/{}_xsel.in'.format(name, name))

    return

def mk_contours(name):

    if not os.path.isfile('./data/{}/{}_events.fits'.format(name, name)):
        return

    # start ds9
    p = subprocess.Popen(['ds9', './data/{}/{}_events.fits'.format(name, name)])

    ds9_cmds = ['xpaset -p ds9 bin factor 8',
                "xpaset -p ds9 bin filter 'pi=50:200'",
                'xpaset -p ds9 smooth 3',
                'xpaset -p ds9 contour yes',
                'xpaset -p ds9 contour smooth 1',
                'xpaset -p ds9 contour color white',
                'xpaset -p ds9 contour convert',
                'xpaset -p ds9 regions save ./data/{}/{}_contours.reg'.format(name,
                                                                         name),
                'xpaset -p ds9 exit']

    # wait for ds9 to start
    print("waiting 5 seconds for ds9 to start")
    sleep(5)

    for cmd in ds9_cmds:
        os.system(cmd)

    p.wait()
    p.kill()

    return

def put_contours(name):

    if not os.path.isfile('./data/{}/{}_PS1stack_i.fits'.format(name, name)):
        return

    if not os.path.isfile('./data/{}/{}_contours.reg'.format(name, name)):
        return

    gc = aplpy.FITSFigure('./data/{}/{}_PS1stack_i.fits'.format(name, name))
    gc.show_rgb('./data/{}/{}irg.tiff'.format(name, name))
    gc.show_regions('./data/{}/{}_contours.reg'.format(name, name))
    plt.tight_layout()
    gc.save('./data/{}/{}_contours.png'.format(name, name))

    plt.close()

    return


# get file data
data = np.genfromtxt('./PSZcatalogs/PSZ2_unconfirmed_catalog - Master.csv',
           delimiter=',', names=True, dtype=None)

for name in data['Name']:
    print(name.decode())
    #combineXRT(name.decode())
    #mk_contours(name.decode())
    put_contours(name.decode())
