import os
import numpy as np
import subprocess
from glob import glob
import aplpy
from time import sleep
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

def combineXRT(name):
    # find the x-ray files
    files = glob('./data/{}/**/*xpcw3po_cl.evt.gz'.format(name), recursive=True)

    if len(files) < 1:
        return

    # remove the old file if it is there
    if os.path.isfile('./data/{}/{}_events.fits'.format(name, name)):
        os.remove('./data/{}/{}_events.fits'.format(name, name))
    if os.path.isfile('./data/{}/{}_img_50-200_bl8.fits'.format(name, name)):
        os.remove('./data/{}/{}_img_50-200_bl8.fits'.format(name, name))

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

        f.writelines('set phaname PI\n')
        f.writelines('filter pha_cutoff 50 200\n')
        f.writelines('set xybinsize 8\n')

        f.writelines('extract image\n')
        f.writelines('save image {}/{}_img_50-200_bl8.fits\n'.format('/'.join(
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
    print("waiting 2 seconds for ds9 to start")
    sleep(2)

    for cmd in ds9_cmds:
        os.system(cmd)

    p.wait()
    p.kill()

    return

def mk_contours2(name):

    if not os.path.isfile('./data/{}/{}_img_50-200_bl8.fits'.format(name,
                                                                    name)):
        return

    if not os.path.isfile('./data/{}/{}_img_50-200_bl8.det'.format(name,
                                                                    name)):
        return

    # start ds9
    p = subprocess.Popen(['ds9',
                          './data/{}/{}_img_50-200_bl8.fits'.format(name,
                                                                    name)])

    # figure out the contour levels
    with open('./data/{}/{}_img_50-200_bl8.det'.format(name, name), 'r') as f:
        for l in f.readlines():
            if 'Back' in l:
                background = float(l.split(':')[-1])
                break

    # now scale the background appropriately -- 5 sigma
    cat = fits.getdata('./data/{}/{}xrt.fits'.format(name, name))
    exp_time = cat['xrt_exposure'].sum()

    # 8x8 binning
    min_pixel = background * 8 * 8 * exp_time * 5

    # now figure out the top level
    img = fits.getdata('./data/{}/{}_img_50-200_bl8.fits'.format(name, name))
    # smooth the image
    kernal = Gaussian2DKernel(stddev=1)
    smoothed_img = convolve(img, kernal)
    max_pixel = smoothed_img.max()

    # generate contour levels using a sqrt scale
    # i got this from the ds9 source code
    nlvls = 5
    lvls = np.linspace(0, 1, nlvls) / (nlvls - 1)
    lvls = lvls**2 * (max_pixel - min_pixel) + min_pixel
    lvls_str = [str(i) for i in lvls]
    lvls_str = ' '.join(lvls_str)

    ds9_cmds = ['xpaset -p ds9 smooth 3',
                'xpaset -p ds9 contour yes',
                'xpaset -p ds9 contour nlevels 5',
                'xpaset -p ds9 contour levels "{{{}}}"'.format(lvls_str),
                'xpaset -p ds9 contour smooth 1',
                'xpaset -p ds9 contour color white',
                'xpaset -p ds9 contour convert',
                'xpaset -p ds9 regions save ./data/{}/{}_contours.reg'.format(name,
                                                                         name),
                'xpaset -p ds9 exit']

    # wait for ds9 to start
    print("waiting 2 seconds for ds9 to start")
    sleep(2)

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

    cat = fits.getdata('./data/{}/{}xrt.fits'.format(name, name))

    exp_time = cat['xrt_exposure'].sum()
    text = 'exp time: {}s'.format(exp_time)

    gc = aplpy.FITSFigure('./data/{}/{}_PS1stack_i.fits'.format(name, name))
    gc.show_rgb('./data/{}/{}irg.tiff'.format(name, name))
    gc.show_regions('./data/{}/{}_contours.reg'.format(name, name))
    plt.tight_layout()

    # add exposure info
    xo, yo = (80, 80)
    plt.text(xo + 2, yo + 2, text, color='black', fontsize=18)
    plt.text(xo, yo, text, color='white', fontsize=18)

    gc.save('./data/{}/{}_contours.png'.format(name, name))

    plt.close()

    return

def put_contours2(name):
    if not os.path.isfile('./data/{}/{}_img_50-200_bl8.fits'.format(name,
                                                                    name)):
        return

    if not os.path.isfile('./data/{}/{}_img_50-200_bl8.det'.format(name,
                                                                    name)):
        return

    # figure out the contour levels
    with open('./data/{}/{}_img_50-200_bl8.det'.format(name, name), 'r') as f:
        for l in f.readlines():
            if 'Back' in l:
                background = float(l.split(':')[-1])
                break

    # now scale the background appropriately -- 5 sigma
    cat = fits.getdata('./data/{}/{}xrt.fits'.format(name, name))
    exp_time = cat['xrt_exposure'].sum()

    # 8x8 binning
    min_pixel = background * 8 * 8 * exp_time * 5

    # now figure out the top level
    img = fits.getdata('./data/{}/{}_img_50-200_bl8.fits'.format(name, name))
    # smooth the image
    kernal = Gaussian2DKernel(stddev=1, x_size=3, y_size=3)
    smoothed_img = convolve(img, kernal)
    max_pixel = smoothed_img.max()

    # generate contour levels using a sqrt scale
    # i got this from the ds9 source code
    nlvls = 5
    lvls = np.linspace(0, 1, nlvls) / (nlvls - 1)
    lvls = lvls**2 * (max_pixel - min_pixel) + min_pixel
    lvls_str = [str(i) for i in lvls]
    lvls_str = ' '.join(lvls_str)

    cat = fits.getdata('./data/{}/{}xrt.fits'.format(name, name))

    exp_time = cat['xrt_exposure'].sum()
    text = 'exp time: {}s'.format(exp_time)

    gc = aplpy.FITSFigure('./data/{}/{}_img_50-200_bl8.fits'.format(name, name))
    gc.show_grayscale()
    gc.show_contour(smooth=3, levels=lvls, cmap='viridis')
    plt.tight_layout()

    # add exposure info
    xo, yo = (80, 80)
    plt.text(xo + 2, yo + 2, text, color='black', fontsize=18)
    plt.text(xo, yo, text, color='white', fontsize=18)

    gc.save('./data/{}/{}_XRT_contours.png'.format(name, name))

    plt.close()

    return

def source_detec(name):
    if not os.path.isfile('./data/{}/{}_img_50-200_bl8.fits'.format(name,
                                                                    name)):
        return

    with open('./data/{}/{}_ximg.in'.format(name, name), 'w') as f:
        f.writelines('read ./data/{}/{}_img_50-200_bl8.fits\n'.format(name,
                                                                      name))

        f.writelines('detect/snr_threshold=3/'
                'fitsdet={{./data/{}/{}_img_50-200_bl8.det.fits}}\n'.format(name,
                                                                        name))
        f.writelines('exit\n')

    # call xselect
    os.system('ximage < ./data/{}/{}_ximg.in'.format(name, name))

    return


if __name__ == "__main__":
    # get file data
    data = np.genfromtxt('./PSZcatalogs/PSZ2_unconfirmed_catalog - proc2.csv',
            delimiter=',', names=True, dtype=None)

    for name in data['Name']:
        print(name.decode())
        #combineXRT(name.decode())
        source_detec(name.decode())
        #mk_contours(name.decode())
        mk_contours2(name.decode())
        #put_contours(name.decode())
        put_contours2(name.decode())
