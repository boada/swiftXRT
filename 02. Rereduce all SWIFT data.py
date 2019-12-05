#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import subprocess
import tempfile

import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

# parallel processor
from utilities import parallel_process

# In[ ]:

# this allows us to run the pipeline not interactively.
os.environ['HEADASNOQUERY'] = ''
os.environ['HEADASPROMPT'] = '/dev/null'

# In[ ]:


def get_immediate_subdirectories(a_dir):
    ''' Get a list of a directorys immediate subdirectories'''
    return [
        os.path.join(a_dir, name) for name in os.listdir(a_dir)
        if os.path.isdir(os.path.join(a_dir, name))
    ]


def get_immediate_subfiles(a_dir):
    ''' Get a list of all the FILES in a directory'''
    return [
        os.path.join(a_dir, name) for name in os.listdir(a_dir)
        if os.path.isfile(os.path.join(a_dir, name))
    ]


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
    if verb:
        # it wasn't found
        print("# ERROR: Couldn't find %s" % exe, file=sys.stderr)
    raise FileNotFoundError(exe)
    return False


# In[ ]:


def run_pipeline(name, outpath, overwrite=False):

    if not os.path.isdir(f'{outpath}/{name}'):
        return

    if not check_exe('xrtpipeline', verb=False):
        return

    # set up enviroment
    if not os.path.isdir(f'/tmp/{os.getpid()}.tmp/pfiles'):
        os.makedirs(f'/tmp/{os.getpid()}.tmp/pfiles')

    # get a list of the XRT obs.
    obs = get_immediate_subdirectories(f'{outpath}/{name}')

    # if there aren't any observations, keep going
    if not len(obs):
        return

    if not os.path.isdir(f'{outpath}/{name}/reduced'):
        os.makedirs(f'{outpath}/{name}/reduced')

    # get a list of the reduced obs.
    reduc = get_immediate_subdirectories(f'{outpath}/{name}/reduced')

    for ob_dir in obs:
        ob_id = ob_dir.split('/')[-1]
        reduc_ids = [r.split('/')[-1] for r in reduc]
        if ob_id == 'reduced':
            continue
        if ob_id in reduc_ids and not overwrite:
            continue
        env_cmd = f'export PFILES="/tmp/{os.getpid()}.tmp/pfiles;{os.environ["HEADAS"]}/syspfiles" \n'
        pipe_cmd = (
            f'xrtpipeline indir={ob_dir} outdir={outpath}/{name}/reduced/{ob_id} steminputs=sw{ob_id} '
            'srcra=OBJECT srcdec=OBJECT datamode=PC cleanup=yes vigflag=yes clobber=yes '
            f' > {outpath}/{name}/reduced/{ob_id}_reduce.log \n')

        fd, path = tempfile.mkstemp()
        try:
            with os.fdopen(fd, 'w') as tmp:
                # do stuff with temp file
                tmp.write(env_cmd)
                tmp.write(pipe_cmd)
            subproc = subprocess.Popen(f'sh {path}', shell=True)
            subproc.wait()
        finally:
            os.remove(path)
    return name


# In[ ]:


def validate():
    PSZs = get_immediate_subdirectories('./data_full')
    for psz in tqdm_notebook(PSZs):
        # get a list of the XRT obs.
        for ob_dir in get_immediate_subdirectories(psz):
            # check for events and exposure map
            evts = False
            expm = False
            if ob_id == 'reduced':
                continue
            else:
                for f in get_immediate_subfiles(f'{psz}/reduced/{ob_id}'):
                    if 'xpcw3po_cl.evt' in f:
                        evts = True
                    elif 'xpcw3po_ex.img' in f:
                        expm = True
                    else:
                        continue

                if not evts and expm:
                    print(f'{psz}/reduced/{ob_id} NOT VALID')


# In[ ]:

# get file data
data = load_PSZcatalog()
data = data.sort_index(axis=1)

outpath = './data_full'

arr = [{
    'name': n.replace(' ', '_'),
    'outpath': outpath,
    'overwrite': True
} for n in data['NAME']]
parallel_process(arr, run_pipeline, use_kwargs=True, n_jobs=6)

# In[ ]:

#outpath = './data_full'
#name = 'PSZ2_G359.07-32.12'
#run_pipeline(name, outpath, True)
