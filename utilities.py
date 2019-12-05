from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from tqdm import tqdm_notebook
import os
import sys
import subprocess
import shlex


def parallel_process(array,
                     function,
                     n_jobs=None,
                     use_kwargs=False,
                     front_num=0):
    """
        A parallel version of the map function with a progress bar.

        Args:
            array (array-like): An array to iterate over.
            function (function): A python function to apply to the elements of
            array n_jobs (int, default=16): The number of cores to use
            use_kwargs (boolean, default=False): Whether to consider the
            elements of array as dictionaries of keyword arguments to function
            front_num (int, default=3): The number of iterations to run
            serially before kicking off the parallel job. This can be useful
            for catching bugs
        Returns:
            [function(array[0]), function(array[1]), ...]

    """

    # We run the first few iterations serially to catch bugs
    if front_num > 0:
        front = [
            function(**a) if use_kwargs else function(a)
            for a in array[:front_num]
        ]

    # If we set n_jobs to 1, just run a list comprehension.
    # This is useful for benchmarking and debugging.
    if n_jobs == 1:
        out = [
            function(**a) if use_kwargs else function(a)
            for a in tqdm_notebook(array[front_num:])
        ]
        if all(v is None for v in out):
            return
        else:
            return out

    #Assemble the workers
    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
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
        # Print out the progress as tasks complete
        out = []
#        for f in tqdm_notebook(as_completed(futures),
#                               desc=f'{function.__name__}',
#                               **kwargs):
        for f in tqdm(as_completed(futures),
                                desc=f'{function.__name__}',
                                **kwargs):
            try:
                out.append(f.result())
            except Exception as e:
                out.append(e)

        if front_num:
            if all(v is None for v in front + out):
                return
            else:
                return front + out
        elif all(v is None for v in out):
            return
        else:
            return out


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


def system_call(cmd, checkexe=False, shlexify=True):
    if shlexify:
        args = shlex.split(cmd)
    else:
        args = cmd
        
    if checkexe:
        check_exe(args[0])

    with subprocess.Popen(args,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True,
                          shell=True) as proc:

        #proc.wait(timeout=60)
        stdout, stderr = proc.communicate()

    return stdout, stderr


def set_env():
    os.environ['PFILES'] = f"/tmp/{os.getpid()}.tmp/pfiles;$HEADAS/syspfiles"
    pfiles = f'/tmp/{os.getpid()}.tmp/pfiles'
    if not os.path.isdir(pfiles):
        os.makedirs(pfiles)
