import os

import dask
from dask.distributed import progress
import xarray as xr


import numpy as np 

# try to force flushing of memory
import gc, ctypes

#---------- cluster tools

def wait_cluster_ready(cluster, nworkers):
    """
    Wait for the client to be ready (all workers started)
    """
    import time, sys
    Nw = 0.
    while(len(cluster.scheduler.workers)<nworkers):
        p = len(cluster.scheduler.workers)/nworkers
        time.sleep(2)
        # print('{} % of the workers started'.format(int(p*100)))
        if p>=Nw:
            Nw = int(p*10)/10
            Nw+=0.1
    print('100 % of the workers started: {} workers'.format(nworkers))
    
def trim_memory() -> int:
    libc = ctypes.CDLL("libc.so.6")
    return libc.malloc_trim(0)

def clear_memory(client):
    client.run(gc.collect)
    client.run(trim_memory)
    return

def dask_compute_batch(computations, client, batch_size=None):
    """ breaks down a list of computations into batches
    """
    from dask.diagnostics import ProgressBar
    ProgressBar().register()

    # compute batch size according to number of workers
    if batch_size is None:
        # batch_size = len(client.scheduler_info()["workers"])
        batch_size = sum(list(client.nthreads().values()))
    # find batch indices
    total_range = range(len(computations))
    splits = max(1, np.ceil(len(total_range)/batch_size))
    batches = np.array_split(total_range, splits)
    # launch computations
    outputs = []
    for b in batches:
        out = dask.compute(*computations[slice(b[0], b[-1]+1)])
        progress(out)

        # try to manually clean up memory
        # https://coiled.io/blog/tackling-unmanaged-memory-with-dask/
        # client.run(gc.collect)
        client.run(trim_memory)  # should not be done systematically
        outputs.append(out)
    return outputs
    
