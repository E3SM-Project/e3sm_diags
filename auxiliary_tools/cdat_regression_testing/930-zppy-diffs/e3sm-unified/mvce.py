import dask
import esmpy as ESMF
import numpy as np
from dask import bag
from xesmf.backend import Grid

config = {"scheduler": "processes", "multiprocessing.context": "fork"}

def grid_code():
    Grid(
        np.array((576,361)),
        staggerloc=ESMF.StaggerLoc.CENTER,
        coord_sys=ESMF.CoordSys.SPH_DEG,
        num_peri_dims=None
    )

run_workers = [1, 2]


for num_workers in run_workers:
    with dask.config.set(config):
            results = bag.map(grid_code).compute(num_workers=num_workers)
