import dask
import esmpy as ESMF
import numpy as np
from dask import bag
from xesmf.backend import Grid

# import debugpy
# debugpy.listen(("0.0.0.0", 5678))
# print("Waiting for debugger attach...")
# debugpy.wait_for_client()

config = {"scheduler": "processes", "multiprocessing.context": "fork"}

def grid_code(_):
    Grid(
        np.array((576, 361)),
        staggerloc=ESMF.StaggerLoc.CENTER,
        coord_sys=ESMF.CoordSys.SPH_DEG,
        num_peri_dims=None
    )

run_workers = [1, 2]

for num_workers in run_workers:
    with dask.config.set(config):
        b = bag.from_sequence(range(num_workers))
        results = b.map(grid_code).compute()
