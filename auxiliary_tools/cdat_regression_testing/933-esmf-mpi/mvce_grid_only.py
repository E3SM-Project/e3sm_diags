import dask
import esmpy as ESMF
import numpy as np
from dask import bag
from xesmf.backend import Grid
import multiprocessing

# Explicitly initialize ESMF in the main process
ESMF.Manager(debug=True)

def xesmf_code(_):
    Grid(
        np.array((576, 361)),
        staggerloc=ESMF.StaggerLoc.CENTER,
        coord_sys=ESMF.CoordSys.SPH_DEG,
        num_peri_dims=None
    )

if __name__ == '__main__':
    multiprocessing.set_start_method('forkserver', force=True)  # Avoid 'fork'

    config = {"scheduler": "processes", "multiprocessing.context": "forkserver"}
    with dask.config.set(config):
        print(f"Running with 2 workers and 'forkserver'")
        print("*"*50)
        try:
            b = bag.from_sequence(range(2))
            results = b.map(xesmf_code).compute(num_workers=2)
        except Exception as e:
            print(f"{type(e).__name__}: Message: {str(e)}")
        else:
            print("Done")