# 1. srun -N 1 -t 01:00:00 --pty bash
# 2. source /lcrc/soft/climate/e3sm-unified/test_e3sm_unified_1.11.0rc6_chrysalis.sh
# 3. ssh -4 -L 5678:localhost:5678 user@compute-node
# 4. Run and Debug View -> Python: Remote Attach -> Attach to Remote

import dask
import esmpy as ESMF
import numpy as np
from dask import bag
from xesmf.backend import Grid

# Initialize ESMF with MPI only once
ESMF.Manager(logkind=ESMF.LogKind.MULTI, debug=True)

def xesmf_code(_):
    Grid(
        np.array((576, 361)),
        staggerloc=ESMF.StaggerLoc.CENTER,
        coord_sys=ESMF.CoordSys.SPH_DEG,
        num_peri_dims=None
    )


#%%
# 1. Try running with fork
# Result: `concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending.`
config = {"scheduler": "processes", "multiprocessing.context": "fork"}
run_workers = [2]

with dask.config.set(config):
    print(f"Running with 2 workers and 'fork'")
    print("*"*50)
    try:
        b = bag.from_sequence(range(2))
        results = b.map(xesmf_code).compute(num_workers=2)
    except Exception as e:
        print(f"{type(e).__name__}: Message: {str(e)}")
    else:
        print("Done")

# # 2. Try running with spawn context
# Result: `concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending.`
# config_new = {"scheduler": "processes", "multiprocessing.context": "spawn"}

# if __name__ == "__main__":
#     with dask.config.set(config_new):
#         print(f"Running with 2 workers and 'spawn'")
#         print("*"*50)
#         try:
#             b = bag.from_sequence(range(2))
#             results = b.map(xesmf_code).compute(num_workers=2)
#         except Exception as e:
#             print(f"{type(e).__name__}: Message: {str(e)}")
#         else:
#             print("Done")
