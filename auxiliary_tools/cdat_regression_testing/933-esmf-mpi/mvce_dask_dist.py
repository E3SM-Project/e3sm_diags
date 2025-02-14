# 1. srun -N 1 -t 01:00:00 --pty bash
# 2. source /lcrc/soft/climate/e3sm-unified/test_e3sm_unified_1.11.0rc6_chrysalis.sh
# 3. ssh -4 -L 5678:localhost:5678 user@compute-node
# 4. Run and Debug View -> Python: Remote Attach -> Attach to Remote

import esmpy as ESMF
import numpy as np
from dask.distributed import Client
from xesmf.backend import Grid

ESMF.Manager(debug=True, mpi=False)


import debugpy
debugpy.listen(("0.0.0.0", 5678))
print("Waiting for debugger attach...")
debugpy.wait_for_client()

def xesmf_code(_):
    Grid(
        np.array((576, 361)),
        staggerloc=ESMF.StaggerLoc.CENTER,
        coord_sys=ESMF.CoordSys.SPH_DEG,
        num_peri_dims=None
    )

# Open MPI has detected that this process has attempted to initialize
# MPI (via MPI_INIT or MPI_INIT_THREAD) more than once.  This is
# # erroneous.
if __name__ == "__main__":
    run_workers = [2]  # Example worker counts
    client = Client(processes=False, n_workers=2)

    for num_workers in run_workers:
        print(f"\nRunning with {num_workers} workers using Dask distributed")
        print("*" * 50)
        try:
            futures = client.map(xesmf_code, range(num_workers))
            results = client.gather(futures)
        except Exception as e:
            print(f"{type(e).__name__}: Message: {str(e)}")
        else:
            print("Done")


    client.close()