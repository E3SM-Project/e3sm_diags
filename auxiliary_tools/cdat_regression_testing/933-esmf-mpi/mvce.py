# srun -N 1 -t 01:00:00 --pty bash
# source /lcrc/soft/climate/e3sm-unified/test_e3sm_unified_1.11.0rc7_chrysalis.sh
# ssh -4 -L 5678:localhost:5678 user@compute-node
# Run and Debug View -> Python: Remote Attach -> Attach to Remote
#
# Compute env (Conda): source /lcrc/soft/climate/e3sm-unified/base/etc/profile.d/conda.sh && conda activate e3sm_unified_1.11.0rc7_chrysalis
# Compute env (Spack): source /lcrc/soft/climate/e3sm-unified/test_e3sm_unified_1.11.0rc7_chrysalis.sh


import dask
import esmpy as ESMF
import numpy as np
from dask import bag
from xesmf.backend import Grid

# Initialize ESMF with MPI only once
ESMF.Manager(debug=True)

def xesmf_code(_):
    Grid(
        np.array((576, 361)),
        staggerloc=ESMF.StaggerLoc.CENTER,
        coord_sys=ESMF.CoordSys.SPH_DEG,
        num_peri_dims=None
    )

config = {"scheduler": "processes", "multiprocessing.context": "fork"}

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

# Finalize ESMF with MPI.
ESMF.ESMP_Finalize()