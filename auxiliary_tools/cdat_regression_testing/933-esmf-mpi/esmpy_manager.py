# srun -N 1 -t 01:00:00 --pty bash
# source /lcrc/soft/climate/e3sm-unified/test_e3sm_unified_1.11.0rc7_chrysalis.sh
# ssh -4 -L 5678:localhost:5678 user@compute-node
# Run and Debug View -> Python: Remote Attach -> Attach to Remote
#
# Compute env (Conda): source /lcrc/soft/climate/e3sm-unified/base/etc/profile.d/conda.sh && conda activate e3sm_unified_1.11.0rc7_chrysalis
# Compute env (Spack): source /lcrc/soft/climate/e3sm-unified/test_e3sm_unified_1.11.0rc7_chrysalis.sh

import os

# ESMF_MPIRUN: Disables the use of mpirun by setting it to "no".
os.environ["ESMF_MPIRUN"] = "no"
# Set the communication method to mpiuni to disable MPI
os.environ["ESMF_COMM"] = "mpiuni"
import esmpy as ESMF

# Initialize ESMF Manager with MPI. -- still crashes
ESMF.Manager(debug=True)

# Finalize ESMF with MPI. -- still crashes
ESMF.ESMP_Finalize()