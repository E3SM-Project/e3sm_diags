"""
# Need to do this before importing MPI.
import mpi4py
mpi4py.rc.initialize = False  # do not initialize MPI automatically
mpi4py.rc.finalize = False    # do not finalize MPI automatically
"""

# import mpi4py
# mpi4py.rc.initialize = False  # do not initialize MPI automatically
# mpi4py.rc.finalize = False    # do not finalize MPI automatically

# PMI2_Init failed to intialize.  Return code: 14
from mpi4py import MPI

