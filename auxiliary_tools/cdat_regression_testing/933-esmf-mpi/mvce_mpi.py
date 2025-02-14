from mpi4py import MPI

# Initialize MPI
if not MPI.Is_initialized():
    print("MPI is not initialized. Initializing.")
    MPI.Init()
else:
    print("MPI is already initialized.")

# Perform your MPI operations here
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(f"MPI process with rank {rank} out of {size} is running.")

# Finalize MPI
if not MPI.Is_finalized():
    print("MPI is not finalized. Finalizing.")
    MPI.Finalize()
else:
    print("MPI is already finalized.")