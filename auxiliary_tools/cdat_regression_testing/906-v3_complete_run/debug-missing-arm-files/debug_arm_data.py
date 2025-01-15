# %%
"""
Issue 1 - Sub-optimal `CLOUD` and `time_bnds chunking
  * Related dataset: "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/CLOUD_twpc3_200001_201412.nc"
  * Dataset shape is (time: 131400, bound: 2, lev: 80)
  * `CLOUD` variable has chunks of (1, 80), resulting in 131400 chunks in 2 graph layers. (very bad, slow loading)
  * `time_bnds` has chunks of (1, 2), resulting in 131400 chunks in 3 graph layers. (very bad, slow loading)
"""

import xcdat as xc

ds = xc.open_mfdataset(
    [
        "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/CLOUD_twpc3_200001_201412.nc"
    ]
)

print(ds.CLOUD.data)
#               Array	        Chunk
# Bytes	        40.10 MiB	    320 B
# Shape	        (131400, 80)	(1, 80)
# Dask graph	131400 chunks in 2 graph layers
# Data type	    float32 numpy.ndarray

print(ds.time_bnds.data)
# Array	Chunk
# Bytes	2.01 MiB	16 B
# Shape	(131400, 2)	(1, 2)
# Dask graph	131400 chunks in 3 graph layers
# Data type	object numpy.ndarray


# %%
"""
Issue 2 - Sub-optimal `time_bnds` chunking
  * Related dataset: "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/FLDS_sgpc1_200001_201412.nc"
  * Dataset shape is (time: 131400, bound: 2, lev: 80)
  * `FLDS` variable has chunks of (1019,), resulting in 129 in 2 graph layers (okay)
  * `time_bnds` has chunks of (1, 2), resulting in 131400 chunks in 3 graph layers (very bad, slow loading)
"""

ds2 = xc.open_mfdataset(
    [
        "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/FLDS_sgpc1_200001_201412.nc"
    ]
)

print(ds2.FLDS.data)
#               Array	Chunk
# Bytes	        513.28 kiB	3.98 kiB
# Shape	        (131400,)	(1019,)
# Dask graph	129 chunks in 2 graph layers
# Data type	    float32 numpy.ndarray

print(ds2.time_bnds.data)
#               Array	Chunk
# Bytes	        2.01 MiB	16 B
# Shape	        (131400, 2)	(1, 2)
# Dask graph	131400 chunks in 3 graph layers
# Data type	    object numpy.ndarray

# %%
"""
Issue 3 - Sub-optimal `time_bnds` chunking
  * Related dataset: "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/PRECT_sgpc1_200001_201412.nc"
  * Dataset shape is (time: 131400, bound: 2, lev: 80)
  * `PRECT` variable has chunks of (1019,), resulting in 129 in 2 graph layers (okay)
  * `time_bnds` has chunks of (1, 2), resulting in 131400 chunks in 3 graph layers (very bad, slow loading)
"""

ds3 = xc.open_mfdataset(
    [
        "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/PRECT_sgpc1_200001_201412.nc"
    ]
)

print(ds3.PRECT.data)
#               Array	Chunk
# Bytes	        513.28 kiB	3.98 kiB
# Shape	        (131400,)	(1019,)
# Dask graph	129 chunks in 2 graph layers
# Data type	    float32 numpy.ndarray

print(ds3.time_bnds.data)
#               Array	Chunk
# Bytes	        2.01 MiB	16 B
# Shape	        (131400, 2)	(1, 2)
# Dask graph	131400 chunks in 3 graph layers
# Data type	    object numpy.ndarray

# %%
"""
Issue 4 - Sub-optimal `num_a1` and `time_bnds` chunking
  * Related dataset: "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/num_a1_enac1_200001_201412.nc"
  * Dataset shape is (time: 131400, bound: 2, lev: 80)
  * `num_a1` variable has chunks of (1, 80), resulting in 131400 chunks in 2 graph layers. (very bad, slow loading)
  * `time_bnds` has chunks of (1, 2), resulting in 131400 chunks in 3 graph layers (very bad, slow loading)
"""

ds4 = xc.open_mfdataset(
    [
        "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/num_a1_enac1_200001_201412.nc"
    ]
)

print(ds4.num_a1.data)
#               Array	Chunk
# Bytes	        40.10 MiB	320 B
# Shape	        (131400, 80)	(1, 80)
# Dask graph	131400 chunks in 2 graph layers
# Data type	    float32 numpy.ndarray

print(ds4.time_bnds.data)
#               Array	Chunk
# Bytes	        2.01 MiB	16 B
# Shape	        (131400, 2)	(1, 2)
# Dask graph	131400 chunks in 3 graph layers
# Data type	    object numpy.ndarray


# %%
"""
Issue 5 - Sub-optimal `time_bnds` chunking
  * Related dataset: "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/PRECT_twpc3_200001_201412.nc"
  * Dataset shape is (time: 131400, bound: 2, lev: 80)
  * `PRECT` variable has chunks of (1019,), resulting in 129 in 2 graph layers (okay)
  * `time_bnds` has chunks of (1, 2), resulting in 131400 chunks in 3 graph layers (very bad, slow loading)
"""

ds5 = xc.open_mfdataset(
    [
        "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/PRECT_twpc3_200001_201412.nc"
    ]
)

print(ds5.time_bnds.data)
#               Array	Chunk
# Bytes	        2.01 MiB	16 B
# Shape	        (131400, 2)	(1, 2)
# Dask graph	131400 chunks in 3 graph layers
# Data type	    object numpy.ndarray
