# %%
import xcdat as xc

filepaths = [
    "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/CLOUD_twpc3_200001_201412.nc"
]

ds = xc.open_mfdataset(
    filepaths,
    # add_bounds=["X", "Y", "T"],
    decode_times=True,
    use_cftime=True,
    coords="minimal",
    compat="override",
)

# %%
ds.load(scheduler="sync")

# %%
print(ds.CLOUD.data)

#               Array	        Chunk
# Bytes	        40.10 MiB	    320 B
# Shape	        (131400, 80)	(1, 80)
# Dask graph	131400 chunks in 2 graph layers
# Data type	    float32 numpy.ndarray

# The ncdump shows _ChunkSizes = 1, 80, which aligns with the Xarray info above.
# ncdump -hs "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/CLOUD_twpc3_200001_201412.nc"


        # float CLOUD(time, lev) ;
        #         CLOUD:mdims = 1 ;
        #         CLOUD:units = "fraction" ;
        #         CLOUD:long_name = "Cloud fraction" ;
        #         CLOUD:cell_methods = "time: point" ;
        #         CLOUD:basename = "CLOUD" ;
        #         CLOUD:missing_value = 1.e+20f ;
        #         CLOUD:_FillValue = 1.e+20f ;
        #         CLOUD:_Storage = "chunked" ;
        #         CLOUD:_ChunkSizes = 1, 80 ;
        #         CLOUD:_Endianness = "little" ;


# %%

# Issue 2 - Time bounds are pre-defined with odd chunking scheme .
# The chunking scheme is (1, 2), resulting in 131400 chunks in 3 graph layers.
# 131400 chunk is a massive number of chunks for a 2D array, resulting in
# inefficient loading.
            Array	Chunk
Bytes	    2.01 MiB	16 B
Shape	    (131400, 2)	(1, 2)
Dask graph	131400 chunks in 3 graph layers
Data type	object numpy.ndarray

filepaths2 = [
    "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/post/atm/site/FLDS_sgpc1_200001_201412.nc"
]

ds = xc.open_mfdataset(
    filepaths2,
    add_bounds=["X", "Y", "T"],
    decode_times=True,
    use_cftime=True,
    coords="minimal",
    compat="override",
)

# %%
