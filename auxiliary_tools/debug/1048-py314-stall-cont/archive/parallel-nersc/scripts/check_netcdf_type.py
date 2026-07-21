import glob
from netCDF4 import Dataset


def check_netcdf_type(file_path):
    try:
        with Dataset(file_path, "r") as ds:
            file_format = ds.file_format
        return file_format
    except Exception as e:
        return f"Error: {e}"


def main():
    pattern = "/global/cfs/projectdirs/e3sm/vo13/1048-py314-stall-cont/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/180x360_aave/clim/30yr/*.nc"
    files = glob.glob(pattern)
    if not files:
        print("No files found.")
        return

    for f in files:
        file_type = check_netcdf_type(f)
        print(f"{f}: {file_type}")


if __name__ == "__main__":
    main()
