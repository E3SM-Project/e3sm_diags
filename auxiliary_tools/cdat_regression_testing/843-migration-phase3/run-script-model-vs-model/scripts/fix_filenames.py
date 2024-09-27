"""
A script to update the filename for annual_cycle_zonal_mean files to align
with the dev branch.

NOTE: Make sure to run this script before regression testing.
"""
import os


def replace_annual_cycle(root_dir):
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if "Annual-Cycle" in filename:
                old_file_path = os.path.join(dirpath, filename)
                new_file_name = filename.replace("Annual-Cycle", "ANNUALCYCLE-global")
                new_file_path = os.path.join(dirpath, new_file_name)

                print(f"Renaming file: {old_file_path} to {new_file_path}")
                os.rename(old_file_path, new_file_path)


# Replace 'your_root_directory' with the path to your root directory
root_directory = "/global/cfs/cdirs/e3sm/www/cdat-migration-fy24/main-model-vs-model/annual_cycle_zonal_mean"
replace_annual_cycle(root_directory)
