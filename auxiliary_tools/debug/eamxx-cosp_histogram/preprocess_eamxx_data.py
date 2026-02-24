#!/usr/bin/env python
"""
Preprocess EAMxx COSP histogram files to add missing coordinate values.

This script copies files from input directory to output directory and adds
default coordinate values to COSP dimensions.
"""

import os
import glob
import shutil
import numpy as np
import xarray as xr

# Input and output directories
INPUT_DIR = "/pscratch/sd/t/terai/EAMxx/ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1/rgr/climo"
OUTPUT_DIR = "/pscratch/sd/c/chengzhu/EAMxx/ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1/rgr/climo"

# Default coordinate values for EAMxx COSP variables
DEFAULT_COORDS = {
    "cosp_prs": np.array([90000, 74000, 62000, 50000, 37500, 24500, 9000]),
    "cosp_tau": np.array([0.15, 0.8, 2.45, 6.5, 16.2, 41.5, 100]),
    "cosp_cth": np.array([0, 250, 750, 1250, 1750, 2250, 2750, 3500, 4500, 6000,
                          8000, 10000, 12000, 14500, 16000, 18000]),
}

def add_cosp_coordinates(ds):
    """
    Add default coordinate values to EAMxx COSP dimensions in dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with EAMxx COSP data

    Returns
    -------
    xr.Dataset
        Dataset with added coordinates
    """
    coords_added = False

    # Check each COSP dimension and add coordinates if missing
    for dim_name, coord_values in DEFAULT_COORDS.items():
        if dim_name in ds.dims:
            dim_size = ds.dims[dim_name]

            # Check if coordinate values exist and are not empty
            has_coords = dim_name in ds.coords and ds.coords[dim_name].size > 0

            if not has_coords:
                # Add default coordinates
                if dim_size == len(coord_values):
                    print(f"    Adding coordinates for dimension: {dim_name} (size={dim_size})")
                    ds = ds.assign_coords({dim_name: coord_values})
                    coords_added = True
                else:
                    print(f"    WARNING: Dimension '{dim_name}' has size {dim_size} "
                          f"but defaults have {len(coord_values)} values. Skipping.")

    return ds, coords_added


def process_file(input_path, output_path):
    """Process a single file."""
    print(f"Processing: {os.path.basename(input_path)}")

    # Open dataset
    ds = xr.open_dataset(input_path)

    # Add coordinates
    ds_modified, coords_added = add_cosp_coordinates(ds)

    if coords_added:
        print(f"    Saving to: {output_path}")
        ds_modified.to_netcdf(output_path)
        print(f"    Success!")
    else:
        print(f"    No COSP coordinates to add, copying file...")
        shutil.copy2(input_path, output_path)

    ds.close()


def main():
    """Main processing function."""
    print("=" * 70)
    print("Preprocessing EAMxx COSP Data")
    print("=" * 70)
    print(f"Input directory:  {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print()

    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Created output directory: {OUTPUT_DIR}")
    print()

    # Get all .nc files from input directory
    input_files = glob.glob(os.path.join(INPUT_DIR, "*.nc"))

    if not input_files:
        print(f"ERROR: No .nc files found in {INPUT_DIR}")
        return

    print(f"Found {len(input_files)} files to process")
    print()

    # Process each file
    for input_path in sorted(input_files):
        filename = os.path.basename(input_path)
        output_path = os.path.join(OUTPUT_DIR, filename)

        try:
            process_file(input_path, output_path)
        except Exception as e:
            print(f"    ERROR processing {filename}: {e}")

        print()

    print("=" * 70)
    print("Processing complete!")
    print(f"Output files saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
