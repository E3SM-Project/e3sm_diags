import os
import filecmp


def get_png_files(directory):
    """Get a list of .png files in the given directory."""
    png_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".png"):
                png_files.append(os.path.join(root, file))
    return png_files


def compare_directories(dir1, dir2):
    """Compare .png files between two directories."""
    png_files_dir1 = get_png_files(dir1)
    png_files_dir2 = get_png_files(dir2)

    # Create a dictionary with relative paths as keys
    png_files_dict1 = {os.path.relpath(file, dir1): file for file in png_files_dir1}
    png_files_dict2 = {os.path.relpath(file, dir2): file for file in png_files_dir2}

    # Find common files
    common_files = set(png_files_dict1.keys()).intersection(set(png_files_dict2.keys()))

    # Compare common files
    for file in common_files:
        if not filecmp.cmp(png_files_dict1[file], png_files_dict2[file], shallow=False):
            print(f"Difference found in file: {file}")

    # Find files only in dir1
    only_in_dir1 = sorted(set(png_files_dict1.keys()) - set(png_files_dict2.keys()))
    if only_in_dir1:
        print(f"Files only in {dir1}:")
        for file in only_in_dir1:
            print(file)

    # Find files only in dir2
    only_in_dir2 = sorted(set(png_files_dict2.keys()) - set(png_files_dict1.keys()))
    if only_in_dir2:
        print(f"Files only in {dir2}:")
        for file in only_in_dir2:
            print(file)


if __name__ == "__main__":
    dir1 = "/global/cfs/cdirs/e3sm/www/e3sm_diags/complete_run/25-01-15-branch-907/arm_diags"
    dir2 = "/global/cfs/cdirs/e3sm/www/e3sm_diags/complete_run/v2.12.1v2/arm_diags"
    compare_directories(dir1, dir2)
