import os

from PIL import Image, ImageChops


def images_are_identical(img1_path, img2_path):
    # Open both images
    with Image.open(img1_path) as img1, Image.open(img2_path) as img2:
        # Check if dimensions are the same
        if img1.size != img2.size:
            return False
        
        # Check if pixel values are the same
        diff = ImageChops.difference(img1, img2)
        return not diff.getbbox()  # getbbox() returns None if images are identical

def compare_png_images(dir1, dir2):
    # List to hold mismatched files, if any
    mismatched_files = []

    # Walk through both directories simultaneously
    for root, _, files in os.walk(dir1):
        for file in files:
            if file.endswith('.png'):
                # Full path for file in dir1
                path1 = os.path.join(root, file)
                
                # Construct corresponding path for file in dir2
                relative_path = os.path.relpath(path1, dir1)
                path2 = os.path.join(dir2, relative_path)
                
                # Check if the file exists in the second directory
                if not os.path.exists(path2):
                    print(f"File missing in second directory: {relative_path}")
                    mismatched_files.append(relative_path)
                    continue
                
                # Compare images
                if not images_are_identical(path1, path2):
                    print(f"Images do not match: {relative_path}")
                    mismatched_files.append(relative_path)

    # Final output
    if not mismatched_files:
        print("All .png images are identical.")
    else:
        print("Some images do not match or are missing in the second directory.")

# Example usage
compare_png_images('/global/cfs/cdirs/e3sm/www/chengzhu/complete_run_11112024/e3sm_diags_extended', '/global/cfs/cdirs/e3sm/www/chengzhu/tutorial2024/e3sm_diags_extended')
