import os
import re
import urllib.request

from tests.integration.config import TEST_DATA_DIR, TEST_IMAGES_DIR, TEST_ROOT_PATH


# https://stackoverflow.com/questions/49113616/how-to-download-file-using-python
def retrieve_file(url, file_path):
    dir_path = os.path.join(*os.path.split(file_path)[:-1])
    # https://stackoverflow.com/questions/12517451/automatically-creating-directories-with-file-output
    if dir_path:
        os.makedirs(dir_path, exist_ok=True)
    print("Downloading {} to {}".format(url, file_path))
    urllib.request.urlretrieve(url, file_path)
    return file_path


def download_files(url_prefix, url_suffix, directory_prefix=None):
    print(f"Downloading {url_suffix}")
    print("url_prefix={}".format(url_prefix))
    print("url_suffix={}".format(url_suffix))
    print("(local) directory_prefix={}".format(directory_prefix))

    url = os.path.join(url_prefix, url_suffix)

    if directory_prefix:
        links_file_path = os.path.join(directory_prefix, url_suffix)
    else:
        links_file_path = url_suffix

    links_file_path = "{}.html".format(links_file_path)
    print(
        "Downloading files from {}; checking for links on {}".format(
            url, links_file_path
        )
    )
    html_path = retrieve_file(url, links_file_path)
    links: list[str] = []

    with open(html_path, "r") as html:
        for line in html:
            match_object = re.search(r'href=[\'"]?([^\'" >]+)', line)
            if match_object:
                link = match_object.group(1)
                # Ignore parent directory and sorting links
                if (
                    ("../" not in link)
                    and (not link.startswith("/"))
                    and ("?" not in link)
                ):
                    print("Found a link: {}".format(link))
                    links.append(link)

    if os.path.exists(links_file_path):
        os.remove(links_file_path)

    files = []
    directories = []

    for link in links:
        if link.endswith("/"):
            # List directories to download.
            directories.append(link)
        else:
            # List '.csv', '.mat', '.nc', and '.png' files to download.
            files.append(link)

    print("\n###Downloading files")

    if directory_prefix:
        new_directory_prefix = os.path.join(directory_prefix, url_suffix)
    else:
        new_directory_prefix = url_suffix
    for f in files:
        url_to_download = os.path.join(url, f)
        file_path = os.path.join(new_directory_prefix, f)
        retrieve_file(url_to_download, file_path)

    print("\n###Downloading directories")
    for d in directories:
        new_directory = d.rstrip("/")
        download_files(url, new_directory, directory_prefix=new_directory_prefix)


def download():
    download_files(
        "https://web.lcrc.anl.gov/public/e3sm/e3sm_diags_test_data/integration",
        TEST_DATA_DIR,
        directory_prefix=TEST_ROOT_PATH,
    )
    download_files(
        "https://web.lcrc.anl.gov/public/e3sm/e3sm_diags_test_data/integration/expected",
        TEST_IMAGES_DIR,
        directory_prefix=TEST_ROOT_PATH,
    )
    print(f"Downloaded {TEST_DATA_DIR} and {TEST_ROOT_PATH}")


if __name__ == "__main__":
    download()
