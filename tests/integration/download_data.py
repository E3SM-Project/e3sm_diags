from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

from tests.integration.config import TEST_DATA_DIR, TEST_IMAGES_DIR, TEST_ROOT_PATH

LOCAL_TEST_DATA_ROOT = Path("/e3sm_diags_downloaded_data")
CONTAINER_TEST_DATA_ROOT = Path("/e3sm_diags_downloaded_data")
CRANE_CLI = "crane"
TEST_DATA_IMAGE = (
    "ghcr.io/e3sm-project/containers-e3sm-diags-test-data:e3sm-diags-test-data-0.0.2"
)


def _copy_tree_contents(source_dir: Path, destination_dir: Path) -> None:
    if not source_dir.is_dir():
        raise FileNotFoundError(f"Integration test data source not found: {source_dir}")

    destination_dir.mkdir(parents=True, exist_ok=True)
    print(f"Copying {source_dir} to {destination_dir}")

    for entry in source_dir.iterdir():
        destination_path = destination_dir / entry.name
        if entry.is_dir():
            shutil.copytree(entry, destination_path, dirs_exist_ok=True)
        else:
            shutil.copy2(entry, destination_path)


def _copy_requested_directory(source_root: Path, directory_name: str) -> None:
    source_dir = source_root / Path(TEST_ROOT_PATH) / directory_name
    destination_dir = Path(TEST_ROOT_PATH) / directory_name
    _copy_tree_contents(source_dir, destination_dir)


def _remove_existing_path(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def _run_command(command: list[str]) -> None:
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Required CLI '{command[0]}' was not found. Install crane "
            "(for example, `conda install conda-forge::crane`)."
        ) from exc
    except subprocess.CalledProcessError as exc:
        error_output = (exc.stderr or exc.stdout).strip()
        raise RuntimeError(
            f"Command failed: {' '.join(command)}\n{error_output}"
        ) from exc


def _destination_path_from_relative_path(
    destination_root: Path,
    relative_path: PurePosixPath,
) -> Path:
    destination_path = destination_root.joinpath(*relative_path.parts)

    destination_root_resolved = destination_root.resolve()
    destination_path_resolved = destination_path.resolve(strict=False)

    if os.path.commonpath(
        [str(destination_root_resolved), str(destination_path_resolved)]
    ) != str(destination_root_resolved):
        raise RuntimeError(
            f"Refusing to extract outside {destination_root}: {relative_path}"
        )

    return destination_path


def _extract_requested_directories_from_tarball(
    tarball_path: Path,
    directory_names: list[str],
) -> None:
    destination_root = Path(TEST_ROOT_PATH)
    source_root = PurePosixPath(CONTAINER_TEST_DATA_ROOT.as_posix().lstrip("/"))
    requested_root = source_root / PurePosixPath(TEST_ROOT_PATH)
    extracted_directories: set[str] = set()

    for directory_name in directory_names:
        _remove_existing_path(destination_root / directory_name)

    with tarfile.open(tarball_path) as archive:
        for member in archive:
            member_path = PurePosixPath(member.name.lstrip("/"))

            try:
                relative_path = member_path.relative_to(requested_root)
            except ValueError:
                continue

            if not relative_path.parts or relative_path.parts[0] not in directory_names:
                continue

            destination_path = _destination_path_from_relative_path(
                destination_root,
                relative_path,
            )
            extracted_directories.add(relative_path.parts[0])

            if member.isdir():
                destination_path.mkdir(parents=True, exist_ok=True)
                continue

            if member.isfile():
                destination_path.parent.mkdir(parents=True, exist_ok=True)
                extracted_file = archive.extractfile(member)
                if extracted_file is None:
                    raise RuntimeError(
                        f"Unable to extract {member.name} from {tarball_path}"
                    )
                with extracted_file, open(destination_path, "wb") as output_file:
                    shutil.copyfileobj(extracted_file, output_file)
                os.chmod(destination_path, member.mode & 0o777)
                continue

            if member.issym() or member.islnk():
                destination_path.parent.mkdir(parents=True, exist_ok=True)
                _remove_existing_path(destination_path)
                os.symlink(member.linkname, destination_path)

    missing_directories = sorted(set(directory_names) - extracted_directories)
    if missing_directories:
        missing_text = ", ".join(missing_directories)
        raise RuntimeError(
            f"Could not find requested test-data directories in {tarball_path}: {missing_text}"
        )


def _copy_requested_directories_from_image(
    container_image: str,
    directory_names: list[str],
) -> None:
    directory_list = ", ".join(directory_names)
    with tempfile.NamedTemporaryFile(suffix=".tar") as temp_tarball:
        print(f"Exporting {container_image} with {CRANE_CLI}")
        _run_command([CRANE_CLI, "export", container_image, temp_tarball.name])
        print(f"Extracting {directory_list} from {container_image}")
        _extract_requested_directories_from_tarball(
            Path(temp_tarball.name),
            directory_names,
        )


def download(
    source_root: Path = LOCAL_TEST_DATA_ROOT,
    include_data: bool = True,
    include_images: bool = True,
    container_image: str = TEST_DATA_IMAGE,
) -> None:
    if not include_data and not include_images:
        raise ValueError("At least one of include_data or include_images must be True")

    copied_directories: list[str] = []

    use_local_source = source_root.exists()

    if not use_local_source:
        print(
            "Local integration data tree not found; falling back to the test-data "
            f"container image {container_image}."
        )
        directories_to_copy = []
        if include_data:
            directories_to_copy.append(TEST_DATA_DIR)
        if include_images:
            directories_to_copy.append(TEST_IMAGES_DIR)
        _copy_requested_directories_from_image(
            container_image,
            directories_to_copy,
        )

    if include_data:
        if use_local_source:
            _copy_requested_directory(source_root, TEST_DATA_DIR)
        copied_directories.append(TEST_DATA_DIR)

    if include_images:
        if use_local_source:
            _copy_requested_directory(source_root, TEST_IMAGES_DIR)
        copied_directories.append(TEST_IMAGES_DIR)

    copied = ", ".join(copied_directories)
    if use_local_source:
        print(f"Copied {copied} from {source_root}")
    else:
        print(f"Copied {copied} from {container_image}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy integration test assets from a local tree or export them from the test-data image with crane."
    )
    parser.add_argument(
        "--source-root",
        default=str(LOCAL_TEST_DATA_ROOT),
        help=(
            "Root directory containing tests/integration/integration_test_data "
            "and tests/integration/integration_test_images "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--image",
        default=TEST_DATA_IMAGE,
        help=(
            "OCI image containing /e3sm_diags_downloaded_data (default: %(default)s)."
        ),
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--data-only",
        action="store_true",
        help="Copy only tests/integration/integration_test_data.",
    )
    group.add_argument(
        "--images-only",
        action="store_true",
        help="Copy only tests/integration/integration_test_images.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    download(
        source_root=Path(args.source_root),
        include_data=not args.images_only,
        include_images=not args.data_only,
        container_image=args.image,
    )
