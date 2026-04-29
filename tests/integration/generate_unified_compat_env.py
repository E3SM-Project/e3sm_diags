from __future__ import annotations

import argparse
import bz2
import hashlib
import json
import re
import urllib.request
from pathlib import Path
from typing import TypedDict

TARGET_DEPENDENCIES = (
    "cartopy",
    "dask",
    "esmf",
    "esmpy",
    "matplotlib-base",
    "numpy",
    "pandas",
    "python",
    "uxarray",
    "xarray",
    "xcdat",
    "xesmf",
    "xgcm",
)

# Temporary CI validation pin for the latest-released compat environment.
# This isolates the known Cartopy 0.25.0 rendering regression hypothesis
# without changing the main CI environment or the selected e3sm-unified release.
COMPAT_DEPENDENCY_OVERRIDES = {
    "cartopy": "cartopy =0.24.0",
}


class PackageRecord(TypedDict):
    name: str
    version: str
    build: str
    build_number: int
    subdir: str
    timestamp: int
    depends: list[str]
    filename: str


def normalize_dependency_spec(spec: str) -> str | None:
    spec = spec.split("#", 1)[0].strip()
    if not spec:
        return None

    parts = spec.split()
    if not parts:
        return None

    if len(parts) == 1:
        return parts[0]

    return " ".join(parts[:2])


def get_package_name(spec: str) -> str:
    return spec.split()[0]


def extract_env_dependencies(env_text: str) -> list[str]:
    in_dependencies = False
    dependencies_indent = 0
    dependencies: list[str] = []

    for line in env_text.splitlines():
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())

        if stripped == "dependencies:":
            in_dependencies = True
            dependencies_indent = indent
            continue

        if in_dependencies and stripped and indent <= dependencies_indent:
            break

        if not in_dependencies or not stripped.startswith("- "):
            continue

        dependency = normalize_dependency_spec(stripped[2:])
        if dependency is not None:
            dependencies.append(dependency)

    return dependencies


def fetch_repodata(repodata_url: str) -> dict[str, object]:
    with urllib.request.urlopen(repodata_url) as response:
        payload = response.read()

    if repodata_url.endswith(".bz2"):
        payload = bz2.decompress(payload)

    return json.loads(payload.decode("utf-8"))


def parse_version_key(version: str) -> tuple[tuple[int, ...], int, int]:
    match = re.fullmatch(r"(\d+(?:\.\d+)*)(?:rc(\d+))?", version)
    if match is None:
        raise ValueError(f"Unsupported version format: {version}")

    release = tuple(int(part) for part in match.group(1).split("."))
    rc_num = int(match.group(2)) if match.group(2) is not None else 0
    is_final = 1 if match.group(2) is None else 0

    return release, is_final, rc_num


def _coerce_package_record(
    filename: object,
    package: object,
    default_subdir: str | None,
) -> PackageRecord | None:
    if not isinstance(filename, str) or not isinstance(package, dict):
        return None

    name = package.get("name")
    version = package.get("version")
    build = package.get("build")
    subdir = package.get("subdir", default_subdir)
    build_number = package.get("build_number", 0)
    timestamp = package.get("timestamp", 0)
    depends = package.get("depends", [])

    if not isinstance(name, str) or not isinstance(version, str):
        return None
    if not isinstance(build, str) or not isinstance(subdir, str):
        return None
    if not isinstance(build_number, int) or not isinstance(timestamp, int):
        return None
    if not isinstance(depends, list) or not all(
        isinstance(item, str) for item in depends
    ):
        return None

    return {
        "name": name,
        "version": version,
        "build": build,
        "build_number": build_number,
        "subdir": subdir,
        "timestamp": timestamp,
        "depends": depends,
        "filename": filename,
    }


def select_latest_nompi_package(repodata: dict[str, object]) -> PackageRecord:
    candidates: list[PackageRecord] = []

    default_subdir: str | None = None
    info = repodata.get("info")
    if isinstance(info, dict):
        info_subdir = info.get("subdir")
        if isinstance(info_subdir, str):
            default_subdir = info_subdir

    for key in ("packages", "packages.conda"):
        package_map = repodata.get(key, {})
        if not isinstance(package_map, dict):
            continue

        for filename, package in package_map.items():
            package_record = _coerce_package_record(filename, package, default_subdir)
            if package_record is None:
                continue
            if package_record["name"] != "e3sm-unified":
                continue
            if "nompi" not in package_record["build"]:
                continue
            if package_record["subdir"] != "linux-64":
                continue

            candidates.append(package_record)

    if not candidates:
        raise ValueError("No linux-64 nompi e3sm-unified packages found in repodata")

    return max(
        candidates,
        key=lambda package: (
            parse_version_key(package["version"]),
            package["build_number"],
            package["timestamp"],
            package["filename"],
        ),
    )


def select_dependency_specs(
    package_depends: list[str], selected_dependencies: tuple[str, ...]
) -> list[str]:
    selected = set(selected_dependencies)
    result: list[str] = []

    for raw_spec in package_depends:
        spec = normalize_dependency_spec(raw_spec)
        if spec is None:
            continue

        package_name = get_package_name(spec)
        if package_name in selected:
            result.append(spec)

    return sorted(result)


def apply_dependency_overrides(
    dependency_specs: list[str], overrides: dict[str, str]
) -> list[str]:
    specs_by_package = {get_package_name(spec): spec for spec in dependency_specs}
    specs_by_package.update(overrides)

    return [specs_by_package[package_name] for package_name in sorted(specs_by_package)]


def _extract_build_suffix(spec: str) -> str:
    parts = spec.split(maxsplit=1)
    if len(parts) != 2:
        return ""

    token = parts[1]
    operators = (">=", "<=", "==", "!=", "~=", ">", "<")
    operator = next((item for item in operators if token.startswith(item)), "")
    remainder = token[len(operator) :]

    if "=" not in remainder:
        return ""

    version_part, build_part = remainder.split("=", 1)
    if not version_part or not build_part:
        return ""

    if "," in version_part:
        return ""

    return f"={build_part}"


def merge_dependency_spec(base_spec: str, override_spec: str) -> str:
    override_parts = override_spec.split(maxsplit=1)
    if len(override_parts) != 2:
        return override_spec

    if _extract_build_suffix(override_spec):
        return override_spec

    build_suffix = _extract_build_suffix(base_spec)
    if not build_suffix:
        return override_spec

    return f"{override_parts[0]} {override_parts[1]}{build_suffix}"


def extract_python_version(package: PackageRecord) -> str:
    match = re.search(r"py(\d{2,3})", package["build"])
    if match is not None:
        digits = match.group(1)
        return f"{digits[0]}.{digits[1:]}"

    for raw_spec in package["depends"]:
        spec = normalize_dependency_spec(raw_spec)
        if spec is None or get_package_name(spec) != "python":
            continue

        match = re.search(r">=([0-9]+\.[0-9]+)", spec)
        if match is None:
            break

        return match.group(1)

    raise ValueError("Could not determine python version from package metadata")


def build_env_text(
    base_dependency_specs: list[str],
    dependency_specs: list[str],
    python_version: str,
    env_name: str,
) -> str:
    overrides = {get_package_name(spec): spec for spec in dependency_specs}
    dependencies: list[str] = []
    seen: set[str] = set()

    for base_spec in base_dependency_specs:
        package_name = get_package_name(base_spec)

        if package_name == "python":
            dependencies.append(f"python ={python_version}")
            seen.add(package_name)
            continue

        if package_name in overrides:
            dependencies.append(
                merge_dependency_spec(base_spec, overrides[package_name])
            )
            seen.add(package_name)
            continue

        dependencies.append(base_spec)
        seen.add(package_name)

    for spec in dependency_specs:
        package_name = get_package_name(spec)
        if package_name == "python" or package_name in seen:
            continue

        dependencies.append(spec)
        seen.add(package_name)

    env_lines = [
        f"name: {env_name}",
        "channels:",
        "  - conda-forge",
        "dependencies:",
    ]
    env_lines.extend(f"  - {spec}" for spec in dependencies)

    return "\n".join(env_lines) + "\n"


def build_metadata(
    repodata_url: str,
    package: PackageRecord,
    dependency_specs: list[str],
    env_text: str,
    python_version: str,
) -> dict[str, object]:
    env_hash = hashlib.sha256(env_text.encode("utf-8")).hexdigest()
    return {
        "dependency_specs": dependency_specs,
        "env_hash": env_hash,
        "package_build": package["build"],
        "package_filename": package["filename"],
        "package_subdir": package["subdir"],
        "python_version": python_version,
        "repodata_url": repodata_url,
        "version": package["version"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repodata-url", required=True)
    parser.add_argument(
        "--base-env-file",
        type=Path,
        default=Path("conda-env/ci.yml"),
        required=False,
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument(
        "--env-name", default="e3sm_unified_latest_released_compat", required=False
    )
    args = parser.parse_args()

    repodata = fetch_repodata(args.repodata_url)
    package = select_latest_nompi_package(repodata)
    python_version = extract_python_version(package)
    base_env_text = args.base_env_file.read_text(encoding="utf-8")
    base_dependency_specs = extract_env_dependencies(base_env_text)
    dependency_specs = select_dependency_specs(package["depends"], TARGET_DEPENDENCIES)
    dependency_specs = apply_dependency_overrides(
        dependency_specs,
        COMPAT_DEPENDENCY_OVERRIDES,
    )
    env_text = build_env_text(
        base_dependency_specs=base_dependency_specs,
        dependency_specs=dependency_specs,
        python_version=python_version,
        env_name=args.env_name,
    )
    metadata = build_metadata(
        repodata_url=args.repodata_url,
        package=package,
        dependency_specs=dependency_specs,
        env_text=env_text,
        python_version=python_version,
    )

    args.output.write_text(env_text, encoding="utf-8")
    args.metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(f"env_hash={metadata['env_hash']}")
    print(f"package_version={metadata['version']}")
    print(f"python_version={metadata['python_version']}")


if __name__ == "__main__":
    main()
