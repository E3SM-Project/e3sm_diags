from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path

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
)
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


def fetch_recipe_text(recipe_url: str) -> str:
    with urllib.request.urlopen(recipe_url) as response:
        return response.read().decode("utf-8")


def parse_context_version(recipe_text: str) -> str | None:
    in_context = False
    context_indent = 0

    for line in recipe_text.splitlines():
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())

        if stripped == "context:":
            in_context = True
            context_indent = indent
            continue

        if in_context and stripped and indent <= context_indent:
            break

        if in_context:
            match = re.match(r'version:\s*"([^"]+)"', stripped)
            if match is not None:
                return match.group(1)

    return None


def extract_run_requirements(recipe_text: str) -> list[str]:
    lines = recipe_text.splitlines()
    in_requirements = False
    in_run = False
    requirements_indent = 0
    run_indent = 0
    item_indent: int | None = None
    requirements: list[str] = []

    for line in lines:
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())

        if stripped == "requirements:":
            in_requirements = True
            requirements_indent = indent
            continue

        if in_requirements and stripped and indent <= requirements_indent:
            break

        if in_requirements and stripped == "run:":
            in_run = True
            run_indent = indent
            continue

        if in_run and stripped and indent <= run_indent:
            break

        if not in_run or not stripped.startswith("- "):
            continue

        if item_indent is None and not stripped.startswith("- if:"):
            item_indent = indent

        if item_indent is None or indent != item_indent:
            continue

        if stripped.startswith("- if:"):
            continue

        requirements.append(stripped[2:].strip())

    return requirements


def select_dependency_specs(
    run_requirements: list[str], selected_dependencies: tuple[str, ...]
) -> list[str]:
    selected = set(selected_dependencies)
    result: list[str] = []

    for raw_spec in run_requirements:
        spec = normalize_dependency_spec(raw_spec)
        if spec is None:
            continue

        package_name = spec.split()[0]
        if package_name in selected:
            result.append(spec)

    return sorted(result)


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
            dependencies.append(merge_dependency_spec(base_spec, overrides[package_name]))
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
    recipe_url: str,
    display_url: str,
    feedstock_ref: str,
    recipe_version: str | None,
    dependency_specs: list[str],
    env_text: str,
) -> dict[str, object]:
    env_hash = hashlib.sha256(env_text.encode("utf-8")).hexdigest()
    return {
        "dependency_specs": dependency_specs,
        "display_url": display_url,
        "env_hash": env_hash,
        "feedstock_ref": feedstock_ref,
        "recipe_url": recipe_url,
        "recipe_version": recipe_version,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recipe-url", required=True)
    parser.add_argument("--display-url", required=True)
    parser.add_argument("--feedstock-ref", required=True)
    parser.add_argument("--python-version", required=True)
    parser.add_argument(
        "--base-env-file",
        type=Path,
        default=Path("conda-env/ci.yml"),
        required=False,
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument(
        "--env-name", default="e3sm_unified_latest_feedstock_compat", required=False
    )
    args = parser.parse_args()

    recipe_text = fetch_recipe_text(args.recipe_url)
    recipe_version = parse_context_version(recipe_text)
    run_requirements = extract_run_requirements(recipe_text)
    base_env_text = args.base_env_file.read_text(encoding="utf-8")
    base_dependency_specs = extract_env_dependencies(base_env_text)
    dependency_specs = select_dependency_specs(run_requirements, TARGET_DEPENDENCIES)
    env_text = build_env_text(
        base_dependency_specs=base_dependency_specs,
        dependency_specs=dependency_specs,
        python_version=args.python_version,
        env_name=args.env_name,
    )
    metadata = build_metadata(
        recipe_url=args.recipe_url,
        display_url=args.display_url,
        feedstock_ref=args.feedstock_ref,
        recipe_version=recipe_version,
        dependency_specs=dependency_specs,
        env_text=env_text,
    )

    args.output.write_text(env_text, encoding="utf-8")
    args.metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    print(f"env_hash={metadata['env_hash']}")
    print(f"recipe_version={recipe_version}")


if __name__ == "__main__":
    main()
