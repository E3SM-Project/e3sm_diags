"""
This script compares the dependencies of two conda environments relevant to the
E3SM project, and isolates overlapping dependencies with different versions.

It exports the package lists from:
1. A named conda environment (e.g., "e3sm_diags_dev_985").
2. The latest E3SM Unified environment, loaded via a shell script.

The script extracts all package names and versions (including pip-installed packages)
from both environments, then prints the overlapping dependencies with different versions.

Requires `conda` and `pyyaml` to be installed and available in the environment.
"""

import subprocess
import tempfile
import os

def get_env_export(env_name=None, shell_source=None):
    """
    Returns the conda environment export as a dict.
    If shell_source is given, sources the shell script before running conda env export.
    If env_name is given, uses 'conda env export -n ENVNAME'.
    If neither is given, uses 'conda env export' for the current environment.
    """
    if shell_source:
        with tempfile.NamedTemporaryFile(delete=False, mode="w+") as tf:
            tf_name = tf.name
        bash_cmd = (
            f"bash -c 'source {shell_source} && conda env export > {tf_name}'"
        )
        ret = subprocess.run(bash_cmd, shell=True)
        if ret.returncode != 0:
            os.unlink(tf_name)
            raise RuntimeError(f"Failed to source {shell_source} and export conda env")
        with open(tf_name) as f:
            env_yaml = f.read()
        os.unlink(tf_name)
    else:
        cmd = ["conda", "env", "export"]
        if env_name:
            cmd += ["-n", env_name]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to export conda env: {result.stderr}")
        env_yaml = result.stdout

    try:
        import yaml
    except ImportError:
        raise ImportError("Please install pyyaml: pip install pyyaml")
    return yaml.safe_load(env_yaml)

def extract_deps_with_versions(env_dict):
    deps = {}
    for dep in env_dict.get("dependencies", []):
        if isinstance(dep, str):
            parts = dep.split("=")
            name = parts[0]
            version = parts[1] if len(parts) > 1 else ""
            deps[name] = version
        elif isinstance(dep, dict) and "pip" in dep:
            for pipdep in dep["pip"]:
                parts = pipdep.split("==")
                name = parts[0]
                version = parts[1] if len(parts) > 1 else ""
                deps[name] = version
    return deps

def main():
    env1 = get_env_export(env_name="e3sm_diags_dev_985")
    env2 = get_env_export(shell_source="/lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh")

    deps1 = extract_deps_with_versions(env1)
    deps2 = extract_deps_with_versions(env2)

    overlapping = set(deps1.keys()) & set(deps2.keys())
    differing_versions = []
    for dep in sorted(overlapping):
        v1 = deps1[dep]
        v2 = deps2[dep]
        if v1 != v2:
            differing_versions.append((dep, v1, v2))

    print("Overlapping dependencies with different versions:")
    for dep, v1, v2 in differing_versions:
        print(f"{dep}: e3sm_diags_dev_985={v1!r}, e3sm_unified={v2!r}")

if __name__ == "__main__":
    main()