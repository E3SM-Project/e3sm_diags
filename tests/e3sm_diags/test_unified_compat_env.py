from __future__ import annotations

import json

from tests.integration.generate_unified_compat_env import (
    build_env_text,
    build_metadata,
    extract_env_dependencies,
    extract_run_requirements,
    merge_dependency_spec,
    normalize_dependency_spec,
    parse_context_version,
    select_dependency_specs,
)

RECIPE_TEXT = """
schema_version: 1

context:
  name: e3sm-unified
  version: "1.13.0rc1"

requirements:
  host:
    - python
  run:
    - python
    - e3sm_diags ==3.1.0
    - uxarray >=2024.12.0
    - xcdat ==0.10.1
    - if: mpi != "hpc"
      then:
        - nco ==5.3.6
    - dask ==2025.9.1
    - esmf ==8.9.0 ${{ mpi_prefix }}_*
    - esmpy ==8.9.0
    - matplotlib-base ==3.10.6
    - numpy >=2.0.0
    - xarray ==2025.9.0
    - xesmf ==0.8.8
    - if: mpi == "hpc"
      then:
        - pandas
"""

CI_ENV_TEXT = """
name: e3sm_diags_ci
channels:
  - conda-forge
dependencies:
  - python >=3.11,<3.14
  - pip
  - setuptools
  - cartopy >=0.17.0,<0.25.0
  - cartopy_offlinedata
  - cf-units
  - dask !=2024.12.0,!=2024.12.1
  - esmf >=8.8.0=nompi_*
  - matplotlib-base >=3.8.2
  - numpy >=2.0.0,<3.0.0
  - pytest
  - pytest-cov
  - uxarray >=2023.3.0
  - xarray >=2024.3.0
  - xcdat >=0.11.1,<1.0.0
  - xesmf >=0.8.7
"""


class TestGenerateUnifiedCompatEnv:
    def test_parses_context_version(self):
        result = parse_context_version(RECIPE_TEXT)

        assert result == "1.13.0rc1"

    def test_extracts_only_top_level_run_requirements(self):
        result = extract_run_requirements(RECIPE_TEXT)

        assert "python" in result
        assert "uxarray >=2024.12.0" in result
        assert 'if: mpi != "hpc"' not in result
        assert "nco ==5.3.6" not in result
        assert "pandas" not in result

    def test_extracts_base_env_dependencies(self):
        result = extract_env_dependencies(CI_ENV_TEXT)

        assert "python >=3.11,<3.14" in result
        assert "cf-units" in result
        assert "esmf >=8.8.0=nompi_*" in result
        assert "xesmf >=0.8.7" in result

    def test_normalizes_dependency_spec_to_solver_safe_form(self):
        assert normalize_dependency_spec("python") == "python"
        assert normalize_dependency_spec("numpy >=2.0.0") == "numpy >=2.0.0"
        assert (
            normalize_dependency_spec("esmf ==8.9.0 ${{ mpi_prefix }}_*")
            == "esmf ==8.9.0"
        )
        assert normalize_dependency_spec("xarray ==2025.9.0  # comment") == (
            "xarray ==2025.9.0"
        )

    def test_merges_base_build_suffix_into_override_spec(self):
        assert (
            merge_dependency_spec("esmf >=8.8.0=nompi_*", "esmf ==8.9.0")
            == "esmf ==8.9.0=nompi_*"
        )
        assert (
            merge_dependency_spec("xcdat >=0.11.1,<1.0.0", "xcdat ==0.10.1")
            == "xcdat ==0.10.1"
        )

    def test_selects_target_dependency_subset(self):
        run_requirements = extract_run_requirements(RECIPE_TEXT)

        result = select_dependency_specs(
            run_requirements,
            (
                "python",
                "uxarray",
                "xcdat",
                "dask",
                "esmf",
                "esmpy",
                "matplotlib-base",
                "numpy",
                "xarray",
                "xesmf",
                "pandas",
            ),
        )

        assert result == [
            "dask ==2025.9.1",
            "esmf ==8.9.0",
            "esmpy ==8.9.0",
            "matplotlib-base ==3.10.6",
            "numpy >=2.0.0",
            "python",
            "uxarray >=2024.12.0",
            "xarray ==2025.9.0",
            "xcdat ==0.10.1",
            "xesmf ==0.8.8",
        ]

    def test_builds_env_text_with_ci_support_packages(self):
        base_dependency_specs = extract_env_dependencies(CI_ENV_TEXT)
        env_text = build_env_text(
            base_dependency_specs=base_dependency_specs,
            dependency_specs=[
                "python",
                "esmf ==8.9.0",
                "matplotlib-base ==3.10.6",
                "esmpy ==8.9.0",
                "xarray ==2025.9.0",
            ],
            python_version="3.13",
            env_name="compat",
        )

        assert "name: compat" in env_text
        assert "  - python =3.13" in env_text
        assert "  - pip" in env_text
        assert "  - setuptools" in env_text
        assert "  - pytest" in env_text
        assert "  - pytest-cov" in env_text
        assert "  - cf-units" in env_text
        assert "  - cartopy_offlinedata" in env_text
        assert "  - esmf ==8.9.0=nompi_*" in env_text
        assert "  - esmpy ==8.9.0" in env_text
        assert "  - matplotlib-base ==3.10.6" in env_text
        assert "  - xarray ==2025.9.0" in env_text

    def test_builds_metadata_with_env_hash(self):
        metadata = build_metadata(
            recipe_url="https://example.com/raw.yaml",
            display_url="https://example.com/blob.yaml",
            feedstock_ref="main",
            recipe_version="1.13.0rc1",
            dependency_specs=["xarray ==2025.9.0"],
            env_text="name: compat\n",
        )

        assert metadata["recipe_version"] == "1.13.0rc1"
        assert metadata["feedstock_ref"] == "main"
        assert metadata["dependency_specs"] == ["xarray ==2025.9.0"]
        assert isinstance(metadata["env_hash"], str)
        assert len(metadata["env_hash"]) == 64

        serialized = json.dumps(metadata, sort_keys=True)
        assert "blob.yaml" in serialized
