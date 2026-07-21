#!/usr/bin/env bash
# Create the Python 3.13.12-3.13.14 and 3.14.0-3.14.6 test environments with
# xarray=2026.7.0. Run from the repository root:
#
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/1_create_python_release_envs_xr2026070.sh
#
# Existing environments are skipped. To remove and recreate them:
#
#   bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/1_create_python_release_envs_xr2026070.sh --force

set -euo pipefail

readonly XARRAY_VERSION="2026.7.0"
readonly PYTHON_VERSIONS=(
  "3.13.12"
  "3.13.13"
  "3.13.14"
  "3.14.0"
  "3.14.1"
  "3.14.2"
  "3.14.3"
  "3.14.4"
  "3.14.5"
  "3.14.6"
)
readonly ENV_NAMES=(
  "ed_1048_xr_2026070_py31312"
  "ed_1048_xr_2026070_py31313"
  "ed_1048_xr_2026070_py31314"
  "ed_1048_xr_2026070_py3140"
  "ed_1048_xr_2026070_py3141"
  "ed_1048_xr_2026070_py3142"
  "ed_1048_xr_2026070_py3143"
  "ed_1048_xr_2026070_py3144"
  "ed_1048_xr_2026070_py3145"
  "ed_1048_xr_2026070_py3146"
)

FORCE=0
SCRIPT_DIR=""
REPO_ROOT=""
ENV_YAML=""
TEMP_DIR=""


usage() {
  cat <<'EOF'
Usage: 1_create_python_release_envs_xr2026070.sh [--force]

Create the Python patch-release comparison environments for the E3SM Diags
#1048 regression test. Every environment uses xarray=2026.7.0:
  - ed_1048_xr_2026070_py31312 (python=3.13.12)
  - ed_1048_xr_2026070_py31313 (python=3.13.13)
  - ed_1048_xr_2026070_py31314 (python=3.13.14)
  - ed_1048_xr_2026070_py3140  (python=3.14.0)
  - ed_1048_xr_2026070_py3141  (python=3.14.1)
  - ed_1048_xr_2026070_py3142  (python=3.14.2)
  - ed_1048_xr_2026070_py3143  (python=3.14.3)
  - ed_1048_xr_2026070_py3144  (python=3.14.4)
  - ed_1048_xr_2026070_py3145  (python=3.14.5)
  - ed_1048_xr_2026070_py3146  (python=3.14.6)

Each environment is created from conda-env/dev_latest.yml with exact Python
and Xarray pins. The current local e3sm_diags checkout is then installed with
`python -m pip install`.

Options:
  --force       Remove and recreate target environments that already exist.
                Without --force, existing environments are skipped.
  -h, --help    Show this help text.
EOF
}


log() {
  printf '[create_python_release_envs_xr2026070] %s\n' "$*"
}


die() {
  printf '[create_python_release_envs_xr2026070] Error: %s\n' "$*" >&2
  exit 1
}


cleanup() {
  if [[ -n "${TEMP_DIR}" && -d "${TEMP_DIR}" ]]; then
    rm -rf "${TEMP_DIR}"
  fi
}


trap cleanup EXIT


require_cmd() {
  local cmd=$1
  command -v "${cmd}" >/dev/null 2>&1 || die "Required command not found: ${cmd}"
}


parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --force)
        FORCE=1
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "Unknown argument: $1"
        ;;
    esac
  done
}


init_paths() {
  SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
  REPO_ROOT=$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)
  ENV_YAML="${REPO_ROOT}/auxiliary_tools/debug/1048-py314-stall-cont/conda-env/dev_latest.yml"
  [[ -f "${ENV_YAML}" ]] || die "Missing environment YAML: ${ENV_YAML}"

  TEMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/e3sm-diags-python-releases.XXXXXX")
}


env_exists() {
  local env_name=$1
  conda env list | awk 'NR > 2 {gsub(/\*/, "", $1); if ($1 != "") print $1}' | grep -Fxq "${env_name}"
}


prepare_existing_env() {
  local env_name=$1

  if ! env_exists "${env_name}"; then
    return
  fi

  if [[ ${FORCE} -eq 1 ]]; then
    log "Removing existing environment: ${env_name}"
    conda env remove -y -n "${env_name}"
  else
    log "Skipping existing environment: ${env_name}"
  fi
}


render_env_yaml() {
  local python_version=$1
  local rendered_yaml=$2

  sed \
    -e "s/^\\([[:space:]]*-[[:space:]]*python\\).*/\\1=${python_version}/" \
    -e "s/^\\([[:space:]]*-[[:space:]]*xarray\\).*/\\1=${XARRAY_VERSION}/" \
    "${ENV_YAML}" > "${rendered_yaml}"

  grep -Fqx "  - python=${python_version}" "${rendered_yaml}" \
    || die "Failed to pin python=${python_version} in ${rendered_yaml}"
  grep -Fqx "  - xarray=${XARRAY_VERSION}" "${rendered_yaml}" \
    || die "Failed to pin xarray=${XARRAY_VERSION} in ${rendered_yaml}"
}


create_env() {
  local env_name=$1
  local python_version=$2
  local rendered_yaml="${TEMP_DIR}/${env_name}.yml"

  render_env_yaml "${python_version}" "${rendered_yaml}"

  log "Creating ${env_name} with python=${python_version}, xarray=${XARRAY_VERSION}"
  conda env create -f "${rendered_yaml}" -n "${env_name}"

  log "Installing local e3sm_diags into ${env_name}"
  conda run -n "${env_name}" python -m pip install \
    --no-deps --no-build-isolation "${REPO_ROOT}"
}


print_summary() {
  local i

  printf '\n%-34s  %-10s  %-12s\n' "Environment" "Python" "Xarray"
  printf '%-34s  %-10s  %-12s\n' \
    "----------------------------------" "----------" "------------"

  for i in "${!ENV_NAMES[@]}"; do
    printf '%-34s  %-10s  %-12s\n' \
      "${ENV_NAMES[$i]}" "${PYTHON_VERSIONS[$i]}" "${XARRAY_VERSION}"
  done

  printf '\nRun the minimum QA matrix with:\n'
  printf '  bash %s/2_run_python_release_qa_lcrc_xr2026070.sh --repro-runs 3 --timeout 60m\n' \
    "${SCRIPT_DIR}"
}


main() {
  local i
  local env_name

  parse_args "$@"
  require_cmd conda
  require_cmd git
  require_cmd mktemp
  init_paths

  for i in "${!ENV_NAMES[@]}"; do
    env_name=${ENV_NAMES[$i]}
    prepare_existing_env "${env_name}"

    if env_exists "${env_name}"; then
      continue
    fi

    create_env "${env_name}" "${PYTHON_VERSIONS[$i]}"
  done

  print_summary
}


main "$@"
