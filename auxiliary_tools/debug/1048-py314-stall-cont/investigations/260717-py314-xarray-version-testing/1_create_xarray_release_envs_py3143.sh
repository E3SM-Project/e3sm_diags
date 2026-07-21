#!/usr/bin/env bash
# bash auxiliary_tools/debug/1048-py314-stall-cont/create_xarray_release_envs_py3143.sh

set -euo pipefail

readonly PYTHON_VERSION="3.14.3"
readonly ENV_NAMES=(
  "ed_1048_xr_2025120_py3143"
  "ed_1048_xr_2026010_py3143"
  "ed_1048_xr_latest_2026070_py3143"
)
readonly XARRAY_VERSIONS=(
  "2025.12.0"
  "2026.01.0"
  "2026.07.0"
)

FORCE=0
SCRIPT_DIR=""
REPO_ROOT=""
ENV_YAML=""
TEMP_DIR=""


usage() {
  cat <<'EOF'
Usage: create_xarray_release_envs_py3143.sh [--force]

Create the Python 3.14.3 conda-release comparison environments for the E3SM
Diags #1048 xarray regression workflow:
  - ed_1048_xr_2025120_py3143
  - ed_1048_xr_2026010_py3143
  - ed_1048_xr_latest_2026070_py3143

Each environment is created from dev_latest.yml with exact conda pins for
python=3.14.3 and the target xarray version, then the current local e3sm_diags
checkout is installed into it with `python -m pip install`.

Options:
  --force       Remove and recreate target envs if they already exist.
                Without --force, existing envs are skipped.
  -h, --help    Show this help text.
EOF
}


log() {
  printf '[create_xarray_release_envs_py3143] %s\n' "$*"
}


die() {
  printf '[create_xarray_release_envs_py3143] Error: %s\n' "$*" >&2
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

  if REPO_ROOT=$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel 2>/dev/null); then
    :
  else
    REPO_ROOT=$(cd -- "${SCRIPT_DIR}/../../.." && pwd)
  fi

  ENV_YAML="${SCRIPT_DIR}/dev_latest.yml"
  [[ -f "${ENV_YAML}" ]] || die "Missing environment YAML: ${ENV_YAML}"

  TEMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/e3sm-diags-xarray-py3143.XXXXXX")
}


env_exists() {
  local env_name=$1
  conda env list | awk 'NR > 2 {gsub(/\*/, "", $1); if ($1 != "") print $1}' | grep -Fxq "${env_name}"
}


remove_env_if_requested() {
  local env_name=$1

  if env_exists "${env_name}"; then
    if [[ ${FORCE} -eq 1 ]]; then
      log "Removing existing environment: ${env_name}"
      conda env remove -y -n "${env_name}"
    else
      log "Skipping existing environment: ${env_name}"
    fi
  fi
}


should_create_env() {
  local env_name=$1

  if env_exists "${env_name}"; then
    if [[ ${FORCE} -eq 1 ]]; then
      printf 'yes\n'
    else
      printf 'no\n'
    fi
  else
    printf 'yes\n'
  fi
}


render_env_yaml() {
  local xarray_version=$1
  local rendered_yaml=$2

  sed \
    -e "s/^\\([[:space:]]*-[[:space:]]*python\\).*/\\1=${PYTHON_VERSION}/" \
    -e "s/^\\([[:space:]]*-[[:space:]]*xarray\\).*/\\1=${xarray_version}/" \
    "${ENV_YAML}" > "${rendered_yaml}"

  grep -Fqx "  - python=${PYTHON_VERSION}" "${rendered_yaml}" \
    || die "Failed to pin python=${PYTHON_VERSION} in ${rendered_yaml}"
  grep -Fqx "  - xarray=${xarray_version}" "${rendered_yaml}" \
    || die "Failed to pin xarray=${xarray_version} in ${rendered_yaml}"
}


create_env() {
  local env_name=$1
  local xarray_version=$2
  local rendered_yaml="${TEMP_DIR}/${env_name}.yml"

  render_env_yaml "${xarray_version}" "${rendered_yaml}"

  log "Creating conda environment ${env_name} with python=${PYTHON_VERSION}, xarray=${xarray_version}"
  conda env create -f "${rendered_yaml}" -n "${env_name}"
}


install_local_e3sm_diags() {
  local env_name=$1

  log "Installing local e3sm_diags into ${env_name}"
  conda run -n "${env_name}" python -m pip install --no-deps --no-build-isolation "${REPO_ROOT}"
}


print_summary() {
  local i

  printf '\n'
  printf '%-34s  %-14s  %-20s\n' "Environment" "Python" "Xarray source"
  printf '%-34s  %-14s  %-20s\n' "----------------------------------" "--------------" "--------------------"

  for i in "${!ENV_NAMES[@]}"; do
    printf '%-34s  %-14s  %-20s\n' \
      "${ENV_NAMES[$i]}" \
      "${PYTHON_VERSION}" \
      "conda xarray=${XARRAY_VERSIONS[$i]}"
  done

  printf '\n'
  log "Ready-to-run examples:"
  printf '  %s\n' "bash auxiliary_tools/debug/1048-py314-stall-cont/run_xarray_release_qa_lcrc_py3143.sh"
  printf '  %s\n' "bash auxiliary_tools/debug/1048-py314-stall-cont/run_xarray_release_qa_lcrc_py3143.sh --repro-runs 5 --timeout 60m"
}


main() {
  local i
  local env_name
  local xarray_version

  parse_args "$@"
  require_cmd bash
  require_cmd git
  require_cmd conda
  require_cmd mktemp

  init_paths

  for i in "${!ENV_NAMES[@]}"; do
    remove_env_if_requested "${ENV_NAMES[$i]}"
  done

  for i in "${!ENV_NAMES[@]}"; do
    env_name=${ENV_NAMES[$i]}
    xarray_version=${XARRAY_VERSIONS[$i]}

    if [[ "$(should_create_env "${env_name}")" != "yes" ]]; then
      continue
    fi

    create_env "${env_name}" "${xarray_version}"
    install_local_e3sm_diags "${env_name}"
  done

  print_summary
}


main "$@"
