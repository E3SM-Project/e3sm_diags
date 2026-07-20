#!/usr/bin/env bash
# bash auxiliary_tools/debug/1048-py314-stall-cont/create_xarray_bisect_envs.sh

set -euo pipefail

readonly XR_COMMIT_018AD08B="018ad08b12e8471b8bcc0135ce59b227f50da54b"
readonly XR_COMMIT_0A2D81C7="0a2d81c7a17aab867aed362b0882d34cb89e1311"

FORCE=0
XARRAY_DIR=""

SCRIPT_DIR=""
REPO_ROOT=""
ENV_YAML=""
WORKTREE_ROOT=""

SUMMARY_ROWS=()
ACTIVE_WORKTREES=()


usage() {
  cat <<'EOF'
Usage: create_xarray_bisect_envs.sh [--force] [--xarray-dir PATH]

Create four conda dev environments for the E3SM Diags #1048 xarray bisect:
  - ed_1048_xr_before_018ad08b
  - ed_1048_xr_after_018ad08b
  - ed_1048_xr_before_0a2d81c7
  - ed_1048_xr_after_0a2d81c7

Options:
  --force             Remove and recreate any of the four target envs if present.
  --xarray-dir PATH   Reuse or create the local xarray clone at PATH.
  -h, --help          Show this help text.

Environment:
  XARRAY_SRC_DIR      Default local xarray clone path when --xarray-dir is omitted.
EOF
}


log() {
  printf '[create_xarray_bisect_envs] %s\n' "$*"
}


die() {
  printf '[create_xarray_bisect_envs] Error: %s\n' "$*" >&2
  exit 1
}


cleanup() {
  local worktree

  for worktree in "${ACTIVE_WORKTREES[@]}"; do
    if [[ -d "${worktree}" ]]; then
      git -C "${XARRAY_DIR}" worktree remove --force "${worktree}" >/dev/null 2>&1 || true
    fi
  done

  if [[ -n "${WORKTREE_ROOT}" && -d "${WORKTREE_ROOT}" ]]; then
    rm -rf "${WORKTREE_ROOT}"
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
      --xarray-dir)
        [[ $# -ge 2 ]] || die "--xarray-dir requires a path argument"
        XARRAY_DIR=$2
        shift 2
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
  local default_xarray_dir

  SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

  if REPO_ROOT=$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel 2>/dev/null); then
    :
  else
    REPO_ROOT=$(cd -- "${SCRIPT_DIR}/../../.." && pwd)
  fi

  ENV_YAML="${SCRIPT_DIR}/dev_xr2026010.yml"
  [[ -f "${ENV_YAML}" ]] || die "Missing environment YAML: ${ENV_YAML}"

  default_xarray_dir="${XARRAY_SRC_DIR:-${TMPDIR:-/tmp}/e3sm-diags-xarray-bisect/xarray}"
  XARRAY_DIR=${XARRAY_DIR:-${default_xarray_dir}}
  XARRAY_DIR=$(mkdir -p "$(dirname "${XARRAY_DIR}")" && cd -- "$(dirname "${XARRAY_DIR}")" && pwd)/$(basename "${XARRAY_DIR}")

  WORKTREE_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/e3sm-diags-xarray-worktrees.XXXXXX")
}


ensure_xarray_clone() {
  if [[ -d "${XARRAY_DIR}/.git" ]]; then
    log "Reusing local xarray clone at ${XARRAY_DIR}"
  elif [[ -e "${XARRAY_DIR}" ]]; then
    die "Path exists but is not a git clone: ${XARRAY_DIR}"
  else
    log "Cloning xarray into ${XARRAY_DIR}"
    git clone https://github.com/pydata/xarray.git "${XARRAY_DIR}"
  fi

  git -C "${XARRAY_DIR}" rev-parse --is-inside-work-tree >/dev/null 2>&1 \
    || die "Not a valid git checkout: ${XARRAY_DIR}"

  log "Fetching xarray history and tags"
  git -C "${XARRAY_DIR}" fetch --all --tags --prune
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
      die "Environment already exists: ${env_name} (rerun with --force to recreate)"
    fi
  fi
}


resolve_commit() {
  local revspec=$1
  git -C "${XARRAY_DIR}" rev-parse "${revspec}^{commit}"
}


create_xarray_worktree() {
  local env_name=$1
  local resolved_sha=$2
  local worktree_path="${WORKTREE_ROOT}/${env_name}"

  git -C "${XARRAY_DIR}" worktree add --quiet --detach "${worktree_path}" "${resolved_sha}"
  ACTIVE_WORKTREES+=("${worktree_path}")
  printf '%s\n' "${worktree_path}"
}


remove_active_worktree() {
  local worktree_path=$1
  local updated=()
  local existing

  git -C "${XARRAY_DIR}" worktree remove --force "${worktree_path}" >/dev/null 2>&1 || true

  for existing in "${ACTIVE_WORKTREES[@]}"; do
    if [[ "${existing}" != "${worktree_path}" ]]; then
      updated+=("${existing}")
    fi
  done

  if [[ ${#updated[@]} -gt 0 ]]; then
    ACTIVE_WORKTREES=("${updated[@]}")
  else
    ACTIVE_WORKTREES=()
  fi
}


create_env() {
  local env_name=$1

  log "Creating conda environment ${env_name} from ${ENV_YAML}"
  conda env create -f "${ENV_YAML}" -n "${env_name}"
}


install_xarray_revision() {
  local env_name=$1
  local resolved_sha=$2
  local worktree_path

  worktree_path=$(create_xarray_worktree "${env_name}" "${resolved_sha}")

  log "Replacing xarray in ${env_name} with ${resolved_sha}"
  conda run -n "${env_name}" python -m pip uninstall -y xarray >/dev/null 2>&1 || true
  conda run -n "${env_name}" python -m pip install --no-deps --force-reinstall --no-build-isolation "${worktree_path}"

  remove_active_worktree "${worktree_path}"
}


install_local_e3sm_diags() {
  local env_name=$1

  log "Installing local e3sm_diags into ${env_name}"
  conda run -n "${env_name}" python -m pip install --no-deps --no-build-isolation "${REPO_ROOT}"
}


record_summary() {
  local env_name=$1
  local label=$2
  local resolved_sha=$3

  SUMMARY_ROWS+=("${env_name}|${label}|${resolved_sha}")
}


print_summary() {
  local row
  local env_name
  local label
  local resolved_sha

  printf '\n'
  printf '%-32s  %-22s  %s\n' "Environment" "Requested label" "Resolved xarray SHA"
  printf '%-32s  %-22s  %s\n' "--------------------------------" "----------------------" "----------------------------------------"

  for row in "${SUMMARY_ROWS[@]}"; do
    IFS='|' read -r env_name label resolved_sha <<< "${row}"
    printf '%-32s  %-22s  %s\n' "${env_name}" "${label}" "${resolved_sha}"
  done

  printf '\n'
  log "Ready-to-run examples:"
  printf '  %s\n' "conda run -n ed_1048_xr_before_018ad08b python auxiliary_tools/debug/1048-py314-stall-cont/xarray-min/xr_mvce.py"
  printf '  %s\n' "conda run -n ed_1048_xr_after_0a2d81c7 python auxiliary_tools/debug/1048-py314-stall-cont/parallel/min-scripts/qa.py"
}


main() {
  local env_names=(
    "ed_1048_xr_before_018ad08b"
    "ed_1048_xr_after_018ad08b"
    "ed_1048_xr_before_0a2d81c7"
    "ed_1048_xr_after_0a2d81c7"
  )
  local labels=(
    "before_018ad08b"
    "after_018ad08b"
    "before_0a2d81c7"
    "after_0a2d81c7"
  )
  local revspecs=(
    "${XR_COMMIT_018AD08B}^"
    "${XR_COMMIT_018AD08B}"
    "${XR_COMMIT_0A2D81C7}^"
    "${XR_COMMIT_0A2D81C7}"
  )
  local i
  local env_name
  local label
  local revspec
  local resolved_sha
  local resolved_shas=()

  parse_args "$@"
  require_cmd bash
  require_cmd git
  require_cmd conda
  require_cmd mktemp

  init_paths
  ensure_xarray_clone

  for revspec in "${revspecs[@]}"; do
    resolved_shas+=("$(resolve_commit "${revspec}")")
  done

  for env_name in "${env_names[@]}"; do
    remove_env_if_requested "${env_name}"
  done

  for i in "${!env_names[@]}"; do
    env_name=${env_names[$i]}
    label=${labels[$i]}
    resolved_sha=${resolved_shas[$i]}

    log "Preparing ${env_name} (${label} -> ${resolved_sha})"
    create_env "${env_name}"
    install_xarray_revision "${env_name}" "${resolved_sha}"
    install_local_e3sm_diags "${env_name}"
    record_summary "${env_name}" "${label}" "${resolved_sha}"
  done

  print_summary
}


main "$@"
