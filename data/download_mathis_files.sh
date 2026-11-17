#!/usr/bin/env bash
#
# Download all data files of a given phase (default: phase_3) listed in
# data/mathis_data_info_final.csv from CC-IN2P3 into data/mathis_files_downloaded via scp.
#
# Usage:
#   ./data/download_mathis_files.sh [phase]            # e.g. phase_3 (default)
#   DRY_RUN=1 ./data/download_mathis_files.sh          # only print what would be copied
#   REMOTE_HOST=cca018.in2p3.fr ./data/download_mathis_files.sh
#
# Files already downloaded completely (local size == remote size) are skipped, so an
# interrupted run can simply be restarted. Incomplete files are re-downloaded.
# Set FORCE=1 to re-download everything.
# A single shared SSH connection is used, so you authenticate only once.

set -euo pipefail

PHASE="${1:-phase_3}"
REMOTE_HOST="${REMOTE_HOST:-cca.in2p3.fr}"
REMOTE_USER="${REMOTE_USER:-mpetro}"
# Project directory on CC; file_path entries in the CSV are relative to it.
REMOTE_PROJECT_DIR="${REMOTE_PROJECT_DIR:-~/sps_mpetro/Sensitivity_training}"
DRY_RUN="${DRY_RUN:-0}"
FORCE="${FORCE:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CSV="$PROJECT_DIR/data/mathis_data_info_final.csv"
OUT_DIR="$PROJECT_DIR/data/mathis_files_downloaded"

[[ -f "$CSV" ]] || { echo "CSV not found: $CSV" >&2; exit 1; }

REMOTE="$REMOTE_USER@$REMOTE_HOST"
CONTROL_PATH="$(mktemp -u /tmp/ssh_mathis_XXXXXX)"
SSH_OPTS=(-o ControlMaster=auto -o ControlPath="$CONTROL_PATH" -o ControlPersist=600
          -o ServerAliveInterval=30 -o ServerAliveCountMax=4 -o ConnectTimeout=30)
REMOTE_SIZES="$(mktemp /tmp/mathis_sizes_XXXXXX)"

cleanup() {
    ssh "${SSH_OPTS[@]}" -O exit "$REMOTE" &>/dev/null || true
    rm -f "$REMOTE_SIZES"
}
trap cleanup EXIT

# Map a CSV file_path to its remote path and to a local path relative to OUT_DIR.
remote_path() {
    case "$1" in
        /*) echo "$1" ;;
        *)  echo "$REMOTE_PROJECT_DIR/$1" ;;
    esac
}
local_path() {
    case "$1" in
        data/mathis_files/*) echo "${1#data/mathis_files/}" ;;
        data/*)              echo "${1#data/}" ;;
        /*)                  echo "${1#/}" ;;
        *)                   echo "$1" ;;
    esac
}

FILES=()
while IFS= read -r f; do FILES+=("$f"); done < <(
    awk -F, -v phase="$PHASE" 'NR > 1 && $1 == phase { print $4 }' "$CSV" | sort -u
)

[[ ${#FILES[@]} -gt 0 ]] || { echo "No files found for $PHASE in $CSV" >&2; exit 1; }
echo "Found ${#FILES[@]} files for $PHASE. Destination: $OUT_DIR"

echo "Opening SSH connection to $REMOTE ..."
ssh "${SSH_OPTS[@]}" -fN "$REMOTE"

# Fetch all remote file sizes in one call ("<size> <path>" per line; missing files are omitted).
echo "Fetching remote file sizes ..."
ssh "${SSH_OPTS[@]}" "$REMOTE" "cd $REMOTE_PROJECT_DIR && stat -L -c '%s %n' ${FILES[*]}" \
    > "$REMOTE_SIZES" 2>/dev/null || true

remote_size() { awk -v f="$1" '$2 == f { print $1; exit }' "$REMOTE_SIZES"; }
local_size()  { [[ -f "$1" ]] && wc -c < "$1" | tr -d ' ' || echo ""; }

n_ok=0; n_skip=0; n_todo=0; failed=()
for f in "${FILES[@]}"; do
    src="$(remote_path "$f")"
    dst="$OUT_DIR/$(local_path "$f")"

    rsize="$(remote_size "$f")"
    lsize="$(local_size "$dst")"
    rm -f "$dst.part"

    if [[ -n "$lsize" && "$FORCE" != "1" ]]; then
        if [[ -n "$rsize" && "$lsize" == "$rsize" ]]; then
            echo "[skip] $dst already downloaded ($lsize bytes)"
            n_skip=$((n_skip + 1))
            continue
        elif [[ -z "$rsize" && "$lsize" -gt 0 ]]; then
            echo "[skip] $dst exists (remote size unknown, cannot verify)"
            n_skip=$((n_skip + 1))
            continue
        fi
        echo "[redo] $dst incomplete ($lsize of ${rsize:-?} bytes)"
    fi

    echo "[get ] $REMOTE:$src -> $dst"
    if [[ "$DRY_RUN" == "1" ]]; then
        n_todo=$((n_todo + 1))
        continue
    fi

    mkdir -p "$(dirname "$dst")"
    if scp "${SSH_OPTS[@]}" -p "$REMOTE:$src" "$dst.part" \
        && [[ -z "$rsize" || "$(local_size "$dst.part")" == "$rsize" ]] \
        && mv "$dst.part" "$dst"; then
        n_ok=$((n_ok + 1))
    else
        rm -f "$dst.part"
        failed+=("$f")
    fi
done

echo
if [[ "$DRY_RUN" == "1" ]]; then
    echo "Dry run: $n_todo would be downloaded, $n_skip already present. Run without DRY_RUN=1 to download."
    exit 0
fi
echo "Done: $n_ok downloaded, $n_skip skipped, ${#failed[@]} failed."
if [[ ${#failed[@]} -gt 0 ]]; then
    printf '  failed: %s\n' "${failed[@]}" >&2
    exit 1
fi
