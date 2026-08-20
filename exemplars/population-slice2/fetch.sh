#!/usr/bin/env bash
#
# fetch.sh -- download the WS-E slice 2 subset of OpenNeuro ds003745.
#
# Pinned to snapshot 2.1.1. Every file, byte count and md5 in manifest.csv was
# read from that snapshot; nothing here is inferred. Files land under data/
# with their BIDS relpaths preserved, so data/ is a valid (partial) BIDS tree
# with a derivatives/fmriprep/ subtree beside it.
#
#   ./fetch.sh --dry-run              # list what would be fetched, print total
#   ./fetch.sh                        # fetch everything in the manifest
#   ./fetch.sh --tier core            # only the task-trust RSA tier
#   ./fetch.sh --role events,confounds --role mask-mni
#   ./fetch.sh --subject sub-104 --subject sub-105
#   ./fetch.sh --verify               # re-check md5 of what is already on disk
#   ./fetch.sh --verify --repair      # ...and delete anything that fails, so a
#                                     #    following plain run re-fetches it
#
# Re-running is cheap: a file whose size already matches the manifest is
# skipped without a network request. Note that a plain run trusts the byte
# count: to catch a file that is corrupt *at the right size*, use --verify.

set -euo pipefail

DATASET="ds003745"
SNAPSHOT="2.1.1"

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${HERE}/manifest.csv"
DEST="${HERE}/data"

# Filters are kept as "|"-delimited strings rather than arrays: macOS still
# ships bash 3.2, where an empty array trips `set -u`.
DRY_RUN=0
VERIFY_ONLY=0
REPAIR=0
RETRIES=5
TIERS=""
ROLES=""
SUBJECTS=""

die() { printf 'fetch.sh: %s\n' "$*" >&2; exit 1; }

usage() {
  # Print the header comment block: everything from line 3 up to the first line
  # that is not a comment. Found by scanning rather than by hard-coded line
  # numbers, so editing the header cannot silently truncate or leak code.
  awk 'NR<3 {next} /^#/ {sub(/^#[[:space:]]?/,""); print; next} {exit}' "${BASH_SOURCE[0]}"
  exit 0
}

# Never leave a half-written .part behind on Ctrl-C.
# Never leave a half-written .part behind.
#
# Two details matter here. First, the explicit `return 0`: in bash the exit
# status of an EXIT trap becomes the script's exit status, so a falsy last
# command would turn every successful run into a failure. Second, INT/TERM get
# their own handler that *exits* -- a bare handler would run and then let the
# script resume, so Ctrl-C would abort only the file in flight and carry on
# downloading the remaining several GB.
CURRENT_PART=""
cleanup() {
  if [[ -n "$CURRENT_PART" ]]; then rm -f "$CURRENT_PART"; fi
  return 0
}
on_signal() {
  printf '\nfetch.sh: interrupted\n' >&2
  cleanup
  trap - EXIT
  exit "$1"
}
trap cleanup EXIT
trap 'on_signal 130' INT
trap 'on_signal 143' TERM

# --- md5 helper: macOS ships `md5`, Linux ships `md5sum` ---------------------
if command -v md5sum >/dev/null 2>&1; then
  md5_of() { md5sum "$1" | awk '{print $1}'; }
elif command -v md5 >/dev/null 2>&1; then
  md5_of() { md5 -q "$1"; }
else
  md5_of() { printf ''; }   # no md5 tool: size check only
fi

file_size() {
  # portable stat; -1 for "cannot stat", so callers compare unequal rather than
  # dying under `set -e` inside a command substitution
  if stat --version >/dev/null 2>&1; then
    stat -c%s "$1" 2>/dev/null || echo -1
  else
    stat -f%z "$1" 2>/dev/null || echo -1
  fi
}

# --- argument parsing -------------------------------------------------------
# add_filter <varname> <comma-separated values>  -> appends "|v1|v2" to var
add_filter() {
  local cur="$1" add="$2" flag="$3"
  [[ -n "$add" ]] || die "$flag needs a non-empty value"
  case "$add" in (*[!A-Za-z0-9,_-]*) die "$flag: invalid value '$add'" ;; esac
  printf '%s|%s' "$cur" "$(printf '%s' "$add" | tr ',' '|')"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1; shift ;;
    --verify)     VERIFY_ONLY=1; shift ;;
    --repair)     REPAIR=1; shift ;;
    --tier)       [[ $# -ge 2 ]] || die "--tier needs a value"
                  TIERS="$(add_filter "$TIERS" "$2" --tier)"; shift 2 ;;
    --role)       [[ $# -ge 2 ]] || die "--role needs a value"
                  ROLES="$(add_filter "$ROLES" "$2" --role)"; shift 2 ;;
    --subject)    [[ $# -ge 2 ]] || die "--subject needs a value"
                  SUBJECTS="$(add_filter "$SUBJECTS" "$2" --subject)"; shift 2 ;;
    --retries)    [[ $# -ge 2 ]] || die "--retries needs a value"
                  case "$2" in (''|*[!0-9]*) die "--retries must be a non-negative integer, got '$2'" ;; esac
                  RETRIES="$2"; shift 2 ;;
    -h|--help)    usage ;;
    *)            die "unknown argument '$1' (try --help)" ;;
  esac
done

if (( REPAIR )) && (( ! VERIFY_ONLY )); then
  die "--repair only makes sense together with --verify"
fi

[[ -f "$MANIFEST" ]] || die "manifest not found: $MANIFEST"
command -v curl >/dev/null 2>&1 || die "curl is required"

in_filter() {  # in_filter <needle> <"|a|b|c" filter> ; empty filter matches all
  local needle="$1" filter="$2"
  [[ -z "$filter" ]] && return 0
  case "${filter}|" in (*"|${needle}|"*) return 0 ;; esac
  return 1
}

human() {  # bytes -> human readable, integer arithmetic only
  local b="$1"
  if   (( b >= 1073741824 )); then printf '%d.%02d GiB' $(( b / 1073741824 )) $(( (b % 1073741824) * 100 / 1073741824 ))
  elif (( b >= 1048576 ));    then printf '%d.%02d MiB' $(( b / 1048576 ))    $(( (b % 1048576) * 100 / 1048576 ))
  elif (( b >= 1024 ));       then printf '%d.%02d KiB' $(( b / 1024 ))       $(( (b % 1024) * 100 / 1024 ))
  else printf '%d B' "$b"; fi
}

printf '%s snapshot %s -> %s\n' "$DATASET" "$SNAPSHOT" "$DEST"
[[ -n "$TIERS"    ]] && printf '  filter tier:    %s\n' "${TIERS#|}"
[[ -n "$ROLES"    ]] && printf '  filter role:    %s\n' "${ROLES#|}"
[[ -n "$SUBJECTS" ]] && printf '  filter subject: %s\n' "${SUBJECTS#|}"
printf '\n'

n_total=0; n_sel=0; n_have=0; n_get=0; n_fail=0; n_badmd5=0
bytes_sel=0; bytes_get=0; bytes_ok=0; n_sel_data=0

# `|| [[ -n "$relpath" ]]` keeps the final row even if the manifest happens to
# lack a trailing newline: without it `read` returns false on the last line and
# a whole file is silently dropped, with a successful exit status.
# `sha256` is read but unused -- it is empty in this manifest by design.
while IFS=, read -r relpath role tier source bytes md5 sha256 url s3_uri || [[ -n "${relpath:-}" ]]; do
  relpath="${relpath%$'\r'}"; s3_uri="${s3_uri%$'\r'}"   # tolerate CRLF manifests
  [[ "$relpath" == "relpath" ]] && continue      # header
  [[ -z "${relpath:-}" ]] && continue
  n_total=$(( n_total + 1 ))

  in_filter "$tier" "$TIERS" || continue
  in_filter "$role" "$ROLES" || continue
  if [[ -n "$SUBJECTS" ]]; then
    keep=0
    # Dataset- and derivative-level metadata is always kept: it is a few KB, and
    # derivatives/fmriprep/dataset_description.json is what makes the fetched
    # derivative tree readable as a BIDS derivative at all.
    if [[ "$role" == "dataset-metadata" || "$role" == "derivative-metadata" ]]; then
      keep=1
    else
      # Match whole path segments only, so --subject sub-1 does not match
      # sub-104/sub-105/... and quietly select the entire dataset.
      old_ifs="$IFS"; IFS='|'
      for s in $SUBJECTS; do
        [[ -z "$s" ]] && continue
        case "/${relpath}/" in (*"/${s}/"*) keep=1; break ;; esac
        case "${relpath##*/}" in ("${s}_"*) keep=1; break ;; esac
      done
      IFS="$old_ifs"
    fi
    (( keep )) || continue
  fi

  n_sel=$(( n_sel + 1 ))
  bytes_sel=$(( bytes_sel + bytes ))
  case "$role" in (dataset-metadata|derivative-metadata) ;; (*) n_sel_data=$(( n_sel_data + 1 )) ;; esac
  out="${DEST}/${relpath}"

  # already present and the right size?
  if [[ -f "$out" ]] && [[ "$(file_size "$out")" == "$bytes" ]]; then
    if (( VERIFY_ONLY )) && [[ -n "$md5" ]]; then
      got="$(md5_of "$out")"
      if [[ -n "$got" && "$got" != "$md5" ]]; then
        printf '  BAD MD5  %s\n     expected %s\n     got      %s\n' "$relpath" "$md5" "$got"
        n_badmd5=$(( n_badmd5 + 1 ))
        if (( REPAIR )); then rm -f "$out"; printf '     removed (re-run without --verify to re-fetch)\n'; fi
      fi
    fi
    n_have=$(( n_have + 1 ))
    continue
  fi

  n_get=$(( n_get + 1 ))
  bytes_get=$(( bytes_get + bytes ))

  if (( VERIFY_ONLY )); then
    if [[ -f "$out" ]]; then
      printf '  WRONG SIZE  %s (expected %s, on disk %s)\n' "$relpath" "$bytes" "$(file_size "$out")"
      if (( REPAIR )); then rm -f "$out"; printf '     removed (re-run without --verify to re-fetch)\n'; fi
    else
      printf '  MISSING  %s\n' "$relpath"
    fi
    continue
  fi
  if (( DRY_RUN )); then
    printf '  FETCH  %10s  %s\n' "$(human "$bytes")" "$relpath"
    continue
  fi

  mkdir -p "$(dirname "$out")"
  tmp="${out}.part"
  CURRENT_PART="$tmp"
  printf '  [%d] %s (%s)\n' "$n_get" "$relpath" "$(human "$bytes")"
  if ! curl -fL --silent --show-error \
            --retry "$RETRIES" --retry-delay 2 --retry-connrefused \
            --connect-timeout 30 --max-time 3600 \
            -o "$tmp" "$url"; then
    printf '     DOWNLOAD FAILED: %s\n' "$url" >&2
    rm -f "$tmp"; CURRENT_PART=""; n_fail=$(( n_fail + 1 )); continue
  fi

  got_size="$(file_size "$tmp")"
  if [[ "$got_size" != "$bytes" ]]; then
    printf '     SIZE MISMATCH: expected %s got %s -- discarding\n' "$bytes" "$got_size" >&2
    rm -f "$tmp"; CURRENT_PART=""; n_fail=$(( n_fail + 1 )); continue
  fi
  if [[ -n "$md5" ]]; then
    got="$(md5_of "$tmp")"
    if [[ -n "$got" && "$got" != "$md5" ]]; then
      printf '     MD5 MISMATCH: expected %s got %s -- discarding\n' "$md5" "$got" >&2
      rm -f "$tmp"; CURRENT_PART=""; n_badmd5=$(( n_badmd5 + 1 )); n_fail=$(( n_fail + 1 )); continue
    fi
  fi
  mv -f "$tmp" "$out"
  CURRENT_PART=""
  bytes_ok=$(( bytes_ok + bytes ))
done < "$MANIFEST"

printf '\n'
printf 'manifest rows        : %d\n' "$n_total"
printf 'selected             : %d  (%s)\n' "$n_sel" "$(human "$bytes_sel")"

# A filter that matches nothing is almost always a typo. Fail loudly rather
# than reporting a successful no-op.
if (( n_sel == 0 )); then
  printf '\n'
  die "no manifest rows matched the given filters -- check --tier/--role/--subject spelling"
fi
# Metadata is kept unconditionally under --subject, so a --subject that matches
# no real subject would otherwise look like a tiny successful fetch rather than
# a typo. Only guard when the user did not explicitly ask for a metadata role.
if [[ -n "$SUBJECTS" && -z "$ROLES" ]] && (( n_sel_data == 0 )); then
  printf '\n'
  die "--subject matched no data files (only dataset metadata) -- check the subject spelling"
fi

printf 'already on disk      : %d\n' "$n_have"
if (( VERIFY_ONLY )); then
  printf 'missing or wrong size: %d\n' "$n_get"
  printf 'md5 mismatches       : %d\n' "$n_badmd5"
  (( n_get == 0 && n_badmd5 == 0 )) || exit 1
elif (( DRY_RUN )); then
  printf 'would fetch          : %d  (%s)\n' "$n_get" "$(human "$bytes_get")"
  printf '\n(dry run -- nothing downloaded)\n'
else
  printf 'fetched this run     : %d  (%s)\n' "$(( n_get - n_fail ))" "$(human "$bytes_ok")"
  printf 'failures             : %d\n' "$n_fail"
  (( n_fail == 0 )) || exit 1
fi
