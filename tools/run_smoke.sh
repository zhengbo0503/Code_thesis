#!/usr/bin/env bash
# run_smoke.sh - run every chapter's smoke test, each in its own process.
#
# Usage:
#   tools/run_smoke.sh [MATLAB_APP]
# e.g.
#   tools/run_smoke.sh /Applications/MATLAB_R2026a.app
#
# Each MATLAB test runs in a separate process on purpose: a mismatched
# OpenBLAS/libomp segfaults MATLAB rather than raising an error, so
# isolation stops one bad MEX file from taking down the whole run.
#
# NOTE ON EXIT CODES: Advanpix currently segfaults MATLAB during shutdown
# (this happens on R2025b too, not just R2026a). All computation and file
# output completes first, so results are valid, but the process exit code
# is unreliable. This script therefore scores each test by parsing its
# printed "N/M passed" line, not by the exit status.

set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATLAB_APP="${1:-/Applications/MATLAB_R2026a.app}"
MATLAB="$MATLAB_APP/bin/matlab"
ADVANPIX="/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac"

if [ ! -x "$MATLAB" ]; then echo "No MATLAB at $MATLAB" >&2; exit 2; fi
echo "Using $($MATLAB -help 2>/dev/null | head -1 || echo "$MATLAB_APP")"

LOG="$(mktemp -d)"
declare -a NAMES RESULTS
fails=0

run_matlab () {         # name  dir  command
    local name="$1" dir="$2" cmd="$3"
    printf '\n########## %s ##########\n' "$name"
    ( cd "$ROOT/$dir" && timeout 3600 "$MATLAB" -nodisplay -nosplash \
        -batch "addpath('$ADVANPIX'); $cmd" ) 2>&1 | tee "$LOG/$name.txt" \
        | grep -Ev '^\s*$' | sed -n '/===/,/----/p'
    local line
    line="$(grep -Eo '[0-9]+/[0-9]+ passed, [0-9]+ failed' "$LOG/$name.txt" | tail -1)"
    if [ -z "$line" ]; then
        # scripts that report by assertion rather than a tally
        if grep -q 'validation passed' "$LOG/$name.txt"; then line="preflight passed"
        else line="NO RESULT (check $LOG/$name.txt)"; fails=$((fails+1)); fi
    else
        local nf; nf="$(sed -E 's/.*passed, ([0-9]+) failed/\1/' <<<"$line")"
        [ "$nf" != "0" ] && fails=$((fails+nf))
    fi
    NAMES+=("$name"); RESULTS+=("$line")
}

run_julia () {
    local name="$1"
    printf '\n########## %s ##########\n' "$name"
    ( cd "$ROOT/ch4_eigenvalue_eigenvector/julia" && \
      timeout 1800 julia --project=JacobiEigen smoke_test.jl ) 2>&1 \
      | tee "$LOG/$name.txt"
    local line
    line="$(grep -Eo '[0-9]+/[0-9]+ passed, [0-9]+ failed' "$LOG/$name.txt" | tail -1)"
    if [ -z "$line" ]; then line="NO RESULT"; fails=$((fails+1));
    else
        local nf; nf="$(sed -E 's/.*passed, ([0-9]+) failed/\1/' <<<"$line")"
        [ "$nf" != "0" ] && fails=$((fails+nf))
    fi
    NAMES+=("$name"); RESULTS+=("$line")
}

# MEX first: everything in chapter 5 depends on it.
run_matlab "ch5-mex"     "ch5_singular_value_singular_vector" "smoke_test_mex"
run_matlab "ch3"         "ch3_preconditioning"                "smoke_test"
run_matlab "ch4"         "ch4_eigenvalue_eigenvector"         "smoke_test"
run_matlab "ch5"         "ch5_singular_value_singular_vector" "smoke_test"
run_matlab "ch5-complex" "ch5_singular_value_singular_vector/complex" "validate_complex_stack"
run_julia  "ch4-julia"

printf '\n================ SUMMARY ================\n'
for i in "${!NAMES[@]}"; do printf '  %-14s %s\n' "${NAMES[$i]}" "${RESULTS[$i]}"; done
printf '  logs: %s\n' "$LOG"
printf '  total failures: %d\n' "$fails"
exit $(( fails > 0 ? 1 : 0 ))
