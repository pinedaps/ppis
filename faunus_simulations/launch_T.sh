#!/usr/bin/env bash

# ----------------------------------------
#   launch_T.sh
#   Counts temperature simulations and submits
#   lunarc_submission_T.sh with the right
#   number of CPUs.
# ----------------------------------------

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]
Same options as T_analysis_faunus_lunarc.sh.

Example:
  $0 --pdb ../pdbs/XXXX --temps 290,300,310 --pH 7.1 --saltcon 0.115 --epsilon 0.8368 --outdir XXXX
  $0 --pdb ../pdbs/XXXX --tmin 280 --tmax 320 --tstep 10 --pH 7.1 --saltcon 0.115 --epsilon 0.8368 --outdir XXXX
EOF
}

# ---- Parse just the temperature-related args to count N ---------------------

TMIN="" TMAX="" TSTEP="" USER_TEMPS=""

_args=("$@")
i=0
while [[ $i -lt ${#_args[@]} ]]; do
    case "${_args[$i]}" in
        --temps)  USER_TEMPS="${_args[$((i+1))]}"; ((i+=2)) ;;
        --tmin)   TMIN="${_args[$((i+1))]}";       ((i+=2)) ;;
        --tmax)   TMAX="${_args[$((i+1))]}";       ((i+=2)) ;;
        --tstep)  TSTEP="${_args[$((i+1))]}";      ((i+=2)) ;;
        -h|--help) usage; exit 0 ;;
        *)        ((i+=1)) ;;
    esac
done

# ---- Count simulations -------------------------------------------------------

if [[ -n "$USER_TEMPS" ]]; then
    N=$(echo "$USER_TEMPS" | tr ',' '\n' | wc -l)
elif [[ -n "$TMIN" && -n "$TMAX" && -n "$TSTEP" ]]; then
    N=$(python3 -c "import math; print(int(math.floor(($TMAX - $TMIN) / $TSTEP)) + 1)")
else
    echo "ERROR: provide --temps or --tmin/--tmax/--tstep" >&2
    usage; exit 1
fi

echo "Submitting $N parallel temperature simulation(s) with --cpus-per-task=$N"
sbatch --cpus-per-task="$N" \
    "$(dirname "$0")/lunarc_submission_T.sh" "$@"
