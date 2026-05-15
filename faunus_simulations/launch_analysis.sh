#!/usr/bin/env bash
# launch_analysis.sh — count simulations and submit lunarc_submission.sh via SLURM
# Accepts the same arguments as analysis_faunus.sh.

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]
Same options as analysis_faunus.sh. Submits to SLURM with --cpus-per-task=N.

Examples:
  $0 --pdb ../pdbs/4LZT.pdb --epsilons 0.5,0.8368,1.0 --T 298.15 --pH 7.1 --saltcon 0.115 --outdir 4LZT
  $0 --pdb ../pdbs/4LZT.pdb --tmin 280 --tmax 320 --tstep 10 --pH 7.1 --saltcon 0.115 --epsilon 0.8368 --outdir 4LZT
EOF
}

#######################################
# Parse just the sweep args to count N
#######################################

EPSMIN="" EPSMAX="" EPSSTEP="" USER_EPSILONS=""
TMIN="" TMAX="" TSTEP="" USER_TEMPS=""

_args=("$@")
i=0
while [[ $i -lt ${#_args[@]} ]]; do
    case "${_args[$i]}" in
        --epsilons) USER_EPSILONS="${_args[$((i+1))]}"; ((i+=2)) ;;
        --epsmin)   EPSMIN="${_args[$((i+1))]}";        ((i+=2)) ;;
        --epsmax)   EPSMAX="${_args[$((i+1))]}";        ((i+=2)) ;;
        --epsstep)  EPSSTEP="${_args[$((i+1))]}";       ((i+=2)) ;;
        --temps)    USER_TEMPS="${_args[$((i+1))]}";    ((i+=2)) ;;
        --tmin)     TMIN="${_args[$((i+1))]}";          ((i+=2)) ;;
        --tmax)     TMAX="${_args[$((i+1))]}";          ((i+=2)) ;;
        --tstep)    TSTEP="${_args[$((i+1))]}";         ((i+=2)) ;;
        -h|--help)  usage; exit 0 ;;
        *)          ((i+=1)) ;;
    esac
done

#######################################
# Count simulations
#######################################

if [[ -n "$USER_EPSILONS" ]]; then
    N=$(echo "$USER_EPSILONS" | tr ',' '\n' | wc -l)
elif [[ -n "$EPSMIN" && -n "$EPSMAX" && -n "$EPSSTEP" ]]; then
    N=$(python3 -c "import math; print(int(math.floor(($EPSMAX - $EPSMIN) / $EPSSTEP)) + 1)")
elif [[ -n "$USER_TEMPS" ]]; then
    N=$(echo "$USER_TEMPS" | tr ',' '\n' | wc -l)
elif [[ -n "$TMIN" && -n "$TMAX" && -n "$TSTEP" ]]; then
    N=$(python3 -c "import math; print(int(math.floor(($TMAX - $TMIN) / $TSTEP)) + 1)")
else
    echo "ERROR: provide epsilon or temperature sweep arguments." >&2
    usage; exit 1
fi

echo "Submitting $N parallel simulation(s) with --cpus-per-task=$N"
sbatch --cpus-per-task="$N" \
    "$(dirname "$0")/lunarc_submission.sh" "$@"
