#!/usr/bin/env bash

# ----------------------------------------
#   launch_eps.sh
#   Counts epsilon simulations and submits
#   lunarc_submission_eps.sh with the right
#   number of CPUs.
# ----------------------------------------

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]
Same options as eps_analysis_faunus_lunarc.sh.

Example:
  $0 --pdb ../pdbs/XXXX --epsilons 0.5,0.8368,1.0 --T 298.15 --pH 7.1 --saltcon 0.115 --outdir XXXX
  $0 --pdb ../pdbs/XXXX --epsmin 0.5 --epsmax 1.0 --epsstep 0.1 --T 298.15 --pH 7.1 --saltcon 0.115 --outdir XXXX
EOF
}

# ---- Parse just the epsilon-related args to count N -------------------------

EPSMIN="" EPSMAX="" EPSSTEP="" USER_EPSILONS=""

_args=("$@")
i=0
while [[ $i -lt ${#_args[@]} ]]; do
    case "${_args[$i]}" in
        --epsilons)  USER_EPSILONS="${_args[$((i+1))]}"; ((i+=2)) ;;
        --epsmin)    EPSMIN="${_args[$((i+1))]}";        ((i+=2)) ;;
        --epsmax)    EPSMAX="${_args[$((i+1))]}";        ((i+=2)) ;;
        --epsstep)   EPSSTEP="${_args[$((i+1))]}";       ((i+=2)) ;;
        -h|--help)   usage; exit 0 ;;
        *)           ((i+=1)) ;;
    esac
done

# ---- Count simulations -------------------------------------------------------

if [[ -n "$USER_EPSILONS" ]]; then
    N=$(echo "$USER_EPSILONS" | tr ',' '\n' | wc -l)
elif [[ -n "$EPSMIN" && -n "$EPSMAX" && -n "$EPSSTEP" ]]; then
    N=$(python3 -c "import math; print(int(math.floor(($EPSMAX - $EPSMIN) / $EPSSTEP)) + 1)")
else
    echo "ERROR: provide --epsilons or --epsmin/--epsmax/--epsstep" >&2
    usage; exit 1
fi

echo "Submitting $N parallel epsilon simulation(s) with --cpus-per-task=$N"
sbatch --cpus-per-task="$N" \
    "$(dirname "$0")/lunarc_submission_eps.sh" "$@"
