#!/usr/bin/env bash

# ----------------------------------------
#   eps_analysis_faunus.sh
#   Epsilon workflow:
#     • Parse epsilon inputs
#     • Run topology generation
#     • Run faunus
#     • Calculate and print execution time
# ----------------------------------------

set -euo pipefail

# Start timing
SCRIPT_START_TIME=$(date +%s)
SCRIPT_START_DATE=$(date)

echo "Script started at: $SCRIPT_START_DATE"

#######################################
# Help message
#######################################

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Epsilon input options (choose one):
  --pdb <path>
        Pdb path
  --epsmin <value> --epsmax <value> --epsstep <value>
        Generate an epsilon array from min to max with given step
  --epsilons <comma-separated list>
        Provide an explicit list, e.g. --epsilons 0.5,0.8368,1.0
  --T <value>
        Temperature in Kelvin
  --pH <value>
        Provide a pH value
  --saltcon <value>
        Salt concentration in mol/L

Other options:
  --outdir <path>     Output directory (default: ./<input_pdb>)
  -h, --help          Show this help message

Example:
  $0 --pdb ../pdbs/XXXX --epsmin 0.5 --epsmax 1.0 --epsstep 0.1 --T 298.15 --pH 7.1 --saltcon 0.115 --outdir XXXX
  $0 --pdb ../pdbs/XXXX --epsilons 0.5,0.8368,1.0 --T 298.15 --pH 7.1 --saltcon 0.115 --outdir XXXX
EOF
}

#######################################
# Parse arguments
#######################################

while [[ $# -gt 0 ]]; do
    case "$1" in
	--pdb)
            PDB="$2"
            shift 2
            ;;
        --epsmin)
            EPSMIN="$2"
            shift 2
            ;;
        --epsmax)
            EPSMAX="$2"
            shift 2
            ;;
        --epsstep)
            EPSSTEP="$2"
            shift 2
            ;;
        --epsilons)
            USER_EPSILONS="$2"
            shift 2
            ;;
        --T)
            TC="$2"
            shift 2
            ;;
        --pH)
            PH="$2"
            shift 2
            ;;
	--saltcon)
            SC="$2"
            shift 2
            ;;
	--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

#######################################
# Default values
#######################################

FILE="${PDB##*/}"
OUTDIR="${OUTDIR:-$FILE}"
FILE="${FILE%.*}"

# Resolve to absolute paths so subshells and cd's don't break relative refs
PDB=$(realpath "$PDB")
mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")

echo "The output directory is $OUTDIR"

#######################################
# Build epsilon array
#######################################

EPS_ARRAY=()
if [[ -n "${USER_EPSILONS:-}" ]]; then
    IFS=',' read -ra EPS_ARRAY <<< "$USER_EPSILONS"
elif [[ -n "${EPSMIN:-}" && -n "${EPSMAX:-}" && -n "${EPSSTEP:-}" ]]; then
    EPS="$EPSMIN"
    while (( $(echo "$EPS <= $EPSMAX" | bc -l) )); do
        EPS_ARRAY+=("$EPS")
        EPS=$(echo "$EPS + $EPSSTEP" | bc -l)
    done
else
    echo "ERROR: You must provide either --epsilons or --epsmin/--epsmax/--epsstep." >&2
    exit 1
fi

#########################################
# Create output directory and filenames
#########################################

OUTDIR_BASE=$(basename "$OUTDIR")
TOPO_DIR="$OUTDIR/topologies"
PLOT_DIR="$OUTDIR/plots"

mkdir -p "$TOPO_DIR" "$PLOT_DIR"

echo "pdb: $FILE"
echo "Epsilons: ${EPS_ARRAY[*]}"
echo "T: $TC"
echo "pH: $PH"
echo "[Salt]: $SC"
echo

#######################################
# Steps 1+2: Topology + Faunus (parallel)
#######################################

PIDS=()
for EPS in "${EPS_ARRAY[@]}"; do
    (
        EPS_SUBDIR="${OUTDIR}/${OUTDIR_BASE}_${EPS}"
        EPS_DAT_DIR="$EPS_SUBDIR/results/dat"
        EPS_YAML_DIR="$EPS_SUBDIR/results/yaml"
        EPS_TRAJ_DIR="$EPS_SUBDIR/results/traj"
        mkdir -p "$EPS_DAT_DIR" "$EPS_YAML_DIR" "$EPS_TRAJ_DIR"
        cd "$EPS_SUBDIR"

        TOPO_OUT="topology_${FILE}_eps${EPS}.yaml"
        echo "  Running topology for pdb = $FILE at epsilon = $EPS → $TOPO_OUT"
        python3 $HOME/projects/ppis/pdb2xyz/__init__AH_Hakan_Lambda_faunus_MB.py \
             -i "$PDB" \
             -o "${FILE}.xyz" \
             -t "$TOPO_OUT" \
             --T "$TC"  \
             --pH "$PH" \
             --saltcon "$SC"  \
             --epsilon "$EPS"
        echo "  Topology generation complete for epsilon = $EPS"

        echo "  Faunus simulation for epsilon = $EPS"
        $HOME/projects/faunus-rs/target/release/faunus run --input "$TOPO_OUT"

        mv "$TOPO_OUT"   "$TOPO_DIR/"
        mv output.yaml   "$EPS_YAML_DIR/"
        mv *.gz          "$EPS_DAT_DIR/"
        mv traj*         "$EPS_TRAJ_DIR/"
    ) &
    PIDS+=($!)
done

echo "Waiting for ${#PIDS[@]} parallel simulations to finish..."
FAILED=0
for PID in "${PIDS[@]}"; do
    wait "$PID" || { echo "ERROR: simulation job $PID failed" >&2; FAILED=1; }
done
[[ $FAILED -eq 0 ]] || exit 1

#######################################
# Step 3: Plot results
#######################################

echo "Generating plots..."
python3 $HOME/projects/ppis/faunus_simulations/plot_scripts/plot_dat.py --dat_dir "$OUTDIR" --plots_dir "$PLOT_DIR"

echo "Faunus simulations complete."
echo

###########################################
# Step 4: Calculate and log execution time
###########################################

SCRIPT_END_TIME=$(date +%s)
SCRIPT_END_DATE=$(date)
ELAPSED=$((SCRIPT_END_TIME - SCRIPT_START_TIME))

HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo "========================================"
echo "Execution Time Summary"
echo "========================================"
echo "Script started:  $SCRIPT_START_DATE"
echo "Script ended:    $SCRIPT_END_DATE"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
echo "======================================="
echo
