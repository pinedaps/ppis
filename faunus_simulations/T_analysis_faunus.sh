#!/usr/bin/env bash

# ----------------------------------------
#   T_analysis_faunus.sh
#   Temperature workflow:
#     • Parse T inputs
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

Temperature input options (choose one):
  --pdb <path>        
        Pdb path
  --tmin <value> --tmax <value> --tstep <value>
        Generate a temperature array from min to max with given step
  --temps <comma-separated list>
        Provide an explicit list, e.g. --temps 300,310,320
  --pH <value>
        Provide a pH value
  --saltcon <value>
        Salt concentration in mol/L  
  --epsilon <value>
        Provide the reference epsilon for the LJ potential at 293 K

Other options:
  --outdir <path>     Output directory (default: ./<input_pdb>)  
  -h, --help          Show this help message

Example:
  $0 --pdb ../pdbs/XXXX --tmin 280 --tmax 320 --tstep 10 --pH 7.1 --saltcon 0.115 --epsilon 0.8368 --outdir XXXX 
  $0 --pdb ../pdbs/XXXX --temps 290,300,310 --pH 7.1 --saltcon 0.115 --epsilon 0.8368 --outdir XXXX
EOF
}

#######################################
# Parse arguments
#######################################

# To support both long options and getopts,
# we manually parse long options:
while [[ $# -gt 0 ]]; do
    case "$1" in
	--pdb)
            PDB="$2"
            shift 2
            ;;
        --tmin)
            TMIN="$2"
            shift 2
            ;;
        --tmax)
            TMAX="$2"
            shift 2
            ;;
        --tstep)
            TSTEP="$2"
            shift 2
            ;;
        --temps)
            USER_TEMPS="$2"
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
	--epsilon)
            EC="$2"
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
# Default values and build faunus-rs
#######################################

FILE="${PDB##*/}"
OUTDIR="${OUTDIR:-$FILE}"
FILE="${FILE%.*}"

PDB=$(realpath "$PDB")
mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")

echo "The output directory is $OUTDIR"

#######################################
# Build temperature array
#######################################

T_ARRAY=()
if [[ -n "${USER_TEMPS:-}" ]]; then
    IFS=',' read -ra T_ARRAY <<< "$USER_TEMPS"
elif [[ -n "${TMIN:-}" && -n "${TMAX:-}" && -n "${TSTEP:-}" ]]; then
    T_ARRAY=()
    T="$TMIN"
    while (( $(echo "$T <= $TMAX" | bc -l) )); do
        T_ARRAY+=("$T")
        T=$(echo "$T + $TSTEP" | bc -l)
    done
else
    echo "ERROR: You must provide either --temps or --tmin/--tmax/--tstep." >&2
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
echo "Temperatures: ${T_ARRAY[*]}"
echo "pH: $PH"
echo "[Salt]: $SC"
echo "epsilon: $EC"
echo

#######################################
# Steps 1+2: Topology + Faunus (parallel)
#######################################

PIDS=()
for T in "${T_ARRAY[@]}"; do
    (
        T_SUBDIR="${OUTDIR}/${OUTDIR_BASE}_${T}"
        T_DAT_DIR="$T_SUBDIR/results/dat"
        T_YAML_DIR="$T_SUBDIR/results/yaml"
        T_TRAJ_DIR="$T_SUBDIR/results/traj"
        mkdir -p "$T_DAT_DIR" "$T_YAML_DIR" "$T_TRAJ_DIR"
        cd "$T_SUBDIR"

        TOPO_OUT="topology_${FILE}_T${T}.yaml"
        echo "  Running topology for pdb = $FILE at T = $T → $TOPO_OUT"
        python3 $HOME/projects/ppis/pdb2xyz/__init__AH_Hakan_Lambda_faunus_MB.py \
             -i "$PDB" \
             -o "${FILE}.xyz" \
             -t "$TOPO_OUT" \
             --T "$T"  \
             --pH "$PH" \
             --saltcon "$SC"  \
             --epsilon "$EC"
        echo "  Topology generation complete for T = $T"

        echo "  Faunus simulation for T = $T"
        $HOME/projects/faunus-rs/target/release/faunus run --input "$TOPO_OUT"

        mv "$TOPO_OUT"   "$TOPO_DIR/"
        mv output.yaml   "$T_YAML_DIR/"
        mv *.gz          "$T_DAT_DIR/"
        mv traj*         "$T_TRAJ_DIR/"
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
# Step 3: Calculate and log execution time
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
