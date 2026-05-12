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

mkdir -p "$OUTDIR"
TOPO_DIR="$OUTDIR/topologies"
DAT_DIR="$OUTDIR/results/dat"
YAML_DIR="$OUTDIR/results/yaml"
TRAJ_DIR="$OUTDIR/results/traj"
PLOT_DIR="$OUTDIR/plots"
XYZ_OUT="${FILE}.xyz"

mkdir -p "$TOPO_DIR" "$DAT_DIR" "$YAML_DIR" "$TRAJ_DIR" "$PLOT_DIR"

echo "pdb: $FILE"
echo "Epsilons: ${EPS_ARRAY[*]}"
echo "T: $TC"
echo "pH: $PH"
echo "[Salt]: $SC"
echo

#######################################
# Step 1: Generate topology files
#######################################

for EPS in "${EPS_ARRAY[@]}"; do
    TOPO_OUT="topology_${FILE}_eps${EPS}.yaml"
    echo "  Running topology for pdb = $FILE at epsilon = $EPS → $TOPO_OUT"
    python3 $HOME/projects/ppis/pdb2xyz/__init__AH_Hakan_Lambda_faunus_MB.py \
         -i "$PDB" \
         -o "$XYZ_OUT" \
         -t "$TOPO_OUT" \
	 --T "$TC"  \
         --pH "$PH" \
         --saltcon "$SC"  \
         --epsilon "$EPS"

    echo "  Topology generation complete."

#######################################
# Step 2: Run Faunus
#######################################

    echo "  Faunus simulation for epsilon = $EPS "
    echo
    $HOME/projects/faunus-rs/target/release/faunus run --input "$TOPO_OUT"
    mv "$TOPO_OUT" "$TOPO_DIR"
    mv "output.yaml" "output_${EPS}.yaml"
done

mv *.gz       	"$DAT_DIR"
mv $XYZ_OUT   	"$TOPO_DIR"
mv *.yaml     	"$YAML_DIR"
mv traj*       "$TRAJ_DIR"

#######################################
# Step 3: Plot results
#######################################

echo "Generating plots..."
python3 $HOME/projects/ppis/faunus_simulations/plot_scripts/plot_dat.py --dat_dir "$DAT_DIR" --plots_dir "$PLOT_DIR"

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
