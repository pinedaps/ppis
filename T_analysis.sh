#!/usr/bin/env bash

# -----------------------------
#   T_analysis.sh
#   Temperature workflow:
#     • Parse T inputs
#     • Run topology generation
#     • Run duello scan
#     • Run plotting
# -----------------------------

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
        Generate a temperature array from min to max with given step.
  --temps <comma-separated list>
        Provide an explicit list, e.g. --temps 300,310,320
  --pH <value>
        Provide a pH value
  --epsilon <value>
        Provide the reference epsilon for the LJ potential at 293 K
  --ionic_strength <value>
        Provide the ionic strength at the provided pH in M
  --constant_T_dependence <value>
	Provide the constant to be used in the Lambda(T) equation
  --prot_radius <value>
        Provide the radius of the protein to calculate B2_HS in Å

Other options:
  --outdir <path>     Output directory (default: ./<input_pdb>)  
  -h, --help          Show this help message

Example:
  $0 --pdb pdbs/XXXX --tmin 280 --tmax 320 --tstep 10 --pH 7.1 --epsilon 0.8368 --ionic_strength 0.1 --constant_T_dependence 0.0252 --prot_radius 17 --outdir XXXX
  $0 --pdb pdbs/XXXX --temps 290,300,310 --pH 7.1 --epsilon 0.8368 --ionic_strength 0.1 --constant_T_dependence 0.0252 --prot_radius 17 --outdir XXXX
EOF
}

#######################################
# Parse arguments
#######################################

# To support both long options and getopts,
# we manually parse long options:
while [[ $# -gt 0 ]]; do
    case "$1" in
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
	--epsilon)
            EC="$2"
            shift 2
            ;;
	--ionic_strength)
            IS="$2"
            shift 2
            ;;
	--constant_T_dependence)
            CT="$2"
            shift 2
            ;;
	--prot_radius)
            PR="$2"
            shift 2
            ;;
        --pdb)
            PDB="$2"
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

mkdir -p "$OUTDIR"
XYZ_OUT="${OUTDIR}/${FILE}.xyz"
TOPO_DIR="$OUTDIR/topologies"
SCAN_DIR="$OUTDIR/scans"
PLOT_DIR="$OUTDIR/plots"

mkdir -p "$TOPO_DIR" "$SCAN_DIR" "$PLOT_DIR"

echo "pH: $PH"
echo "epsilon_c: $EC"
echo "Temperatures: ${T_ARRAY[*]}"
echo "Ionic strength: $IS"
echo "Coordinates files: $XYZ_OUT"
echo

#######################################
# Step 1: Generate topology files
#######################################

echo "=== Generating topology files ==="

for T in "${T_ARRAY[@]}"; do
    TOPO_OUT="${TOPO_DIR}/topology_${FILE}_T${T}.yaml"
    echo "  Running topology for pdb = $FILE at T = $T → $TOPO_OUT"
    #python3 pdb2xyz/__init__AH_Hakan_Lambda_duello_buried_aa.py \
    python3 pdb2xyz/__init__AH_Hakan_Epsilon.py \
        -i "$PDB" \
    	-o "$XYZ_OUT" \
    	-t "$TOPO_OUT" \
    	--pH "$PH" \
    	--T "$T"  \
    	--epsilon "$EC" \
	--constant_T_dependence "$CT" \
	--sidechains
done

echo "Topology generation complete."
echo

#######################################
# Step 2: Run duello scan
#######################################

echo "=== Running duello scan ==="

for T in "${T_ARRAY[@]}"; do
    TOPO_IN="${TOPO_DIR}/topology_${FILE}_T${T}.yaml"
    SCAN_OUT="${SCAN_DIR}/scan_T${T}.dat"
    echo "  duello scan for T = $T → $SCAN_OUT"
    
     duello scan --mol1 "$XYZ_OUT" \
                 --mol2 "$XYZ_OUT" \
                 --rmin 23 \
                 --rmax 60 \
                 --dr 0.5 \
	    	 --resolution 0.7 \
	  	 --cutoff 100  \
                 --top "$TOPO_IN"  \
		 --molarity "$IS"  \
		 --temperature "$T" \
		 --pmf "$SCAN_OUT" \
		 --backend reference  #\
                 #--grid "type=invr2,size=200,shift=false,energy_cap=50"
done

echo "Duello scans complete."
echo

#######################################
# Step 3: Plot results
#######################################

echo "=== Plotting results ==="

python3 plot_scripts/plot_potential.py -p "${SCAN_DIR}/" -pe "./experiments/${OUTDIR}.dat" -s $PR

echo "Plots generated in: $PLOT_DIR"
echo
echo "=== Done! ==="

#######################################
# Calculate and log execution time
#######################################

SCRIPT_END_TIME=$(date +%s)
SCRIPT_END_DATE=$(date)
ELAPSED=$((SCRIPT_END_TIME - SCRIPT_START_TIME))

HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo ""
echo "========================================"
echo "Execution Time Summary"
echo "========================================"
echo "Script started:  $SCRIPT_START_DATE"
echo "Script ended:    $SCRIPT_END_DATE"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
echo "========================================"
