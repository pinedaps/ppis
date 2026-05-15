#!/usr/bin/env bash
# analysis_faunus.sh — unified epsilon/temperature sweep with equilibration + production
# Auto-detects LUNARC ($HOME/pinedaps) vs local ($HOME/projects) environment.

set -euo pipefail

SCRIPT_START_TIME=$(date +%s)
SCRIPT_START_DATE=$(date)
echo "Script started at: $SCRIPT_START_DATE"

#######################################
# Environment auto-detection
#######################################

if [[ -d "$HOME/pinedaps" ]]; then
    PROJ_ROOT="$HOME/pinedaps"
else
    PROJ_ROOT="$HOME/projects"
fi
FAUNUS="$PROJ_ROOT/faunus-rs/target/release/faunus"
PDB2XYZ="$PROJ_ROOT/ppis/pdb2xyz/__init__AH_Hakan_Lambda_faunus_MB.py"
PLOT_SCRIPT="$PROJ_ROOT/ppis/faunus_simulations/plot_scripts/plot_dat.py"

#######################################
# Help
#######################################

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Sweep mode — provide exactly one of:
  Epsilon sweep (--T required):
    --epsilons <v1,v2,...>              e.g. 0.5,0.8368,1.0
    --epsmin <v> --epsmax <v> --epsstep <v>
    --T <value>                         Temperature in K

  Temperature sweep (--epsilon required):
    --temps <v1,v2,...>                 e.g. 290,300,310
    --tmin <v> --tmax <v> --tstep <v>
    --epsilon <value>                   Reference epsilon at 293 K

Common:
  --pdb <path>       PDB file (required)
  --pH <value>       pH
  --saltcon <value>  Salt concentration mol/L
  --outdir <path>    Output directory (default: basename of pdb)
  -h, --help

Examples:
  $0 --pdb ../pdbs/4LZT.pdb --epsilons 0.5,0.8368,1.0 --T 298.15 --pH 7.1 --saltcon 0.115 --outdir 4LZT
  $0 --pdb ../pdbs/4LZT.pdb --tmin 280 --tmax 320 --tstep 10 --pH 7.1 --saltcon 0.115 --epsilon 0.8368 --outdir 4LZT
EOF
}

#######################################
# Parse arguments
#######################################

PDB="" PH="" SC="" OUTDIR=""
TC="" EC=""
EPSMIN="" EPSMAX="" EPSSTEP="" USER_EPSILONS=""
TMIN="" TMAX="" TSTEP="" USER_TEMPS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pdb)       PDB="$2";           shift 2 ;;
        --pH)        PH="$2";            shift 2 ;;
        --saltcon)   SC="$2";            shift 2 ;;
        --outdir)    OUTDIR="$2";        shift 2 ;;
        --T)         TC="$2";            shift 2 ;;
        --epsilon)   EC="$2";            shift 2 ;;
        --epsilons)  USER_EPSILONS="$2"; shift 2 ;;
        --epsmin)    EPSMIN="$2";        shift 2 ;;
        --epsmax)    EPSMAX="$2";        shift 2 ;;
        --epsstep)   EPSSTEP="$2";       shift 2 ;;
        --temps)     USER_TEMPS="$2";    shift 2 ;;
        --tmin)      TMIN="$2";          shift 2 ;;
        --tmax)      TMAX="$2";          shift 2 ;;
        --tstep)     TSTEP="$2";         shift 2 ;;
        -h|--help)   usage; exit 0 ;;
        *)           echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

#######################################
# Resolve paths
#######################################

FILE="${PDB##*/}"
OUTDIR="${OUTDIR:-$FILE}"
FILE="${FILE%.*}"

PDB=$(realpath "$PDB")
mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
echo "Output directory: $OUTDIR"

#######################################
# Build sweep arrays
#   SWEEP_TAGS  — topo filename suffix (eps0.8368 / T298.15)
#   SWEEP_VALS  — directory suffix     (0.8368    / 298.15)
#   SWEEP_T     — temperature for each point
#   SWEEP_EPS   — epsilon for each point
#######################################

SWEEP_TAGS=()
SWEEP_VALS=()
SWEEP_T=()
SWEEP_EPS=()

if [[ -n "${USER_EPSILONS:-}" || ( -n "${EPSMIN:-}" && -n "${EPSMAX:-}" && -n "${EPSSTEP:-}" ) ]]; then
    # --- Epsilon sweep ---
    [[ -n "${TC:-}" ]] || { echo "ERROR: --T is required for epsilon sweep." >&2; exit 1; }
    RAW=()
    if [[ -n "${USER_EPSILONS:-}" ]]; then
        IFS=',' read -ra RAW <<< "$USER_EPSILONS"
    else
        V="$EPSMIN"
        while (( $(echo "$V <= $EPSMAX" | bc -l) )); do
            RAW+=("$V"); V=$(echo "$V + $EPSSTEP" | bc -l)
        done
    fi
    for EPS in "${RAW[@]}"; do
        SWEEP_TAGS+=("eps${EPS}"); SWEEP_VALS+=("$EPS")
        SWEEP_T+=("$TC");          SWEEP_EPS+=("$EPS")
    done

elif [[ -n "${USER_TEMPS:-}" || ( -n "${TMIN:-}" && -n "${TMAX:-}" && -n "${TSTEP:-}" ) ]]; then
    # --- Temperature sweep ---
    [[ -n "${EC:-}" ]] || { echo "ERROR: --epsilon is required for temperature sweep." >&2; exit 1; }
    RAW=()
    if [[ -n "${USER_TEMPS:-}" ]]; then
        IFS=',' read -ra RAW <<< "$USER_TEMPS"
    else
        V="$TMIN"
        while (( $(echo "$V <= $TMAX" | bc -l) )); do
            RAW+=("$V"); V=$(echo "$V + $TSTEP" | bc -l)
        done
    fi
    for T in "${RAW[@]}"; do
        SWEEP_TAGS+=("T${T}"); SWEEP_VALS+=("$T")
        SWEEP_T+=("$T");       SWEEP_EPS+=("$EC")
    done

else
    echo "ERROR: provide an epsilon sweep (--epsilons or --epsmin/--epsmax/--epsstep)" >&2
    echo "       or a temperature sweep (--temps or --tmin/--tmax/--tstep)." >&2
    exit 1
fi

#########################################
# Shared output directories
#########################################

OUTDIR_BASE=$(basename "$OUTDIR")
EQUIL_OUTDIR="${OUTDIR}_equilibration"
TOPO_DIR="$OUTDIR/topologies"
PLOT_DIR="$OUTDIR/plots"

mkdir -p "$EQUIL_OUTDIR" "$TOPO_DIR" "$PLOT_DIR"

echo "pdb: $FILE   pH: $PH   [Salt]: $SC"
echo "Sweep (${#SWEEP_TAGS[@]} points): ${SWEEP_TAGS[*]}"
echo

#######################################
# Per-simulation function
#   $1 TAG  — topo filename suffix (eps0.8368 / T298.15)
#   $2 VAL  — directory suffix
#   $3 T    — temperature
#   $4 EPS  — epsilon
# Globals: FILE PDB PH SC OUTDIR OUTDIR_BASE EQUIL_OUTDIR TOPO_DIR FAUNUS PDB2XYZ
#######################################

run_simulation() {
    local TAG="$1" VAL="$2" T="$3" EPS="$4"
    local TOPO="topology_${FILE}_${TAG}.yaml"
    local TOPO_EQUIL="topology_${FILE}_${TAG}_equil.yaml"

    # ---- EQUILIBRATION ----
    local EQUIL_SUBDIR="${EQUIL_OUTDIR}/${OUTDIR_BASE}_${VAL}"
    mkdir -p "$EQUIL_SUBDIR/results/"{dat,yaml,traj}
    cd "$EQUIL_SUBDIR"

    echo "  [${TAG}] Generating topology"
    python3 "$PDB2XYZ" -i "$PDB" -o "${FILE}.xyz" -t "$TOPO" \
        --T "$T" --pH "$PH" --saltcon "$SC" --epsilon "$EPS"

    sed 's/repeat: 50000/repeat: 5000/' "$TOPO" > "$TOPO_EQUIL"

    echo "  [${TAG}] Equilibration (5000 repeats)"
    "$FAUNUS" run --input "$TOPO_EQUIL" -s state.yaml

    cp "$TOPO"     "$TOPO_DIR/"
    mv output.yaml results/yaml/
    mv *.gz        results/dat/
    mv traj*       results/traj/

    # ---- PRODUCTION ----
    local PROD_SUBDIR="${OUTDIR}/${OUTDIR_BASE}_${VAL}"
    mkdir -p "$PROD_SUBDIR/results/"{dat,yaml,traj}

    cp "$TOPO"       "$PROD_SUBDIR/"
    cp "${FILE}.xyz" "$PROD_SUBDIR/"
    cp state.yaml    "$PROD_SUBDIR/"
    cd "$PROD_SUBDIR"

    echo "  [${TAG}] Production (50000 repeats)"
    "$FAUNUS" run --input "$TOPO" -s state.yaml

    mv "$TOPO" "$TOPO_DIR/${TOPO%.yaml}_prod.yaml"
    mv output.yaml results/yaml/
    mv *.gz        results/dat/
    mv traj*       results/traj/
}

#######################################
# Launch all simulations in parallel
#######################################

PIDS=()
for i in "${!SWEEP_TAGS[@]}"; do
    ( run_simulation "${SWEEP_TAGS[$i]}" "${SWEEP_VALS[$i]}" "${SWEEP_T[$i]}" "${SWEEP_EPS[$i]}" ) &
    PIDS+=($!)
done

echo "Waiting for ${#PIDS[@]} parallel simulations to finish..."
FAILED=0
for PID in "${PIDS[@]}"; do
    wait "$PID" || { echo "ERROR: simulation job $PID failed" >&2; FAILED=1; }
done
[[ $FAILED -eq 0 ]] || exit 1

#######################################
# Plot
#######################################

echo "Generating plots..."
python3 "$PLOT_SCRIPT" --dat_dir "$OUTDIR" --plots_dir "$PLOT_DIR"
echo "Faunus simulations complete."
echo

#######################################
# Execution time
#######################################

SCRIPT_END_TIME=$(date +%s)
ELAPSED=$((SCRIPT_END_TIME - SCRIPT_START_TIME))

echo "========================================"
echo "Execution Time Summary"
echo "========================================"
echo "Script started:  $SCRIPT_START_DATE"
echo "Script ended:    $(date)"
printf "Total elapsed time: %dh %dm %ds (%d seconds)\n" \
    $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60)) "$ELAPSED"
echo "======================================="
