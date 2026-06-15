#!/bin/bash
# Usage: ./scripts/builder.sh <input.json>
# Example: ./scripts/builder.sh examples/test_case.json
#
# This script sets up the system using the provided JSON parameter file,
# runs the python builder, and stages the simulation.

if [ "$1" == "--help" ] || [ "$1" == "-h" ] || [ "$#" -ne 1 ]; then
    echo "========================================================================="
    echo " Membrane Builder Pipeline "
    echo "========================================================================="
    echo "Usage: ./scripts/builder.sh <input.json>"
    echo ""
    echo "This script automates the entire system generation, including:"
    echo "  - Creating arbitrary polymer/additive systems based on config."
    echo "  - Adding ions, solvating, and cleaning core-penetrating waters."
    echo "  - Perfectly generating complex GROMACS index.ndx files."
    echo "  - Staging cluster scripts (job.sh) and all MDPs for equilibration."
    echo ""
    echo "Examples:"
    echo "  ./scripts/builder.sh examples/test_case.json"
    echo ""
    echo "To build a PURE polymer system without additives:"
    echo "  Set \"additive_type\": \"none\" and \"percent_label\": 0 in your JSON."
    echo "  Set \"polymer_type\" to your custom topology prefix (e.g. \"pb51peo27\")."
    echo "========================================================================="
    exit 1
fi

INPUT_JSON=$1

if [ ! -f "$INPUT_JSON" ]; then
    echo "Error: Parameter file $INPUT_JSON not found!"
    exit 1
fi

PERCENT=$(python3 -c "import json; print(json.load(open('$INPUT_JSON')).get('percent_label', 'unknown'))")
ADDITIVE=$(python3 -c "import json; print(json.load(open('$INPUT_JSON')).get('additive_type', 'unknown'))")

# Parse the JSON file to get percentage and polymer type
percent_label=$(grep '"percent_label"' "$INPUT_JSON" | sed 's/[^0-9]*//g')
polymer_type=$(grep '"polymer_type"' "$INPUT_JSON" | awk -F'"' '{print $4}')

if [ -z "$percent_label" ]; then
    echo "Error: percent_label not found in $INPUT_JSON"
    exit 1
fi

if [ -z "$polymer_type" ]; then
    polymer_type="pb22peo14"
fi

target_dir="output_systems/${percent_label}_${polymer_type}"

MAX_RETRIES=5
attempt=1
SUCCESS=0

while [ $attempt -le $MAX_RETRIES ]; do
    echo "=========================================================="
    echo " Building system (Attempt $attempt) with ${PERCENT}% ${ADDITIVE}"
    echo "=========================================================="

    # We don't need to mkdir because build_system.py already creates it
    python3 scripts/build_system.py "$INPUT_JSON"
    
    cd "$target_dir" || exit 1

    # Copy mdps and scripts from the github repo structure
    echo "-> Copying MDP templates and job scripts..."
    cp -r ../../mdp_templates/* .
    cp ../../scripts/delete_waters.py .

    # Update job names
    sed -i "s/#SBATCH -J job_name/#SBATCH -J ${PERCENT}${ADDITIVE}/g" job*.sh

    # Generate index.ndx
    echo "-> Generating index.ndx..."
    gmx_mpi make_ndx -f system_300.gro -o index.ndx << INDEX_EOF > /dev/null 2>&1
keep 0
r TIP3* | r NA* | r CL* | r SOL* | r ION*
name 1 SOL_ION
! "SOL_ION"
name 2 MEMB
q
INDEX_EOF

    echo "-> Running Automatic System Minimization..."
    echo "-> 1) Running Steepest Descent to relax gross clashes..."
    gmx_mpi grompp -f step6.0_minimization.mdp -c system_300.gro -r system_300.gro -p system_300.top -n index.ndx -o step6.0.tpr -maxwarn 2 > grompp_min.log 2>&1
    gmx_mpi mdrun -deffnm step6.0 -v > mdrun_min.log 2>&1

    echo "-> 2) Running Conjugate Gradient to achieve negative energy..."
    gmx_mpi grompp -f step6.0b_minimization_cg.mdp -c step6.0.gro -r system_300.gro -p system_300.top -n index.ndx -o step6.0b.tpr -maxwarn 2 > grompp_cg.log 2>&1
    gmx_mpi mdrun -deffnm step6.0b -v > mdrun_cg.log 2>&1

    echo "-> Minimization complete! Overwriting system_300.gro with negative-energy coordinates..."
    mv step6.0b.gro system_300.gro

    # Perform final validation check on the Potential Energy
    FINAL_ENERGY=$(grep "Potential Energy" mdrun_cg.log | tail -n 1 | awk '{print $4}')
    echo "=========================================================="
    if [[ "$FINAL_ENERGY" == "-"* ]]; then
        echo " [SUCCESS] System is fully minimized!"
        echo " Final Potential Energy: $FINAL_ENERGY kJ/mol (Strictly Negative)"
        SUCCESS=1
        break
    else
        echo " [WARNING] Final Potential Energy: $FINAL_ENERGY kJ/mol"
        echo " Energy is NOT negative! Rebuilding the system with a new random seed..."
        cd ../../
        rm -rf "$target_dir"
        attempt=$((attempt + 1))
    fi
done

if [ $SUCCESS -eq 0 ]; then
    echo "=========================================================="
    echo " [FATAL ERROR] Could not achieve negative energy after $MAX_RETRIES attempts!"
    echo "=========================================================="
    exit 1
fi

echo "=========================================================="
echo " The relaxed coordinate file is: system_300.gro"
echo " "
echo " To run the equilibrations on the cluster:"
echo " cd ${target_dir}/"
echo " sbatch job.sh"
echo "=========================================================="
