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

echo "=========================================================="
echo " Building system with ${PERCENT}% ${ADDITIVE} using ${INPUT_JSON}"
echo "=========================================================="

# Run the unified builder, passing the json config
echo "-> Generating Coordinates and Topology..."
python3 scripts/build_system.py "$INPUT_JSON"

# Navigate to the newly created directory (which the python script names after percent_label)
cd ${PERCENT} || exit 1

# Copy mdps and scripts from the github repo structure
echo "-> Copying MDP templates and job scripts..."
cp -r ../mdp_templates/* .
cp ../scripts/delete_waters.py .

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

echo "=========================================================="
echo " Success! The generated structure and topology are ready in ${PERCENT}/"
echo " To run the minimizations and equilibrations on the cluster:"
echo " cd ${PERCENT}/"
echo " sbatch job.sh"
echo "=========================================================="
