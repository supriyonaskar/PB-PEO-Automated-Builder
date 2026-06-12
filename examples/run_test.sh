#!/bin/bash
# Small test script to verify the pipeline

echo "========================================="
echo " Running Automated Pipeline Test Case "
echo "========================================="

cd /home/naskars/Project_Landfester/github

# 1. Run the builder using the test JSON
./scripts/builder.sh examples/test_case.json

# 2. Navigate to the created folder (10/)
cd 10/

# 3. Modify minimization mdp to run a very short minimization (e.g. 50 steps) just to test it works
sed -i 's/nsteps                  = 5000/nsteps                  = 50/g' step6.0_minimization.mdp

# 4. Run grompp to verify topology matches coordinate file
echo "-> Testing grompp for minimization..."
gmx_mpi grompp -f step6.0_minimization.mdp -o min_test -c system_300.gro -p system_300.top -r system_300.gro -maxwarn 2

# 5. Run a short minimization without srun for local testing
echo "-> Running short 50-step minimization test..."
gmx_mpi mdrun -ntomp 4 -deffnm min_test -pin on

echo "========================================="
echo " Test Complete! Check for any errors above."
echo "========================================="
