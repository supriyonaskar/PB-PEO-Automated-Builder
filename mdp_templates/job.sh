#!/bin/bash -l
# Standard output and error:
#SBATCH -D ./
# Job name
#SBATCH -J job_name
#
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#
# --- default case: use a single GPU on a shared node ---
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=110000
#
#SBATCH --time=24:00:00

# Load compiler and MPI modules
module purge
module load gcc/15 openmpi/5.0 rocm/7.2 openmpi_gpu/5.0 gromacs/2026.1

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2025.06/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2025.06/etc/profile.d/conda.sh" ]; then
        . "/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2025.06/etc/profile.d/conda.sh"
    else
        export PATH="/mpcdf/soft/RHEL_9/packages/x86_64/python-waterboa/2025.06/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate /ptmp/naskars/software/openmm

echo "Starting Equilibration Protocol at $(date)"

# Minimization is now handled automatically by builder.sh!
# The output is perfectly relaxed and saved directly as system_300.gro.

# 6 Equilibration steps WITH restraints
echo "==> Running 6 Equilibration steps WITH restraints"
gmx_mpi grompp -f step6.1_equilibration.mdp -o step6.1 -c system_300.gro -r system_300.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm step6.1 -pin on

gmx_mpi grompp -f step6.2_equilibration.mdp -o step6.2 -c step6.1.gro -r system_300.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm step6.2 -pin on

python3 merge_gap.py step6.2.gro step6.2_merged.gro
mv step6.2_merged.gro step6.2.gro
sleep 20s

gmx_mpi grompp -f step6.3_equilibration.mdp -o step6.3 -c step6.2.gro -r step6.2.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm step6.3 -pin on

gmx_mpi grompp -f step6.4_equilibration.mdp -o step6.4 -c step6.3.gro -r step6.2.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm step6.4 -pin on

gmx_mpi grompp -f step6.5_equilibration.mdp -o step6.5 -c step6.4.gro -r step6.2.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm step6.5 -pin on

gmx_mpi grompp -f step6.6_equilibration.mdp -o step6.6 -c step6.5.gro -r step6.2.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm step6.6 -pin on

echo "==> Running 6 Equilibration steps WITHOUT restraints"
gmx_mpi grompp -f free_step6.1_equilibration.mdp -o free_step6.1 -c step6.6.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm free_step6.1 -pin on

gmx_mpi grompp -f free_step6.2_equilibration.mdp -o free_step6.2 -c free_step6.1.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm free_step6.2 -pin on

gmx_mpi grompp -f free_step6.3_equilibration.mdp -o free_step6.3 -c free_step6.2.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm free_step6.3 -pin on

gmx_mpi grompp -f free_step6.4_equilibration.mdp -o free_step6.4 -c free_step6.3.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm free_step6.4 -pin on

gmx_mpi grompp -f free_step6.5_equilibration.mdp -o free_step6.5 -c free_step6.4.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm free_step6.5 -pin on



gmx_mpi grompp -f free_step6.6_equilibration.mdp -o free_step6.6 -c free_step6.5.gro -p system_300.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm free_step6.6 -pin on

# -------------------------------------------------------------
# Custom Step: Delete penetrating waters and create new topology
# -------------------------------------------------------------
echo "==> Deleting penetrating waters from free_step6.6.gro..."
python3 delete_waters.py free_step6.6.gro system_300.top free_step6.6_clean.gro system_300_new.top

echo "==> Regenerating index.ndx for cleaned system..."
gmx_mpi make_ndx -f free_step6.6_clean.gro -o index.ndx << EOF > /dev/null 2>&1
keep 0
r TIP3* | r NA* | r CL* | r SOL* | r ION*
name 1 SOL_ION
! "SOL_ION"
name 2 MEMB
q
EOF

# 10 ns Production
echo "==> Waiting 30s for filesystem sync..."
sleep 30
echo "==> Running 10 ns Production"
gmx_mpi grompp -f production_10ns.mdp -o production_10ns -c free_step6.6_clean.gro -p system_300_new.top -n index.ndx -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm production_10ns -pin on
echo "Protocol Complete at $(date)"

gmx_mpi grompp -f prod.mdp -c production_10ns.gro -t production_10ns.cpt -p system_300_new.top -n index.ndx -o prod.tpr -maxwarn 20
srun gmx_mpi mdrun -ntomp 24 -deffnm prod -pin on

