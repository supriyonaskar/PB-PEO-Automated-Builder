import MDAnalysis as mda
import numpy as np
import sys
import shutil

def delete_waters(input_gro, input_top, output_gro, output_top):
    print(f"Loading {input_gro} to check for penetrating waters...")
    u = mda.Universe(input_gro)
    
    # 1. Identify Membrane and Core Z-region
    membrane = u.select_atoms("resname mol_01")
    if len(membrane) == 0:
        membrane = u.select_atoms("not (resname SOL TIP3 NA CL HOH WAT ION CHL1)")
        
    z_coords = membrane.positions[:, 2]
    midplane_z = np.mean(z_coords)
    
    # Define hydrophobic core as the middle 50% of the total membrane thickness
    z_span = np.max(z_coords) - np.min(z_coords)
    core_min = midplane_z - (z_span * 0.25)
    core_max = midplane_z + (z_span * 0.25)
    
    print(f"Membrane midplane: {midplane_z/10:.2f} nm. Core region: {core_min/10:.2f} to {core_max/10:.2f} nm.")
    
    # 2. Select Waters
    waters = u.select_atoms("resname SOL TIP3 HOH WAT")
    water_oxygens = waters.select_atoms("name OW OH2 O")
    
    # 3. Find penetrating waters
    inside_oxygens = water_oxygens.select_atoms(f"prop z >= {core_min} and prop z <= {core_max}")
    num_deleted = len(inside_oxygens)
    
    if num_deleted == 0:
        print("No waters found inside the hydrophobic core. All good!")
        shutil.copy(input_gro, output_gro)
        shutil.copy(input_top, output_top)
    else:
        print(f"Found {num_deleted} waters in the core! Deleting them...")
        bad_residues = inside_oxygens.residues
        good_atoms = u.select_atoms("all and not (group bad)", bad=bad_residues.atoms)
        
        # Write output gro
        with mda.Writer(output_gro, len(good_atoms)) as W:
            W.write(good_atoms)
        print(f"Saved cleaned coordinates to {output_gro}.")
        
        # Update Topology
        with open(input_top, 'r') as f:
            lines = f.readlines()
            
        with open(output_top, 'w') as f:
            for line in lines:
                if line.strip().startswith("SOL") or line.strip().startswith("TIP3"):
                    parts = line.split()
                    if len(parts) >= 2:
                        old_sol = int(parts[1])
                        new_sol = old_sol - num_deleted
                        f.write(f"{parts[0]:<10} {new_sol}\n")
                        print(f"Updated topology {parts[0]} count from {old_sol} to {new_sol}")
                        continue
                f.write(line)
        print(f"Saved updated topology to {output_top}.")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python delete_waters.py <input.gro> <input.top> <output.gro> <output.top>")
        sys.exit(1)
    delete_waters(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
