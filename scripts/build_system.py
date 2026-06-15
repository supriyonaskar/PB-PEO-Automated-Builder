import numpy as np
import sys
import random
import json
from scipy.spatial import cKDTree
import os

import sys

config_file = sys.argv[1] if len(sys.argv) > 1 else "config.json"
with open(config_file, "r") as f:
    config = json.load(f)


box_x = config.get("box_x", 11.39)
box_y = config.get("box_y", 11.39)
import math
chains_per_leaflet = config.get("chains_per_leaflet", 150)
z_padding = config.get("z_padding", 3.5)
water_spacing = config.get("water_spacing", 0.33)
salt_conc = config.get("salt_conc", 0.150)

# Cholesterol settings
percent_label = config.get("percent_label", 20)

additive_type = config.get("additive_type", "miglyol").lower()
if additive_type == "none" or percent_label == 0:
    additive_resname = "NONE"
    additive_gro = None
    additive_itp = None
elif additive_type == "miglyol":
    additive_resname = "MIG"
    additive_gro = "building_blocks/miglyol.gro"
    additive_itp = "MIG.itp"
elif additive_type in ["cholesterol", "chol"]:
    additive_resname = "CHL1"
    additive_gro = "building_blocks/chol_gmx.gro"
    additive_itp = "CHL1.itp"
else:
    # Generic additive handling
    additive_gro = f"building_blocks/{additive_type}.gro"
    if not os.path.exists(additive_gro):
        print(f"Error: Could not find {additive_gro}. Please place your additive structure in building_blocks/")
        sys.exit(1)
        
    # Dynamically extract residue name from the first atom in the .gro file
    with open(additive_gro, "r") as f:
        lines = f.readlines()
        if len(lines) > 2:
            first_atom_line = lines[2]
            # In GROMACS .gro format, chars 5-10 contain the residue name
            additive_resname = first_atom_line[5:10].strip()
        else:
            additive_resname = "UNK"
            
    # Assume the topology file matches the residue name
    additive_itp = f"{additive_resname}.itp"
    
    # Antechamber feature: if additive itp is not in toppar, build it!
    additive_itp_path = f"forcefield/toppar/{additive_itp}"
    if not os.path.exists(additive_itp_path) and additive_type != "none":
        print(f"Topology {additive_itp} not found in toppar/! Attempting to auto-generate using Antechamber/ACPYPE...")
        # Convert gro to pdb for antechamber using obabel or gmx trjconv
        # Here we just write the system call logic
        os.makedirs("scratch_ff", exist_ok=True)
        ret = os.system(f"gmx_mpi editconf -f {additive_gro} -o scratch_ff/temp.pdb")
        if ret == 0:
            print("Running Antechamber...")
            ret2 = os.system(f"antechamber -i scratch_ff/temp.pdb -fi pdb -o scratch_ff/temp.mol2 -fo mol2 -c bcc -s 2")
            if ret2 == 0:
                print("Running ACPYPE to generate GROMACS topology...")
                os.system(f"acpype -i scratch_ff/temp.mol2 -o gmx -b {additive_resname}")
                # Move generated itp
                os.system(f"cp {additive_resname}.GMX/{additive_resname}_GMX.itp {additive_itp_path}")
                print(f"Successfully generated {additive_itp}!")
            else:
                print("Antechamber failed. Ensure AmberTools is installed.")
                sys.exit(1)
        else:
            print("Failed to convert .gro to .pdb for Antechamber.")
            sys.exit(1)

if additive_type == "none" or percent_label == 0:
    n_additive_total = 0
else:
    n_additive_total = int((chains_per_leaflet * 2) * (percent_label / 100.0))
    # Ensure even number to split between leaflets equally
    if n_additive_total % 2 != 0:
        n_additive_total += 1
n_additive_per_leaflet = n_additive_total // 2

polymer_type = config.get("polymer_type", "pb22peo14").lower()
polymer_itp_path = f"forcefield/toppar/{polymer_type}.itp"

import re
# Auto-generate custom polymer topology if missing
if not os.path.exists(polymer_itp_path):
    print(f"Custom topology {polymer_itp_path} not found. Auto-generating using build_top.py...")
    match = re.match(r"pb(\d+)peo(\d+)", polymer_type)
    if match:
        target_pb = int(match.group(1))
        target_peo = int(match.group(2))
        print(f"Detected target: PB={target_pb}, PEO={target_peo}")
        
        ret = os.system(f"python3 scripts/build_top.py {target_pb} {target_peo} {polymer_itp_path}")
        if ret != 0:
            print("Error: Topology auto-generation failed.")
            sys.exit(1)
    else:
        print(f"Error: Polymer topology {polymer_itp_path} not found and naming format not recognized.")
        sys.exit(1)
additive_z_var = config.get("additive_z_variance_nm", 1.0)

n_additive_per_leaflet = n_additive_total // 2
additive_z_offset = config.get("additive_z_offset", 0.0)
additive_z_random = config.get("additive_z_random", 1.0)

work_dir = f"output_systems/{percent_label}_{polymer_type}"
if not os.path.exists(work_dir):
    os.makedirs(work_dir, exist_ok=True)

import shutil
print(f"Creating output directory {work_dir}/ and copying forcefield files...")
shutil.copytree('forcefield/charmm36-mar2019.ff', os.path.join(work_dir, 'charmm36-mar2019.ff'), dirs_exist_ok=True)
shutil.copytree('forcefield/toppar', os.path.join(work_dir, 'toppar'), dirs_exist_ok=True)
if os.path.exists('restrains'):
    shutil.copytree('restrains', os.path.join(work_dir, 'restrains'), dirs_exist_ok=True)

out_gro_name = os.path.join(work_dir, config.get("output_gro", "system_300.gro"))
out_top_name = os.path.join(work_dir, config.get("output_top", "system_300.top"))

print(f"Building membrane with {chains_per_leaflet*2} PB-PEO chains and {n_additive_total} cholesterol molecules...")

# 1. Parse itp_mol_01.itp to get atom order and bonds for polymer
atoms = []
bonds = []
with open(polymer_itp_path, "r") as f:
    in_atoms = False
    in_bonds = False
    for line in f:
        line = line.split(";")[0].strip()
        if not line: continue
        if line.startswith("[ atoms ]"):
            in_atoms = True
            in_bonds = False
            continue
        elif line.startswith("[ bonds ]"):
            in_atoms = False
            in_bonds = True
            continue
        elif line.startswith("["):
            in_atoms = False
            in_bonds = False
            
        if in_atoms:
            parts = line.split()
            if len(parts) >= 6:
                atom_id = int(parts[0])
                atom_type = parts[1]
                resid = int(parts[2])
                resname = parts[3]
                atomname = parts[4]
                atoms.append((atom_id, atom_type, resid, resname, atomname))
        
        if in_bonds:
            parts = line.split()
            if len(parts) >= 2:
                bonds.append((int(parts[0]), int(parts[1])))

# 2. Build a linear chain using graph traversal
coords = np.zeros((len(atoms), 3))
placed = set()
atom_dict = {a[0]: {"idx": i, "name": a[4], "resid": a[2]} for i, a in enumerate(atoms)}

adj = {a[0]: [] for a in atoms}
for a, b in bonds:
    if a in adj and b in adj:
        adj[a].append(b)
        adj[b].append(a)

backbone_names = {"C10", "C5", "C8", "C1", "C2", "O1"}
start_atom = 4
coords[atom_dict[start_atom]["idx"]] = [0.0, 0.0, 0.0]
placed.add(start_atom)
queue = [start_atom]
child_counts = {a[0]: 0 for a in atoms}
bb_counts = {a[0]: 0 for a in atoms}

h_offsets = [
    np.array([ 0.11,  0.0,   0.0]),
    np.array([-0.11,  0.0,   0.0]),
    np.array([ 0.0,   0.11,  0.0]),
    np.array([ 0.0,  -0.11,  0.0]),
    np.array([ 0.08,  0.08,  0.0]),
    np.array([-0.08, -0.08,  0.0])
]

while queue:
    curr = queue.pop(0)
    curr_idx = atom_dict[curr]["idx"]
    curr_pos = coords[curr_idx]
    
    for nbr in adj[curr]:
        if nbr not in placed:
            placed.add(nbr)
            nbr_idx = atom_dict[nbr]["idx"]
            nbr_name = atom_dict[nbr]["name"]
            
            offset = np.array([0.0, 0.0, 0.0])
            if "H" in nbr_name:
                c = child_counts[curr]
                offset = h_offsets[c % len(h_offsets)]
                child_counts[curr] += 1
            elif nbr_name in backbone_names:
                base_x = 0.10 if curr_idx % 2 == 0 else -0.10
                c = bb_counts[curr]
                if c == 0:
                    offset = np.array([base_x, 0.0, 0.12])
                elif c == 1:
                    offset = np.array([-base_x, 0.0, 0.12])
                elif c == 2:
                    offset = np.array([0.0, 0.10, 0.12])
                else:
                    offset = np.array([0.0, -0.10, 0.12])
                bb_counts[curr] += 1
            else:
                offset = np.array([0.0, 0.12, 0.05])
                if nbr_name == "C7":
                    offset = np.array([0.0, 0.12, -0.05])
            
            coords[nbr_idx] = curr_pos + offset
            queue.append(nbr)

for a in atoms:
    if a[0] not in placed:
        coords[atom_dict[a[0]]["idx"]] = [0, 0, -1.0]

# 3. Read Cholesterol Coordinates
mig_atoms = []
mig_coords = []
if n_additive_total > 0:
    with open(additive_gro, "r") as f:
        lines = f.readlines()[2:-1] # skip header and box
        for line in lines:
            resname = line[5:10].strip()
            aname = line[10:15].strip()
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            mig_atoms.append((resname, aname))
            mig_coords.append([x, y, z])
    mig_coords = np.array(mig_coords)
    
    # Center cholesterol
    mig_coords -= np.mean(mig_coords, axis=0)
    
    mig_coords_centered = mig_coords - np.mean(mig_coords, axis=0)


# 4. Build Bilayer
import math
# Use box_x and box_y loaded from config directly
# Calculate min_dist dynamically based on available area per chain
min_dist = math.sqrt((box_x * box_y) / chains_per_leaflet) * 0.70

def generate_random_positions(num_pts, bx, by, md):
    pts = []
    max_attempts = 500000
    attempts = 0
    while len(pts) < num_pts and attempts < max_attempts:
        x = random.uniform(0, bx)
        y = random.uniform(0, by)
        valid = True
        for px, py in pts:
            # Periodic boundary condition distance
            dx = abs(x - px)
            dy = abs(y - py)
            dx = min(dx, bx - dx)
            dy = min(dy, by - dy)
            if (dx*dx + dy*dy) < md*md:
                valid = False
                break
        if valid:
            pts.append((x, y))
        attempts += 1
        
    if len(pts) < num_pts:
        print(f"Error: Could only place {len(pts)} out of {num_pts} points. Reduce spacing or min_dist.")
        sys.exit(1)
    return pts

bilayer_coords = []
bilayer_resnames = []
bilayer_resids = []
bilayer_atomnames = []

coords[:, 0] -= np.mean(coords[:, 0])
coords[:, 1] -= np.mean(coords[:, 1])

single_chain_coords = coords
res_offset = 1

print(f"Generating random coordinates for box {box_x:.2f} x {box_y:.2f} nm...")
random.seed(42) # for reproducibility
selected_pts_lower = generate_random_positions(chains_per_leaflet, box_x, box_y, min_dist)
selected_pts_upper = generate_random_positions(chains_per_leaflet, box_x, box_y, min_dist)

print("Building lower leaflet (Polymers)...")
for x, y in selected_pts_lower:
    for i in range(len(atoms)):
        bilayer_coords.append([single_chain_coords[i, 0] + x, single_chain_coords[i, 1] + y, single_chain_coords[i, 2]])
        bilayer_resnames.append(atoms[i][3])
        bilayer_resids.append(res_offset)
        bilayer_atomnames.append(atoms[i][4])
    res_offset += 1

min_z_chain = np.min(single_chain_coords[:, 2])
z_shift_upper = 2 * min_z_chain - 0.1

print("Building upper leaflet (Polymers)...")
for x, y in selected_pts_upper:
    for i in range(len(atoms)):
        bilayer_coords.append([-single_chain_coords[i, 0] + x, single_chain_coords[i, 1] + y, -single_chain_coords[i, 2] + z_shift_upper])
        bilayer_resnames.append(atoms[i][3])
        bilayer_resids.append(res_offset)
        bilayer_atomnames.append(atoms[i][4])
    res_offset += 1

# 5. Insert Cholesterol Randomly
if n_additive_total > 0:
    print(f"Inserting {n_additive_total} Cholesterol molecules...")
    bilayer_arr = np.array(bilayer_coords)
    # box_x and box_y are already set from area calculation
    
    # Find Z ranges for the hydrophobic core and hydrophilic outer edge
    pb_indices_single = [i for i, a in enumerate(atoms) if a[3] == "BDE"]
    z_interface_lower = np.max(single_chain_coords[pb_indices_single, 2]) # true PB length
    z_interface_upper = z_shift_upper - z_interface_lower
    
    z_outer_lower = np.max(single_chain_coords[:, 2]) # PEO outer edge
    z_outer_upper = z_shift_upper - z_outer_lower
    
    placed_migs = 0
    max_attempts = 20000
    attempts = 0
    
    # Determine exact target depth from surface
    depth_config = config.get("additive_depth_from_surface_nm", "middle")
    thickness = z_interface_lower - min_z_chain
    if depth_config == "middle":
        D = thickness / 2.0
    else:
        D = float(depth_config)
    
    # Lower leaflet PB spans [min_z_chain, z_interface_lower] (surface is z_interface_lower)
    target_lower = z_interface_lower - D
    
    # Upper leaflet PB spans [-z_interface_lower + z_shift_upper, min_z_chain - 0.1] (surface is -z_interface_lower + z_shift_upper)
    target_upper = (-z_interface_lower + z_shift_upper) + D
    
    while placed_migs < n_additive_total and attempts < max_attempts:
        attempts += 1
        x = random.uniform(0, box_x)
        y = random.uniform(0, box_y)
        
        # Half in lower leaflet, half in upper leaflet
        if placed_migs < n_additive_total / 2.0:
            z = random.uniform(target_lower - additive_z_var, target_lower + additive_z_var)
        else:
            z = random.uniform(target_upper - additive_z_var, target_upper + additive_z_var)
            
        mig = mig_coords_centered.copy()
        
        # Random 3D rotation
        from scipy.spatial.transform import Rotation
        rot = Rotation.random().as_matrix()
        mig = np.dot(mig, rot.T)
        
        mig[:, 0] += x
        mig[:, 1] += y
        mig[:, 2] += z
        
        # Check clash
        tree = cKDTree(bilayer_arr)
        dists, _ = tree.query(mig, k=1)
        
        # 0.25 nm prevents LJ potentials from reaching floating point infinity
        if np.min(dists) > 0.25:
            bilayer_arr = np.vstack([bilayer_arr, mig])
            for i in range(len(mig)):
                bilayer_resnames.append(mig_atoms[i][0])
                bilayer_resids.append(res_offset)
                bilayer_atomnames.append(mig_atoms[i][1])
            res_offset += 1
            placed_migs += 1
            
    if placed_migs < n_additive_total:
        print(f"WARNING: Could only place {placed_migs}/{n_additive_total} miglyol molecules.")
            
    bilayer_coords = bilayer_arr.tolist()

bilayer_coords = np.array(bilayer_coords)
min_z = np.min(bilayer_coords[:, 2])
bilayer_coords[:, 2] -= min_z
bilayer_coords[:, 2] += z_padding

min_poly_z = np.min(bilayer_coords[:, 2])
max_poly_z = np.max(bilayer_coords[:, 2])
box_z = max_poly_z + z_padding

# 6. Solvate with water grid
print("Solvating using cKDTree...")
grid_pts = []
for x in np.arange(0, box_x, water_spacing):
    for y in np.arange(0, box_y, water_spacing):
        for z in np.arange(0, box_z, water_spacing):
            grid_pts.append([x, y, z])
grid_pts = np.array(grid_pts)

wrapped_coords = bilayer_coords.copy()
wrapped_coords[:, 0] %= box_x
wrapped_coords[:, 1] %= box_y
wrapped_coords[:, 2] %= box_z
tree = cKDTree(wrapped_coords, boxsize=[box_x, box_y, box_z])
dists, _ = tree.query(grid_pts, k=1)
valid_mask = dists >= 0.3

# Find Z-range of the ENTIRE polymer slab to exclude water from the bilayer completely
poly_z = bilayer_coords[:300*len(atoms), 2]
min_poly_z = np.min(poly_z) - 0.2
max_poly_z = np.max(poly_z) + 0.2

# Mask out water inside the entire polymer slab
core_mask = (grid_pts[:, 2] < min_poly_z) | (grid_pts[:, 2] > max_poly_z)
valid_mask = valid_mask & core_mask

water_coords = grid_pts[valid_mask]
num_water = len(water_coords)
print(f"Added {num_water} waters (excluded from entire bilayer Z range: {min_poly_z:.2f} to {max_poly_z:.2f}).")

# 7. Calculate 0.150 M NaCl and replace water
nacl_ratio = salt_conc / 55.5
num_na = int(round(num_water * nacl_ratio))
num_cl = num_na
total_ions = num_na + num_cl

print(f"Adding {num_na} NA and {num_cl} CL ions.")

replace_idx = random.sample(range(num_water), total_ions)
na_idx = set(replace_idx[:num_na])
cl_idx = set(replace_idx[num_na:])

# 8. Write .gro file
print(f"Writing {out_gro_name}...")
with open(out_gro_name, "w") as f:
    f.write(f"PB22-PEO14 Bilayer 300 Chains {n_additive_total} Chol with NaCl\n")
    f.write(f"{len(bilayer_coords) + (num_water - total_ions)*3 + total_ions}\n")
    
    for i in range(len(bilayer_coords)):
        resid = (bilayer_resids[i] % 100000)
        resname = bilayer_resnames[i]
        aname = bilayer_atomnames[i]
        x, y, z = bilayer_coords[i]
        f.write(f"{resid:5d}{resname:<5s}{aname:>5s}{(i+1)%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
    
    atom_id = len(bilayer_coords) + 1
    w_resid = 1
    actual_water_count = 0
    
    for w in na_idx:
        x, y, z = water_coords[w]
        f.write(f"{(w_resid)%100000:5d}{'NA':<5s}{'NA':>5s}{atom_id%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        atom_id += 1
        w_resid += 1
        
    for w in cl_idx:
        x, y, z = water_coords[w]
        f.write(f"{(w_resid)%100000:5d}{'CL':<5s}{'CL':>5s}{atom_id%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        atom_id += 1
        w_resid += 1
        
    for w in range(num_water):
        if w not in na_idx and w not in cl_idx:
            x, y, z = water_coords[w]
            actual_water_count += 1
            f.write(f"{(w_resid)%100000:5d}{'TIP3':<5s}{'OH2':>5s}{atom_id%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
            f.write(f"{(w_resid)%100000:5d}{'TIP3':<5s}{'H1':>5s}{(atom_id+1)%100000:5d}{x+0.05:8.3f}{y+0.05:8.3f}{z+0.05:8.3f}\n")
            f.write(f"{(w_resid)%100000:5d}{'TIP3':<5s}{'H2':>5s}{(atom_id+2)%100000:5d}{x-0.05:8.3f}{y-0.05:8.3f}{z+0.05:8.3f}\n")
            atom_id += 3
            w_resid += 1
            
    f.write(f"{box_x:10.5f} {box_y:10.5f} {box_z:10.5f}\n")

# 9. Write .top file
with open(out_top_name, "w") as f:
    f.write('#include "charmm36-mar2019.ff/forcefield.itp"\n')
    f.write('#include "toppar/block_charm.itp"\n')
    f.write(f'#include "toppar/{polymer_type}.itp"\n')
    f.write('#include "charmm36-mar2019.ff/ions.itp"\n')
    if n_additive_total > 0:
        f.write(f'#include "toppar/{additive_itp}"\n')
    f.write('#include "toppar/TIP3.itp"\n\n')
    f.write('[system]\n')
    f.write(f'PB22-PEO14 Bilayer 300 Chains {n_additive_total} Chol\n\n')
    f.write('[molecules]\n')
    f.write(f'mol_01 300\n')
    if n_additive_total > 0:
        f.write(f'{additive_resname} {n_additive_total}\n')
    f.write(f'NA {num_na}\n')
    f.write(f'CL {num_cl}\n')
    f.write(f'TIP3 {actual_water_count}\n')

print(f"Generated {out_gro_name} and {out_top_name}.")
