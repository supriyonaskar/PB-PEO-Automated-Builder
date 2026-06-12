# Membrane Building Methodology

This repository contains a unified, production-level pipeline for generating and simulating complex polymer-additive membrane mixtures (e.g. PB-PEO membranes with varying percentages of Cholesterol or Miglyol).

## Overview

The methodology leverages a highly robust and generic python/bash pipeline that allows you to dynamically inject additives into the hydrophobic core of the bilayer, clear penetrating water molecules, and seamlessly set up the entire simulation for GROMACS.

## Features
- **Dynamic Additive Placement**: Automatically splits and places your target additive (e.g. Miglyol or Cholesterol) evenly between leaflets based on configured surface depth and core limits.
- **Accurate Core Clearance**: A post-solvation script (`delete_waters.py`) specifically targets the hydrophobic core defined by the polymer boundaries and securely deletes all penetrating `SOL` or `TIP3` molecules.
- **Robust Group Generation**: Uses advanced index group regeneration (`keep 0`) and wildcards to seamlessly patch the system topologies and prevent GROMACS `make_ndx` crashing, even across systems with varying salt concentrations.
- **Parametrized Configuration**: Almost all settings (chains per leaflet, percentage, box size, depth) can be modified in a single `config.json`.

## Repository Structure

- `scripts/`: Contains the core automation logic.
  - `builder.sh`: The master bash script to initiate the system build.
  - `build_system.py`: Python script that reads `config.json`, generates random coordinates, packs the polymer, and inserts the additive (Cholesterol/Miglyol) using `cKDTree`.
  - `delete_waters.py`: Scans the system post-solvation and deletes any waters penetrating the membrane core.
- `mdp_templates/`: Contains all GROMACS `.mdp` files required for minimization, equilibration (steps 6.1 through 6.6), and 10ns production. Also includes standard cluster `job.sh` templates.
- `building_blocks/`: Stores the single-molecule coordinate files (`.gro`) for your additives (e.g., `miglyol.gro`, `chol_gmx.gro`).
- `forcefield/`: Stores the global `charmm36-mar2019.ff` and `toppar` parameters needed for system generation.
- `examples/`: Example configurations to test run the pipeline.

## Usage

You can use the unified `builder.sh` to generate systems on the fly using a configuration JSON file.

**Check the help menu at any time:**
```bash
./scripts/builder.sh --help
```

**Run the builder:**
```bash
./scripts/builder.sh examples/test_case.json
```

### JSON Configuration (`examples/test_case.json`)
```json
{
    "percent_label": 10,
    "box_x": 13.0,
    "box_y": 13.0,
    "chains_per_leaflet": 150,
    "z_padding": 3.5,
    "water_spacing": 0.33,
    "salt_conc": 0.150,
    "additive_type": "miglyol",
    "polymer_type": "pb22peo14",
    "additive_z_variance_nm": 1.0,
    "additive_depth_from_surface_nm": "middle"
}
```

### What `builder.sh` does:
1. Reads your input JSON file.
2. Runs `build_system.py` which:
   - Copies the necessary forcefield parameters.
   - Packs the required number of polymer chains.
   - Disperses the appropriate number of additive molecules into the hydrophobic core.
   - Solvates the system and adds neutralizing ions.
3. Stages the folder with the `.mdp` files and the GROMACS execution pipeline (`job.sh`).
4. Generates a perfectly robust `index.ndx` using `gmx make_ndx`.

## Adding Custom Additives
The pipeline is fully generalized and can accept **any** arbitrary additive automatically, without changing any code.

To add a new additive (e.g., `squalene`):
1. **Provide the Coordinate File**: Save a single-molecule `.gro` file inside the `building_blocks/` directory, named exactly after your additive type (e.g., `building_blocks/squalene.gro`).
2. **Provide the Topology File**: Ensure the `.itp` file for your molecule is placed in the `forcefield/toppar/` directory. The `.itp` file must be named *exactly* after the 3-5 letter residue name found inside your `.gro` file (e.g., if your `.gro` file uses the residue name `SQUA`, then your topology file must be `forcefield/toppar/SQUA.itp`).

Then, simply run the builder with your custom additive name:
```bash
./scripts/builder.sh 15 squalene pb22peo14
```
The python builder will dynamically parse your `.gro` file, extract the correct residue name, map the corresponding `.itp` file, and inject it perfectly into the core of your membrane!

## Building Pure Polymer Systems (No Additives)
You can easily use this pipeline to build a pure `pb_x_peo_y` bilayer without any additives. 

To do this, specify `none` as your additive type in your JSON config:
```json
{
    "percent_label": 0,
    "additive_type": "none",
    "polymer_type": "pb51peo27",
    ...
}
```
**Important Note:** For custom polymer sizes (like `pb51peo27`), simply place your `pb51peo27.itp` file into `forcefield/toppar/`. The builder will dynamically read the connectivity graph from your `.itp` file and fold the 3D chain automatically!

To run the simulation, simply navigate to the newly created percentage directory and execute `sbatch job.sh`!
