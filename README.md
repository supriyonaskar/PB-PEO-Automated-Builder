# PB-PEO Block Copolymer Membrane Builder 🚀

Welcome to the fully automated, highly-scalable **PB-PEO Membrane Builder Pipeline**!

This repository contains a unified, mathematically generalized pipeline designed to automatically construct, solvate, and equilibrate massive block-copolymer membranes using GROMACS.

## ✨ Key Features
- **Dynamic Topology Generation (`pb_X_peo_Y`):** Automatically extrapolate and build custom polymer topologies (e.g. `pb100peo50`) from a template on the fly.
- **Arbitrary Additives:** Automatically packs any concentration of custom additives (e.g., Miglyol, Cholesterol) perfectly distributed in the membrane core.
- **Automated Geometry:** Folds polymers into 3D using Breadth-First-Search (BFS) and packs the bilayer.
- **Flawless Solvation:** Uses high-performance K-D Trees to insert water *around* the membrane without penetrating the hydrophobic core.
- **Ready-to-Run Pipelines:** Spits out perfect GROMACS index (`.ndx`) files, aligns constraints, and provides SLURM `job.sh` execution scripts.

---

## 🛠️ Usage & JSON Configuration

The entire simulation is controlled via a single JSON configuration file. 
You can run the builder by pointing it to your config:

```bash
./scripts/builder.sh examples/test_case.json
```

**Example JSON (`examples/test_case.json`):**
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

---

## 1. Building Arbitrary Polymer Sizes (`pb_X_peo_Y`)

This repository is powered by a custom topology extrapolator (`build_top.py`) that can mathematically scale up the standard `pb22peo14` forcefield to generate the topology for **any** custom block copolymer length!

If you want to simulate a membrane with a massive `pb51peo27` polymer:
1. Simply set `"polymer_type": "pb51peo27"` in your JSON config.
2. If the file `forcefield/toppar/pb51peo27.itp` does not exist, the pipeline will **automatically trigger `build_top.py`**.
3. It will extract your target numbers (PB=51, PEO=27), duplicate the internal repeating structural units from the `pb22peo14` template, mathematically shift all the atom numbers, and generate your new `.itp` file on the fly!

*To build a pure membrane without any additives, just set `"additive_type": "none"` and `"percent_label": 0`.*

---

## 2. Adding Custom Additives

Once your polymer is set, you can easily dope the membrane with custom molecules! The pipeline is completely generalized and will mathematically space the molecules inside the hydrophobic core.

To add a brand new additive (e.g., `squalene`):
1. **Coordinate File:** Place your single-molecule `.gro` file into `building_blocks/squalene.gro`.
2. **Topology File:** Place your molecule's `.itp` file into `forcefield/toppar/SQUA.itp` (make sure the filename matches the 3-5 letter residue name defined inside your `.gro` file).
3. **JSON Config:** Set `"additive_type": "squalene"` and set your desired percentage (e.g. `"percent_label": 15`).

The Python builder will dynamically parse your `.gro` file, match it to the `.itp` forcefield, and seamlessly inject it into the hydrophobic core of the bilayer, exactly centered and distributed!

---

## 🚀 Running the Simulation

When `./scripts/builder.sh` finishes, it will generate a self-contained directory named after your percentage (e.g. `10/` or `0/`).

To launch the 6-step equilibration and production run on your HPC cluster, simply:
```bash
cd 10/
sbatch job.sh
```

Enjoy your automated membrane building! 🎉
