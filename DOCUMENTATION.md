# HEROINE: High-Efficiency Resolution of Integrated Nanoparticle Encapsulation

Welcome to the documentation for **HEROINE**, a specialized Molecular Dynamics pipeline designed for the pharmaceutical research and cheminformatics community.

## 📖 Table of Contents
1. [Introduction](#introduction)
2. [Installation & Setup](#installation--setup)
3. [Core Concepts](#core-concepts)
4. [Using the Web UI](#using-the-web-ui)
5. [Using the Command Line](#using-the-command-line)
6. [Technical Pipeline](#technical-pipeline)
7. [Scientific Metrics](#scientific-metrics)
8. [Troubleshooting](#troubleshooting)

---

## 1. Introduction
HEROINE is an automated "Virtual Laboratory" for simulating drug encapsulation within polymer nanoparticles. By bridging RDKit's cheminformatics capabilities with the LAMMPS physics engine, it allows researchers to predict the compatibility of drug-polymer pairs without the immediate need for expensive wet-lab experimentation.

## 2. Installation & Setup
The project uses `uv` for modern, fast Python dependency management.

### Prerequisites
- **Python 3.12+**
- **LAMMPS**: The `lmp` executable must be in your system PATH.
- **Ovito (Optional)**: Required for rendering trajectory GIFs.

### Installation
```bash
# Sync dependencies
uv sync
```

---

## 3. Core Concepts
- **SMILES**: HEROINE uses SMILES strings as the primary chemical identifier.
- **Langevin Dynamics**: Simulations use an implicit solvent model to represent the viscosity and thermal kicks of water without the computational overhead of individual water molecules.
- **Topology**: The system automatically assigns bond, angle, and dihedral parameters based on a heuristic OPLS-AA-like force field for CHONS atoms.

---

## 4. Using the Web UI
Launch the interactive interface:
```bash
uv run ui.py
```
### Workflow:
1. **Search**: Enter a molecule name (e.g., "aspirin") or SMILES.
2. **Preview**: Click "🔍 Preview" to verify the structure via RDKit.
3. **Build**: Add polymers and drugs to your "System Composition." You can add multiple different drug types.
4. **Simulate**: Adjust steps and temperature, then hit "🔥 Run Simulation."
5. **Analyze**: Download the generated **Lab Report PDF** and view stability plots directly in the browser.

---

## 5. Using the Command Line
For batch processing or remote servers:
```bash
uv run main.py --polymer "PEG" --drug "aspirin" --drug_count 10 --steps 50000 --name "my_simulation"
```
**Key Flags:**
- `--polymer` / `--drug`: Names or SMILES.
- `--polymer_count` / `--drug_count`: Quantity of molecules.
- `--steps`: Simulation duration (default 50,000).
- `--report`: Generate the PDF lab report.

---

## 6. Technical Pipeline
1. **Resolution**: `utils.py` queries PubChem API for names.
2. **Topology**: `topology.py` builds the 3D system and LAMMPS data file.
3. **MD Execution**: `simulation.py` generates the `.in` script and manages the `lmp` subprocess.
4. **Analytics**: `analysis.py` parses the trajectory and log for physical metrics.
5. **Reporting**: `report.py` assembles the final PDF.

---

## 7. Scientific Metrics
### Encapsulation Efficiency (%)
Calculated by verifying if a drug's center of mass is within $5.0 \text{\AA}$ of any polymer atom in the final frame. This proximity-based approach is robust for non-spherical nanoparticles.

### Radius of Gyration ($R_g$)
Tracks the compactness of the polymer. A decreasing $R_g$ indicates successful nanoparticle formation (folding) around the payload.

### Thermodynamic Stability
Monitored via Potential Energy and Temperature convergence plots to ensure the simulation reached physical equilibrium.

---

## 8. Troubleshooting
- **"lmp not found"**: Ensure LAMMPS is installed. Try `lmp -h` in your terminal to verify.
- **"Invalid SMILES"**: If a name search fails, verify the name on [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
- **"Lost Atoms"**: This usually means the simulation timestep is too large or the system energy is too high. Try reducing the `--temp` or increasing `--padding`.

---
*HEROINE is part of the Proteus Project Suite.*
