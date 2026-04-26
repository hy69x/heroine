# Heroine Simulation Variables

This document serves as the single source of truth for all configurable variables in the **Heroine** drug encapsulation pipeline.

## Input & Composition
| Flag | Variable | Default | Description |
| :--- | :--- | :--- | :--- |
| `--polymer` | `polymer` | **Required** | The SMILES string or common name of the polymer chain. Names are fetched from PubChem automatically. |
| `--drug` | `drug` | **Required** | The SMILES string or common name of the drug/payload molecule. Names are fetched from PubChem automatically. |
| `--name` | `name` | `encapsulation_run` | Name of the job (creates `output/<name>` directory). |
| `--polymer_count` | `polymer_count` | `1` | Number of polymer chains to simulate. |
| `--drug_count` | `drug_count` | `5` | Number of drug molecules to inject. |
| `--report` | `report` | `False` | Generate a professional PDF lab report (`lab_report.pdf`). |
| `--render` | `render` | `False` | Render a GIF animation of the simulation trajectory. |
| `--plot` | `plot` | `True` | Generates a stability plot (`stability.png`). |

## Physics & Environment
| Flag | Variable | Default | Description |
| :--- | :--- | :--- | :--- |
| `--steps` | `steps` | `50000` | **Total Simulation Time**. Total number of integration steps to run. |
| `--temp` | `temp` | `300.0` | **Temperature (K)**. Controls thermal energy in the system. |
| `--damp` | `damp` | `20.0` | **Damping (fs)**. Langevin thermostat parameter (viscosity). |
| `--timestep` | `timestep` | `1.0` | **Time Step (fs)**. Resolution of the simulation integration. |
| `--padding` | `padding` | `30.0` | **Padding (Å)**. Extra space around molecules to determine Simulation Box size. |

> **Note on Duration:** The total physical time simulated is `steps * timestep`. For example, 50,000 steps at 1.0 fs = 50 picoseconds (ps).

## Force Field (Lennard-Jones)
*Note: CHONS atoms (C, H, O, N, S) use specific OPLS-AA parameters by default.*

- **Epsilon**: Interaction strength (kcal/mol).
- **Sigma**: Particle size (Å).
- **Scaling**: 1-4 interactions are scaled by 0.5 (standard molecular mechanics).

## Standard Outputs
All files are saved in `output/<name>/`:

| File | Description |
| :--- | :--- |
| `system.data` | LAMMPS topology file (Atoms, Bonds, Angles, Dihedrals). |
| `simulation.in` | LAMMPS input script. |
| `simulation.log` | Raw thermodynamic data log. |
| `trajectory.dump` | Atom positions over time (coordinates). |
| `stability.png` | Graph showing Temperature and Potential Energy equilibrium. |
| `encapsulation.gif` | Trajectory animation (if `--render` is used). |
| `lab_report.pdf` | Professional PDF summary (if `--report` is used). |
