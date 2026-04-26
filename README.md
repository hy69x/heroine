# HEROINE: Drug Encapsulation Simulator

**HEROINE** (Extract from the Proteus Project) is an automated pipeline for simulating the encapsulation of drug molecules within polymer nanoparticles. It converts chemical SMILES into physical topologies, executes LAMMPS molecular dynamics, and calculates encapsulation efficiency.

## 🌟 Features

- **Automated SMILES Fetching**: Provide a drug or polymer name (e.g., "aspirin"), and HEROINE will fetch the SMILES from PubChem.
- **Physical Topology**: Uses RDKit to generate 3D coordinates and LAMMPS topology with heuristic OPLS-AA parameters.
- **Molecular Dynamics**: High-performance simulation using LAMMPS (supports GPU acceleration if available).
- **Robust Analytics**: Calculates Radius of Gyration and encapsulation efficiency using a proximity-based algorithm.
- **Smart Pipeline**: Built-in checkpointing skips completed phases, and improved error diagnostics parse LAMMPS logs for failures.

## 🚀 Quick Start (Web UI)

HEROINE now includes a modern, easy-to-use Gradio Web UI for configuring simulations directly from your browser without using the command line.

To launch the Web UI:

```bash
uv run ui.py
```
Then open your browser to `http://127.0.0.1:7860`. The UI allows you to search for molecules by name, adjust simulation counts, and view the generated PDF report, plots, and GIFs directly!

## 🚀 Quick Start (CLI)

To run a drug encapsulation simulation:

```bash
uv run main.py --polymer "C=CC(=O)N" --drug "CC(C)C1=C(C(=C(C=C1)O)C)O" --name "my_run" --polymer_count 1 --drug_count 10
```

### Key Arguments:
- `--polymer`: SMILES string of the polymer chain.
- `--drug`: SMILES string of the drug/payload molecule.
- `--polymer_count`: Number of polymer chains to simulate (default: 1).
- `--drug_count`: Number of drug molecules to simulate (default: 5).
- `--steps`: Number of simulation steps (default: 50,000).
- `--render`: Generate a GIF animation of the trajectory (requires Ovito).
- `--report`: Generate a professional PDF lab report summarizing the results.

## 📦 Requirements

- **Python**: 3.12+ (managed by `uv`).
- **LAMMPS**: The simulation engine `lmp` must be in your system PATH.
- **RDKit**: Used for chemical topology and geometry optimization.
- **Ovito (Optional)**: Used for 3D trajectory rendering.

## 📁 Output Structure

All results are saved in the `output/<run_name>/` directory:
- `system.data`: LAMMPS topology file.
- `simulation.in`: LAMMPS input script.
- `trajectory.dump`: Raw coordinate data.
- `stability.png`: Stability plot (Temperature & Potential Energy).
- `encapsulation.gif`: 3D animation (if `--render` is used).
- `lab_report.pdf`: Professional summary (if `--report` is used).

## 🛠 Tech Stack
- **Dependency Management**: [uv](https://github.com/astral-sh/uv)
- **Molecular Dynamics**: [LAMMPS](https://www.lammps.org/)
- **Cheminformatics**: [RDKit](https://www.rdkit.org/)
- **Visualization**: [Ovito](https://www.ovito.org/)
- **Analytics**: NumPy, Matplotlib, ReportLab
