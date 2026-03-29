"""
Heroine: Drug Encapsulation Simulator
Extract from Proteus Project
"""

import argparse
import sys
import os
from pathlib import Path
from rdkit import Chem
from datetime import datetime

# Add src to python path
sys.path.append(str(Path(__file__).parent / "src"))

from src import topology, simulation, analysis, visualization, report

def validate_smiles(smiles: str, label: str = "Molecule"):
    """Validates a SMILES string using RDKit."""
    if not smiles: return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Invalid {label} SMILES: {smiles}")
        return None
    return smiles

def run_encapsulation_pipeline(args):
    """
    Orchestrates the drug encapsulation pipeline.
    """
    # Validation
    if not validate_smiles(args.polymer, "Polymer"):
        raise ValueError(f"Invalid Polymer SMILES: {args.polymer}")
    if not validate_smiles(args.drug, "Drug"):
        raise ValueError(f"Invalid Drug SMILES: {args.drug}")
    
    # Construct system SMILES
    polymer_list = [args.polymer] * args.polymer_count
    drug_list = [args.drug] * args.drug_count
    system_smiles = ".".join(polymer_list + drug_list)
    
    # Setup Paths
    output_dir = Path(os.getcwd()) / "output" / args.name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    paths = {
        "data": output_dir / "system.data",
        "input": output_dir / "simulation.in",
        "log": output_dir / "simulation.log",
        "dump": output_dir / "trajectory.dump",
        "gif": output_dir / "encapsulation.gif",
        "plot": output_dir / "stability.png" if args.plot else None
    }
    
    print("=" * 50)
    print(f"HEROINE: Drug Encapsulation Pipeline - {args.name}")
    print(f"System: {args.polymer_count}x Polymer, {args.drug_count}x Drug")
    print("=" * 50)
    
    # 1. Topology
    print("[*] Phase 1: Topology Generation (Combining Polymer & Drug)")
    bond_p, angle_p, dihedral_p = topology.generate_topology(system_smiles, paths["data"], padding=args.padding)
    
    # 2. Simulation
    print("[*] Phase 2: Molecular Dynamics Simulation (LAMMPS)")
    simulation.generate_input_file(
        paths["data"], paths["input"], paths["dump"], 
        steps=args.steps, temp=args.temp, damp=args.damp,
        timestep=args.timestep,
        bond_params=bond_p, angle_params=angle_p, dihedral_params=dihedral_p
    )
    simulation.run_simulation(paths["input"], paths["log"])
    
    # 3. Analysis
    print("[*] Phase 3: Analytics (Calculating Encapsulation Efficiency)")
    results = analysis.analyze_results(
        paths["log"], 
        output_plot=paths["plot"],
        polymer_count=args.polymer_count,
        payload_count=args.drug_count,
        dump_path=paths["dump"]
    )

    # 4. Visualization (Optional)
    if args.render:
        print("[*] Phase 4: Visualization (Generating Trajectory GIF)")
        visualization.render_trajectory(dump_path=paths["dump"], output_gif=paths["gif"])

    # 5. Automated Lab Notebook (Optional)
    if args.report:
        print("[*] Phase 5: Automated Lab Notebook (Generating PDF Report)")
        report.generate_report(
            output_dir,
            args.name,
            args.polymer,
            args.steps,
            args.temp,
            results["rg"],
            efficiency=results.get("efficiency"),
            plot_path=paths["plot"]
        )
    
    print("=" * 50)
    print(f"HEROINE: Simulation Finished Successfully for: {args.name}")
    print(f"Encapsulation Efficiency: {results.get('efficiency'):.2f}%")
    print(f"Output Directory: {output_dir}")
    print("=" * 50)

def main():
    parser = argparse.ArgumentParser(description="Heroine: Drug Encapsulation Simulator")
    parser.add_argument("--polymer", type=str, required=True, help="SMILES string of the polymer")
    parser.add_argument("--drug", type=str, required=True, help="SMILES string of the drug/payload")
    parser.add_argument("--name", type=str, default="encapsulation_run", help="Name of the run")
    parser.add_argument("--polymer_count", type=int, default=1, help="Number of polymer chains")
    parser.add_argument("--drug_count", type=int, default=5, help="Number of drug molecules")
    parser.add_argument("--steps", type=int, default=50000, help="Number of simulation steps")
    parser.add_argument("--temp", type=float, default=300.0, help="Temperature (K)")
    parser.add_argument("--damp", type=float, default=20.0, help="Langevin damping parameter")
    parser.add_argument("--timestep", type=float, default=1.0, help="Timestep (fs)")
    parser.add_argument("--padding", type=float, default=30.0, help="Box padding (Angstroms)")
    parser.add_argument("--render", action="store_true", help="Render a GIF")
    parser.add_argument("--plot", action="store_true", default=True, help="Generate stability plot")
    parser.add_argument("--report", action="store_true", help="Generate PDF report")
    
    args = parser.parse_args()
    
    try:
        run_encapsulation_pipeline(args)
    except Exception as e:
        print(f"\n[!] Heroine Failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
