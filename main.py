"""
HEROINE: Drug Encapsulation Simulator
Part of the Proteus Project Suite
"""

import argparse
import sys
import os
from pathlib import Path
from rdkit import Chem
from datetime import datetime

# Add src to python path
sys.path.append(str(Path(__file__).parent / "src"))

from src import topology, simulation, analysis, visualization, report, utils

def resolve_molecule(identifier: str, label: str = "Molecule"):
    """
    Validates identifier as SMILES, or attempts to fetch it from PubChem if it's a name.
    """
    if not identifier: return None
    
    # Try RDKit validation first
    mol = Chem.MolFromSmiles(identifier)
    if mol is not None:
        return identifier
    
    # If invalid SMILES, try fetching as a name
    print(f"[*] '{identifier}' is not a valid SMILES. Attempting to fetch from PubChem...")
    fetched = utils.fetch_smiles_by_name(identifier)
    if fetched:
        return fetched
        
    print(f"Error: Could not resolve {label}: {identifier}")
    return None

def run_encapsulation_pipeline(args):
    """
    Orchestrates the drug encapsulation pipeline.
    """
    # Validation and Resolution
    polymer_smiles = resolve_molecule(args.polymer, "Polymer")
    if not polymer_smiles:
        raise ValueError(f"Invalid Polymer identifier: {args.polymer}")
        
    drug_smiles = resolve_molecule(args.drug, "Drug")
    if not drug_smiles:
        raise ValueError(f"Invalid Drug identifier: {args.drug}")
    
    # Construct system SMILES
    # If UI passed total_poly_count, it means it already pre-expanded the strings.
    if hasattr(args, "total_poly_count"):
        system_smiles = f"{polymer_smiles}.{drug_smiles}"
        p_count_for_analysis = args.total_poly_count
        d_count_for_analysis = args.total_drug_count
    else:
        polymer_list = [polymer_smiles] * args.polymer_count
        drug_list = [drug_smiles] * args.drug_count
        system_smiles = ".".join(polymer_list + drug_list)
        p_count_for_analysis = args.polymer_count
        d_count_for_analysis = args.drug_count
    
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
    
    # Checkpointing: Check if simulation is already done
    skip_simulation = False
    if paths["log"].exists():
        with open(paths["log"], "r") as f:
            content = f.read()
            if "Total wall time:" in content or "Loop time" in content:
                print(f"[*] Simulation log found and appears complete. Skipping MD Phase.")
                skip_simulation = True

    # 1. Topology
    if not skip_simulation or not paths["data"].exists():
        print("[*] Phase 1: Topology Generation (Combining Polymer & Drug)")
        bond_p, angle_p, dihedral_p = topology.generate_topology(system_smiles, paths["data"], padding=args.padding)
    else:
        print("[*] Phase 1: Skipping Topology Generation (Using existing data file)")
        # We need these params if we were to rerun simulation, but if we skip simulation we don't strictly need them
        # unless they are needed for analysis (they aren't currently).
        bond_p, angle_p, dihedral_p = {}, {}, {}
    
    # 2. Simulation
    if not skip_simulation:
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
        polymer_count=p_count_for_analysis,
        payload_count=d_count_for_analysis,
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
            polymer_smiles,
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
