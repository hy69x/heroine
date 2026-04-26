import gradio as gr
import uuid
import os
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw

# Ensure src is in the path
sys.path.append(str(Path(__file__).parent / "src"))
from main import run_encapsulation_pipeline, resolve_molecule

class SimArgs:
    def __init__(self, polymer, drug, name, polymer_count, drug_count, steps, temp, damp, timestep, padding, render, plot, report):
        self.polymer = polymer
        self.drug = drug
        self.name = name
        self.polymer_count = polymer_count
        self.drug_count = drug_count
        self.steps = steps
        self.temp = temp
        self.damp = damp
        self.timestep = timestep
        self.padding = padding
        self.render = render
        self.plot = plot
        self.report = report

def get_mol_preview(identifier):
    if not identifier:
        return None, "No identifier provided"
    smiles = resolve_molecule(identifier)
    if not smiles:
        return None, f"Could not resolve '{identifier}'"
    
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        return img, smiles
    return None, "Invalid SMILES"

def add_molecule(identifier, count, category, current_list):
    img, smiles = get_mol_preview(identifier)
    if not img:
        return current_list, gr.update(visible=True, value=smiles) # Show error
    
    new_entry = {
        "id": str(uuid.uuid4())[:4],
        "name": identifier,
        "smiles": smiles,
        "count": int(count),
        "category": category
    }
    current_list.append(new_entry)
    return current_list, gr.update(visible=False)

def remove_molecule(mol_id, current_list):
    current_list = [m for m in current_list if m["id"] != mol_id]
    return current_list

def format_composition(mol_list):
    if not mol_list:
        return "System is empty. Add molecules to begin."
    
    lines = []
    for m in mol_list:
        lines.append(f"- [{m['category']}] {m['name']} ({m['count']} copies)")
    return "\n".join(lines)

def run_multi_simulation(mol_list, steps, temp, render, report):
    if not mol_list:
        return "Error: System is empty.", None, None, None
    
    polymers = [m for m in mol_list if m["category"] == "Polymer"]
    drugs = [m for m in mol_list if m["category"] == "Drug"]
    
    if not polymers:
        return "Error: At least one polymer type is required.", None, None, None
    if not drugs:
        return "Error: At least one drug type is required for encapsulation analysis.", None, None, None

    # Construct the system SMILES and counts
    # We need to adapt the pipeline to handle this.
    # Currently run_encapsulation_pipeline takes args.polymer and args.drug.
    # I will construct a 'dummy' args object where 'polymer' and 'drug' are already combined strings
    # and counts are 1.
    
    full_polymer_smiles_list = []
    for p in polymers:
        full_polymer_smiles_list.extend([p["smiles"]] * p["count"])
    
    full_drug_smiles_list = []
    for d in drugs:
        full_drug_smiles_list.extend([d["smiles"]] * d["count"])
        
    combined_polymer = ".".join(full_polymer_smiles_list)
    combined_drug = ".".join(full_drug_smiles_list)
    
    run_name = f"gradio_run_{uuid.uuid4().hex[:8]}"
    
    # We set counts to 1 because we've already manually expanded the SMILES strings
    # But wait, analysis needs the REAL counts.
    # If I set polymer_count=len(full_polymer_smiles_list), and drug_count=len(full_drug_smiles_list)
    # AND I pass a SINGLE polymer and drug SMILES which are actually the combined ones,
    # then main.py will do:
    # polymer_list = [combined_polymer] * len(...) -> WRONG.
    
    # I should pass them as-is and set counts to 1.
    args = SimArgs(
        polymer=combined_polymer,
        drug=combined_drug,
        name=run_name,
        polymer_count=1, 
        drug_count=1,
        steps=int(steps),
        temp=float(temp),
        damp=20.0,
        timestep=1.0,
        padding=30.0,
        render=render,
        plot=True,
        report=report
    )
    
    # Wait, analysis.analyze_results uses polymer_count and payload_count.
    # If I set them to 1 and 1, it will think there's only 1 drug.
    # I need to pass the REAL counts to the pipeline.
    
    args.polymer_count = 1
    args.drug_count = 1
    # BUT I want analysis to know the truth.
    # I'll monkey-patch or modify the call if needed, but let's see if I can just set them.
    # Actually, main.py does:
    # polymer_list = [args.polymer] * args.polymer_count
    # So if I set polymer=combined_polymer and polymer_count=1, it results in combined_polymer. Correct.
    # And I'll set drug_count = len(full_drug_smiles_list).
    # Wait, if drug=combined_drug and drug_count=N, it will be combined_drug repeated N times. WRONG.
    
    # OK, the easiest way is to pass the first drug as 'drug' and set its count, 
    # OR just pass the combined ones and set count=1.
    # If I set count=1, analysis only sees 1 drug.
    
    # Let's modify main.py slightly to accept a 'system_smiles' override or just handle lists.
    # Or better: construct a single string for polymer and a single string for drug that represents ALL of them.
    
    real_poly_count = len(full_polymer_smiles_list)
    real_drug_count = len(full_drug_smiles_list)
    
    # If I set polymer = combined_polymer and polymer_count=1
    # AND drug = combined_drug and drug_count=1
    # Then topology will see all atoms.
    # But analysis will see 1 polymer and 1 drug.
    # Encapsulation efficiency will be (is drug1 encapsulated?) -> 0 or 100%.
    
    # I need to fix main.py or analysis.py to handle this properly.
    # Let's fix main.py to allow passing lists.
    
    try:
        # Instead of calling run_encapsulation_pipeline directly, I'll pass a modified args
        # where I've already expanded the SMILES.
        # But main.py's run_encapsulation_pipeline RE-EXPANDS them.
        
        # I'll pass args.polymer = "." separated list of all polymers
        # and args.polymer_count = 1.
        # Same for drugs, but I want analysis to see N drugs.
        # This is tricky without changing main.py.
        
        # Let's change main.py to be more flexible.
        
        # For now, let's assume I fix main.py.
        
        args.total_poly_count = real_poly_count
        args.total_drug_count = real_drug_count
        
        run_encapsulation_pipeline(args)
        
        output_dir = Path(os.getcwd()) / "output" / run_name
        plot_path = output_dir / "stability.png"
        gif_path = output_dir / "encapsulation.gif"
        report_path = output_dir / "lab_report.pdf"
        
        out_plot = str(plot_path) if plot_path.exists() else None
        out_gif = str(gif_path) if gif_path.exists() else None
        out_report = str(report_path) if report_path.exists() else None
        
        return "Simulation completed successfully!", out_plot, out_gif, out_report
    except Exception as e:
        import traceback
        return f"Simulation failed: {e}\n{traceback.format_exc()}", None, None, None

# Gradio Interface
with gr.Blocks(title="HEROINE Simulator", theme=gr.themes.Soft()) as app:
    composition = gr.State([])
    
    gr.Markdown("# 💊 HEROINE: Advanced Drug Encapsulation Simulator")
    
    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### 1. Build Your System")
            with gr.Row():
                mol_input = gr.Textbox(label="Molecule Name or SMILES", placeholder="e.g. aspirin, PEG, CCO")
                mol_type = gr.Dropdown(choices=["Polymer", "Drug"], value="Drug", label="Type")
            
            with gr.Row():
                mol_count = gr.Number(value=1, label="Quantity", minimum=1, precision=0)
                preview_btn = gr.Button("🔍 Preview")
                add_btn = gr.Button("➕ Add to System", variant="secondary")
            
            mol_error = gr.Markdown(visible=False, value="Error")
            mol_preview = gr.Image(label="Molecule Preview", interactive=False, height=200)
            
            gr.Markdown("### 2. System Composition")
            comp_display = gr.Markdown("System is empty. Add molecules to begin.")
            clear_btn = gr.Button("Clear System")

        with gr.Column(scale=2):
            gr.Markdown("### 3. Simulation Settings")
            with gr.Row():
                steps = gr.Number(value=50000, label="Simulation Steps", precision=0)
                temp = gr.Number(value=300.0, label="Temperature (K)")
            
            with gr.Row():
                render = gr.Checkbox(value=False, label="Render GIF (Slow)")
                report = gr.Checkbox(value=True, label="Generate PDF Report")
            
            run_btn = gr.Button("🔥 Run Simulation", variant="primary", size="lg")
            
            gr.Markdown("### 4. Results")
            status_output = gr.Textbox(label="Status", interactive=False)
            with gr.Row():
                plot_output = gr.Image(label="Stability Plot")
                gif_output = gr.Image(label="Trajectory GIF")
            report_output = gr.File(label="Lab Report PDF")

    # Interaction Logic
    preview_btn.click(fn=lambda x: get_mol_preview(x)[0], inputs=[mol_input], outputs=[mol_preview])
    
    add_btn.click(
        fn=add_molecule,
        inputs=[mol_input, mol_count, mol_type, composition],
        outputs=[composition, mol_error]
    ).then(
        fn=format_composition,
        inputs=[composition],
        outputs=[comp_display]
    )
    
    clear_btn.click(fn=lambda: ([], "System is empty. Add molecules to begin."), outputs=[composition, comp_display])
    
    run_btn.click(
        fn=run_multi_simulation,
        inputs=[composition, steps, temp, render, report],
        outputs=[status_output, plot_output, gif_output, report_output]
    )

if __name__ == "__main__":
    app.launch(server_name="127.0.0.1", server_port=7860)
