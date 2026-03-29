"""
Module VI: Automated Lab Notebook
Generates a professional PDF report summarizing the simulation results.
"""

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from datetime import datetime

def generate_2d_molecule(smiles: str, output_path: Path):
    """Generates a clean 2D representation of the molecule."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Draw.MolToFile(mol, str(output_path), size=(400, 400))
            return True
    except Exception as e:
        print(f"Error generating 2D molecule image: {e}")
    return False

def generate_report(
    output_dir: Path,
    name: str,
    smiles: str,
    steps: int,
    temp: float,
    rg: float,
    efficiency: float = None,
    plot_path: Path = None
):
    """
    Creates a professional PDF report in the output directory.
    """
    report_path = output_dir / "lab_report.pdf"
    mol_img_path = output_dir / "structure_2d.png"
    
    print(f"[*] Generating Lab Report: {report_path.name}")
    
    # 1. Generate 2D Image
    generate_2d_molecule(smiles, mol_img_path)
    
    # 2. Setup PDF Document
    doc = SimpleDocTemplate(str(report_path), pagesize=A4)
    styles = getSampleStyleSheet()
    
    # Custom Styles
    title_style = ParagraphStyle(
        'TitleStyle',
        parent=styles['Heading1'],
        fontSize=24,
        spaceAfter=20,
        textColor=colors.HexColor("#2C3E50")
    )
    
    header_style = ParagraphStyle(
        'HeaderStyle',
        parent=styles['Heading2'],
        fontSize=14,
        spaceBefore=10,
        spaceAfter=10,
        textColor=colors.HexColor("#2980B9")
    )

    elements = []
    
    # Title & Metadata
    elements.append(Paragraph(f"Proteus Simulation Report: {name}", title_style))
    elements.append(Paragraph(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles['Normal']))
    elements.append(Spacer(1, 20))
    
    # Section: Molecular Structure
    elements.append(Paragraph("1. Molecular Composition", header_style))
    elements.append(Paragraph(f"<b>SMILES String:</b> {smiles}", styles['Normal']))
    
    if mol_img_path.exists():
        img = Image(str(mol_img_path), width=200, height=200)
        img.hAlign = 'CENTER'
        elements.append(img)
        elements.append(Paragraph("<para align='center'><i>2D Representation of the Monomer/System</i></para>", styles['Normal']))

    elements.append(Spacer(1, 20))
    
    # Section: Simulation Parameters
    elements.append(Paragraph("2. Simulation Environment", header_style))
    data = [
        ["Parameter", "Value"],
        ["Total Steps", f"{steps:,}"],
        ["Temperature", f"{temp} K"],
        ["Integrator", "Langevin Dynamics (NVE + Fix)"],
        ["Force Field", "Lennard-Jones (Generic CHONS)"],
    ]
    t = Table(data, colWidths=[150, 200])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#ECF0F1")),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.HexColor("#2C3E50")),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
    ]))
    elements.append(t)
    elements.append(Spacer(1, 20))

    # Section: Physical Metrics
    elements.append(Paragraph("3. Physical Analytics", header_style))
    metrics_data = [
        ["Metric", "Result"],
        ["Final Radius of Gyration (Rg)", f"{rg:.4f} \u00c5"]
    ]
    if efficiency is not None:
        metrics_data.append(["Encapsulation Efficiency", f"{efficiency:.2f} %"])
        
    t_metrics = Table(metrics_data, colWidths=[200, 150])
    t_metrics.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#D5F5E3")),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
    ]))
    elements.append(t_metrics)
    elements.append(Spacer(1, 20))

    # Section: Stability Plots
    if plot_path and plot_path.exists():
        elements.append(Paragraph("4. System Equilibrium Stability", header_style))
        plot_img = Image(str(plot_path), width=450, height=270)
        plot_img.hAlign = 'CENTER'
        elements.append(plot_img)
        elements.append(Paragraph("<para align='center'><i>Thermodynamic Stability: Temperature & Potential Energy</i></para>", styles['Normal']))

    # Build PDF
    doc.build(elements)
    print(f"[*] Report saved successfully: {report_path}")

if __name__ == "__main__":
    pass
