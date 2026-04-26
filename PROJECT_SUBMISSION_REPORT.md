# HEROINE: Project Submission Report
**Project Title:** HEROINE (High-Efficiency Resolution of Integrated Nanoparticle Encapsulation)  
**Date:** April 2026  
**Subject:** Computational Molecular Modeling / Software Engineering  

---

## 1. Project Overview
HEROINE is an automated pipeline for the simulation and analysis of drug-polymer encapsulation. The project aims to provide a high-fidelity, user-friendly tool for predicting how effectively a specific polymer chain can trap drug molecules within a nanoparticle structure.

## 2. Methodology
The software follows a five-stage pipeline:
1. **Input Resolution**: Resolving names/SMILES via PubChem API.
2. **Topological Generation**: Converting chemical graphs to 3D physical coordinates and assigning force-field parameters.
3. **Molecular Dynamics (MD)**: Simulating the system using LAMMPS with a Langevin thermostat.
4. **Analytics**: Calculating the Radius of Gyration ($R_g$) and Encapsulation Efficiency.
5. **Visualization & Reporting**: Rendering 3D trajectories and generating a consolidated PDF lab report.

## 3. Technical Implementation
- **Backend**: Python 3.12, RDKit, NumPy.
- **Simulation Engine**: LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator).
- **Frontend**: Gradio Web UI.
- **Reporting**: ReportLab for PDF generation and Matplotlib for thermodynamic stability plotting.

## 4. Key Improvements & Innovations
- **Multi-Drug Support**: Ability to simulate complex mixtures of different drug types simultaneously.
- **Proximity-Based Analytics**: A robust algorithm that checks for drug-polymer contact across the entire polymer surface rather than a simple center-of-mass heuristic.
- **Error Diagnostics**: Automated parsing of simulation logs to detect and report physical instabilities.
- **Checkpointing**: Intelligent skipping of completed simulation phases to optimize resource usage.

## 5. Conclusion
HEROINE demonstrates a successful integration of high-performance physics engines with modern web technologies, providing a bridge between theoretical chemistry and practical pharmaceutical application.
