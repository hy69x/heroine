# HEROINE Web UI

This directory contains the Gradio-based Web UI for the HEROINE Drug Encapsulation Simulator.

## Features

- **Interactive Molecule Search**: Input chemical names (e.g., "aspirin", "PEG") or valid SMILES strings directly. The app automatically fetches SMILES codes from PubChem for you.
- **Dynamic Configuration**: Adjust the number of polymer chains, drug molecules, simulation steps, and temperature via sliders and numeric inputs.
- **Real-Time Results**: The UI streams back the simulation results, displaying:
  - The final stability plot showing potential energy and temperature convergence.
  - An animated GIF of the encapsulation process (if Ovito is available and enabled).
  - A downloadable, dynamically generated PDF Lab Report.

## Installation

The UI depends on `gradio`, which is included in the project's dependencies. If you haven't already, install the dependencies using `uv`:

```bash
uv sync
```

## Running the UI

Start the Gradio server by executing:

```bash
uv run ui.py
```

After starting, it will provide a local URL (typically `http://127.0.0.1:7860`). Open this link in your web browser to access the graphical interface.
