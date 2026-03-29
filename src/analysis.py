"""
Module IV: Analytics
Parses simulation logs to calculate physical metrics (Radius of Gyration) and generates stability plots.
"""

from pathlib import Path
import sys

try:
    import matplotlib
    matplotlib.use('Agg') # Use non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
    PLOT_AVAILABLE = True
except ImportError:
    PLOT_AVAILABLE = False

def moving_average(data, window_size):
    if window_size <= 1:
        return data
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def analyze_results(log_path: Path, output_plot: Path = None, polymer_count: int = 1, payload_count: int = 0, dump_path: Path = None):
    """
    Parses the LAMMPS log file to find physical metrics and optionally generates a stability plot.
    Memory-efficient: reads line-by-line.
    """
    print(f"[*] Analyzing results from {log_path}")
    
    results = {"rg": 0.0, "efficiency": None}

    if not log_path.exists():
        print(f"Error: Log file {log_path} not found.")
        return results
        
    print(f"[*] Parsing log file: {log_path.name}")

    steps, rg_values, temp_values, energy_values = [], [], [], []
    
    try:
        # Step 1: Find the start of the LAST run section and its headers
        last_run_start_byte = 0
        headers = []
        col_map = {}

        with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
            while True:
                line = f.readline()
                if not line: break
                if line.strip().startswith("Step") and "c_rg" in line:
                    last_run_start_byte = f.tell() - len(line)
                    headers = line.strip().split()
                    col_map = {h: idx for idx, h in enumerate(headers)}
        
        if not col_map:
             print("Error: Could not find thermo output with c_rg in log.")
             return results

        rg_col = col_map.get("c_rg")
        temp_col = col_map.get("Temp")
        energy_col = col_map.get("PotEng") or col_map.get("epair") or col_map.get("E_pair")

        # Step 2: Parse data from the last run section only
        with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
            f.seek(last_run_start_byte)
            f.readline() # Skip header
            
            for line in f:
                tokens = line.strip().split()
                if not tokens: continue
                
                # Check if first token is an integer (Step)
                if tokens[0].isdigit():
                    try:
                        steps.append(int(tokens[0]))
                        rg_values.append(float(tokens[rg_col]))
                        temp_values.append(float(tokens[temp_col]))
                        energy_values.append(float(tokens[energy_col]))
                    except (ValueError, IndexError):
                        pass
                elif "Loop time" in line:
                    break
                
        if not rg_values:
            print("Warning: No data extraction.")
        else:
            results["rg"] = rg_values[-1]
            print(f"[*] Final Radius of Gyration (Rg): {results['rg']:.4f} Angstroms")

        # Generate stability plot if requested and possible
        if output_plot and steps and PLOT_AVAILABLE:
            print(f"[*] Generating stability plot: {output_plot.name}")
            window = max(1, len(steps) // 50)
            s_steps = steps[window-1:]
            s_temp = moving_average(temp_values, window)
            s_energy = moving_average(energy_values, window)

            fig, ax1 = plt.subplots(figsize=(10, 6))
            color = 'tab:red'
            ax1.set_xlabel('Step')
            ax1.set_ylabel('Temperature (K)', color=color)
            ax1.plot(steps, temp_values, color=color, alpha=0.2)
            ax1.plot(s_steps, s_temp, color=color, linewidth=2, label='Temp (MA)')
            ax1.tick_params(axis='y', labelcolor=color)

            ax2 = ax1.twinx()
            color = 'tab:blue'
            ax2.set_ylabel('Potential Energy (kcal/mol)', color=color)
            ax2.plot(steps, energy_values, color=color, alpha=0.2)
            ax2.plot(s_steps, s_energy, color=color, linewidth=2, label='Energy (MA)')
            ax2.tick_params(axis='y', labelcolor=color)

            plt.title(f'System Equilibrium Stability (Smoothing Window: {window})')
            fig.tight_layout()
            plt.savefig(output_plot)
            plt.close()

        # Phase 3: Encapsulation Efficiency
        if payload_count > 0 and dump_path and dump_path.exists():
            print(f"[*] Calculating Encapsulation Efficiency...")
            eff = calculate_encapsulation_efficiency(dump_path, polymer_count, payload_count, results["rg"])
            results["efficiency"] = eff
            print(f"[*] Encapsulation Efficiency: {eff:.2f}%")

        return results

    except Exception as e:
        print(f"Error analyzing log file: {e}")
        import traceback
        traceback.print_exc()
        return results

def calculate_encapsulation_efficiency(dump_path: Path, polymer_count: int, payload_count: int, rg: float):
    """
    Parses the LAST frame of the trajectory dump without reading the whole file.
    Dynamically calculates buffer size based on the number of atoms in the system.
    """
    try:
        filesize = dump_path.stat().st_size
        if filesize == 0: return 0.0

        # 1. Read the first ~2KB to find the total number of atoms
        num_atoms = 0
        with open(dump_path, 'r', encoding='utf-8', errors='replace') as f:
            for _ in range(20): # Atom count is usually in the first 10 lines
                line = f.readline()
                if not line: break
                if "ITEM: NUMBER OF ATOMS" in line:
                    num_atoms = int(f.readline().strip())
                    break
        
        if num_atoms == 0:
            print("[!] Warning: Could not determine atom count from dump header.")
            num_atoms = 10000 # Fallback
            
        # 2. Estimate buffer size: ~100 bytes per atom line + 500 bytes header
        # We read 1.5x the estimated frame size to ensure we catch the last 'ITEM: ATOMS'
        estimated_frame_size = num_atoms * 100 + 500
        buffer_size = min(filesize, int(estimated_frame_size * 1.5))
        
        with open(dump_path, 'rb') as f:
            f.seek(filesize - buffer_size)
            chunk = f.read(buffer_size).decode('utf-8', errors='replace')
            
            lines = chunk.splitlines()
            last_atoms_idx = -1
            # Search backwards for the LAST frame start
            for i in range(len(lines)-1, -1, -1):
                if "ITEM: ATOMS" in lines[i]:
                    last_atoms_idx = i
                    break
            
            if last_atoms_idx == -1:
                print(f"[!] Warning: Could not find last frame in trailing {buffer_size/1024:.1f}KB.")
                return 0.0
            
            header = lines[last_atoms_idx].strip().split()
            # ITEM: ATOMS id mol type x y z ...
            col_map = {h: idx-2 for idx, h in enumerate(header)}
            
            mol_col = col_map.get("mol")
            x_col = col_map.get("x")
            y_col = col_map.get("y")
            z_col = col_map.get("z")
            
            polymer_atoms = []
            payload_mols = {m: [] for m in range(polymer_count + 1, polymer_count + payload_count + 1)}
            
            for line in lines[last_atoms_idx+1:]:
                tokens = line.strip().split()
                if not tokens or len(tokens) < 5: continue
                if "ITEM:" in tokens[0]: break # Next frame or section
                
                try:
                    mol_id = int(tokens[mol_col])
                    pos = [float(tokens[x_col]), float(tokens[y_col]), float(tokens[z_col])]
                    
                    if mol_id <= polymer_count:
                        polymer_atoms.append(pos)
                    elif mol_id in payload_mols:
                        payload_mols[mol_id].append(pos)
                except (ValueError, IndexError):
                    continue
                    
            if not polymer_atoms:
                print("[!] Warning: No polymer atoms found in last frame.")
                return 0.0
                
            p_atoms_np = np.array(polymer_atoms)
            polymer_com = np.mean(p_atoms_np, axis=0)
            
            encapsulated_count = 0
            threshold = 1.5 * rg
            
            for mol_id, atoms in payload_mols.items():
                if not atoms: continue
                payload_com = np.mean(np.array(atoms), axis=0)
                dist = np.linalg.norm(payload_com - polymer_com)
                if dist <= threshold:
                    encapsulated_count += 1
                    
            return (encapsulated_count / payload_count) * 100.0 if payload_count > 0 else 0.0
            
    except Exception as e:
        print(f"Error calculating efficiency: {e}")
        return 0.0


if __name__ == "__main__":
    pass