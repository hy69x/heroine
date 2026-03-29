"""
Module V: Visualization
Uses Ovito to render the simulation trajectory into an animation.
"""

import sys
from pathlib import Path
import warnings

# Suppress OVITO PyPI warning as requested by the library
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

# Debugging Import
try:
    from ovito.io import import_file
    from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer
    from ovito.modifiers import AmbientOcclusionModifier, ColorCodingModifier, SliceModifier
    OVITO_AVAILABLE = True
    IMPORT_ERROR = None
except ImportError as e:
    OVITO_AVAILABLE = False
    IMPORT_ERROR = e

def render_trajectory(dump_path: Path, output_gif: Path, max_frames: int = 50, cut: bool = False):
    """
    Renders the trajectory.dump file into a GIF animation using Ovito.
    Includes element-specific radii and frame sampling for efficiency.
    """
    print(f"[*] Visualization Module Invoked for: {dump_path}")
    
    if not OVITO_AVAILABLE:
        print(f"Warning: Ovito python module not found. Skipping visualization.")
        return

    if not dump_path.exists():
        print(f"Error: Dump file not found at {dump_path}")
        return

    print(f"[*] Rendering visualization to {output_gif.name}...")

    try:
        # 1. Load Data
        pipeline = import_file(str(dump_path))
        pipeline.add_to_scene()
        
        num_frames = pipeline.source.num_frames
        print(f"[*] Loaded Trajectory: {num_frames} frames found.")

        if num_frames == 0:
            print("Warning: Trajectory has 0 frames. Nothing to render.")
            return

        # 2. Visual settings & Element Radii (VdW-ish)
        # 1:C, 2:H, 3:O, 4:N, 5:S
        radii_map = {1: 1.7, 2: 1.2, 3: 1.5, 4: 1.55, 5: 1.8}
        
        def setup_particles(frame, data):
            # Request mutable access to particle types using the '_' notation
            types = data.particles_.particle_types_
            for t_id, radius in radii_map.items():
                # Check if this type ID exists in the data
                try:
                    ptype = types.type_by_id(t_id)
                    ptype.radius = radius
                except (IndexError, KeyError, RuntimeError):
                    # Skip if type doesn't exist or isn't accessible
                    continue
        
        pipeline.modifiers.append(setup_particles)

        # Optional: Cross-section (Slice)
        if cut:
            print("[*] Applying SliceModifier for cross-section view.")
            # Slice along the Z-axis at the center
            pipeline.modifiers.append(SliceModifier(normal=(0, 0, 1), distance=0.0))

        # Ambient Occlusion for depth
        pipeline.modifiers.append(AmbientOcclusionModifier(intensity=0.4))
        
        # Color by Molecule ID to distinguish polymer from payload
        pipeline.modifiers.append(ColorCodingModifier(property='Molecule Identifier', start_value=1, end_value=10))

        # 3. Setup Camera/Viewport
        vp = Viewport()
        vp.type = Viewport.Type.Perspective
        vp.zoom_all()
        
        # 4. Frame Sampling Logic
        # If we have 500 frames but max_frames is 50, we render every 10th frame
        every_nth = max(1, num_frames // max_frames)
        render_frames = list(range(0, num_frames, every_nth))
        
        print(f"[*] Rendering {len(render_frames)} sampled frames (every {every_nth}th frame).")

        # 5. Render Animation
        render_args = {
            "filename": str(output_gif),
            "size": (800, 600),
            "range": (render_frames[0], render_frames[-1]),
            "every_nth": every_nth,
            "fps": 10
        }
        
        try:
            # Try GPU Renderer first
            print("[*] Attempting GPU (OpenGL) rendering...")
            vp.render_anim(renderer=OpenGLRenderer(), **render_args)
            print("[*] GPU Rendering finished.")
        except Exception as e:
            print(f"Warning: GPU Rendering failed ({e}). Falling back to CPU Raytracing...")
            # Fallback to CPU Tachyon
            vp.render_anim(renderer=TachyonRenderer(antialiasing_samples=1), **render_args)
            print("[*] CPU Rendering finished.")
        
        print(f"[*] Visualization saved successfully: {output_gif}")
        
    except Exception as e:
        print(f"Error rendering visualization: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    pass