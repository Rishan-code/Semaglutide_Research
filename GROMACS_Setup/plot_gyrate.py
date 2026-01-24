import matplotlib.pyplot as plt
import numpy as np
import os

def read_xvg(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            # expected format: time, Rg, Rg_x, Rg_y, Rg_z
            data.append([float(x) for x in parts])
    return np.array(data)

def plot_rg(xvg_file, output_image, title='Radius of Gyration'):
    if not os.path.exists(xvg_file):
        print(f"Error: File {xvg_file} not found.")
        return

    data = read_xvg(xvg_file)
    time = data[:, 0]  # Time (ps)
    rg = data[:, 1]    # Radius of Gyration (nm)

    plt.figure(figsize=(10, 6))
    plt.plot(time, rg, label='Rg', color='blue', linewidth=1.5)
    
    plt.title(title, fontsize=16)
    plt.xlabel('Time (ps)', fontsize=14)
    plt.ylabel('Rg (nm)', fontsize=14)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_image, dpi=300)
    print(f"Plot saved to {output_image}")
    plt.close()

if __name__ == "__main__":
    # Standard Rg
    plot_rg(os.path.join("analysis", "gyrate.xvg"), 
            os.path.join("analysis", "gyrate_plot.png"),
            "Radius of Gyration (Aggregate)")
            
    # Monomer Rg (if exists)
    monomer_xvg = os.path.join("analysis", "gyrate_monomer.xvg")
    if os.path.exists(monomer_xvg):
        plot_rg(monomer_xvg,
                os.path.join("analysis", "gyrate_monomer_plot.png"),
                "Radius of Gyration (Single Monomer)")
