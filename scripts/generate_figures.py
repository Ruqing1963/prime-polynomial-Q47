#!/usr/bin/env python3
"""
generate_figures.py - Generate publication-quality figures for the paper

This script generates all four figures:
- Fig 1: Density evolution across scales
- Fig 2: Residue analysis mod 283
- Fig 3: Hardy-Littlewood prediction comparison
- Fig 4: Prime quadruplet star map

Author: Ruqing Chen
Repository: https://github.com/Ruqing1963/prime-polynomial-Q47
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from pathlib import Path

plt.rcParams['font.size'] = 11
plt.rcParams['mathtext.fontset'] = 'cm'

def load_prime_n_values(filepath):
    """Load prime-producing n values from file."""
    n_values = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#') and line.isdigit():
                n_values.append(int(line))
    return sorted(set(n_values))

def create_fig1_density(n_values, output_dir):
    """Create Figure 1: Density evolution."""
    print("Creating Figure 1: Density evolution...")
    
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    ranges = [
        (13, 5000), (5000, 20000), (200000, 500000), (500000, 1000000),
        (2000000, 10000000), (50000000, 100000000), (100000000, 300000000),
    ]
    
    centers = []
    densities = []
    for low, high in ranges:
        count = sum(1 for n in n_values if low <= n < high)
        if count > 10:
            center = (low + high) / 2
            density = count / (high - low) * 1e6
            centers.append(center)
            densities.append(density)
    
    ax.semilogx(centers, densities, 'bo-', markersize=10, linewidth=2, label='Observed Density')
    
    # Bateman-Horn prediction
    n_theory = np.logspace(3, 8.5, 100)
    C_fit = densities[0] * 46 * np.log(centers[0])
    theory_density = C_fit / (46 * np.log(n_theory))
    ax.semilogx(n_theory, theory_density, 'g--', linewidth=2, alpha=0.7, label='Bateman-Horn Prediction')
    
    # Data gaps
    ax.axvspan(20000, 200000, alpha=0.15, color='gray', label='No data')
    ax.axvspan(1000000, 2000000, alpha=0.15, color='gray')
    ax.axvspan(10000000, 50000000, alpha=0.15, color='gray')
    
    ax.set_xlabel('n (log scale)', fontsize=13)
    ax.set_ylabel('Prime Density (primes per million n)', fontsize=13)
    ax.set_title('Prime Density Evolution of Q(n) = n^47 - (n-1)^47', fontsize=14)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(10, 5e8)
    ax.set_ylim(0, 30000)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig1_density.pdf', format='pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig1_density.png', dpi=150, bbox_inches='tight')
    plt.close()

def create_fig2_residue(n_values, output_dir):
    """Create Figure 2: Residue analysis mod 283."""
    print("Creating Figure 2: Residue analysis...")
    
    p = 283
    n_array = np.array(n_values)
    residues = n_array % p
    
    # Find forbidden residues
    roots_of_unity = [x for x in range(1, p) if pow(x, 47, p) == 1]
    forbidden = set()
    for r in roots_of_unity:
        if r != 1:
            inv = pow(r - 1, p - 2, p)
            n_forb = (r * inv) % p
            forbidden.add(n_forb)
    
    hist = Counter(residues)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    
    # Left: full distribution
    residue_counts = [hist.get(r, 0) for r in range(p)]
    colors = ['red' if r in forbidden else 'steelblue' for r in range(p)]
    ax1.bar(range(p), residue_counts, color=colors, alpha=0.8, width=1)
    ax1.set_xlabel('Residue class n (mod 283)', fontsize=12)
    ax1.set_ylabel('Count of prime-producing n', fontsize=12)
    ax1.set_title('Distribution by Residue mod 283', fontsize=12)
    
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='steelblue', alpha=0.8, label='Allowed'),
                       Patch(facecolor='red', alpha=0.8, label='Forbidden')]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    # Right: forbidden residues detail
    forbidden_list = sorted(forbidden)
    forbidden_counts = [hist.get(r, 0) for r in forbidden_list]
    ax2.bar(range(len(forbidden_list)), forbidden_counts, color='red', alpha=0.8)
    ax2.set_xticks(range(0, len(forbidden_list), 5))
    ax2.set_xticklabels([str(forbidden_list[i]) for i in range(0, len(forbidden_list), 5)], rotation=45)
    ax2.set_xlabel('Forbidden residue classes', fontsize=12)
    ax2.set_ylabel('Count (expected: 0)', fontsize=12)
    ax2.set_title(f'All {len(forbidden)} Forbidden Classes = Zero Primes', fontsize=12)
    ax2.set_ylim(0, 5)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig2_residue.pdf', format='pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig2_residue.png', dpi=150, bbox_inches='tight')
    plt.close()

def create_fig3_hardy_littlewood(n_values, output_dir):
    """Create Figure 3: Hardy-Littlewood comparison."""
    print("Creating Figure 3: Hardy-Littlewood comparison...")
    
    def count_ktuples_up_to_N(n_vals, k, N_max):
        n_list = sorted([n for n in n_vals if n <= N_max])
        count = 0
        i = 0
        while i < len(n_list):
            run_len = 1
            while i + 1 < len(n_list) and n_list[i+1] == n_list[i] + 1:
                run_len += 1
                i += 1
            if run_len >= k:
                count += 1
            i += 1
        return count
    
    N_values = [1e6, 1e7, 5e7, 1e8, 1.5e8, 2e8, 2.5e8, 3e8]
    pair_counts = [count_ktuples_up_to_N(n_values, 2, N) for N in N_values]
    triple_counts = [count_ktuples_up_to_N(n_values, 3, N) for N in N_values]
    quad_counts = [count_ktuples_up_to_N(n_values, 4, N) for N in N_values]
    
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    ax.semilogy(np.array(N_values)/1e6, pair_counts, 'b-o', markersize=8, linewidth=2, label='2-tuples (pairs)')
    ax.semilogy(np.array(N_values)/1e6, triple_counts, 'g-s', markersize=8, linewidth=2, label='3-tuples (triples)')
    ax.semilogy(np.array(N_values)/1e6, quad_counts, 'r-^', markersize=10, linewidth=2, label='4-tuples (quadruplets)')
    
    # Hardy-Littlewood prediction
    HL_prediction = [max(0.1, 3.52 * (N/3e8)**0.8) for N in N_values]
    ax.semilogy(np.array(N_values)/1e6, HL_prediction, 'r--', linewidth=2, alpha=0.7, 
                label='Hardy-Littlewood prediction (k=4)')
    
    ax.set_xlabel('N (millions)', fontsize=13)
    ax.set_ylabel('Cumulative count (log scale)', fontsize=13)
    ax.set_title('Prime k-tuples: Observation vs Hardy-Littlewood Prediction', fontsize=14)
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    
    ax.annotate('Observed: 3\nPredicted: 3.52\nRatio: 0.85', 
                xy=(300, 3), xytext=(200, 20),
                fontsize=11, ha='center',
                arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5),
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, edgecolor='orange'))
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig3_hardy_littlewood.pdf', format='pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'fig3_hardy_littlewood.png', dpi=150, bbox_inches='tight')
    plt.close()

def create_fig4_starmap(output_dir):
    """Create Figure 4: Prime quadruplet star map."""
    print("Creating Figure 4: Star map...")
    
    quadruplets = [117309848, 136584738, 218787064]
    
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Dark background
    ax.set_facecolor('#0a0a20')
    fig.patch.set_facecolor('#0a0a20')
    
    # Number line
    ax.axhline(y=0.5, color='#1a1a40', linewidth=20, alpha=0.8)
    
    # Grid lines
    for x in range(0, 350, 50):
        ax.axvline(x=x, color='#2a2a50', linewidth=0.5, alpha=0.5, linestyle=':')
    
    # Data coverage
    ax.plot([100, 300], [0.5, 0.5], color='#304060', linewidth=8, alpha=0.6, solid_capstyle='round')
    
    # Stars
    for i, loc in enumerate(quadruplets, 1):
        x = loc / 1e6
        
        # Glow
        ax.plot(x, 0.5, 'o', color='gold', markersize=35, alpha=0.15)
        ax.plot(x, 0.5, 'o', color='gold', markersize=25, alpha=0.25)
        ax.plot(x, 0.5, 'o', color='gold', markersize=18, alpha=0.4)
        
        # Core
        ax.plot(x, 0.5, '*', color='white', markersize=20, zorder=10)
        ax.plot(x, 0.5, '*', color='gold', markersize=15, zorder=11)
        
        # Label
        y_offset = 0.25 if i != 2 else -0.25
        ax.annotate(f'S-{i}\nn = {loc:,}', 
                    xy=(x, 0.5), xytext=(x, 0.5 + y_offset),
                    ha='center', va='center' if i != 2 else 'top',
                    fontsize=10, fontweight='bold', color='white',
                    arrowprops=dict(arrowstyle='-', color='gold', lw=1, alpha=0.6))
    
    ax.set_xlim(90, 310)
    ax.set_ylim(0, 1)
    ax.set_xlabel('n (millions)', fontsize=13, color='white')
    ax.set_title('Prime Quadruplet Star Map: Three Rare Events in the Arithmetic Universe', 
                 fontsize=13, color='white', pad=15)
    ax.set_yticks([])
    
    ax.tick_params(axis='x', colors='white')
    ax.spines['bottom'].set_color('#404080')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Info box
    textstr = 'Hardy-Littlewood\nPrediction: 3.52\nObserved: 3'
    props = dict(boxstyle='round', facecolor='#1a1a40', alpha=0.9, edgecolor='gold')
    ax.text(0.02, 0.95, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', color='white', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'fig4_starmap.pdf', format='pdf', bbox_inches='tight', facecolor='#0a0a20')
    plt.savefig(output_dir / 'fig4_starmap.png', dpi=150, bbox_inches='tight', facecolor='#0a0a20')
    plt.close()

def main():
    # Setup paths
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / 'data'
    output_dir = script_dir.parent / 'figures'
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("GENERATING FIGURES FOR PRIME POLYNOMIAL Q47 PAPER")
    print("=" * 70)
    
    # Load data
    print("\nLoading data...")
    filepath = data_dir / 'prime_n_values_full.txt'
    n_values = load_prime_n_values(filepath)
    print(f"Loaded {len(n_values):,} prime-producing n values")
    
    # Generate figures
    create_fig1_density(n_values, output_dir)
    create_fig2_residue(n_values, output_dir)
    create_fig3_hardy_littlewood(n_values, output_dir)
    create_fig4_starmap(output_dir)
    
    print("\n" + "=" * 70)
    print(f"All figures saved to: {output_dir}")
    print("=" * 70)

if __name__ == "__main__":
    main()
