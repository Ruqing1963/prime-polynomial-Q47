# Prime Polynomial Q47

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Large-scale computational study of prime values of Q(n) = n⁴⁷ - (n-1)⁴⁷**

## Key Results

| Metric | Value |
|--------|-------|
| Total primes found | **2,597,698** |
| Search range | n ≤ 3×10⁸ |
| Value size | 53 - 392 digits |
| Prime quadruplets | **3** (observed) vs **3.52** (Hardy-Littlewood prediction) |

### 🌟 Highlights

- **Perfect Hardy-Littlewood Match**: Observed 3 quadruplets vs predicted 3.52 (ratio 0.85)
- **Small-Prime Immunity Verified**: Zero primes in all 46 forbidden residue classes mod 283
- **Bateman-Horn Consistency**: Density decay from ~25,000/M to ~9,800/M across 6 orders of magnitude

## Prime Quadruplet Locations

| ID | Starting n | Status |
|----|------------|--------|
| S-1 | 117,309,848 | ✅ Verified |
| S-2 | 136,584,738 | ✅ Verified |
| S-3 | 218,787,064 | ✅ Verified |

## Repository Structure

```
prime-polynomial-Q47/
├── README.md
├── LICENSE
├── data/
│   ├── prime_n_values_full.txt      # All 2,597,698 prime-producing n values
│   ├── quadruplets.txt              # 4-tuple locations
│   ├── triples.txt                  # 3-tuple locations
│   └── density_by_range.csv         # Density statistics
├── paper/
│   ├── main.tex                     # LaTeX source
│   └── main.pdf                     # Compiled paper
├── figures/
│   ├── fig1_density.pdf             # Density evolution
│   ├── fig2_residue.pdf             # Residue analysis mod 283
│   ├── fig3_hardy_littlewood.pdf    # H-L prediction comparison
│   └── fig4_starmap.pdf             # Quadruplet star map
└── scripts/
    ├── analyze_primes.py            # Main analysis script
    ├── generate_figures.py          # Figure generation
    ├── verify_quadruplets.py        # Quadruplet verification
    └── residue_analysis.py          # Mod 283 analysis
```

## Quick Start

```bash
# Clone the repository
git clone https://github.com/Ruqing1963/prime-polynomial-Q47.git
cd prime-polynomial-Q47

# Install dependencies
pip install numpy matplotlib scipy

# Run analysis
python scripts/analyze_primes.py

# Generate figures
python scripts/generate_figures.py
```

## Mathematical Background

### The Polynomial

$$Q(n) = n^{47} - (n-1)^{47}$$

This is a degree-46 polynomial (the discrete derivative of n⁴⁷).

### Small-Prime Immunity Theorem

**Theorem:** For all primes p < 283 with p ≢ 1 (mod 47), Q(n) is never divisible by p (except when p | n(n-1)).

**Proof:** If Q(n) ≡ 0 (mod p) with p ∤ n(n-1), then (n/(n-1))⁴⁷ ≡ 1 (mod p). The order divides gcd(47, p-1) = 1 when p ≢ 1 (mod 47), contradiction.

### Hardy-Littlewood Prediction

The expected number of prime quadruplets is:

$$E[\text{quadruplets}] = C_{\text{sys}} \cdot \int_1^{N} \frac{dn}{(\log Q(n))^4} \approx 3.52$$

where C_sys ≈ 6,475 due to small-prime immunity.

## Data Format

### prime_n_values_full.txt
```
13
53
83
...
```
One n value per line, sorted in ascending order.

### density_by_range.csv
```csv
range_start,range_end,prime_count,density_per_million
13,5000,127,25466.2
5000,20000,302,20133.3
...
```

## Citation

```bibtex
@dataset{chen2026prime,
  author    = {Chen, Ruqing},
  title     = {Structure Formation in Arithmetic: Prime Values of n^47-(n-1)^47},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX},
  url       = {https://github.com/Ruqing1963/prime-polynomial-Q47}
}
```

## License

This project is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

- **Author:** Ruqing Chen
- **Email:** ruqing@hotmail.com
- **Affiliation:** GUT Geoservice Inc., Montreal, Canada

## References

1. Bateman, P.T. & Horn, R.A. (1962). "A heuristic asymptotic formula concerning the distribution of prime numbers." Math. Comp. 16, 363-367.

2. Hardy, G.H. & Littlewood, J.E. (1923). "Some problems of 'Partitio Numerorum' III." Acta Math. 44, 1-70.
