# Logistic growth regression models applied to HIV transmitted founder and chronic control isolates

[![DOI](https://zenodo.org/badge/508808726.svg)](https://zenodo.org/badge/latestdoi/508808726)

Code to fit (and illustrate) the logistic growth models used in Sugrue *et al*. (2022). *The apparent interferon resistance of transmitted HIV-1 is possibly a consequence of enhanced replicative fitness*. PLOS Pathogens 18[11]: e1010973. [doi:10.1371/journal.ppat.1010973](https://doi.org/10.1371/journal.ppat.1010973).

## Requirements

- Install the [conda package manager](https://conda.io/)
- All other dependencies can then be installed using:

```
conda env create -f environment.yml
```

## Usage

Re-run all analyses using:

```
conda activate hiv_growth_rate
make
```
