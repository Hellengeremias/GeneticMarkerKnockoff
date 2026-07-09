# GeneticMarkerKnockoff
Code and analysis scripts for evaluating knockoff generators in high-dimensional genetic marker selection.

## Knockoff Filter for Genetic Marker Selection
This repository contains the R code used in the manuscript:

### “Knockoff filter for genetic marker selection: performance under different knockoff generators”

The study evaluates knockoff-based variable selection for high-dimensional genetic marker data. Three knockoff generators are investigated:

1. a generator based on the Haldane mapping function (Haldane-k);
2. a generator based on empirical Markov transition probabilities (EM-k); and
3. a generator based on classification trees (Tree-k).

The proposed methods are compared with existing knockoff generation approaches under simulation scenarios with different dependence structures among genetic markers. The knockoff filter is also applied to genomic data from the ISA-Nutrição 2015 study.

## Repository contents

simulation_study.R

This script contains the code used in the simulation study, including response simulation, knockoff generation, computation of LASSO-based variable importance statistics, variable selection using the knockoff filter, and evaluation of selection performance.

real_data_analysis_HMM.R

This script contains the analysis of the ISA-Nutrição 2015 genomic data using knockoffs generated from a hidden Markov model.

real_data_analysis_other_generators.R

This script contains the analysis of the ISA-Nutrição 2015 genomic data using the other knockoff generators evaluated in the study, including Model-X knockoffs and the proposed generators (Haldane-k, EM-k and Tree-k).

## Citation

If you use the code or methods provided in this repository, please cite the associated manuscript:

Geremias, H. G. G. S. and Zuanetti, D. A. Knockoff filter for genetic marker selection: performance under different knockoff generators. Manuscript submitted for publication.

The citation information will be updated after publication.

## Data availability

The genotype data used in the simulations and real-data application are available under controlled access from the European Genome-phenome Archive (EGA) under accession number EGAD00010002678.

## Software

The analyses were performed in R. The required R packages and additional software dependencies are specified within the corresponding scripts.

## Contact

For questions related to the code or analyses, please open an issue in this repository or contact the corresponding author.
