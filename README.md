# GeneticMarkerKnockoff
Code and analysis scripts for evaluating knockoff generators in high-dimensional genetic marker selection.

## Knockoff Filter for Genetic Marker Selection
This repository contains the R code used in the manuscript:

### “Knockoff filter for genetic marker selection: performance under different knockoff generators”

The study evaluates knockoff-based variable selection for high-dimensional genetic marker data. Three knockoff generators are investigated:

1. a generator based on the Haldane mapping function;
2. a generator based on empirical Markov transition probabilities; and
3. a generator based on classification trees.

The proposed methods are compared with existing knockoff generation approaches under simulation scenarios with different dependence structures among genetic markers. The knockoff filter is also applied to genomic data from the ISA-Nutrição 2015 study.

## Knockoff generators

The repository includes implementations of three knockoff generators for genetic marker data.

### Haldane-based generator (Haldane-k)

This approach uses genetic distances and the Haldane mapping function to characterize dependence between adjacent genetic markers and generate knockoff variables.

### Empirical Markov generator (EM-k)

This approach estimates transition probabilities directly from the observed genotype data and generates knockoffs using an empirical Markov dependence structure.

### Classification-tree generator (Tree-k)

This approach models the conditional distribution of each genetic marker using classification trees and generates knockoff variables from the estimated conditional probabilities.

## Variable importance statistics

For the knockoff analyses, the original and knockoff variables are combined into an augmented design matrix. Variable importance statistics are obtained using LASSO regression, and knockoff statistics are constructed by contrasting the importance assigned to each original variable with that assigned to its knockoff counterpart.

## Simulation studies

The simulation scripts evaluate the knockoff generators under scenarios with different dependence structures among genetic markers.

The main simulation workflow consists of:

1. loading the genotype matrix;
2. simulating the response variable;
3. generating knockoff variables;
4. computing variable importance statistics;
5. applying the knockoff threshold;
6. evaluating variable selection performance; and
7. summarizing the results across simulation replicates.

Performance measures include False Discovery Rate (FDR), power,  F1 score, and the Matthews correlation coefficient.

## Real data analysis

The real-data application uses genomic data from the ISA-Nutrição 2015 study.

Two outcomes are considered:

- the logarithm of body mass index (BMI), as a quantitative outcome; and obesity status defined according to BMI, as a binary outcome.

- The genotype data are available under controlled access from the European Genome-phenome Archive (EGA) under accession number: EGAD00010002678

Because the genotype data are subject to controlled-access restrictions, they are not distributed through this repository.

## Reproducibility

The analyses were performed in R. Before reproducing the analyses, users should install the required packages and external software described in the corresponding analysis directories.

## Citation

If you use the code or methods provided in this repository, please cite the associated manuscript:

Geremias, H. G. G. S. and Zuanetti, D. A. Knockoff filter for genetic marker selection: performance under different knockoff generators. Manuscript submitted for publication.

The citation information will be updated after publication.

## Data availability

Analysis scripts are available in this repository.

The genotype data used in the real-data application are available under controlled access from the European Genome-phenome Archive (EGA) under accession number EGAD00010002678.

## Contact

For questions related to the code or analyses, please open an issue in this repository or contact the corresponding author.
