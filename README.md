# disease-via-community-assembly

This repository contains the code and data used in the analysis for the paper:

Community assembly reveals environmental controls over rodent competition drive deer mouse density and hantavirus infection.  
Angela D. Luis and Dean E. Pearson
Submitted to Ecology Letters, 2025.

---

## Repository Contents
- `Code.R` – R script used to perform the analyses.  

- `WebSpeciesData.csv` – Dataset of average species abundances (columns) per trapping web (rows). Two-letter codes corresponding to species names are provided in "SpeciesTraitData.csv" dataset and SI Appendix Table S1. Additionally contains number of months the web was trapped ("n_months"), average rodent diversity by inverse Simpson's D diversity index ("invSimpson"), average SNV prevalence among deer mice ("prevalence"), and whether or not the site was included in Luis et al. (2018) modeling study ("in2018pub").
  
- `EnvironmentalData.csv` – Average environmental conditions per web, including daily temperature minima ("tmin") and maxima ("tmax") (°C), daily total precipitation in mm/day (sum of all forms converted to water, "precip"), snow-water equivalent (km/m2) – a measurement of the amount of water contained within the snowpack ("swe"), and elevation. Additionally, from the Rangeland Analysis Platform (rangelands.app), new herbaceous above ground biomass over each 16-day period ("biomass") which is then summed to an annual total and averaged, average percent tree cover, and average percent bare ground.
    
- `SpeciesTraitData.csv` – Rodent species 2-letter codes and species trait data from EltonTraits 1.0 database (Wilman et al. 2014) including average adult body mass, whether nocturnal, and the proportion of their diet belonging to the following categories: invertebrates, endotherms, ectotherms, fish, vertebrates (general), scavenging, fruit, nectar, seed, and plant.

- `PhyloTree.txt` - Phylogenetic tree of rodent species in the dataset, taken from the Bininda-Emonds et al. (2007) mammalian supertree.

---

## Citation

If you use this code or data, please cite:

### Data and Code Archive
Angela D. Luis. Disease via Community Assembly. Zenodo. [DOI]

### Associated Paper 
Luis, A.D. and D.E. Pearson. 2025. Community assembly reveals environmental controls over rodent competition drive deer mouse density and hantavirus infection. Ecology Letters, in Review. 

---

## Acknowledgements

This research was supported by the National Science Foundation under EEID grant 2109828.
