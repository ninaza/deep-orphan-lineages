# DEEPOL analysis

Analysis code for the manuscript “Three novel deep-branching protist lineages from the marine rare biosphere” (Pohl et al., in preparation).

## Project overview

This repository contains the code behind the identification and biogeography of three deep-branching marine protist lineages, Deep Orphan Lineages 1 to 3 (DEEPOL1 to DEEPOL3). The lineages were found among orphan sequences from long-read metabarcoding of marine environmental DNA. Their ecological niches, community structure and phylogenetic specialization across the global ocean were then studied using EukBank short-read data.

The workflow is modular, so each lineage can be configured and run independently.

## Repository structure

```
DEEPOL-analysis/
├── README.md
├── 01_clade_extraction/
│   └── extract_orphan_clades.R
├── 02_biogeography/
│   ├── DEEPOLs_biogeography.Rmd          main analysis notebook
│   ├── data_preparation.R                prepares metadata and per-lineage abundance tables
│   ├── sequencing_depth_boxplots.R       sequencing depth per sample group
│   ├── size_fraction_upset_plots.R       ASV overlap across size fractions
│   ├── eukbank_percent_identity.R        mean percent identity of EukBank hits
│   └── scripts/                          helper functions (niche traits, reporting)
├── data/                                 input data (see "Data" below)
└── results/                              generated automatically, one folder per lineage
```

File names in the tree are illustrative; they should match the files in the repository.

## 01 Clade extraction

The sequence similarity network (SSN) analysis, which identifies maximal cliques of orphan sequences, was performed with script from Miguel Méndez Sandín and can be found under (https://github.com/MiguelMSandin/SSNetworks). Its output (the clique assignments) is the input to this step.

`extract_orphan_clades.R` identifies candidate deep-branching orphan clades in a phylogeny and combines them with these SSN clique assignments.

- An orphan clade is a monophyletic group in which every tip is an orphan sequence and which contains at least two sequences. By default only the largest (maximal) such clades are reported.
- Inputs: orphan sequence IDs, the clique assignments from the SSN, and a rooted tree in Newick format.
- Output: a table with the SSN clique and the tree clade of each orphan sequence.

The settings (file paths, minimum clade size, treatment of nested clades) are at the top of the script, and its header describes the input and output formats.

## 02 Biogeography

The main notebook (`DEEPOLs_biogeography.Rmd`) is organised into the following modules:

- **Data integration:** merging ASV presence/absence matrices with global environmental metadata and size fraction data. Ocean basin and polar region are assigned from sample coordinates using IHO sea area polygons (see "Data").
- **Coastal vs open-ocean affinity:** distance to the nearest coastline for every sample, followed by
  - logistic regression of detection and linear regression of relative abundance on distance to coast,
  - PERMANOVA and PERMDISP of ASV composition between coastal and open-ocean samples,
  - indicator ASV analysis (IndVal) for coastal and open-ocean waters,
  - confound and robustness checks (association with ocean layer, polar region, size fraction and study, leave-one-study-out, stratified models, interaction tests),
  - a forest plot of the effects of distance to the coast.
- **Community structure:**
  - UpSet plots: ASV intersections across size fractions.
  - PCoA (Jaccard): community similarity with environmental vector fitting (envfit).
- **Drivers of distribution:**
  - dbRDA: distance-based redundancy analysis (Jaccard) of environmental variables, conditioned on spatial coordinates, with total and adjusted R².
  - Variation partitioning: relative influence of environment and geography (Jaccard distances).
- **Statistical significance:**
  - PERMANOVA: differences between ocean basins and layers.
  - PERMDISP: multivariate dispersion, to check that community shifts are not caused by differences in spread.
- **Phylogenetic niche traits:** habitat specialization (polarity, ocean layer), coastal vs open-ocean indicator ASVs and niche breadth (occupancy) mapped onto lineage-specific phylogenies with `ggtree`.
- **Summary:** a master table with the main statistics for each lineage.

### Additional scripts

- `data_preparation.R`: prepares the metadata (ocean basin, polar group, sequencing depth, size fraction) and the per-lineage abundance tables used by the notebook. Run it before knitting the notebook.
- `sequencing_depth_boxplots.R`: sequencing depth across sample groups.
- `size_fraction_upset_plots.R`: UpSet plots of ASV overlap across size fractions.
- `eukbank_percent_identity.R`: mean percent identity of the EukBank hits to the DEEPOL reference sequences.

## Data

The `data/` folder is expected to contain:

- lineage-specific filtered abundance matrices and the global metadata table,
- the orphan sequence IDs, SSN clique assignments and the tree used in `01_clade_extraction`,
- the IHO Sea Areas polygons (`World_Seas_IHO_v3.shp`, Flanders Marine Institute) used for the ocean assignment.

The coastline used for the distance to coast is obtained through the `rnaturalearth` package.

Public data used:

EukBank (ASV tables and metadata): https://zenodo.org/records/7804946
Reference tree for the clade extraction: Jamy, M. et al. (2022). Data for "Global patterns and rates of habitat transitions across the eukaryotic tree of life". figshare. https://doi.org/10.6084/m9.figshare.15164772.v3

## Running the analysis

1. Place the input data in `data/` and adjust the paths in the settings at the top of each script.
2. Run the SSN analysis with the external code (see "01 Clade extraction"), then `01_clade_extraction/extract_orphan_clades.R`.
3. Run `02_biogeography/data_preparation.R`.
4. Knit `02_biogeography/DEEPOLs_biogeography.Rmd`. Outputs (PDFs and CSVs) are written to `results/`, in one folder per lineage.

## Requirements

The analysis was performed in R (version [x.y.z]). The SSN analysis (not included in this repository) has its own requirements, described in the repository linked above.

- **Data wrangling:** `tidyverse`, `tibble`, `reshape2`, `here`, `conflicted`
- **Statistics:** `vegan`, `ape`, `indicspecies`, `broom`
- **Spatial:** `sf`, `rnaturalearth`, `rnaturalearthdata`
- **Visualization:** `ggplot2`, `patchwork`, `ggforce`, `ggtree`, `ggnewscale`
  
