# Bioconductor workflow for single-cell RNA sequencing: Normalization, dimensionality reduction, clustering, and lineage inference

This repository is designed to provide a tutorial for the analysis of scRNA-seq data in R. It covers four main steps: (1) dimensionality reduction accounting for zero inflation and over-dispersion and adjusting for gene and cell-level covariates; (2) robust and stable cell clustering using resampling-based sequential ensemble clustering; (3) inference of cell lineages and ordering of the cells by developmental progression along lineages; and (4) DE analysis along lineages. The workflow is general and flexible, allowing the user to sustitute the statistical method used in each step by a different method. We hope our proposed workflow will ease technical aspects of scRNA-seq data analysis and help with the discovery of novel biological insights.


## Dependencies

To be able to run [workflow.Rmd](https://github.com/fperraudeau/singlecellworkflow/blob/master/workflow/workflow.Rmd), you need

### Bioconductor

- BiocParallel
- clusterExperiment
- scone
- zinbwave

### GitHub
- slingshot (https://github.com/kstreet13/slingshot)

### CRAN
- doParallel
- gam
- RColorBrewer


Note that you need the devel versions
of the Bioconductor packages `scone (>=1.1.2)`, `zinbwave (>=0.99.6)`, and `clusterExperiment (>=1.3.2)`. We recommend running Bioconductor 3.6 (currently the devel version; see https://www.bioconductor.org/developers/how-to/useDevel/).


