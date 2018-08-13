# Summary of changes

* We make full use of the new Bioconductor class `SingleCellExperiment` and we now apply all the steps of the workflow on the same object.
* We now provide detailed instructions on how to install the packages.
* We provide a better explanation of the biological system used to illustrate the workflow.
* We provide a brief explanation of data containers, "objects", and classes in Bioconductor.
* We provide a more thorough discussion on why it may be a beneficial step to focus on the most variable genes and we discuss variance stabilizing transformations that may lead to better variance estimates.
* We introduce to new functions, `makeFilterStats` and `filterData`, to make it easier to select the most variable genes.
* We provide a brief description of batch effects in the context of scRNA-seq.
* We use a more consistent color-scheme throughout the manuscript.
* We provide a brief explanation of Generalized Additive Models (GAMs).
* We have added guidance on how to choose the value of `K`.
* We have added a section that describes alternative approaches.
* We have fixed the url to GEO that were broken in the F1000 Research version of the workflow.
