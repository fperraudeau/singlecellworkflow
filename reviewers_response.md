For now, this are just notes, we should make this a polished response.

## Open questions

1. Do we need PCA? I think that we don't and I have removed it, but if you feel strongly for including it, we can put it back.
2. Do we need to include all the convoluted preprocessing to go from GEO to a SingleCellExperiment? I actually have a SingleCellExperiment object saved on Github. I see arguments both for and against... it may actually be useful for people to see how to create a SCE from scratch.
3. Should we include spike-ins in the SingleCellExperiment object? I think so because one reviewer explicitly mentions it, but it makes it harder to carry out all the steps.
4. Discussion about choice of K. AIC/BIC don't seem to provide a good answer. What should we do?

## Reviewer 1 (Stephanie Hicks)

> 1. In this workflow, the authors start with a count table. However, the majority of researchers will start with raw reads (e.g. a FASTQ file). It would be great if the author discussed current best practices for the quantification step of scRNA-seq data. Alternatively, the authors could point to other references that have already been developed.

We agree that giving guidance to the readers on how to quantify gene expression in scRNA-seq data is an important topic. However, we feel that it is outside the scope of this paper. We have included a reference to the `scPipe` Bioconductor package, which provides functions to quantify gene expression starting from raw data.
 
> 2. I would like to see the authors take advantage of the rich functionality and data exploration tools for cell- and gene-specific quality control (QC) introduced in low-level analysis workflows such as the one from Lun et al. (2016). Also, in this workflow, the authors create multiple SummarizedExperiment objects (e.g. one with only the top 1000 highly variable genes (HVGs), one with all genes, etc). This doesn’t seem efficient, especially with large single cell data sets such as the 1.3 million cells from embryonic mouse brains. I think both of these concerns can now be addressed with efforts such as the recently developed SingleCellExperiment Bioconductor object (https://github.com/drisso/SingleCellExperiment). For example, the authors could add a “USE” column in the gene- or cell-specific meta table to represent whether or not a particular gene in a particular cell met the filtering criteria applied. The authors could store W in the reduceDim assay of the SingleCellExperiment object.

We agree with the reviewer that the previous version of the workflow was not very efficient in terms of data representation. This was due to the fact that the tools employed in the different steps expected slightly different inputs. We have now harmonized all the tools to expect as input and produce as output an object of class `SingleCellExperiment`, including storing the dimensionality reduced matrix `W` in the `reducedDim` slot and storing the gene-wise variances as a column in the `rowData` for better subsetting and filtering.

> 3. In ZINB-WaVE, the authors specify the number of dimensions for the low-dimensional space (K) to be K=50. Could the authors add more details for the reader explaining why they picked K=50 and describe situations in which a user would want to specify a higher or lower K? In particular, it would be useful to discuss computational time in terms of number of genes and cells. Also, it would be useful to note that if you only wanted to use ZINB-WaVE to remove known covariates for normalization, you can use K=0. 

The question of how to choose K is a very important one and a more comprehensive discussion on the role of K, both on accuracy and computational time, is in the original ZINB-WaVE paper. However, we agree that giving more guidance on the choice of K is a needed addition to this workflow. We have added a section, called "Choice of K", that discusses two functions, `zinbAIC` and `zinbBIC`, that can be used to decide the best value of `K` for a given dataset.

> Minor comments:
> 1. When selecting the top 1000 HVGs, why do the authors not take into account the overall mean-variance relationship and only select genes based on the variance?

This is mostly for illustration purposes, but we do find in our experience that naively selecting most variable genes gives a good set of informative genes, perhaps because it selects on/off genes. We have expanded this section to add a discussion about different strategies that account for the structure of the data can be employed.

> 2. It would be great if the authors referenced other tools available for similar analyses currently available. For example there are several available packages for normalization of scRNA-seq data, such as calculating global scaling factors can be done with scran (https://bioconductor.org/packages/release/bioc/html/scran.html) or gene and cell-specific scaling factors using SCnorm (https://github.com/rhondabacher/SCnorm). Alternatively, users might want to try using relative transcript counts using Census (https://bioconductor.org/packages/release/bioc/html/monocle.html).

This is a good suggestion. We have added a new section at the end of the workflow, named "Alternative approaches", that discusses these and other strategies for dimensionality reduction, normalization, clustering, and batch correction.

## Reviewer 2 (Andrew McDavid)

Major comments:

> 1. There is something wrong with the data download link in the F1000 version so that I am unable to download these files and actually reproduce the workflow. I experimented a bit to see if I could figure out how to download the data anyways, but will reserve further evaluation of this submission until this issue can be resolved by the authors.

```{r}
urls = c("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95601&format=file&file=GSE95601%5FoeHBCdiff% "https://raw.githubusercontent.com/rufletch/p63-HBC-diff/master/ref/oeHBCdiff_clusterLabels.)
```

We apologized with the reviewer. Something must have happened during the formatting of the article for production on the F1000 Research website that omitted some of the characters from the url. The links work with no issues in the R markdown version of the workflow, available at https://github.com/fperraudeau/singlecellworkflow. We will make sure that the new version of the article does not have this problem.

> 2a. This workflow will likely be out-of-date when the underlying packages transition to use SingleCellExperiment. This is actually a positive thing because many of the more opaque lines of code (involving subsetting ERCC genes, etc) will be more streamlined.

The reviewer is correct. We now make full use of the features of SingleCellExperiment and as a result the code is much cleaner and easier to understand.

> 2b. It requires installation of the development branch of bioconductor, which impacts the usefulness of the workflow to the average user. I expect the authors will revise this tutorial when Bioconductor 3.6 is released and use of the devel branch is no longer necessary. Additionally `slingshot` is an requirement, but currently only exists on github and no SHA1 provided. I hope that `slingshot` will be added as a bioconductor package shortly. In the meantime, a tag must be added to the git repo for the release being used in this workflow and instructions provided for how to install this tag. Additionally, the authors may wish to note that installation instructions for the packages will be provided at the end of the workflow so that someone proceeding sequentially will not be tripped up.

Slingshot is now a Bioconductor package, albeit in the devel branch for now. Since the new release is only a few months away, we have decided to require the use of Bioconductor 3.8 (currently devel, but soon release). We think that this is the best choice to avoid the workflow to be outdated later this year. We have added the code at the beginning to show how to install the required packages, using the new recommended `BiocManager` package.

> 2c. Opaque code is presented in order to generate plots, e.g.
```{r}
palDF <- ceObj@clusterLegend[[1]]
pal <- palDF[, "color"]
names(pal) <- palDF[, "name"]
pal["-1"] = "transparent"
plot(fit$points, col = pal[primaryClusterNamed(ceObj)], main = "", pch = 20, xlab = "Component1", ylab = "Component2") legend(x = "topleft", legend = names(pal), cex = .5, fill = pal, title = "Sample")
```

> While this complexity may be necessary, perhaps some of it could be encapsulated as accessor functions in the package? Too much complexity here may cause users to miss the forest for the trees.  

We have created the `plotReduceDims()` function in `clusterExperiment` to address this issue.
 
> The authors could better motivate (or at least explain the impact of) some of the default parameters and procedures.
> Why do we set a zcut threshold of 3 for the `scone` filtering? 
> Why K=50 for zinbwave?
> RSEC parameters

> How should a user decide on a value for these parameters?

We agree with the reviewer that more details on the choice of these parameters are needed.
As the choice of `K` for the `zinbwave` model is critical, we have added a new section that describes a way to compare different values of `K` and select the best.
As for the `zcut` threshold and the `RSEC` parameters, we have added a paragraph in the appropriate section describing in more details what the parameters control and what are good rules of thumb to set them to an appropriate value. [TODO]

## Reviewer 3 (Mike Love)

> It would be useful to put a link to the source code (or to the section where the link to the source code exists) near the top of the document.

We have modified the section "Package version" to include a link to the source code as well as instructions on how to install the packages needed by the workflow.

> I was confused a bit by "the first major bifurcation in the HBC lineage trajectory occurs prior to cell division". Can you be more specific about what you are referring to here by cell division, as without knowledge of the system, I'm not sure where the cell division you refer to should appear.

Thank you, we agree that the paragraph taken out of context is unclear. The first bifurcation identified by slingshot (occurring right after the clusters marked by DeltaHBC1 and DeltaHBC2 in Figure 2) produces two lineage trajectories: the first one (purple curve in Figure 2) generates sustentacular cells through a differentiation process characterized by the absence of cell cycle and division. To simplify, each HBC stem cells gives rise to one sustentacular cell. The other lineage trajectory (blue and orange curves in Figure 2) generates proliferative GBCs, which are characterized by cell cycle and division. These, in turn, generate olfactory sensory neurons and microvillous cells. We hope that this more detailed explanation will be useful to better understand our system. We have updated the text to add this more detailed explanation.

> "within a single object": It may be good to explain what an "object" here is. You could, for example, refer to Figure 2 of the Bioconductor Nature Methods paper.

We have added a paragraph discussing what a SingleCellExperiment object is and what are the advantages of using a single object throughout a data analysis workflow. 

> Misspelling: "reasonnable"

Fixed, thanks.

> On filtering for most variable genes, I understand this decision, and I also recommend it during workshops before making ordination plots. I know that students are not always certain why we care about variance (unsupervised). I like to mention that these are the genes where the "action" is. A side point, the log(x+1) is not variance stabilizing for RNA-seq counts in general. This filter can give higher priority to low count genes than to genes where there is interesting biological variability (though I do not doubt that the very high biological variance genes will be preserved). It might be useful to show a vsn::meanSdPlot() for the matrix log1p(assay(se))? 

Thank you, this is a very good point. See also comment by Reviewer 1 about modeling the mean-variance relationship. We have expanded this section to add a discussion about why we care about the most variable genes and how different strategies that account for the structure of the data can be employed.

> "correcting for batch effects": What are batch effects? (Of course, I know what they are, but a reader may not, and you could cite some of the single cell literature here.)

We have expanded our discussion on batch effects to include an explanation of what we mean by this term in this specific example and more generally.

> "Note that, in this case, the low-dimensional matrix W is not included in the computation of residuals to avoid the removal of the biological signal of interest.": I understood this sentence only on a second pass through. One problem is that you haven't defined W in the text yet (only in Figure 4). I would only reference the matrix W if you have defined it.

We agree, that the reference to the model is too specific. We have modified the sentence. [TODO -- still not sure how to best describe this step]

> Figure 6: Can you change the figure width so that PC1 is not squished?

We actually dropped Figure 6 as using the low-dimensional signal from ZINB-WaVE (rather than PCA on the normalized counts) is the recommended way to visualize the results, and we thought that it was confusing for the reader to show PCA at this stage of the workflow.

> Is there a circularity to the recovery of published clusters in Figure 6? Was ZINB-WaVE used in Fletcher (2017)?

ZINB-WaVE was not used in the original paper, which used a more complex normalization strategy, based on regressing out a set of QC measures. One key result in favor of the ZINB-WaVE projection is that with this method slingshot is able to estimate the correct lineages with minimal supervision (only the starting cluster), while in the original paper, one of the terminal states had to be specified.

DR: Make sure that with the original analysis we had to use a semi-supervised approach.

> Can you say what the meaning of the color white is in Figure 8 (in the text or caption near this figure)?

We have added a better description of the Figure in the caption.

> Figure 15 refers back to Figure 2 but does not use the same color scheme for the known cell types, so the reader cannot verify if you've recovered the lineages from the publication. It would be good therefore to have a legend for these figures (Fig 15 and following) indicating which cell types the colors refer to (this information is in the unlabeled table above, but should be included as a legend here).

Thanks for the suggestion, we agree that having a similar color scheme will help the reader make the connection between the results of Figure 15 and Figure 2. We have changed the colors to match as close as possible those of the published figure. Note that we do this right after the clustering step, so that all figures (from Figure 8 to Figure 17) have consistent colors.

> Can you briefly describe what a GAM is ahead of Figure 17?

We have added a sentence describing GAMs ahead of Figure 17.
