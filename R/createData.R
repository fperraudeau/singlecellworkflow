library(GEOquery)
library(scone)
if(!file.exists("../data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda")) {
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95601&format=file&file=GSE95601%5FoeHBCdiff%5FCufflinks%5FeSet%2ERda%2Egz",
                "../data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda.gz")
  gunzip("../data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda.gz")
}
load("../data/GSE95601_oeHBCdiff_Cufflinks_eSet.Rda")

if(!file.exists("../data/oeHBCdiff_clusterLabels.txt")) {
  download.file("https://raw.githubusercontent.com/rufletch/p63-HBC-diff/master/ref/oeHBCdiff_clusterLabels.txt",
                "../data/oeHBCdiff_clusterLabels.txt")
}

## count matrix
E <- assayData(Cufflinks_eSet)$counts_table

# Remove undetected genes
E <- na.omit(E)
E <- E[rowSums(E)>0,]

# Remove ERCC and CreER
cre <- E["CreER",]
ercc <- E[grep("^ERCC-", rownames(E)),]
E <- E[grep("^ERCC-", rownames(E), invert = TRUE), ]
E <- E[-which(rownames(E)=="CreER"), ]

# Extract QC metrics
qc <- as.matrix(protocolData(Cufflinks_eSet)@data)[,c(1:5, 10:18)]
qc <- cbind(qc, CreER = cre, ERCC_reads = colSums(ercc))

# Extract metadata
batch <- droplevels(pData(Cufflinks_eSet)$MD_c1_run_id)
bio <- droplevels(pData(Cufflinks_eSet)$MD_expt_condition)

clusterLabels<-read.table("../data/oeHBCdiff_clusterLabels.txt",sep="\t",stringsAsFactors=FALSE)
m<-match(colnames(E) ,clusterLabels[,1])
metadata<-data.frame("Experiment"=bio,"Batch"=batch,"clusterLabels"=clusterLabels[m,2],qc)

# symbol for samples missing from original clustering
metadata$clusterLabels[is.na(metadata$clusterLabels)] <- -2

se <- SummarizedExperiment(assays = list(counts = E),
                           colData = metadata)


# Initial Gene Filtering
num_reads = quantile(assay(se)[assay(se) > 0])[4]
num_cells = 0.25*ncol(assay(se))
is_common = rowSums(assay(se) >= num_reads ) >= num_cells
table(is_common)

# Metric-based Filtering sample-filtering
data("housekeeping")
hk = rownames(se)[toupper(rownames(se)) %in% housekeeping$V1]

mfilt <- metric_sample_filter(assay(se), nreads = colData(se)$NREADS,
                             ralign = colData(se)$RALIGN,
                             pos_controls = rownames(se) %in% hk,
                             zcut = 3, mixture = FALSE,
                             plot = TRUE)

# Simplify to a single logical
mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
table(mfilt)

se <- se[,mfilt]

# filtering to top 1000 most variable:
library(matrixStats)
vars <- rowVars(log1p(assay(se)))
names(vars) <- rownames(se)
vars <- sort(vars, decreasing = TRUE)
core <- se[names(vars)[1:1000],]

save(core, file="../data/oe_se_1000Var.rda")
