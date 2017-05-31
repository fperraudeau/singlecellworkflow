load("../data/GSE95601_oeHBCdiff_Cufflinks_eSet")
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


metaData<-data.frame("Experiment"=bio,"Batch"=batch,qc)

# filtering to top 1000 most variable:
library(matrixStats)
vars <- rowVars(log1p(E))
names(vars) <- rownames(E)
vars <- sort(vars, decreasing = TRUE)
core <- E[names(vars)[1:1000],]


write.table(core,file="../data/oeCufflinkCountData_1000Var.txt",sep="\t",row.names=TRUE,col.names=TRUE)
write.table(metaData,file="../data/oeMetadata.txt",sep="\t",row.names=TRUE,col.names=TRUE)