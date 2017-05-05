require(clusterExperiment)
require(SummarizedExperiment)
require(limma)

load('data/Expt4c2b_filtdata.Rda')
load('data/E4c2b_slingshot_wsforkelly.RData')
load('R/olfactory.rda')

names(batch) <- colnames(counts)
counts <- counts[,names(clus.labels)]
batch <- droplevels(batch[names(clus.labels)])
qc <- qc[names(clus.labels),]

se = SummarizedExperiment(assays = list(counts = counts),
                           colData = qc)

assay(se)[1:2,1:5]
colData(se)[,1:5]
pass_filter = apply(assay(se), 1, function(x) length(x[x >= 10]) >= 10)
dim(se)
sum(!pass_filter)
se <- se[pass_filter,]
dim(se)
fq <- round(limma::normalizeQuantiles(assay(se)))
assays(se) <- list(normalized_counts=fq)


# clusterMany
# Parameters that seem relevant when input is zinbwave W matrix
# isCount = FALSE
# no dimReduce
ce <- clusterMany(se, clusterFunction = "pam", ks = 5:10,
                  isCount = TRUE, dimReduce = c("PCA", "var"),
                  nVarDims = c(100, 500, 1000), nPCADims = c(5, 15, 50),
                  run = TRUE)

defaultMar <- par("mar")
plotCMar <- c(1.1, 16.1, 4.1, 1.1)
par(mar = plotCMar)
plotClusters(ce, main = "Clusters from clusterMany",  axisLine = -1)

# remove 'features' from the names of the labels
cl<-clusterLabels(ce)
cl<-gsub("Features","",cl)
clusterLabels(ce)<-cl

# different ordering of the clusterings
cl <- clusterLabels(ce)
ndims = sapply(cl, function(x) strsplit(x, '=|,')[[1]][2])
ord <- order(as.numeric(ndims))
par(mar = plotCMar)
plotClusters(ce, main = "Clusters from clusterMany", whichClusters = ord, axisLine=-1)

# combineMany
clusterMatrix(ce)[1:5,1:10]
ce<-combineMany(ce)
par(mar=plotCMar)
plotClusters(ce,whichClusters="workflow")
wh<-which(clusterLabels(ce)=="combineMany")
if(length(wh)!=1) stop() else clusterLabels(ce)[wh]<-"combineMany,default"
# change proportion to 0.7
ce<-combineMany(ce,proportion=0.7,clusterLabel="combineMany,0.7")
par(mar=plotCMar)
plotClusters(ce,whichClusters="workflow")
# chamge minSize to 3 and proportion to 0.7
ce<-combineMany(ce,proportion=0.7,minSize=3,clusterLabel="combineMany,final")
par(mar=plotCMar)
plotClusters(ce,whichClusters="workflow", main="Min. Size=3\nProportion=0.7")

# plot dendogram
# plotCoClustering(ce)
# error
# Error in gridPLT() : Figure region too small and/or viewport too large

# makeDendrogram
ce<-makeDendrogram(ce,dimReduce="PCA",ndims=500)
plotDendrogram(ce)
ce

#mergeClusters
mergeClusters(ce,mergeMethod="adjP",plot="mergeMethod")
ce<-mergeClusters(ce,mergeMethod="adjP",cutoff=0.1)

par(mar=plotCMar)
plotClusters(ce,whichClusters="workflow")
#plotCoClustering(ce,whichClusters=c("mergeClusters","combineMany"), annLegend=FALSE)
#plotHeatmap(ce,clusterSamplesData="dendrogramValue",breaks=.99)


# Step 4: Finding Features related to the clusters
pairsAll<-getBestFeatures(ce,contrastType="Pairs",p.value=0.05,
                          number=nrow(ce),isCount=TRUE)
head(pairsAll)
length(pairsAll$Feature)==length(unique(pairsAll$Feature))

plotHeatmap(ce, clusterSamplesData="dendrogramValue",
            clusterFeaturesData=unique(pairsAll[,"IndexInOriginal"]),
            main="Heatmap of features w/ significant pairwise differences",
            breaks=.99)



