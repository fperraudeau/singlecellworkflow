#!/usr/bin/env Rscript
library(zinbwave)
library(matrixStats)
load("../data/Expt4c2b_filtdata.Rda")
load("../data/E4c2b_slingshot_wsforkelly.RData")
set.seed(20)

names(batch) <- colnames(counts)

counts <- counts[,names(clus.labels)]
batch <- droplevels(batch[names(clus.labels)])
qc <- qc[names(clus.labels),]

vars <- rowVars(log1p(counts))
names(vars) <- rownames(counts)
vars <- sort(vars, decreasing = TRUE)
core <- counts[names(vars)[1:1000],]

NCORES = 7
time_no_batch = list()
time_batch = list()
K = c(0, 100)
for (i in 1:length(K)){
  print(i)
  time_no_batch[[i]] = system.time(zinb0 <- zinbFit(core, K = K[i], ncores = NCORES))
}


mod = model.matrix( ~ batch)
for (i in 1:length(K)){
  print(i)
  time_batch[[i]] = system.time(zinb0 <- zinbFit(core, K = K[i], X = mod, ncores = NCORES))
}
print(time_batch)
print(time_no_batch)

save(time_batch, time_no_batch, file = 'time_zinb_k0.rda')
