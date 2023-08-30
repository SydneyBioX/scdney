## ------------------------------------------------------------------------
# load library
suppressPackageStartupMessages({
  library(scdney)
  library(mclust)
})

# load dataset
data("GSE82187.sample")
dat <- GSE82187

dat <- log2(dat+1)

# show the data
dat[1:5,1:5]

# set number of clusters (classes defined in colnames)
nCs <- length(table(colnames(dat)))

# cell types
cellTypes = colnames(dat)

dat.processed = dat

## ----warning=F, error=F--------------------------------------------------
# since we have already filtered genes, we will set 'geneFilter=1'. Alternatively, you can log transform your data without filtering and set 'geneFilter=0.8'. This will return same result.
dat.scClust <- dat.processed
colnames(dat.scClust) <- cellTypes
simlr.result <- scClust(dat.scClust, nCs, method = "simlr", similarity = "pearson", geneFilter = 1, seed = 1, cores.ratio = 0)

# Evaluate
adjustedRandIndex(cellTypes, simlr.result$y$cluster)

## ----warning=F, error=F--------------------------------------------------
km.result <- scClust(dat.scClust, nCs, method = "kmeans", similarity = "pearson", geneFilter = 1, seed = 1)

# Evaluate
adjustedRandIndex(cellTypes, km.result$cluster)

## ----warning=F, error=F--------------------------------------------------
simlr.bench <- scClustBench(dat.scClust, nCs, method = "simlr", similarity = c("euclidean", "pearson"), cores = 1, rep = 2, cores.ratio = 0)

simlr.bench.eval <- evalScClustBench(simlr.bench, method = "simlr")

p <- plotSimlrEval(simlr.bench.eval)
p

## ----warning=F, error=F--------------------------------------------------
start = proc.time()
km.bench <- scClustBench(dat.scClust, nCs, method = "kmeans", similarity = NULL, rep = 2, cores = 1)
end = proc.time()
end- start
km.bench.eval <- evalScClustBench(km.bench, method = "kmeans")

p <- plotKmeansEval(km.bench.eval)
p


## ------------------------------------------------------------------------
# load dataset
data(GSE87795_liver.development.data)
dat <- GSE87795_liver.development.data$data
cellTypes <- GSE87795_liver.development.data$cellTypes

# Dimensions of the dataset
dim(dat)

# Show dataset
dat[1:5,1:5]

# Cell types
table(cellTypes)

# number of clusters
nCs <- length(table(cellTypes))

## ------------------------------------------------------------------------
# # Filter low expressed genes
# del <- which(rowSums(dat == 0) / ncol(dat) >= 0.8)
# dat.filtered <- dat[-del,]
#
# # log2 transformation
# dat.processed <- log2(dat.filtered + 1)
dat.processed <- dat

## ------------------------------------------------------------------------
dat.selected <- matPCs(dat.processed, 0.7)

## ------------------------------------------------------------------------
lab <- cellTypes

set.seed(1)
noisyCls <- function(dat, rho, cls.truth){
  cls.noisy <- cls.truth
  names(cls.noisy) <- colnames(dat)
  for(i in 1:length(table(cls.noisy))) {
    # class label starts from 0
    if (i != length(table(cls.noisy))) {
      cls.noisy[sample(which(cls.truth == names(table(cls.noisy))[i]), floor(sum(cls.truth == names(table(cls.noisy))[i]) * rho))] <- names(table(cls.noisy))[i+1]
    } else {
      cls.noisy[sample(which(cls.truth == names(table(cls.noisy))[i]), floor(sum(cls.truth == names(table(cls.noisy))[i]) * rho))] <- names(table(cls.noisy))[1]
    }
  }

  print(sum(cls.truth != cls.noisy))
  return(cls.noisy)
}

cls.noisy01 <- noisyCls(dat.selected, rho=0.1, lab)
cls.noisy02 <- noisyCls(dat.selected, rho=0.2, lab)
cls.noisy03 <- noisyCls(dat.selected, rho=0.3, lab)
cls.noisy04 <- noisyCls(dat.selected, rho=0.4, lab)
cls.noisy05 <- noisyCls(dat.selected, rho=0.5, lab)

## ---- fig.height=10, fig.width=10----------------------------------------
# names(lab) <- colnames(dat.selected)

###################################
# SVM
###################################
acc01 <- acc02 <- acc03 <- acc04 <- acc05 <- c()
ari01 <- ari02 <- ari03 <- ari04 <- ari05 <- c()
base <- "svm"

for(j in 1:10) {
  final <- multiAdaSampling(dat.selected, cls.noisy01, seed=j, classifier=base, percent=1, L=10)$final
  ari01 <- c(ari01, mclust::adjustedRandIndex(lab, final))
  acc01 <- c(acc01, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy02, seed=j, classifier=base, percent=1, L=10)$final
  ari02 <- c(ari02, mclust::adjustedRandIndex(lab, final))
  acc02 <- c(acc02, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy03, seed=j, classifier=base, percent=1, L=10)$final
  ari03 <- c(ari03, mclust::adjustedRandIndex(lab, final))
  acc03 <- c(acc03, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy04, seed=j, classifier=base, percent=1, L=10)$final
  ari04 <- c(ari04, mclust::adjustedRandIndex(lab, final))
  acc04 <- c(acc04, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy05, seed=j, classifier=base, percent=1, L=10)$final
  ari05 <- c(ari05, mclust::adjustedRandIndex(lab, final))
  acc05 <- c(acc05, bAccuracy(lab, final))
}


###################################
# RF
###################################
rfacc01 <- rfacc02 <- rfacc03 <- rfacc04 <- rfacc05 <- c()
rfari01 <- rfari02 <- rfari03 <- rfari04 <- rfari05 <- c()
base <- "rf"
for(j in 1:10) {
  final <- multiAdaSampling(dat.selected, cls.noisy01, seed=j, classifier=base, percent=1, L=10)$final
  rfari01 <- c(rfari01, mclust::adjustedRandIndex(lab, final))
  rfacc01 <- c(rfacc01, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy02, seed=j, classifier=base, percent=1, L=10)$final
  rfari02 <- c(rfari02, mclust::adjustedRandIndex(lab, final))
  rfacc02 <- c(rfacc02, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy03, seed=j, classifier=base, percent=1, L=10)$final
  rfari03 <- c(rfari03, mclust::adjustedRandIndex(lab, final))
  rfacc03 <- c(rfacc03, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy04, seed=j, classifier=base, percent=1, L=10)$final
  rfari04 <- c(rfari04, mclust::adjustedRandIndex(lab, final))
  rfacc04 <- c(rfacc04, bAccuracy(lab, final))
  final <- multiAdaSampling(dat.selected, cls.noisy05, seed=j, classifier=base, percent=1, L=10)$final
  rfari05 <- c(rfari05, mclust::adjustedRandIndex(lab, final))
  rfacc05 <- c(rfacc05, bAccuracy(lab, final))
}

result = list(
  rfacc01 = rfacc01,
  rfacc02 = rfacc02,
  rfacc03 = rfacc03,
  rfacc04 = rfacc04,
  rfacc05 = rfacc05,
  acc01 = acc01,
  acc02 = acc02,
  acc03 = acc03,
  acc04 = acc04,
  acc05 = acc05,
  rfari01 = rfari01,
  rfari02 = rfari02,
  rfari03 = rfari03,
  rfari04 = rfari04,
  rfari05 = rfari05,
  ari01 = ari01,
  ari02 = ari02,
  ari03 = ari03,
  ari04 = ari04,
  ari05 = ari05

)


plot.new()
par(mfrow = c(2,2))
boxplot(rfacc01, rfacc02, rfacc03, rfacc04, rfacc05, col="orange", main="RF Acc", ylim=c(0.45, 1))
points(x=1:5, y=c(bAccuracy(lab, cls.noisy01), bAccuracy(lab, cls.noisy02),
                  bAccuracy(lab, cls.noisy03), bAccuracy(lab, cls.noisy04),
                  bAccuracy(lab, cls.noisy05)), col="red3", pch=c(2,3,4,5,6), cex=1)
boxplot(acc01, acc02, acc03, acc04, acc05, col="lightblue", main="SVM Acc", ylim=c(0.45, 1))
points(x=1:5, y=c(bAccuracy(lab, cls.noisy01), bAccuracy(lab, cls.noisy02),
                  bAccuracy(lab, cls.noisy03), bAccuracy(lab, cls.noisy04),
                  bAccuracy(lab, cls.noisy05)), col="red3", pch=c(2,3,4,5,6), cex=1)
boxplot(rfari01, rfari02, rfari03, rfari04, rfari05, col="orange", main="RF ARI", ylim=c(0.25, 1))
points(x=1:5, y=c(mclust::adjustedRandIndex(lab, cls.noisy01), mclust::adjustedRandIndex(lab, cls.noisy02),
                  mclust::adjustedRandIndex(lab, cls.noisy03), mclust::adjustedRandIndex(lab, cls.noisy04),
                  mclust::adjustedRandIndex(lab, cls.noisy05)), col="red3", pch=c(2,3,4,5,6), cex=1)
boxplot(ari01, ari02, ari03, ari04, ari05, col="lightblue", main="SVM ARI", ylim=c(0.25, 1))
points(x=1:5, y=c(mclust::adjustedRandIndex(lab, cls.noisy01), mclust::adjustedRandIndex(lab, cls.noisy02),
                  mclust::adjustedRandIndex(lab, cls.noisy03), mclust::adjustedRandIndex(lab, cls.noisy04),
                  mclust::adjustedRandIndex(lab, cls.noisy05)), col="red3", pch=c(2,3,4,5,6), cex=1)

## ------------------------------------------------------------------------
# PCA procedure
print("1")
dat.pc <- matPCs(dat.processed, 0.7)
dim(dat.pc)
print("2")

# run scReClassify
cellTypes.reclassify <- multiAdaSampling(dat.pc, cellTypes, seed = 1, classifier = "svm", percent = 1, L = 10)
print("3")

# Verification by marker genes
End <- c("KDR", "LYVE1")
Meg <- c("ITGA2B", "ITGB3")
Mes <- c("MEST", "MMP2")
Ery <- c("HBA-A1", "HBB-BT")

# check examples
idx <- which(cellTypes.reclassify$final != cellTypes)
library(dplyr)
cbind(original=cellTypes[idx], reclassify=cellTypes.reclassify$final[idx]) %>%
  DT::datatable()
print("4")

c1 <- dat.processed[, which(cellTypes=="Endothelial Cell")]
c2 <- dat.processed[, which(cellTypes=="Erythrocyte")]
c3 <- dat.processed[, which(cellTypes=="Hepatoblast")]
c4 <- dat.processed[, which(cellTypes=="Macrophage")]
c5 <- dat.processed[, which(cellTypes=="Megakaryocyte")]
c6 <- dat.processed[, which(cellTypes=="Mesenchymal Cell")]
cs <- rainbow(length(table(cellTypes)))
print("5")

# (example 1 E13.5_C20)
#####
par(mfrow=c(1,2))
marker <- End[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(1, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C20")], pch=16, col="red", cex=2)
marker <- End[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(1, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C20")], pch=16, col="red", cex=2)
#####
print("6")

# (example 2 E13.5_C14)
#####
par(mfrow=c(1,2))
marker <- Meg[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(5, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C14")], pch=16, col="red", cex=2)
marker <- Meg[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(5, dat.processed[marker, which(colnames(dat.processed) %in% "E13.5_C14")], pch=16, col="red", cex=2)
#####
print("7")

# (example 3 E16.5_2_C65)
#####
par(mfrow=c(1,2))
marker <- Mes[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(6, dat.processed[marker, which(colnames(dat.processed) %in% "E16.5_2_C65")], pch=16, col="red", cex=2)
marker <- Mes[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(6, dat.processed[marker, which(colnames(dat.processed) %in% "E16.5_2_C65")], pch=16, col="red", cex=2)
######
print("8")

# (example 4 E16.5_C17)
#####
par(mfrow=c(1,2))
marker <- Mes[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(6, dat.processed[marker, which(colnames(dat.processed) %in% "E16.5_C17")], pch=16, col="red", cex=2)
marker <- Mes[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(6, dat.processed[marker, which(colnames(dat.processed) %in% "E16.5_C17")], pch=16, col="red", cex=2)
#####
print("9")

# (example 5 E12.5_C72)
#####
par(mfrow=c(1,2))
marker <- Ery[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(2, dat.processed[marker, which(colnames(dat.processed) %in% "E12.5_C72")], pch=16, col="red", cex=2)
marker <- Ery[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(2, dat.processed[marker, which(colnames(dat.processed) %in% "E12.5_C72")], pch=16, col="red", cex=2)
#####
print("10")

# (example 6 E12.5_C07)
#####
par(mfrow=c(1,2))
marker <- Ery[1]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(2, dat.processed[marker, which(colnames(dat.processed) %in% "E12.5_C07")], pch=16, col="red", cex=2)
marker <- Ery[2]
boxplot(c1[marker,], c2[marker,], c3[marker,], c4[marker,], c5[marker,], c6[marker,], col=cs, main=marker)
points(2, dat.processed[marker, which(colnames(dat.processed) %in% "E12.5_C07")], pch=16, col="red", cex=2)
#####
print("11")
