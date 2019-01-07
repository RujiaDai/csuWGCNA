#simu
options(stringsAsFactors = F)
library(WGCNA);library(doParallel)

nSamples<-50
nGenes<-2000
set.seed(6)
v1 = rnorm(nSamples);
v2 = rnorm(nSamples);
v3 = rnorm(nSamples);
v4 = rnorm(nSamples);
v5 = rnorm(nSamples);
eigengenes = cbind(v1,v2,v3,v4,v5)
#modProps = c(0.35,0.18,0.1, 0.05, 0.02, 0.3)
modProps = c(0.4,0.25,0.2, 0.1, 0.05, 0)
corprop<-c(seq(0.01,0.1,0.01),0.15,0.2,0.25,0.3)
data<-list()
for(i in 1:length(corprop)){data[[i]]<-simulateDatExpr(eigengenes, nGenes, modProps,signed = FALSE, propNegativeCor = corprop[[i]])}

datExpr = list()
for(i in 1:length(corprop)){datExpr[[i]]<-data[[i]]$datExpr}

 source('hpickSoftThreshold.r')

powers<-seq(2,20,2)
sft<-list()
for(i in 1:length(datExpr)){

  dat<-datExpr[[i]]

 sft[[i]]<-hpickSoftThreshold(t(dat),powerVector = powers, corFnc = bicor, networkType = "hybrid2", verbose = 5)


}



source('Hadjacency.r')
library(reshape2)

stom<-list()
utom<-list()
htom<-list()

for(i in 1:length(datExpr)){

  dat<-datExpr[[i]]
  uadj<-adjacency(datExpr=dat,type='unsigned',power=12,corFnc='bicor')
  utom[[i]]<-TOMsimilarity(uadj, TOMType = "unsigned", verbose = 1)
  sadj<-adjacency(datExpr=dat,type='signed',power=12,corFnc='bicor')
  stom[[i]]<-TOMsimilarity(sadj, TOMType = "unsigned", verbose = 1)
  hadj<-Hadjacency(datExpr=dat,type='hybrid2',power=12,corFnc='bicor')
  htom[[i]]<-TOMsimilarity(hadj, TOMType = "unsigned", verbose = 1)

}

source('Hadjacency.r')
library(reshape2)

scolor<-list()
ucolor<-list()
hcolor<-list()
ds = 2; minModSize = 30; dthresh = 0.1; pam = FALSE;
for(i in 1:length(stom)){


    tree = hclust(1-as.dist(stom[[i]]), method="average")
    cut = cutreeHybrid(dendro = tree, pamStage=pam, minClusterSize= minModSize, cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-stom[[i]]))
    merged = mergeCloseModules(exprData=datExpr[[i]], colors =  cut$labels, cutHeight=dthresh)
    scolor[[i]]<-labels2colors(merged$colors)
    tree = hclust(1-as.dist(utom[[i]]), method="average")
    cut = cutreeHybrid(dendro = tree, pamStage=pam, minClusterSize= minModSize, cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-utom[[i]]))
    merged = mergeCloseModules(exprData=datExpr[[i]], colors =  cut$labels, cutHeight=dthresh)
    ucolor[[i]]<-labels2colors(merged$colors)
    tree = hclust(1-as.dist(htom[[i]]), method="average")
    cut = cutreeHybrid(dendro = tree, pamStage=pam, minClusterSize= minModSize, cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-htom[[i]]))
    merged = mergeCloseModules(exprData=datExpr[[i]], colors =  cut$labels, cutHeight=dthresh)
    hcolor[[i]]<-labels2colors(merged$colors)

}


