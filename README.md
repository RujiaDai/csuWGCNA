# csuWGCNA
csuWGCNA is a modified WGCNA precedure which is intended to capture negative correlations in the expression data.



## Usage:
library(WGCNA); 
library(doParallel);
source('hpickSoftThreshold.r');
source('Hadjacency.r');
options(stringsAsFactors = F)

## power selection
sft<-hpickSoftThreshold(expr,powerVector = powers, corFnc = bicor, networkType = "hybrid2", verbose = 5)

## adjacency matrix construction
adj<-Hadjacency(datExpr=expr,type='hybrid2',power=power, corFnc='bicor')

