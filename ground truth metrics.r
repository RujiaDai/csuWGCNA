load('gold20181101.RData')

#functions 

##relevence&recovery
jaccard<-function (x,y) 
{
    if (!(is.list(x)|is.list(y))) 
        stop("Input x, y must be list\n")

    ny<-length(y)
	nx<-length(x)
    Mat = matrix(NA, nx, ny)
    rownames(Mat) = names(x)
	colnames(Mat) = names(y)

    for (i in 1:nx) {
        for (j in 1:ny) 
		Mat[i, j]  = length(intersect(x[[i]],y[[j]]))/length(union(x[[i]], y[[j]]))
    }
    Mat
}

relevence<-c()
recovery<-c()
for(i in 1:8){

df<-jaccard(cm[[i]],km[[i]])
relevence[[i]]<-median(apply(df,1,max))
recovery[[i]]<-median(apply(df,2,max))

}

cscore<-cbind(relevence,recovery)

##specificity & sensitivity
specificity<-function(test,gold,id){
  
  mat<-matrix(NA,2,2)
  mat[1,1]<-length(intersect(test,gold))
  mat[1,2]<-length(setdiff(test,gold))
  mat[2,1]<-length(setdiff(gold,test))
  mat[2,2]<-length(ctable[[id]])-length(union(test,gold))
  specificity<-mat[2,2]/(mat[1,2]+mat[2,2])
  specificity
}



speciMatrix<-function (test,ref,dataset) 
{
    if (!(is.list(test)|is.list(ref))) 
        stop("Input test, ref must be list\n")

    nref<-length(ref)
	ntest<-length(test)
    Mat = matrix(NA, ntest, nref)
    rownames(Mat) = names(test)
	colnames(Mat) = names(ref)

    for (i in 1:ntest) {
        for (j in 1:nref) 
		Mat[i, j]  = specificity(test[[i]],ref[[j]],dataset)
    }
    Mat
}

sensitivity<-function(test,gold,id){
  
  mat<-matrix(NA,2,2)
  mat[1,1]<-length(intersect(test,gold))
  mat[1,2]<-length(setdiff(test,gold))
  mat[2,1]<-length(setdiff(gold,test))
  mat[2,2]<-length(ctable[[id]])-length(union(test,gold))
  sensitivity<-mat[1,1]/(mat[1,1]+mat[2,1])
  sensitivity
}




sensiMatrix<-function (test,ref,dataset) 
{
    if (!(is.list(test)|is.list(ref))) 
        stop("Input test, ref must be list\n")

    nref<-length(ref)
	ntest<-length(test)
    Mat = matrix(NA, ntest, nref)
    rownames(Mat) = names(test)
	colnames(Mat) = names(ref)

    for (i in 1:ntest) {
        for (j in 1:nref) 
		Mat[i, j]  = sensitivity(test[[i]],ref[[j]],dataset)
    }
    Mat
}

##NPV && PPV
precision<-function(test,gold,id){
  
  mat<-matrix(NA,2,2)
  mat[1,1]<-length(intersect(test,gold))
  mat[1,2]<-length(setdiff(test,gold))
  mat[2,1]<-length(setdiff(gold,test))
  mat[2,2]<-length(ctable[[id]])-length(union(test,gold))
  precision<-mat[1,1]/(mat[1,1]+mat[1,2])
  precision
}



preciMatrix<-function (test,ref,dataset) 
{
    if (!(is.list(test)|is.list(ref))) 
        stop("Input test, ref must be list\n")

    nref<-length(ref)
	ntest<-length(test)
    Mat = matrix(NA, ntest, nref)
    rownames(Mat) = names(test)
	colnames(Mat) = names(ref)

    for (i in 1:ntest) {
        for (j in 1:nref) 
		Mat[i, j]  = precision(test[[i]],ref[[j]],dataset)
    }
    Mat
}


npv<-function(test,gold,id){
  
  mat<-matrix(NA,2,2)
  mat[1,1]<-length(intersect(test,gold))
  mat[1,2]<-length(setdiff(test,gold))
  mat[2,1]<-length(setdiff(gold,test))
  mat[2,2]<-length(ctable[[id]])-length(union(test,gold))
  npv<-mat[2,2]/(mat[2,1]+mat[2,2])
  npv
}


npvMatrix<-function (test,ref,dataset) 
{
    if (!(is.list(test)|is.list(ref))) 
        stop("Input test, ref must be list\n")

    nref<-length(ref)
	ntest<-length(test)
    Mat = matrix(NA, ntest, nref)
    rownames(Mat) = names(test)
	colnames(Mat) = names(ref)

    for (i in 1:ntest) {
        for (j in 1:nref) 
		Mat[i, j]  = npv(test[[i]],ref[[j]],dataset)
    }
    Mat
}



cmat<-list()
smat<-list()
umat<-list()

for(i in 1:8){

cmat[[i]]<-npvMatrix(cm[[i]],km[[i]],i)
smat[[i]]<-npvMatrix(sm[[i]],km[[i]],i)
umat[[i]]<-npvMatrix(um[[i]],km[[i]],i)

}



cmat2<-lapply(cmat,function(x){apply(x,1,max)})
smat2<-lapply(smat,function(x){apply(x,1,max)})
umat2<-lapply(umat,function(x){apply(x,1,max)})

csunpv<-unlist(lapply(cmat2,mean))
signednpv<-unlist(lapply(smat2,mean))
unsignednpv<-unlist(lapply(umat2,mean))

signed_metrics<-cbind(signedsensitivity,signedspecificity,signedprecision,signednpv,sscore)
unsigned_metrics<-cbind(unsignedsensitivity,unsignedspecificity,unsignedprecision,unsignednpv,uscore)
csu_metrics<-cbind(csusensitivity,csuspecificity,csuprecision,csunpv,cscore)

colnames(unsigned_metrics)<-colnames(csu_metrics)<-colnames(signed_metrics)<-c('sensitivity','specificity','precision','npv','relevence','recovery')
rownames(unsigned_metrics)<-rownames(csu_metrics)<-rownames(signed_metrics)<-id2

apply(d[,2:7],1,function(x){2/(1/x[1]+1/x[2])+2/(1/x[3]+1/x[4])+x[5]+x[6]})
