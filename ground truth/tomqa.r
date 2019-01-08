# function for counting the anti-recorcing triplets in TOM
antCounting<-function(df, n){
n<-nrow(df)
time0<-0
time1<-0
for (i in 1:(n-2)){
  for(j in 2:(n-1)){
    for(k in 3:n){
	if(i<j&j<k){
      time0=time0+1
	  if(df[i,j]<0&df[i,k]<0&df[j,k]<0){time1=time1+1}
      if(df[i,j]<0&df[i,k]>0&df[j,k]>0){time1=time1+1}
	  if(df[i,j]>0&df[i,k]<0&df[j,k]>0){time1=time1+1}
	  if(df[i,j]>0&df[i,k]>0&df[j,k]<0){time1=time1+1}
	  else {time1=time1+0}
	}
    }
 }
}
time1/time0
}




antpro<-list()

id=4
testdata<-qndata[[id]]

time<-c()
for (i in 1:length(km[[id]])){
  
  def<-testdata[,match(km[[id]][[i]],colnames(testdata))]
  corr<-bicor(def)
  time<-c(time,antCounting(corr))
  
  }
antpro[[id]]<-time

dataset<-rep(id2[[8]],length(antpro[[8]]))
dff8<-as.matrix(antpro[[8]])
dff8<-cbind(dff8,dataset)
colnames(dff8)[1]<-'antpro'
