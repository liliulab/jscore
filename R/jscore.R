#' @title jscore
#' @description Calculates the j-score between two clustering assignment.
#' @param truth A numeric vector of truth class labels.
#' @param pred A numeric vector of predicted class labels.
#' @return Returns the j-score of the clustering assignment.
#'
#' @example
#' truth=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)
#' pred= c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)
#' j= jscore(truth, pred)
#' @export
jscore<- function(truth, pred){
  if(length(truth)==length(pred)){
  est.mat<-c()
  for(i in unique(pred)){est.mat<- rbind(est.mat, pred==i)}
  est.num<- apply(est.mat,1, sum)

  truth.mat<-c()
  for(i in unique(truth)){truth.mat<- rbind(truth.mat, truth==i)}
  truth.num<- apply(truth.mat,1, sum)

  est.mat.acc<-matrix(nrow = dim(est.mat)[1], ncol=dim(truth.mat)[1])
  for(i in 1:dim(est.mat)[1]){
    for(j in 1:dim(truth.mat)[1]){
      est.mat.acc[i,j]<-sum(est.mat[i,]&truth.mat[j,])/
        sum(est.mat[i,]|truth.mat[j,])
    }
  }

  M1<-sum(apply(est.mat.acc,1,max)*est.num)/length(pred)
  M1.1<-sum(apply(est.mat.acc,2,max)*truth.num)/length(pred)
  M2<-2*M1*M1.1/(M1+M1.1)
  return('jscore'=M2)
  }
  else{return('Truth and Pred have different lengths.')}
}
