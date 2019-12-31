#' @title A CLuster sampler using R
#' @description Further processing of data based on Fdbscan method
#' @param data a data matrix or a dist object
#' @param m the result from Fdbscan(data,eps,M)
#' @param M number of minimum points in the eps region (for core points) 
#' @return a list of Classification of data
#' @examples
#' \dontrun{
#' data<-iris[,1:2]
#' m<-Fdbscan(data,0.25,3)
#' z<-CLuster(data,m,3)
#' z
#' }
#' @export
CLuster<-function(data,m,M){
  colnames(data)<-NULL
  rownames(data)<-NULL
  table1<-table(m)
  N_1<-as.numeric(rownames(table1))
  rownames(table1)<-NULL
  Ind1<-which(table1>=M)
  K<-length(Ind1)
  CLuster<-list()
  for (j in 1:K) {
    CLuster[[j]]<-data[which(m==N_1[Ind1[j]],arr.ind = TRUE),]
  }
  return(CLuster)
}