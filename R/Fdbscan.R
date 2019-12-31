#' @title A Fdbscan sampler using R
#' @description Fdbscan is a method for preliminary clustering of given data based on DBSCAN algorithm
#' @param data a data matrix or a dist object
#' @param eps size of the epsilon neighborhood
#' @param M number of minimum points in the eps region (for core points) 
#' @return An object of class 'dbscan_fast'
#' @examples
#' \dontrun{
#' data<-iris[,1:2]
#' m<-Fdbscan(data,0.25,3)
#' m
#' }
#' @export
Fdbscan<-function(data,eps,M){
  N<-length(data[,1])
  Eu<-as.matrix(dist(data,diag = TRUE))
  colnames(Eu)<-NULL
  rownames(Eu)<-NULL
  I<-c(1:N)
  k<-1
  m<-numeric(N)
  for (i in 1:N) { m[i]<-0
  }
  while(length(I) > 0){
    i <- sample(I,1)
    I <- I[- which( I==i) ]
    T <- c(which(Eu[i,] <= eps ))
    if(m[i] == 0)
      if(length(T) < M)
        m[i] <- -1
    else
      m[i] <- k
    while(length(T) > 0 ){
      j <- sample(T,1)
      T <- T[ - which(T == j)]
      if(m[j]!= 0 & m[j] != -1)
        m[j] <- m[j]
      else{ m[j]<- k
      if(length(c(which(Eu[j,] > eps ))) >= M)
        T <-union(T, c(which(Eu[j,] <= eps )))
      else
        m[j]<- m[j] }
    }
    k<-k+1
  }
  return(m) }