## -----------------------------------------------------------------------------
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

##__Examle__

data1<-iris[,1:2]
m1<-Fdbscan(data1,0.25,3)
m1
table(m1)
data2<-iris[,1:3]
m2<-Fdbscan(data2,0.3,3)
m2
table(m2)

## -----------------------------------------------------------------------------
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

##__Examle__

data1<-iris[,1:2]
m1<-Fdbscan(data1,0.25,3)
z1<-CLuster(data1,m1,3)
z1
data2<-iris[,1:3]
m2<-Fdbscan(data2,0.3,3)
z2<-CLuster(data2,m2,3)
z2

