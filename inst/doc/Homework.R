## -----------------------------------------------------------------------------
n <- 10000
u <-runif(n)
for(sigma in c(1,5,10,100,1000)){
x<-sqrt(-2*(sigma^2)*log(1-u))#F(x)=1-exp(-x^2/(2*sigma^2)),x>=0
hist(x, prob = TRUE, main= expression(f(x)==x*exp(-x^2/(2*sigma^2))/sigma^2))
}

## -----------------------------------------------------------------------------
n<-1000
u<-runif(n)
for(p1 in c(0.1,0.25,0.5,0.75,0.9,0.95)){
for (i in 1:n) {
  if(u[i]<p1)#we we accept x
    x[i]<-rnorm(1,0,1)else
      x[i]<-rnorm(1,3,1)
    x<-x[1:n]
  }
hist(x,prob=TRUE)
}



## -----------------------------------------------------------------------------
n<-6
d<-3
T<- matrix(rnorm(n*d), nrow = n, ncol = d)
sigma<-cor(T)
Q<-chol(sigma)
Z<-matrix(rnorm(n*d), nrow = n, ncol = d)
mu<-c(0,0,0)
X <- Z %*% Q + matrix(mu, n, d, byrow = TRUE)
#X is an n×d data matrix of a random sample from a Nd(μ,Σ)
M<-t(X)%*%X#Denote M ∼ Wd(Σ,n)
M

## -----------------------------------------------------------------------------
m <- 1e4
x <- runif(m, min=0, max=pi/3)
theat.hat<-mean(sin(x))*pi/3
print(c(theat.hat,-cos(pi/3)+1))
m <- 1e6
x <- runif(m, min=0, max=pi/3)
theat.hat<-mean(sin(x))*pi/3
print(c(theat.hat,-cos(pi/3)+1))


## -----------------------------------------------------------------------------
n<-1e4
x<-runif(n,0,1)
theta.hat<-mean(exp(-x)/(1+x^2))
theta.hat
thetas<-exp(-x)/(1+x^2)
sd(thetas)
se<-sd(thetas)/sqrt(n)
se
round(c(theta.hat-1.96*se,theta.hat+1.96*se),3)#The 95% CI of θ (CLT)
#Antithetic variables
f<-function(u){exp(-u)/(1+u^2)
  }
u <- runif(n)
    MC1<-f(u)
    MC2<-f(1-u) # MC1 and MC2 are negatively correlated
    for (i in 1:n) {
   MC<-numeric(n)   
    
    MC[i]<-(MC1[i]+MC2[i])/2
    
    }
  c(sd(MC),sd(thetas),(sd(thetas)-sd(MC))/sd(thetas))


## -----------------------------------------------------------------------------
M <- 1e6 #number of replicates
k <- 50 #number of strata
r <- M / k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
T1<-numeric(k)

estimates <- matrix(0, N, 3)
g<-function(x){
   exp(-x)/(1+x^2)
}
f<-function(x){
  exp(-x)/(1-exp(-1)) 
}
for (i in 1:N) {
  estimates[i, 1] <- mean(g(runif(M)))
  for (j in 1:k) {
  T2[j] <- mean(g(runif(M/k, (j-1)/k, j/k))) 
  # T2 is the stratified sampling
  
  
  #Stratified Importance Sampling
  u<-runif(M/k,(j-1)/k, j/k) 
  t<- -log(1-(u+j-1)*(1-exp(-1))/k)
  T1[j]<-(mean(g(t)/f(t)/k))
  }
  estimates[i, 2] <- mean(T2)
  estimates[i, 3] <- sum(T1)
}


#Example 5.10
u <- runif(M) #f inverse transform method 
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat <- mean(fg)
se<-var(fg)





apply(estimates, 2, mean)
apply(estimates, 2, var)
c(se,apply(estimates, 2, var))


## -----------------------------------------------------------------------------
n<-20
alpha <- 0.05
m<-10000
p<-numeric(m)
for (i in 1:m) {
x<-rchisq(n,2)
UCL<-t.test(x, alternative="two.sided", mu = 2)
p[i] <-UCL$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m) 
print(c(p.hat, se.hat))
1-p.hat


## -----------------------------------------------------------------------------
LN<-c(10,20,50,100,300,500,1000)
quant<-c(0.025,0.05,0.95,0.975)
m<-1e4
sk<-numeric(m)#computes the sample skewness coeff.

f<-function(x){#Normal distribution density
  exp(-(x^2)/2)/sqrt(2*pi)
}
est1<-matrix(0,length(LN),4,dimnames=list(c(LN),c(quant)))
est2<-matrix(0,length(LN),4,dimnames=list(c(LN),c(quant)))
est3<-matrix(0,length(LN),4,dimnames=list(c(LN),c(quant)))

for (j in 1:length(LN) ) {
  n<-LN[j]
  

for (i in 1:m) {
  x<-rnorm(n)
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<- mean((x-xbar)^2)
  sk[i]<-m3/m2^1.5
}

  cv<- qnorm(c(0.025,0.05,0.95,0.975), 0, 6/n)
 est1[j,] <- c(quantile(sk,c(0.025,0.05,0.95,0.975)))
 est2[j,] <- c(cv)
 
}



for (i in 1:length(LN) ){#Compute the standard error of the estimates
  n<-LN[i]
  for (j in 1:4) {
 quant<-c(0.025,0.05,0.95,0.975)   
  est3[i,j] <-  quant[j]*(1-quant[j])/n/(f(est1[i,j])^2)
    
  }
}

est1 #The estimates of quantiles by Monte Carlo experiment.
est2#The quantiles from  √b1 ≈ N (0, 6/n)
est3#  Standard error of the estimates

## -----------------------------------------------------------------------------

sk <- function(x) {
#computes the sample skewness coeff. 
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
alpha <- 0.05
n <- 100
m <- 2500
epsi <- c(seq(0, .15, .01), seq(.15, 1, .04))
N <- length(epsi)
power1 <-numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon 
  pe <- epsi[j]
skt <-numeric(m)
for (i in 1:m) { #for each replicate
sigm <- sample(c(1, 20), replace = TRUE, size = n, prob = c(1-pe, pe))
x <- rbeta(n, sigm, sigm)

skt[i] <- as.integer(abs(sk(x)) >= cv)
}
       
        power1[j] <- mean(skt)
}
#plot power vs epsi 
plot(epsi, power1, type = "b",
xlab = bquote(epsi), ylim = c(0,1)) 
abline(h =0.1, lty = 3)
se <- sqrt(power1 * (1-power1) / m) 
lines(epsi, power1+se, lty = 3)
lines(epsi, power1-se, lty = 3)


##t(v)
pwr2<-numeric(N)
for (j in 1:N) { #for each epsi
 pe <- epsi[j]
skt1 <-numeric(m)
for (i in 1:m) { #for each replicate
tfre <- sample(c(3, 20), replace = TRUE, size = n, prob = c(1-pe, pe))
x <- rt(n, tfre)

skt1[i] <- as.integer(abs(sk(x)) >= cv)
}
       
        pwr2[j] <- mean(skt1)
}

plot(epsi, pwr2, type = "b",
xlab = bquote(epsi), ylim = c(0,1)) 
abline(h = 0.1, lty = 3)
se <- sqrt(pwr2 * (1-pwr2) / m) #add standard errors
lines(epsi, pwr2+se, lty = 3)
lines(epsi, pwr2-se, lty = 3)


## -----------------------------------------------------------------------------
#sample from χ2(1)
#  H0:  u=1
u<-1
n<-10000
m<-100
alpha<-0.05  ##	confidence level of the interval.
qt1<-qt(1-alpha/2,m-1)
pow<-numeric(n)
for (i in 1:n) {
 x<-rchisq(m,1)
 Test<-sqrt(m)*(mean(x)-u)/sd(x)
 pow[i]<-as.integer(abs(Test)>qt1)

}
mean(pow)
#sample from U(0,2)
#  H0 :  u=1

pow1<-numeric(n)
for (i in 1:n) {
  x<-runif(m,0,2)
 Test<-sqrt(m)*(mean(x)-u)/sd(x)
 pow1[i]<-as.integer(abs(Test)>qt1)

  
}
mean(pow1)

#sample from Exponential(1)
## H0  :  u=1
pow2<-numeric(n)
for (i in 1:n) {
  x<-rexp(m,1)
  Test<-sqrt(m)*(mean(x)-u)/sd(x)
 pow2[i]<-as.integer(abs(Test)>qt1)


}
mean(pow2)


## -----------------------------------------------------------------------------
## Z-test
n<-10000
m<-100
alpha<-0.05
p1<-0.651
p2<-0.676
qt1<-qnorm(1-alpha/2,0,1)
pow1<-numeric(n)
for (i in 1:n){
  x<-rbinom(m,1,p1)
  pow1[i]<-as.integer(mean(x)<(p2-qt1*sqrt(p2*(1-p2)/m)))+as.integer(mean(x)>(p2+qt1*sqrt(p2*(1-p2)/m)))
}  
mean(pow1)# approximately equal to0.09

## Two sample test
pow3<-numeric(n)
tq1<-qt(1-alpha/2,2*m-1)
for (i in 1:n) {
  x<-rbinom(m,1,p1)
  y<-rbinom(m,1,p2)
  sw<-(m-1)*(var(x)+var(y))/(m+m-2)
  Test<-(mean(x)-mean(y))/sqrt(sw*(1/m+1/m))
  pow3[i]<-as.integer(Test>tq1)
}
mean(pow3)
m1<-100
m2<-50
pow31<-numeric(n)
for (i in 1:n) {
  x<-rbinom(m1,1,p1)
  y<-rbinom(m2,1,p2)
  sw<-((m1-1)*var(x)+(m2-1)*var(y))/(m1+m2-2)
  Test<-(mean(x)-mean(y))/sqrt(sw*(1/m1+1/m2))
  pow31[i]<-as.integer(Test>tq1)
}
mean(pow31)
## Paired-test
tq2<-qt(1-alpha/2,2*m-1)
pow4<-numeric(n)
for (i in 1:n) {
  x<-rbinom(m,1,p1)
  y<-rbinom(m,1,p2)
  Test<-(mean(x)-mean(y))/sqrt((var(x)+var(y))/2/(n/2))
  pow4[i]<-as.integer(Test>tq2)
}
mean(pow4)
##  Z-test  compare with above
n<-10000
m<-100
alpha<-0.05
p1<-0.676
p2<-0.651
qt1<-qnorm(1-alpha/2,0,1)
pow2<-numeric(n)
for (i in 1:n){
  x<-rbinom(m,1,p1)
  pow2[i]<-as.integer(mean(x)<(p2-qt1*sqrt(p2*(1-p2)/m)))+as.integer(mean(x)>(p2+qt1*sqrt(p2*(1-p2)/m)))
}  
mean(pow2)# approximately equal to0.07


## -----------------------------------------------------------------------------
library(bootstrap)
data(scor)
pairs(scor,panel = points,main="panel display",col = c("red", "green3", "blue","yellow","dark red"))

cor(scor)## Correlation matrix
n<-length(scor$mec)
p12head<-p34head<-p35head<-p45head<-numeric(n)
B<-1000

for (b in 1:B) {
  MEC<-sample(scor$mec,n,replace=TRUE)
  VEC<-sample(scor$vec,n,replace=TRUE)
  ALG<-sample(scor$alg,n,replace=TRUE)
  ANA<-sample(scor$ana,n,replace=TRUE)
  STA<-sample(scor$sta,n,replace=TRUE)
  
  SCOR<-cbind(MEC,VEC,ALG,ANA,STA)
  p12head[b]<-cor(SCOR)[1,2]
  p34head[b]<-cor(SCOR)[3,4]
  p35head[b]<-cor(SCOR)[3,5]
  p45head[b]<-cor(SCOR)[4,5]
}
se12<-sd(p12head)
se34<-sd(p34head)
se35<-sd(p35head)
se45<-sd(p45head)
c(se12,se34,se35,se45)



## -----------------------------------------------------------------------------
library(boot)
sk <- function(data,ind) {
  x<-data[ind]
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

n<-50

mu1<-0
mu2<-4/sqrt(10)#THe skewness of χ2(5) distributions
B<-1000
ci.norm1<-ci.basic1<-ci.perc1<-matrix(0,B,2)
ci.norm2<-ci.basic2<-ci.perc2<-matrix(0,B,2)
for (i in 1:B) {
  X<-rnorm(n)
  Y<-rchisq(n,5)
  boot.obj1<-boot(X,sk,R=2000)
  boot.obj2<-boot(Y,sk,R=2000)
  CI1<-boot.ci(boot.obj1,conf = 0.95,type =c ("norm","basic","perc"))
  CI2<-boot.ci(boot.obj2,conf = 0.95,type = c("norm","basic","perc"))
ci.norm1[i,]<-CI1$norm[2:3]
ci.basic1[i,]<-CI1$basic[4:5]
ci.perc1[i,]<-CI1$percent[4:5]
  
ci.norm2[i,]<-CI2$norm[2:3]
ci.basic2[i,]<-CI2$basic[4:5]
ci.perc2[i,]<-CI2$percent[4:5]


}

Ncix1<-mean(ci.norm1[,2]<=mu1)
Bcix1<-mean(ci.basic1[,2]<=mu1 )
Pcix1<-mean(ci.perc1[,2]<=mu1)
Ncix2<-mean(ci.norm1[,1]>=mu1)
Bcix2<-mean(ci.basic1[,1]>=mu1 )
Pcix2<-mean(ci.perc1[,1]>=mu1)
Ncix<-Ncix1+Ncix2
Bcix<-Bcix1+Bcix2
Pcix<-Pcix1+Pcix2

Nciy1<-mean(ci.norm2[,2]<=mu2)
Bciy1<-mean(ci.basic2[,2]<=mu2)
Pciy1<-mean(ci.perc2[,2]<=mu2)
Nciy2<-mean(ci.norm2[,1]>=mu2)
Bciy2<-mean(ci.basic2[,1]>=mu2 )
Pciy2<-mean(ci.perc2[,1]>=mu2)
Nciy<-Nciy1+Nciy2
Bciy<-Bciy1+Bciy2
Pciy<-Pciy1+Pciy2



rnames<-c("right","left","sum")
cnames<-c("Norm","Basci","Perc")

matrix(c(Ncix1,Ncix2,Ncix,Bcix1,Bcix2,Bcix,Pcix1,Pcix2,Pcix),3,3,dimnames=list(rnames,cnames))
##"Right"is meaning the confidence intervals miss on the right.
matrix(c(Nciy1,Nciy2,Nciy,Bciy1,Bciy2,Bciy,Pciy1,Pciy2,Pciy),3,3,dimnames=list(rnames,cnames))
cp1<-c(1-Ncix,1-Bcix,1-Pcix)
cp2<-c(1-Nciy,1-Bciy,1-Pciy)

list(cp1,cp2)#Coverage probabilities.


## -----------------------------------------------------------------------------
library("bootstrap")
data("scor")
n<-length(scor$mec)
values<-eigen(cov(scor))$values
theats<-values[1]/sum(values)
theta.jack<-numeric(n)
for (i in 1:n) {
  Scor<-scor[-i,]
  val<-eigen(cov(Scor))$values
  theta.jack[i]<-val[1]/sum(val)
}
bias <- (n - 1) * (mean(theta.jack) - theats)
se.jack <- sqrt((n-1)*mean((theta.jack-theats)^2))
c(bias,se.jack)

## -----------------------------------------------------------------------------
library(DAAG)
data(ironslag)
n <- length(ironslag$magnetic)
magnetic<-ironslag$magnetic
chemical<-ironslag$chemical

a <- seq(10, 40, 0.01)
L4<-lm(magnetic ~ chemical+I(chemical^2)+I(chemical^3))$coef

plot(chemical, magnetic, main="Cubic polynomial", pch=16)
yhat<-L4[1] +L4[2]*a+L4[3]*a^2+L4[4]*a^3 
lines(a, yhat,lwd=2)

e1 <- e2 <- e3 <- e4 <- numeric(n)

for (i in 1:n) {
        y <- magnetic[-i]
        x <- chemical[-i]
        
        #For linear model.
        J1 <- lm(y ~ x)
        yhat1 <- J1$coef[1] + J1$coef[2] * chemical[i] 
        e1[i] <- magnetic[i] - yhat1

        ##For quadratic model.
        J2 <- lm(y ~ x + I(x^2))
        yhat2 <- J2$coef[1] + J2$coef[2] * chemical[i] +J2$coef[3] * chemical[i]^2
        e2[i] <- magnetic[i] -yhat2
         
        ##For exponential  model. 
        J3<-lm(log(y)~x)
        logyhat<-J3$coef[1]+J3$coeff[2]*chemical[i]
        yhat3<-exp(logyhat)
        e3[i]<-magnetic[i]-yhat3

        
         ##For cubic polynomial model.
        J4<-lm(y~ x+I(x^2)+I(x^3))$coef
        yhat4 <- J4[1] + J4[2] * chemical[i] +
 J4[3] * chemical[i]^2+J4[4]*chemical[i]^3
        e4[i]<-magnetic[i]-yhat4
}
SStot<-sum((magnetic-mean(magnetic))^2)
SSres1<-sum(e1^2)
SSres2<-sum(e2^2)
SSres3<-sum(e3^2)
SSres4<-sum(e4^2)
R1<-1-SSres1/SStot
R2<-1-SSres2/SStot
R3<-1-SSres3/SStot
R4<-1-SSres4/SStot
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
c(R1,R2,R3,R4)

L2 <- lm(magnetic ~ chemical + I(chemical^2))

plot(chemical, magnetic, main="Quadratic", pch=16) 
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2

lines(a, yhat2, lwd=2)
L2
##The fitted regression equation for Model 2 is 
#Yˆ = 24.49262 − 1.39334X + 0.05452X 
plot(L2)

## -----------------------------------------------------------------------------
## Fixed sample size with different the maximum number of extreme points.
n<-50
m<-60
N<-999
x<-rnorm(n)
y<-rnorm(m)
ep<-5:20
power1<-numeric(16)
count5test <- function(x, y,ep) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y)) 
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > ep))
}
for (j in 1:16) {
expoint<-ep[j]
t0<-count5test(x,y,expoint)
z<-c(x,y)
L<-length(z)
K<-1:L
p<-numeric(N)
for (i in 1:N) {
  k<-sample(K,n,replace = FALSE)
  x1<-z[k]
  y1<-z[-k]
  p[i]<-count5test(x1,y1,expoint)
}

power1[j]<-mean(c(t0,p))
}


## Fixed maximum number of extreme points with different sample size
n1<-c(seq(40,300,20))
m1<-c(seq(60,320,20))
N1<-length(n1)
power2<-numeric(N1)
for (j in 1:N1) {
n2<-n1[j]
m2<-m1[j]
x<-rnorm(n2)
y<-rnorm(m2)
t0<-count5test(x,y,5)
z<-c(x,y)
L<-length(z)
K<-1:L
p<-numeric(N)
for (i in 1:N) {
  k<-sample(K,n,replace = FALSE)
  x1<-z[k]
  y1<-z[-k]
  p[i]<-count5test(x1,y1,5)
  
}
power2[j]<-mean(c(t0,p))
}
plot(ep,power1,xlab = "maximum number of extreme points")

plot(n1,power2,main = "Fixed maximum number of extreme points with different sample size",)



## -----------------------------------------------------------------------------
#First,generate sample from Model1 and Model2
library("ks")
library("boot")
library(Ball)



Akl <- function(x) {
d <- as.matrix(dist(x)) 
m <- rowMeans(d)
M <- mean(d)
a <- sweep(d, 1, m)
b <- sweep(a, 2, m) 
return(b + M)
}

dCov <- function(x, y) { 
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)
m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")

A <- Akl(x)
B <- Akl(y)
dCov <- sqrt(mean(A * B))
dCov
}

ndCov2 <- function(z, ix, dims) {
#dims contains dimensions of x and y
p <- dims[1]
q1 <- dims[2] + 1
d <- p + dims[2]
x <- z[ , 1:p] #leave x as is
y <- z[ix, q1:d] #permute rows of y 
return(nrow(z) * dCov(x, y)^2)
}
alpha<-0.05
B<-30##The larger the sample value is, the better the effect will be.
n<-seq(30,100,10)##The larger the sample value is, the better the effect will be.
K<-length(n)

Bpb1<-Bpb2<-Bpb3<-Bpd1<-Bpd2<-Bpd3<-numeric(B)
Kpb1<-Kpb2<-Kpb3<-Kpd1<-Kpd2<-Kpd3<-numeric(K)
for (k in 1:K) {
  n1<-n[k]

for (b in 1:B) {
x1<-rmvnorm.mixt(n1,mus = c(0,0),Sigmas = diag(2))
e1<-rmvnorm.mixt(n1,mus = c(0,0),Sigmas = diag(2))
y1<-x1/4+e1
y2<-(x1/4) * e1
z1<-cbind(x1,y1)
z2<-cbind(x1,y2)
z3<-cbind(y1,y2)
boot.obj1 <- boot(data = z1, statistic = ndCov2, R = 100,
sim = "permutation", dims = c(2, 2))##The larger the sample value is, the better the effect will be.
tb1 <- c(boot.obj1$t0, boot.obj1$t)
Bpd1[b]<-mean(tb1>= boot.obj1$t0)
p.ball1 <- bcov.test(z1[,1:2],z1[,3:4],R=100)$p.value
Bpb1[b]<-p.ball1
boot.obj2 <- boot(data = z2, statistic = ndCov2, R = 100,
sim = "permutation", dims = c(2, 2))
tb2 <- c(boot.obj2$t0, boot.obj2$t)
Bpd2[b]<-mean(tb2>= boot.obj2$t0)
p.ball2 <- bcov.test(z2[,1:2],z2[,3:4],R=100)$p.value
Bpb2[b]<-p.ball2
boot.obj3 <- boot(data = z3, statistic = ndCov2, R = 100,
sim = "permutation", dims = c(2, 2))
tb3 <- c(boot.obj3$t0, boot.obj3$t)
Bpd3[b]<-mean(tb3>= boot.obj3$t0)
p.ball3 <- bcov.test(z3[,1:2],z3[,3:4],R=100)$p.value
Bpb3[b]<-p.ball3
}
 Kpd1[k]<- mean(Bpd1<alpha)
 Kpd2[k]<- mean(Bpd2<alpha)
 Kpd3[k]<- mean(Bpd3<alpha)
 Kpb1[k]<- mean(Bpb1<alpha)
 Kpb2[k]<- mean(Bpb1<alpha)
 Kpb3[k]<- mean(Bpb1<alpha)
}
##Test the independant of x1 and y1 .
plot(n,Kpd1,type = "l")
lines(n,Kpb1,col=2,type = "l")

##Test the independant of x1 and y2 .
plot(n,Kpd2,type = "l")
lines(n,Kpb2,type="l" ,col=2)

##Test the independant of y1 and y1 .
plot(n,Kpd3,type = "l")
lines(n,Kpb3,type="l",col=2)
##The larger the sample value is, the better the effect will be.

## -----------------------------------------------------------------------------
f <- function(x) {
  1/2*exp(-abs(x))
}
#g <- function(x,sigma) {dnorm(x,sigma)}
rw.Metropolis<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
   k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= ( f(y) / f(x[i-1]) ))
     x[i] <- y 
    else {
     x[i] <- x[i-1]
     k <- k + 1
   }
  }
return(list(x=x, k=k))
}
sigma<-c(0.05,0.1,0.5,2,5,10,20)
B<-length(sigma)
N<-3000
x0<-20
par(mfrow=c(3,2))
for (i in 1:B) {
Rw<-rw.Metropolis(sigma[i],x0,N)
k<-Rw$k
par(mar=c(1,1,1,1))
plot(1:N,Rw$x,type="l",xlab = paste("sigma=",sigma[i]),main=paste("acceptance rates=",(N-k)/N))
   abline(h=c(log(.05/2),-log(.05/2)))
   
}



## -----------------------------------------------------------------------------
x<-1:50
y1<-exp(x)
y2<-log(x)
isTRUE(all.equal(log(y1), exp(y2)))

isTRUE(all.equal(log(y1),x ))
isTRUE(all.equal(x, exp(y2)))
identical(log(y1),exp(y2))
isTRUE(all.equal(log(y1) -exp(y2),0))
isTRUE(all.equal(log(y1)-x ,0))
isTRUE(all.equal(x-exp(y2),0))

## -----------------------------------------------------------------------------
K1<-c(4:25,100,500,1000)
B1<-length(K1)
eps <- .Machine$double.eps^0.25
root1<-numeric(B1)

for (b in 1:B1) {
  k<-K1[b]
  b0<-eps
b1 <- sqrt(k) - eps



 root1[b]<-uniroot(function(a){
  z1<-sqrt(a^2*(k-1)/(k-a^2))
  z2<-sqrt(a^2*k/(k+1-a^2))
  gx<-pt(z1,k-1,lower.tail = FALSE)-pt(z2,k,lower.tail = FALSE)
  return(gx)
},interval = c(b0, b1))$root
}
root1



####################################################################################################################
K<-c(4:25,100)
B<-length(K)
root2<-numeric(B)

for (b in 1:B) {
  k<-K[b]
  f<-function(a){
  z1<-sqrt(a^2*(k-1)/(k-a^2))
  z2<-sqrt(a^2*k/(k+1-a^2))
  f1<-function(x){
    gamma(k/2)*(1+(z1*x)^2 /(k-1))^(-k/2) *z1/gamma((k-1)/2)/sqrt((k-1)*pi)-gamma((k+1)/2)*(1+(z2*x)^2 /k)^(-(k+1)/2) *z2/gamma(k/2)/sqrt(k*pi) }
  return(integrate(f1,eps,1)$value)
  }
  
 it <- 0
 eps <- .Machine$double.eps^0.25
 b0 <- 0
 b1 <- sqrt(k) - eps
 r <- seq(b0, b1, length=3)
 y <- c(f(r[1]), f(r[2]), f(r[3]))
 if (y[1] * y[3] > 0)
stop("f does not have opposite sign at endpoints")
 while(it < 1000 && abs(y[2]) > eps) { 
   it <- it + 1
if (y[1]*y[2] < 0) { 
  r[3] <- r[2] 
  y[3] <- y[2]
        } 
 else {
   r[1] <- r[2]
  y[1] <- y[2] 
}
r[2] <- (r[1] + r[3]) / 2
y[2] <- f(r[2]) 

}



 root2[b]<-r[2] 
}
root2



####################################################################################################################


## -----------------------------------------------------------------------------
N<-30
n<-163 
p<-0.1
q<-0.1
pm<-numeric(N)
qm<-numeric(N)
for (i in 1:N) {
## E-step
  Naa<-as.integer(n*p^2)## Naa and Nbb are integer.
  Nbb<-as.integer(n*q^2)
  
## M-step
  ## By MLE ,we have
  ## (232-$n_{BB}$)p + (98 + $n_{BB}$) = 98+$n_{AA}$
  ## (94+$n_{BB}$)p + (228-$n_{AA}$)q = 94+$n_{BB}$
  A<-matrix(c(232-Nbb,94+Naa,98+Naa,228-Naa),2,2)
  Ap<-matrix(c(98+Naa,94+Nbb,98+Naa,228-Naa),2,2)
  Aq<-matrix(c(232-Nbb,94+Naa,98+Naa,94+Nbb),2,2)
  p<-det(Ap)/det(A)##Cramer's Law
  q<-det(Aq)/det(A)
  pm[i]<-p
  qm[i]<-q
}
p
q
plot(1:N,pm,main = "p value in M-steps")
plot(1:N,qm,main = "q value in M-steps")




## -----------------------------------------------------------------------------
data("mtcars")
mpg<-mtcars$mpg
disp<-mtcars$disp
wt<-mtcars$wt
formulas <- list(

mpg ~ disp,

mpg ~ I(1 / disp),

mpg~disp+wt, 

mpg~I(1/disp)+wt
)
A1<-lapply(formulas,lm)
A1

## -----------------------------------------------------------------------------
data("mtcars")
bootstraps <- lapply(1:10, function(i) { 

rows <- sample(1:nrow(mtcars), rep = TRUE) 

mtcars[rows, ]

})
A2<-lapply(bootstraps,lm, formula=mpg~disp)
A2

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
lapply(A1, rsq)
lapply(A2, rsq)

## -----------------------------------------------------------------------------
 trials <- replicate( 100,
t.test(rpois(10, 10), rpois(7, 10)),
         simplify = FALSE
       )
sapply(trials,function(x) {
  unlist(x)[3]
})
##Extra challenge: get rid of the anonymous function by using [[ directly.

sapply(trials, '[[', 3)



## -----------------------------------------------------------------------------
mcsapply<-function(x,fun){
library(parallel)
  cl<- makeCluster(3)
f1<-  parSapply(cl,x,fun)
stopCluster(cl)
return(f1)
}


trials <- replicate( 10000,
t.test(rnorm(20, 10,11), rpois(7, 10)),
         simplify = FALSE
       )
list(system.time(mcsapply(trials,function(x) {
  unlist(x)[3]
})),system.time(sapply(trials,function(x) {
  unlist(x)[3]
})))
## The increased computing speed is very obvious



## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

####################################################################
# C++ function.

cppFunction(
  'NumericVector rw_C (double sigma,double x0,int N){
  NumericVector x(N); 
  x[1] = x0; 
  NumericVector u = runif(N); 
  double k = 0;
  for (int i = 2; i < N+1 ;i++){ 
    double y = as<double>(rnorm(1, x[i-1], sigma));
    double t = exp(-abs(y))/exp(-abs(x[i-1]));
    if (u[i] <= t){
      x[i] = y; 
      k = k + 1;
    }
    else{ 
       x[i] = x[i-1];} 
  };
  return x;}')


######################################################################################################################################
f <- function(x) {
  1/2*exp(-abs(x))
}
#g <- function(x,sigma) {dnorm(x,sigma)}
rw_R<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
   k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= ( f(y) / f(x[i-1]) ))
     x[i] <- y 
    else {
     x[i] <- x[i-1]
     k <- k + 1
   }
  }
return(list(x=x, k=k))
}

#######################################################################################################################################

sigma<-c(0.05,0.1,0.5,2,5,10,20,50,100)
B<-length(sigma)
N<-10000
x0<-20

for (i in 1:B) {
  
Rw<-rw_R(sigma[i],x0,N)
rcpp<-rw_C(sigma[i],x0,N)
ts <- microbenchmark(Rw=rw_R(sigma[i],x0,N), rcpp=rw_C(sigma[i],x0,N))

qqplot(Rw$x,rcpp,ylab="C++ function",xlab = "R function",main= paste("sigma=",sigma[i]))
  abline(a = 0 , b = 1, col = "red") 
print(summary(ts)[,c(1,3,5,6)])
}


