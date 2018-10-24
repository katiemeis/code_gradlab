
setwd("~/applied_bayesian/homework2/")

####### QUESTION 1 ########
dat <- read.table("HW2Q1.txt", header=T, sep="\t")
x_Beta <- sum((dat$X - 173)^2)/2
y_Beta <- sum((dat$Y - 173)^2)/2

set.seed(10)
ystar = 180
post_dist <- vector()
for(i in 1:5000){
  sigmaX <- 1/(rgamma(1,15/2,x_Beta))
  sigmaY <- 1/(rgamma(1,15/2,y_Beta))
  x_post <- rnorm(1, 173+(sqrt(sigmaX)/sqrt(sigmaY))*0.5*(ystar-173), sqrt(0.75*sigmaX))
  post_dist <- c(post_dist, x_post)
}

plot(density(post_dist))

summary(post_dist)
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
estimate_mode(post_dist)

quantile(post_dist,c(0.025,0.975))

hist(post_dist, main = "Histogram of x*")
hist(dat$X, main = "Histogram of observed")




####### QUESTION 2 ########
#to get alpha we want summation over y plus 1
set.seed(105)
p0 <- rbeta(2000, 6, 76)
p1 <- rbeta(2000, 4, 78)
post_or <- (p1*(1-p0))/(p0*(1-p1))
plot(density(post_or))
summary(post_or)
quantile(post_or,c(0.025,0.975))
#an odds ratio of < 1 means exposure (drug in this case) 
#is associated with lower odds of death




####### QUESTION 3 ########
set.seed(21295)
samples <- runif(2000,0,1)

x=(-log(1-samples))^2
hist(x)
plot(density(x))





####### QUESTION 4 ########
########################## REJECTION SAMPLING ##################################
dat4 <- read.table("HW2Q4.txt", header=T, sep=" ")
set.seed(100)
xi <- dat4$x
ni <- dat4$n
yi <- dat4$y
beta_vec <- seq(-10, 10, 0.01)

#posterior
fx <- function(y, n, x, beta){
  pi_val = 1/(1+exp(-beta*x))
  lik = (pi_val^y)*((1-pi_val)^(n-y))*choose(n, y)
  prior = ((gamma(5/2)/(sqrt(4*pi)))*((1+((beta^2)/4))^(-5/2)))
  target = prior*prod(lik)
  return(target)
}

#apply to a range of betas
fxx <- vector()
for(i in 1:length(beta_vec)){
  fxx_i <- fx(y=yi,n=ni,x=xi, beta_vec[i])
  fxx <- c(fxx, fxx_i)
}

#plot posterior
plot(beta_vec, fxx)

# proposale envolpe distribution: g(x)= N(mu,sigma2) 
gx= function(x,mu,sigma2) 
{
  target = exp(-1/2*(x-mu)^2/(2*sigma2))/sqrt(2*pi*2*sigma2)
  return(target)
}   

#find where the peak is to get mean for gx
beta_vec[which.max(fxx)] #-.37
gxx= gx(x=beta_vec,-.37,1)
M= max(fxx)/max(gxx)

plot(beta_vec, fxx, type="l")
lines(beta_vec, gxx*M, col='blue')

reject = function(x, M, S=5000)
{
  samples=rep(0,S)
  i=0      # counter for the accepted draws
  count=0  # counter for the NO of proposals/trieds
  while(i<S)
  {
    #draw from proposal
    z = rnorm(1,-.37,1)
    # calculate M*g(x) at z
    Mgx = M*gx(z,-.37,1)
    #calculate f(x) at z
    fx = fx(y=yi,n=ni,x=xi, z)
    u = runif(1)
    if(u<(fx/Mgx)) 
    {
      i=i+1
      samples[i]=z 
    }
    count=count+1
  }
  list(acceptrate=S/count,samples=samples)
}
out <- reject(x=beta_vec, M=M) 
out$acceptrate

hist(out$samples)
plot(density(out$samples))
lines(beta_vec, fxx, type="l", col='red')
mean(out$samples)
sd(out$samples)

########################## IMPORTANCE SAMPLING ##################################
set.seed(10)
xi <- dat4$x
ni <- dat4$n
yi <- dat4$y
beta_vec <- seq(-10, 10, 0.01)

#proposal is standard normal
Up = rnorm(5000,0,1)
px_vals <- dnorm(Up,0,1)

# get the same for our target posterior
fx <- function(y, n, x, beta){
  pi_val = 1/(1+exp(-beta*x))
  lik = (pi_val^y)*((1-pi_val)^(n-y))*choose(n, y)
  prior = ((gamma(5/2)/(sqrt(4*pi)))*((1+((beta^2)/4))^(-5/2)))
  target = prior*prod(lik)
  return(target)
}

fx_vals <- vector()
for(i in 1:length(Up)){
  fxx_i <- fx(y=yi,n=ni,x=xi, Up[i])
  fx_vals <- c(fx_vals, fxx_i)
}

#now we can get the importance ratio
w <- fx_vals/px_vals
mean(w)  
# E(x)
mean(w*Up)/mean(w)

# Sampling IR
rs = sample(Up, 1000, w, replace=FALSE)
hist(rs, breaks = 100)
plot(density(rs), col='blue')
mean(rs)
sd(rs)

