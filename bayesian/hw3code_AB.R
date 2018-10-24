
################# CODE FOR HOMEWORK 3 #####################
setwd("~/applied_bayesian/homework3/")


################## Problem 2 ####################
############### Gibbs Sampler ###################
set.seed(10)
N = 5000
#vectors to hold ys
y1 = numeric(N+1)
y2 = numeric(N+1)
y3 = numeric(N+1)
#initial values for ys
y1[1] = 0
y2[1] = 0
y3[1] = 0
#rates
y1_rate = 1-.5*y2[1]-2*y3[1]
y2_rate = 1-.5*y1[1]-3.5*y3[1]
y3_rate = 1-2*y1[1]-3.5*y2[1]
#loop for gibbs sampler
for(i in 1:N){
  check_const1 = T
  while(check_const1){
    #draw y1 candidate first
    y1_rate = 1-.5*y2[i]-2*y3[i]
    y1[i+1] = rexp(1, y1_rate)
    #check the other rates to make sure they meet constraint
    #y2 and y3 not yet updated
    y2_rate = 1-.5*y1[i+1]-3.5*y3[i]; y3_rate = 1-2*y1[i+1]-3.5*y2[i]
    #if all greater than 0 then good to move on, if not stay in loop until false
    if(y1_rate > 0 & y2_rate > 0 & y3_rate > 0){check_const1 = F}
  }
  
  check_const2 = T
  while(check_const2){
    #draw y2 candidate
    y2[i+1] = rexp(1, y2_rate)
    #check other rates vs constraints
    #y1 and y2 updated, y3 not yet updated
    y1_rate = 1-.5*y2[i+1]-2*y3[i]; y3_rate = 1-2*y1[i+1]-3.5*y2[i+1]
    #if all > 0 then good to move on
    if(y1_rate > 0 & y2_rate > 0 & y3_rate > 0){check_const2 = F}
  }
  
  check_const3 = T
  while(check_const3){
    #draw y3 candidate
    y3[i+1] = rexp(1, y3_rate)
    #Check other rates vs constraints
    #now all updated ys
    y1_rate = 1-.5*y2[i+1]-2*y3[i+1]; y2_rate = 1-.5*y1[i+1]-3.5*y3[i+1]
    #if all > 0 then good to move on, done with this iteration
    if(y1_rate > 0 & y2_rate > 0 & y3_rate > 0){check_const3 = F}
  }
}


#trace
plot(1:length(y1),y1, type="l", xlab= "Iteration")
plot(1:length(y2),y2, type="l", xlab= "Iteration")
plot(1:length(y3),y3, type="l", xlab= "Iteration")
#autocorrelation
acf(y1)
acf(y2)
acf(y3)

#lets remove the first 100 just to be safe
plot(100:5000, y1[100:5000],type="l", ylab = "y1", xlab= "Iteration")
plot(100:5000, y2[100:5000],type="l", ylab= "y2", xlab = "Iteration")
plot(100:5000, y3[100:5000],type="l", ylab= "y3", xlab = "Iteration")
#based on autocorrelation thin lag 4
pick = seq(1,N-100,4)
new_y1 = y1[101:N][pick]
new_y2 = y2[101:N][pick]
new_y3 = y3[101:N][pick]
#lets now check the acf again to be sure we fixed the autocorrelation 
acf(new_y1)
acf(new_y2)
acf(new_y3)

mean(new_y1)
median(new_y1)
quantile(new_y1, c(0.025, 0.975))

mean(new_y2)
median(new_y2)
quantile(new_y2, c(0.025, 0.975))

mean(new_y3)
median(new_y3)
quantile(new_y3, c(0.025, 0.975))



################## Problem 3 ####################
###### Metropolis Hastings algorithm ############
library(tmvtnorm)
set.seed(10)
N=5000
p1 = numeric(N+1)
p2 = numeric(N+1)
p1[1] = 0.65
p2[1] = 0.65
#jumping distribution g(x)
Sigma = matrix(c(1,0.2,0.2, 1), 2,2)
accept_count = 0

for(i in 1:N){
  #select candidate from jumping distn
  gxpred = rtmvnorm(n=1, mean=c(p1[i],p2[i]), sigma=(0.2^2)*Sigma, lower=c(0,0), upper=c(1,1))
  p1c = gxpred[1,1]; p2c = gxpred[1,2]
  #compute A ratio
  fc = dbeta(p1c,0.5,0.5)*dbeta(p2c,0.5,0.5)*dbinom(10, 15, p1c)*dbinom(8, 12, p2c)
  fcurrent = dbeta(p1[i],0.5,0.5)*dbeta(p2[i],0.5,0.5)*dbinom(10, 15, p1[i])*dbinom(8, 12, p2[i])
  gx_current = dtmvnorm(c(p1[i], p2[i]), mean=c(p1c,p2c), sigma=(0.1^2)*Sigma, lower=c(0,0), upper=c(1,1))
  gx_candidate = dtmvnorm(c(p1c,p2c), mean=c(p1[i],p2[i]), sigma=(0.1^2)*Sigma, lower=c(0,0), upper=c(1,1))
  A = (fc/fcurrent)*(gx_current/gx_candidate)
  #selection of current x or candidate
  if(A > 1){p1[i+1] = p1c; p2[i+1] = p2c; accept_count=accept_count+1}
  if(A < 1){
    prob = rbinom(1, 1, A)
    if(prob == 1){
      p1[i+1] = p1c
      p2[i+1] = p2c
      accept_count = accept_count + 1
    }
    if(prob == 0){
      p1[i+1]=p1[i]
      p2[i+1]=p2[i]
    }
  }
}

#acceptance ratio
accept_count/N

#trace
plot(1:length(p1),p1, type="l", xlab = "Iteration")
plot(1:length(p2),p2, type="l", xlab = "Iteration")
#autocorrelation
acf(p1)
acf(p2)

#burnin 100 for 0.1, thinning 18 for 0.1
#burnin 200 for 0.2, thinning 12 for 0.2
pick = seq(1,N-200,12)
new_p1 = p1[201:N][pick]
new_p2 = p2[201:N][pick]

#autocorrelation double check
acf(new_p1)
acf(new_p2)


mean(new_p1-new_p2)
median(new_p1-new_p2)
sd(new_p1-new_p2)
quantile(new_p1-new_p2, c(0.025, 0.975))
