
###################################################################################

############################### HOMEWORK 1 CODE ###################################

###################################################################################


############### Problem 1 #################

##### Part b #####
#draw 5000 samples from gamma and plot
set.seed(0)
first_set <- rgamma(n=5000, shape=0.01, rate=0.01)
hist(first_set, main = "Histogram of Gamma(0.01,0.01)", xlab = "Samples")
mean(first_set)
sd(first_set)

second_set <- rgamma(n=5000, shape=0.1, rate=0.1) 
hist(second_set, main = "Histogram of Gamma(0.1,0.1)", xlab = "Samples")
mean(second_set)
sd(second_set)

third_set <- rgamma(n=5000, shape=1, rate=1)
hist(third_set, main = "Histogram of Gamma(1,1)", xlab = "Samples")
mean(third_set)
sd(third_set)




############### Problem 2 #################

y = seq(from=-5, to=8, by=0.001)
m1_vec <- vector()
m2_vec <- vector()
m3_vec <- vector()

for(i in y){
  m1_vec <- c(m1_vec, dnorm(i, 0.5, 0.8))
  m2_vec <- c(m2_vec, dnorm(i, 2.5, 1))
  m3_vec <- c(m3_vec, dnorm(i, 3, 0.3))
  
}

mixture <- 0.3*m1_vec + 0.3*m2_vec + 0.4*m3_vec
plot(y, mixture, type="l", ylab="Density", main = "Marginal density of y")




############### Problem 5 #################
yi = c(23, 22, 22.5, 21.5, 24)
mu = seq(from=0, to=50, by=0.01)

#the likelihood is given in the problem: 1/(1+(y-mu)^2)
#we're told mu ~ unif[0,50] so this is 1/50
#posterior is likelihood times prior
posterior_func <- function(yi, mu){return((prod((1/(1+(yi-mu)^2))))/50)}

#loop over mu vecotr to get the posterior at each
post_vec <- vector()
for(mui in mu){
  post_vec <- c(post_vec, posterior_func(yi, mui))
}
plot(mu,post_vec)

#normalized means it sums to 1
#this means we need to divide by the current sum
norm_constant <- sum(post_vec)
norm_post <- post_vec/norm_constant
#lets check that this sums to 1
sum(norm_post)
plot(mu, norm_post)

#plot both
plot(mu, norm_post, type = "l", lty=1, lwd = 2, ylab = "Density")
lines(mu, post_vec, type="l", lty=2, lwd =2)
legend("topright", c("Unnormalized Posterior", "Normalized Posterior"), lty=c(1,2), lwd=c(2,2))


