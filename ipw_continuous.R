
library(tidyverse)
library(nnet)

set.seed(1)
n <- 250
W <- matrix(rnorm(n*3), ncol=3)
A <- rnorm(n,-4*W[,1] - 4*W[,2] + 4*W[,3],1)
Y <- A + 4*W[,1] + 4*W[,2] + 4*W[,3] + 4*W[,3]^2 + rnorm(n)

data_1 <- tibble(Y,A,W1=W[,1],W2=W[,2],W3=W[,3])

plot_dat <- tibble(A,Y)

ggplot(plot_dat) + geom_point(aes(A,Y))


## crude model

lm(Y~A,data=data_1)

## correct model

lm(Y~A+W1+W2+W3+I(W3^2),data=data_1)

## construct the IPWs via normal density:
A_mod <- glm(A ~ W1+W2+W3,data=data_1,family=gaussian("identity"))
mu_den_lm <- A_mod$fitted.values
sd_den_lm <- A_mod$deviance/A_mod$df.residual
dens_den <- dnorm(data_1$A,mu_den_lm,sd_den_lm)

A_mod <- glm(A ~ 1,data=data_1,family=gaussian("identity"))
mu_num_lm <- A_mod$fitted.values
sd_num_lm <- A_mod$deviance/A_mod$df.residual
dens_num <- dnorm(data_1$A,mu_num_lm,sd_num_lm)

data_1$sw <- dens_num/dens_den
summary(data_1$sw)

## weighted model, normal pdf

lm(Y~A,data=data_1,weights=sw)

## construct the IPWs via quantile binning:
data_1 <- data_1 %>% mutate(quantileA = as.numeric(cut(A, quantile(A, probs=0:10/10), include.lowest=T)))
modA <- multinom(quantileA ~ W1+W2+W3, data = data_1)

propensA <- predict(modA,type="probs") ## gives a 250 x 10 matrix, have to select observed propens

pA <- NULL
for(i in 1:nrow(data_1)){
    pA <- rbind(pA,propensA[i,data_1$quantileA[i]])
}

data_1$swqb <- as.numeric(.1/pA)

## weighted model, quantile binning
lm(Y~A,data=data_1,weights=swqb)


