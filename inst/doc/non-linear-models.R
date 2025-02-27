## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
library(fido)
library(dplyr)
library(tidyr)
library(ggplot2)
data(mallard_family)

# Just take vessel 1
sample.ids <- mallard_family$sample_data[mallard_family$sample_data$Vessel == 1,]
# Just take hourly samples
sample.ids <- sample.ids[(sample.ids$time > "2015-11-20 15:00:00 UTC") & (sample.ids$time < "2015-11-25 16:00:00 UTC"),]

# Subsetting the sample data and OTU data
subset.sample_data <- mallard_family$sample_data[mallard_family$sample_data$X.SampleID %in% sample.ids$X.SampleID,]

subset.otu_table <- mallard_family$otu_table[rownames(mallard_family$otu_table) %in% sample.ids$X.SampleID,]

# Order samples - to make plotting easy later
o <- order(subset.sample_data$time)
subset.otu_table <- subset.otu_table[o,]
subset.sample_data <- subset.sample_data[o,]

# Extract Data / dimensions from Phyloseq object
Y <- t(as(subset.otu_table, "matrix"))
D <- nrow(Y)
N <- nrow(subset.sample_data)

# X in hours
X <- as.numeric(subset.sample_data$time)
X <- t((X-min(X)) / 3600)

## -----------------------------------------------------------------------------
# Specify Priors
Gamma <- function(X) SE(X, sigma=5, rho=10) # Create partial function 
Theta <- function(X) matrix(0, D-1, ncol(X))
upsilon <- D-1+3
Xi <- matrix(.4, D-1, D-1)
diag(Xi) <- 1

# Now fit the model
fit <- fido::basset(Y, X, upsilon, Theta, Gamma, Xi)

## ----fig.height=4, fig.width=6------------------------------------------------
fit.clr <- to_clr(fit)

# Plot Sigma in CLR
plot(fit.clr, par="Sigma", focus.coord=c("clr_seq_6", "clr_seq_5", "clr_seq_2"))

## -----------------------------------------------------------------------------
# predict not just missing days but also forecast into future
X_predict <- t(1:(max(X)))
predicted <- predict(fit.clr, X_predict, jitter=1) 

## ----fig.height=5, fig.width=7------------------------------------------------
family_names <- as(mallard_family$tax_table$Family, "vector")
Y_clr_tidy <- clr_array(Y+0.65, parts = 1) %>%
  gather_array(mean, coord, sample) %>%
  mutate(time = X[1,sample],
         coord = paste0("CLR(", family_names[coord],")"))

predicted_tidy <- gather_array(predicted, val, coord, sample, iter) %>%
  mutate(time = X_predict[1,sample]) %>%
  filter(!is.na(val)) %>%
  group_by(time, coord) %>%
  summarise_posterior(val, na.rm=TRUE) %>%
  ungroup() %>%
  mutate(coord = paste0("CLR(", family_names[coord],")"))

ggplot(predicted_tidy, aes(x = time, y=mean)) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
  geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9)+
  geom_line(color="blue") +
  geom_point(data = Y_clr_tidy, alpha=0.5) +
  facet_wrap(~coord, scales="free_y") +
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45))


## ----message = FALSE, warning = FALSE-----------------------------------------
library(fido)
library(MCMCpack)
library(LaplacesDemon)

set.seed(2024)
D <- 3 # number of taxa
N <- 100 #number of observation
X <- seq(0,10, length.out = N) # specifying the values of X

## -----------------------------------------------------------------------------
# RBF kernel function
rbf <- function(x,l,sigma) {
  K <- outer(x,x, FUN = function(x,y, sigma,l){sigma^2 * exp(- (x - y)^2 / (2 * l^2))}, sigma = sigma, l = l)
  identity_matrix <- diag(ncol(K))
  K  <-  K  + 1e-8 * identity_matrix
  return(K)
}

# Linear kernel function
linear <- function(x,sigma,c) {
  K <- outer(x,x, FUN = function(x,y, sigma,c){ (sigma^2 + (x-c) * (y-c))}, sigma = sigma, c=c)
  identity_matrix <- diag(ncol(K))
  K  <-  K  + 1e-8 * identity_matrix
  return(K)
}

## -----------------------------------------------------------------------------
## Step 1: Sigma ~ IW(upsilon,Xi)
Sigma <- rinvwishart(D-1, 0.5*diag(D-1))

## Step 2a: Lambda_1 ~ GP(0,Sigma,K)
K <- rbf(X,1,1)
Theta <- matrix(0, D-1, N)
x <- matrix(rnorm(N*(D-1)),nrow = D-1)
l1 <- Theta+ t(chol(Sigma))%*%x%*%(chol(K))

## Step 2b: Lambda_2 ~ GP(0, Sigma, K2)
K2 <- linear(X,.1,0)
Theta <- matrix(0, D-1, N)
x <- matrix(rnorm(N*(D-1)),nrow = D-1)
l2 <- Theta+ t(chol(Sigma))%*%x%*%(chol(K2))

l <- l1 + l2

## Step 3: eta ~ N(Lambda,Sigma, I)
eta <- matrix(rnorm(N*(D-1)),nrow = D-1, byrow = TRUE)
eta <- l + t(chol(Sigma))%*%eta

## Step 4: transform to proportions
pai <- t(alrInv(t(eta)))

## Step 5: Simulate Y.
Y <- matrix(0, D, N)
for (i in 1:N) Y[,i] <- rmultinom(1, 5000, prob = pai[,i])

## -----------------------------------------------------------------------------
rbf <- function(x,l,sigma) {
  x <- x[1,]
  K <- outer(x,x, FUN = function(x,y, sigma,l){sigma^2 * exp(- (x - y)^2 / (2 * l^2))}, sigma = sigma, l = l)
  identity_matrix <- diag(ncol(K))
  K  <-  K  + 1e-8 * identity_matrix
  return(K)
}

# Linear kernel function for Gaussian process
linear <- function(x,sigma,c) {
  x <- x[1,]
  K <- outer(x,x, FUN = function(x,y, sigma,c){ (sigma^2 + (x-c) * (y-c))}, sigma = sigma, c=c)
  identity_matrix <- diag(ncol(K))
  K  <-  K  + 1e-8 * identity_matrix
  return(K)
}

Gamma <- list(function(X) rbf(X,1,1), function(X) linear(X,.1,0))

## -----------------------------------------------------------------------------
Theta <- list(function(X) matrix(0, D-1, ncol(X)), function(X) matrix(0, D-1, ncol(X)))

## -----------------------------------------------------------------------------
mmX <- model.matrix(~X-1)
mod <- basset(Y, t(mmX), Theta = Theta, Gamma = Gamma)

## -----------------------------------------------------------------------------
typeof(mod$Lambda)

## -----------------------------------------------------------------------------
## predictions
preds <- predict(mod)

## summarizing predictions over posterior samples
summary_preds <- apply(preds,MARGIN = c(1,2), FUN = mean)

## ----fig.height=4, fig.width=6------------------------------------------------
## plotting
data.frame("points" = c(eta[1,], summary_preds[1,]), "group" = c(rep("Truth", length(eta[1,])), rep("Predicted", length(summary_preds[1,]))), "Index" = c(1:N, 1:N)) %>%
  ggplot(aes(x=Index, y= points, color = group)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = c("red", "black")) +
  theme(legend.title=element_blank()) +
  ylab("Eta")


## -----------------------------------------------------------------------------
Theta <- list(matrix(0, D-1, ncol(mmX)), function(X) matrix(0, D-1, ncol(X)))
Gamma <- list(diag(ncol(mmX)), function(X) rbf(X,1,1))

mmX <- model.matrix(~X-1)
mod <- basset(Y, t(mmX), Theta = Theta, Gamma = Gamma)

## ----fig.height=4, fig.width=6------------------------------------------------
## predictions
preds <- predict(mod)

## summarizing predictions over posterior samples
summary_preds <- apply(preds,MARGIN = c(1,2), FUN = mean)

## plotting
data.frame("points" = c(eta[1,], summary_preds[1,]), "group" = c(rep("Truth", length(eta[1,])), rep("Predicted", length(summary_preds[1,]))), "Index" = c(1:N, 1:N)) %>%
  ggplot(aes(x=Index, y= points, color = group)) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = c("red", "black")) +
  theme(legend.title=element_blank()) +
  ylab("Eta")

