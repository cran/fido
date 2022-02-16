## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
library(fido)
library(phyloseq)
library(dplyr)
library(ggplot2)
data(mallard_family)

# Just take vessel 1
mallard_family <- prune_samples(sample_data(mallard_family)$Vessel==1, mallard_family)

# Just take hourly samples
mallard_family <- prune_samples((sample_data(mallard_family)$time > "2015-11-20 15:00:00 UTC") &
                      (sample_data(mallard_family)$time < "2015-11-25 16:00:00 UTC"), mallard_family)

# Order samples - to make plotting easy later
o <- order(sample_data(mallard_family)$time)
otu_table(mallard_family) <- otu_table(mallard_family)[o,]
sample_data(mallard_family) <- sample_data(mallard_family)[o,]

# Extract Data / dimensions from Phyloseq object
Y <- t(as(otu_table(mallard_family), "matrix"))
rownames(Y) <- taxa_names(mallard_family)
D <- ntaxa(mallard_family)
N <- nrow(sample_data(mallard_family))

# X in hours
X <- as.numeric(sample_data(mallard_family)$time)
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
family_names <- as(tax_table(mallard_family)[,"Family"], "vector")
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


