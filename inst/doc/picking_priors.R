## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # To put prior on ALR_j coordinates for some j in (1,...,D-1)
# Xi <- clrvar2alrvar(S, j)
# # To put prior in a particular ILR coordinate defined by contrast matrix V
# Xi <- clrvar2ilrvar(S, V)
# # To put prior in CLR coordinates (this one needs two transforms)
# foo <- clrvar2alrvar(S, D)
# Xi <- alrvar2clrvar(foo, D)

## ----eval=FALSE---------------------------------------------------------------
# # Transform from log-absolute-abundance effects to effects on absolute-abundances
# foo <- exp(A)
# # To put prior on ALR_j coordinates for some j in (1,...,D-1)
# Theta <- driver::alr_array(foo, j, parts=1)
# # To put prior in a particular ILR coordinate defined by contrast matrix V
# Theta <- driver::ilr_array(foo, V, parts=1)
# # To put prior in CLR coordinates
# Theta <- driver::clr_array(foo, parts=1)

