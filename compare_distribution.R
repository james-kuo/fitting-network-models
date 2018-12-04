### Compare distributions
library(poweRlaw); library(GLDEX)
m <- displ$new(fun.zero.omit(degree(bit_graph, mode = "in"))) 
min <- estimate_xmin(m)
m$setXmin(min)
estimate_pars(m)
#bs_p = bootstrap_p(m) dislnorm disexp

m_ln <- dislnorm$new(fun.zero.omit(degree(bit_graph, mode = "in"))) 
m_ln$setXmin(m$getXmin())
est_ln = estimate_pars(m_ln)
m_ln$setPars(est_ln)
comp = compare_distributions(m, m_ln)

m_exp <- disexp$new(fun.zero.omit(degree(bit_graph, mode = "in"))) 
m_exp$setXmin(m$getXmin())
est_exp = estimate_pars(m_exp)
m_exp$setPars(est_exp)
comp2 <- compare_distributions(m, m_exp)
