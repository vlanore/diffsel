library(rjags)

prefix <- "tmp_test"
n.chains <- 3

chains <- lapply(1:n.chains, function(i) {
  mcmc(read.table(paste0(prefix,i,".chain"), header=T, sep="\t")[,c(1:184)])
})

gelman.diag(chains)
