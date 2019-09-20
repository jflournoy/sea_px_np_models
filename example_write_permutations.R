library(permute)

#This presumes you've arleady set up the target model
load('../chavg.fear/debug.Rdata')
attach(designmat)
nperm <- 2000 #Make sure this number matches input to NeuroPointillist or is bigger

set.seed(622019*37)
ctrl.series <- how(within = Within(type = 'series'), nperm = nperm, blocks = ID)
perm_set.series <- shuffleSet(n = ID, control = ctrl.series)
saveRDS(perm_set.series, 'permutation_set-series.RDS')


ctrl.free <- how(within = Within(type = 'free'), nperm = nperm, blocks = ID)
perm_set.free <- shuffleSet(n = ID, control = ctrl.free)
saveRDS(perm_set.free, 'permutation_set-free.RDS')
