# # set up
# library(HiCcompare2)
# library(BiocParallel)
# register(MulticoreParam(workers = 3), default = TRUE)
# 
# # no parallel
# start.time = proc.time()
# q = cyclic_loess(hicexp, iterations = 3, span = NA, verbose = F, Plot = F, parallel = F)
# end.time = proc.time()[3]-start.time[3]
# end.time
# 
# # parallel
# start.time = proc.time()
# q = cyclic_loess(hicexp, iterations = 3, span = NA, verbose = F, Plot = F, parallel = T)
# end.time = proc.time()[3]-start.time[3]
# end.time