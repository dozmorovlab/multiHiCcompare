# # figure out progressive pooling by distance
# tab <- hicexp@hic_table
# 
# D <- tab$D %>% unique() %>% sort()
# maxD <- max(D)
# N <- length(D)
# 
# 
# 
# size <- 0
# for (i in 1:length(D)) {
#   chunks[i] <- 1 + size
#   size <- size + 1
# }
# 
# 
# chunks <- vector()
# for(i in 1:floor(maxD/2)) {
#   chunks <- c(chunks, rep(i, i))
# }
# 
# 
# pools <- rep(D, 1:length(D))[1:length(D)]
# object.size(pools, units = "auto", standard = "SI")
# 
# # use triangular number to solve for number of pools
# x <- length(D)
# n <- (sqrt((8 * x) + 1) - 1) / 2 
# n <- ceiling(n)
# pools <- rep(D[1:n], 1:n)[1:x]
