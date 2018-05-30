# # read test data
# r1 <- read.table("data/HIC001.NONE.chr22.100000.txt", header = FALSE)
# r2 <- read.table("data/HIC002.NONE.chr22.100000.txt", header = FALSE)
# r3 <- read.table("data/HIC003.NONE.chr22.100000.txt", header = FALSE)
# r4 <- read.table("data/HIC004.NONE.chr22.100000.txt", header = FALSE)
# r5 <- read.table("data/HIC005.NONE.chr22.100000.txt", header = FALSE)
# r6 <- read.table("data/HIC006.NONE.chr22.100000.txt", header = FALSE)
# r7 <- read.table("data/HIC007.NONE.chr22.100000.txt", header = FALSE)
# #
# r1 <- cbind('22', r1)
# r2 <- cbind('22', r2)
# r3 <- cbind('22', r3)
# r4 <- cbind('22', r4)
# r5 <- cbind('22', r5)
# r6 <- cbind('22', r6)
# r7 <- cbind('22', r7)

# devtools::use_data(r1, r2, r3, r4, r5, r6, r7, compress = 'xz', overwrite = TRUE)
# devtools::use_data(r1, r2, r3, r4, compress = 'xz', overwrite = TRUE)
#
# # test out class
# experiment <- new("hicexp", hic_matrices = list(r1, r2, r3, r4, r5, r6, r7), groups = c(1, 1, 1, 2, 2, 2, 2))
# hicexp <- make_hicexp(r1, r2, r3, r4, r5, r6, r7, groups = c(1, 1, 1, 2, 2, 2, 2))

