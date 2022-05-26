### Create a version of 'sample' that avoids the problem when length(x)==1
#
#     output = sample.vec(x, ...)
#
#     x: vector of any length
#     ...: other arguments to 'sample'
#     -------------------------------------------------------
#     output: 
#     vector of length 'size' with elements drawn from x
#
#     notes:    forces 'sample' to always treat x as a vector
#               (see examples on the ?sample page and also here: http://stackoverflow.com/questions/13990125/sampling-in-r-from-vector-of-varying-length)
# 
#     Written by Ben Koch & Natasja van Gestel


sample.vec <- function(x, ...){
  x[sample(length(x), ...)]
}
