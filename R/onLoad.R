# .Call(setSizes)
# .Call(initDTthreads)
init_collapse <- function() cat(.Call(collapse_init, "Welcome to collapse! See ?collapse"))
#
