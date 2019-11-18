# .Call(setSizes)
# .Call(initDTthreads)
initDT <- function() cat(.Call(collapse_init, "Welcome to collapse! See ?collapse"))
# initDT()
