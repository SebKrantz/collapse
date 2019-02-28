# Helper functions

# convert factor to numeric values

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

# Mode function

Mode <- function(x, na.rm = FALSE) {
  ax <- attributes(x)

  if(na.rm) x <- x[!is.na(x)]
  # Calculate the mode
  ux <- unique(x)
  y <- ux[which.max(tabulate(match(x, ux)))]

  attributes(y) <- ax
  return(y)
}

.datatable.aware = TRUE
