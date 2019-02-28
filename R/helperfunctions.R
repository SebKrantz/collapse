# Helper functions
# ----------------
# None of these functions do the necessary checks. Make sure you
# do them yourself if necessary before using these functions!

# convert factor to numeric values
# --------------------------------
as.numeric.factor <- function(x) {
  as.numeric(levels(x))[x]
}

# Mode function
#--------------
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

# Function NaN to NA
#-------------------
nantona <- function(x) {
  x[is.nan(x)] <- NA
  return(x)
  }

# Quick construction of a data.frame
# ----------------------------------
quickdf <- function(l) {
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  return(l)
}

# Create unique identifiers
# -------------------------
# Reworked to avoid the conversion to and from factor.
ident <- function(x) {

  if(is.factor(x)){
    # Work is done already
    l <- length(attr(x,"levels"))
    x <- as.numeric(x)
  } else{
    # Construct l, levs, x
    levs <- sort.default(unique.default(x))
    l <- length(levs)
    x <- match(x, levs)
  }

  s <- as.character(seq_len(l))
  n <- nchar(s)
  ylev <- paste0(strrep("0", n[l] - n), s)
  return(ylev[x])
}
