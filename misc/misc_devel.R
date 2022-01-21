# Uncommit:
git reset --soft HEAD^

list.files(R.home("include")) # https://www.r-bloggers.com/2012/11/using-r-callhello/
file.edit(file.path(R.home("include"), "Rinternals.h"))
file.edit(file.path(R.home("include"), "Rdefines.h"))



dyn.unload("C:/Users/Sebastian Krantz/Documents/R/C++/truelength.dll")
dyn.load("C:/Users/Sebastian Krantz/Documents/R/C++/truelength.dll")
tl <- function(x) .Call("TL", x)
ltl <- function(x) c(length(x), tl(x))


devtools::check(document = FALSE, args = c('--as-cran', '--use-gct'), build_args = c('--no-docs'))

vignettes/collapse_intro
vignettes/collapse_and_plm
vignettes/collapse_and_dplyr
vignettes/collapse_and_data.table
vignettes/collapse_and_sf
tests


# TODO: Eliminate this function, do it at C-level!! Make sure it works in all cases though... even if passing `+`, not "+"
TtI <- function(x)
  switch(if(is.character(x)) x else x + 1L,
         replace_NA = 0L, replace_fill = 1L, replace = 2L, `-` = 3L, `-+` = 4L, `/` = 5L, `%` = 6L, `+` = 7L, `*` = 8L, `%%` = 9L, `-%%` = 10L,
         stop("Unknown transformation!"))



