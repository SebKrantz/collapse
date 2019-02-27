# collapse R package

The creation of the `collapse` package was inspired by the like-named command in the STATA statistical software,
but its main function `collap` is more than simply a reproduction of STATA's `collapse` for R.
It can be used to perform more advanced aggregation tasks more flexibly and faster, particularly when working with multi-level (panel) data structures and multi-type data. Next to that,
the function `qsu` (shorthand for quick summary) provides a 
way of calculating arbitrary summary statistics for cross-sectional and multi-level (panel) data, and provides a quick way to obtain within-transformed data.


## Installation

The package can be installed in R using the following code:

remotes::install_github("SebKrantz/collapse")

## Warning

The package is currently in development. No warranties whatsoever.
