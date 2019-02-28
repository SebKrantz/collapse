# collapse R package

The creation of the `collapse` package was inspired by the like-named command in the STATA statistical software,
and its main function `collap` can aggregate multi-type data exactly the way `collapse` does in STATA. But `collap` it is more than simply a reproduction of STATA's `collapse` for R: The function provides additional flexibility to perform advanced aggregation tasks more user friendly and faster, particularly when working with multi-level (panel) data structures and multi-type data. Next to that,
the function `qsu` (shorthand for quick summary) provides a way of calculating arbitrary summary statistics for cross-sectional and multi-level (panel) data, and a quick method to obtain within-transformed data. More specifically `qsu` provides the functionality of `summary`, `summaryBy` and STATA's `xtsummarize` and `xtsummarize` by groups in one command, in addition to some further features.


## Installation

The package can be installed in R using the following code:

remotes::install_github("SebKrantz/collapse")

## Helping out

Development takes place in the `devel` branch. All help is welcome ofcourse. If you want to contribute, please fork and create a pull request for merging with the `devel` branch.

## Warning

The package is currently in development. No warranties whatsoever.
