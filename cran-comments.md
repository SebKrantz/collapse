## Resubmission
This is a resubmission. In this version I have:

* Added all data.table and Rcpp authors as 'ctb' since some data.table source code was used and a tiny fraction of Rcpp source code. Everybody who has potentially contributed even one line of code is now a 'ctb'.  

* Adopted Rcpp's GPL (>= 2) license which is the most restrictive and not in conflict with data.tables MPL 2.0. 

* Removed a options(warn = -1), and properly implemented par resetting in plot methods using on.exit.

* Removed unnecessary \dontrun{} calls in documentation.

## Test environments
* local Windows 8.1 install, R 3.6.1
* win-builder (devel and release)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on Rhub)
* macOS 10.11 El Capitan, R-release (on Rhub)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

  * checking installed package size ... NOTE
    installed size is  6.7Mb
    sub-directories of 1Mb or more:
      libs   4.6Mb

This has to do with compiled files (.dll's). Pre-compilation, 
the size of all .R, .c, .cpp, .h, .man, .rda and .Rmd files 
together is about 2.4 Mb, of which 0.5 Mb is data (.rda).
