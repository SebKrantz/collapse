## Resubmission
This is a resubmission. In this version I have:

* Only listed the main authors ('aut'/citation) of data.table and Rcpp as 'ctb' to my package as suggested by Uwe.  

* Properly implemented an options(warn = -1) using on.exit.

* Removed a commented-out \dontrun{} call in documentation (GGDS10S.Rd). (The first example in GGDS10S.Rd is still wrapped in \dontrun{} because it uses non-suggested packages and generates a
runtime note on R CMD check)

* Placed some examples in collapse-package.Rd which were in \dontrun{} before in \donttest{}, to avoid a example runtime note. 

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
