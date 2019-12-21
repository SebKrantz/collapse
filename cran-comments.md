# This is the first CRAN submission of collapse

## License and Authorship
Since C-code from data.table and stats was copied/modified into src, data.table's
Mozilla Public License 2.0 was adopted and the authors listed as 'ctb'. 
Additionally, Simon Gaure, the author of lfe from which I import 
the demeanlist function (written in C), and Dirk Eddelbuettel, whose package
Rcpp I have thoroughly exploited in the creation of this one, were added as 'ctb'.
In addition, all authors have been mentioned on the collapse-package.Rd page.

(An alternative would be listing us all as 'aut'. My current understanding though 
 is that an 'aut' is someone with substantial stakes in the design, purpose or functionality 
 of a package, which is not the case here (none of these people know about or have actively 
 contributed to this package and I have not simply duplicated functionality of packages from which
 C-code was taken without substantial modifications).)

(The C-functions imported can be seen in src/ExportSymbols.cpp lines 5-15, 
 the rest is C++ of which I am the sole author).

## Test environments
* local Windows 8.1 install, R 3.6.1
* ubuntu 12.04 (on travis-ci), R 3.1.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 4 NOTEs:

  * checking installed package size ... NOTE
    installed size is  6.7Mb
    sub-directories of 1Mb or more:
      libs   4.6Mb

I guess this has to do with compiled files. Pre-compilation, 
the size of all .R, .c, .cpp, .h, .man, .rda and .Rmd files 
together is about 2.4 Mb, of which 0.5 Mb is data (.rda).

* checking S3 generic/method consistency ... NOTE
  Found the following apparent S3 methods exported but not registered:
    all.identical

* checking Rd \usage sections ... NOTE
  S3 methods shown with full name in documentation object 'small-helpers':
    'all.identical'
  
all.identical is not an S3 method. I used '.' to align with base functions like
all.equal, all.vars, all.names etc., and all_identical already exists with different
functionalities in 2 other CRAN packages.

* checking pragmas in C/C++ headers and code ... NOTE
  Files which contain pragma(s) suppressing diagnostics:
    'src/data.table_forder.c' 'src/data.table_subset.c'
    
I am currently not using a Makevars to enable/compile openMP parallelism which could 
(but need not) be used by data.table's forder and subsetDT. The reason for this is that the
package is compiled by Rcpp and one inevitably runs into this problem: https://stackoverflow.com/questions/54056594/cran-acceptable-way-of-linking-to-openmp-some-c-code-called-from-rcpp. Therefore I have disabled -Wunknown-pragmas warnings in data.table_forder.c and
data.table_subset.c, which naturally occur when the package is compiled without openMP and with -Wall.
I note that these warnings do not show up in R CMD check (so turning them on again would remove this note), but they are annoying running across the screen of the user installing the package. 
