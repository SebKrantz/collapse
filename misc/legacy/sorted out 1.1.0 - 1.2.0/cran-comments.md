## Resubmission collapse 1.1.0
In this version I have:

* Fixed valgrind issues.

## Test environments
* local Windows 8.1 install, R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

  * checking installed package size ... NOTE
    installed size is  7.9Mb
    sub-directories of 1Mb or more:
      doc    1.6Mb
      libs   4.7Mb

This has to do with compiled files. Data is 0.5 Mb. 
