## Test environments

* local OS X install, R 3.4.1
* ubuntu 12.04 (on travis-ci): oldrel/release/devel
* OSX (on travis-ci): oldrel/release
* win (on appveyor): devel/release
* win (on r-hub): oldrel/release/devel/patched

## R CMD check --as-cran results

Status: OK

There were no ERRORs or WARNINGs or NOTEs, except a WARNING on win (r-hub) due to the fact that pandoc was not installed on that machine.

## Downstream dependencies

There are currently no downstream dependencies for this package.
