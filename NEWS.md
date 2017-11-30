# Version 1.0.1 [2017-11-28]

* Bug fixed to pass R CMD check under Solaris

# Version 1.0.0 [2017-11-27]

* CRAN submission.

# Version 0.8.2-9000 [2017-11-19]

* Minor updates to C code.

# Version 0.8.1-9000 [2017-07-25]

* Minor updates to C code to pass R CMD check cleanly under R 3.4.1
* 'acnr' moved from Depends to Imports

# Version 0.8.0-9000 [2017-04-11]

* Added more tests
* Fixed Issue #4

# Version 0.8.0 [2017-03-22]

* Fixed issue #5: the package does not rely on cghseg:::segmeanCO anymore!

# Version 0.7.4 [2017-03-22]

* Added pruned DPA algorithm [GR]
* Added a test that pDPA and classical DP yield identical results for (small) univariate signals

# Version 0.7.3 [2016-12-16]

* Data generation vignette is now .Rmd instead of .Rnw.

# Version 0.7.2 [2016-12-15]

* Moved documentation from inlinedocs to roxygen2 (as per Issue #2)
* Minor updates to vignette following changes to 'acnr'.
* PSSeg vignette is now .Rmd instead of .Rnw

# Version 0.7.1 [2016-11-24]

* Minor bug fix in 'plotSeg' in the default value of parameter binExclPattern.

# Version 0.7.0 [2015-07-08]

* In 'PSSeg': added argument 'rankTransform' to allow rank-based
segmentation to be performed.

# Version 0.6.4 [2014-12-19]

*  Bug fixed in vignette (remove GFLars segmentation)

# Version 0.6.3 [2014-10-27]

* Bug fixes for submission to CRAN

# Version 0.6.2 [2014-09-09]

* Updates to vignette.
* Updates to CITATION.

# Version 0.6.1 [2014-07-16]

* In 'estimateSd': SPEED UP for method "Hall": vectorization.
* In 'jointSeg', 'PSseg': it is now possible to run 'modelSelection'
 on initial segmentation.  Default behavior of the functions is
 unchanged.
* In 'getCopyNumberDataByResampling', if 'regAnnot' is NULL (the
default), frequencies of regions (0,1), (0,2), (1,1) and (1,2) (the
most common alterations) are set to represent 90% of the regions.

# Version 0.6.0 [2014-06-18]

* Cleanups in DESCRIPTION/NAMESPACE.
* Moved 'doCnaStruct' to 'otherMethods/'.
* Packages passes R CMD check --as-cran.

# Version 0.5.4 [2014-06-06]

* Replaced 'jointSeg' by 'jointseg' in 'inst' scripts.

# Version 0.5.3 [2014-06-04]

* Some changes to the main vignette (now using 'knitr').
* Updates in scripts to perform ASCAT segmentation.

# Version 0.5.2 [2014-05-30]

* Added a vignette to illustrate how to add a new method.

# Version 0.5.1 [2014-05-29]

* Changed Illumina to GSE11976 and Affymetrix to GSE29172
* Bug fixed for PSCBS, DynamicProgramming and CBS methods.
* Updated scripts

# Version 0.5.0 [2014-05-26]

* Major changes in the interface of the main functions 'PSSeg' and
'jointseg'.
* Now 'doGFLars' is a wrapper calling the lower-level 'segmentByGFLars'.
* Now 'doRBS' is a wrapper calling the lower-level 'segmentByRBS'.
* Added 'binMissingValues' so that it may be used not only for GFLars.
* Cleaned up documentation and examples.


# Version 0.4.7 [2014-05-14]

* in 'jointSeg': added argument 'segFUN' and flavor "other" to enable
a user-defined segmentation method to be used.


# Version 0.4.6 [2014-05-09]

* BUG FIX in 'doCBS': replaced start bkp given by CBS by end bkp as we
consider ends of segments for all methods.  Only makes a difference
when there are missing values in the signal.

# Version 0.4.5 [2014-05-06]

* Bug fixed in scripts to build table of the paper.
* Updated mapping for d|het segmentations (median).

# Version 0.4.4 [2014-05-05]

* Added scripts to 'inst/otherMethods' and 'inst/syntheticProfiles'
  to address referees' comments.

# Version 0.4.3 [2014-02-20]

* Added 'inst/figures' to reproduce the figures of the paper.

# Version 0.4.2 [2014-02-14]

* 'plotSeg': Calculation of segment means would fail when a breakpoint
  was detected between the first and second position.

# Version 0.4.1 [2014-02-13]

* Renamed function 'jointSeg' to 'jointseg' to avoid name clash with
  package name.
* Now all documented functions are exported.
* Updates in vignette.

# Version 0.4.0 [2014-02-13]

* Added scripts to reproduce methods comparison of the jointSeg paper
  to dedicated subdirectories.
* Added scripts to reproduce figures of the jointSeg paper to 'figures'.
* Removed defunct package PSCN from 'Suggests' field, and moved function
  'doPSCN' to (newly created) directory 'zzz.defunct'.
* Removed 'R.utils' from 'Depends' field.

# Version 0.3.17 [2014-02-12]

* Updates in 'plotSeg'.

# Version 0.3.16 [2013-12-09]

* 'segmentByNnn' replaced by 'doNnn'.
* Replaced 'segmentByCghseg' by 'doDynamicProgramming'.

# Version 0.3.15 [2013-12-06]

* Wrappers for PSCN and CnaStruct now operate systematically on
  log-scaled copy numbers (ie LRRs).  Thanks to report by Hugues
  Sicotte.
* Segmentation methods can now be run on numeric vectors or
 data.frames (when appropriate).

# Version 0.3.14 [2013-11-29]

* Added flavor 'DP'.

# Version 0.3.13 [2013-11-08]

* Bumped version number to trigger build by R-forge.

# Version 0.3.12 [2013-08-04]

* Added function 'getTpFp' so that performance can be evaluated in a
  scale which does not depend on the signal length.  Updated examples
  and vignette accordingly.

# Version 0.3.11 [2013-05-30]

* Added argument 'regNames' to 'plotSeg' so that region labels are
  plotted if available.

# Version 0.3.10 [2013-05-16]

* Some example code now embedded in 'require()' statements to avoid
  problems in the R CMD check mechanism of R-forge.

# Version 0.3.9 [2013-04-09]

* Added argument 'K' to 'pruneByDP'.

# Version 0.3.8 [2013-03-27]

* Modification in 'segmentByCnaStruct'

# Version 0.3.7 [2013-03-27]

* Added segmentation methods PELT and CnaStruct

# Version 0.3.6 [2013-03-19]

* Added flavor "Hall" to estimate standard deviation in 'estimateSd'.
* Added argument "DP" in 'PSSeg' and 'jointSeg'.
* Added argument "relax" in 'getTprTnr' and definition of true negatives.
* Added return value "rse" in 'segmentByRBS'.
* Added argument "connex" in 'getCopyNumberByResampling'.

# Version 0.3.5 [2013-02-26]

* In 'PSSeg' and 'jointSeg': added argument 'jitter'.
* In getCopyNumberDataByResampling': added argument 'connex' to forces
adjacent regions to be connex.
* In 'getTprTnr': added argument 'relax'.

# Version 0.3.4 [2013-02-18]

* 'segmentByGFLars' now handles missing values.
* In 'PSSeg': Flavor "GFLars" can now be run at full resolution

# Version 0.3.3 [2013-02-18]

* BUG FIX in 'segmentByPSCBS'.
* flavor = 'PSBCS' added.
* Removed segmentByPairedPSCBS.

# Version 0.3.2 [2013-02-14]

* BUG FIX in 'modelSelection': if the dimension is too small, in 'jointSeg' if flavor = 'CBS' or 'PSCN' no using of 'modelSelection'.

# Version 0.3.1 [2013-02-07]

* Bumped acnr dependency (>= 0.1.1) to fix minor case sensitivity issue.
* Reduced graphics size in vignette.

# Version 0.3.0 [2013-01-31]

* Updated vignette.
* Passes R CMD check with no errors.

# Version 0.2.5 [2013-01-31]

* Updated NAMESPACE.
* Moved dilution data sets to new package 'acnr'.

# Version 0.2.4 [2013-01-25]

* Cleanups in doc and return values.
* 'jointSeg' now returns 'bestBkp'
* Removed 'position' from I/O in 'getCopyNumberDataByResampling' and 'PSSeg'.
* In 'plotSeg': added arguments 'exclNames', 'ylabs', and
  'binExclPattern' so that 'plotSeg' can handle not only copy number
  signals.
* Improved doc and examples.
* Now using 'ERMajustment' for model selection.
* BUG FIX in 'segmentByRBS': Empty candidate list would give an
  error. Now returning early when 'minRegionSize' is too large for 'K'.

# Version 0.2.3 [2013-01-16]

* Added 'loadCnRegionData'.
* Added Affymetrix dilution data set (from GEO:GSE29172 and GSE26302).
* Removed 'setNormalContamination' as it is not needed anymore.
* In 'getCopyNumberDataByResampling':
  * Added argument 'regAnnot', through which theoretical frequencies
  for each CN regions can be specified.
  * Added arguments 'bkp' and 'regions' to allow for bypassing the
  breakpoint generation step.
  * Made the constraints on CN transitions more generic.
* Added 'loadCnRegionData'.

# Version 0.2.2 [2013-01-15]

* Replaced all 'jumps' by 'bkp'.
* Added argument 'platform' to 'segmentByPSCN'.
* Added Illumina dilution data set (from GEO:GSE11976).

# Version 0.2.1 [2013-01-06]

* BUG FIX in 'segmentByCghseg': index shift when reshaping
  results.
* Now 'PSSeg' can also run 'CBS' segmentation.

# Version 0.2.0 [2013-01-03]

* Now 'PSSeg' can also run 'PSCN' and 'cghseg' segmentations.
* Updated doc and vignette.

# Version 0.1.21 [2012-12-31]

* Now 'segmentByRBS' handles missing values.  Therefore, 'PSSeg'
  and 'jointSeg' with the corresponding flavor also do.
* Added Argument 'statistic' to 'PSSeg'.
* Now 'segmentByRBS' scales each dimension to unit variance using
  'estimateSd'.

# Version 0.1.20 [2012-12-30]

* Some code and doc cleanups.
* Renamed 'binSeg' to 'segmentByRBS'.
* Renamed 'dpSeg' to 'pruneByDP'.
* Renamed 'segmentByGflars' to 'segmentByGFLars'.

# Version 0.1.19 [2012-12-23]

* SPEEDUP: removed redundant calls to 'getRSE' in 'binSeg'.

# Version 0.1.18 [2012-12-17]

* Added function 'runPSCN' as a convenient wrapper to PSCN + DP.

# Version 0.1.17 [2012-12-15]

* Added function 'getTprTnr' to evaluate the performance of
segmentation methods.
* Added function 'prof' for optional reporting of CPU and memory usage.

# Version 0.1.16 [2012-12-12]

* Added 'setNormalContamination'.

# Version 0.1.15 [2012-12-07]

* BUG FIX in 'binSeg': index shift in correspondence b/w breakpoint position and interval.

# Version 0.1.14 [2012-12-06]

* Added 'binSeg' for binary segmentation.
* Added wrappers 'jointSeg', 'plotSeg', 'PSGFL'.
* Added argument 'flavor' to 'jointSeg' so that one can run GFLARS or
binary segmentation before dynamic programming.

# Version 0.1.13 [2012-12-01]

* Added example data files and scripts based on public data set GSE19539.
* Added scripts to generate these data files from raw CEL files in GSE19539.

# Version 0.1.12 [2012-11-27]

* Added model selection part to 'modelSelection'.
* Removed model selection part in 'dpseg'.

# Version 0.1.11 [2012-11-25]

* Renamed 'randomProfileTCGA' to 'getCopyNumberDataByResampling'.

# Version 0.1.10 [2012-11-16]

* Added 'randomProfileTCGA'.

# Version 0.1.9 [2012-10-19]

* BUG FIX in randomProfile: 'minLength' would not work as expected.

# Version 0.1.8 [2012-09-13]

* Tentative bug fix in gflars: indices in gammaTemp.
* SPEEDUP: removed unnecessary calls to 'complex'.
* Updated documentation

# Version 0.1.7 [2012-09-13]

* Some code cleanups.

# Version 0.1.6 [2012-08-23]

* SPEEDUP: removed explicit calls to 'matrix()' and 'apply(..., rev)'.
* Package documentation now done using 'inlinedocs'.

# Version 0.1.5 [2012-08-22]

* Improved scripts for assessing replication of the MATLAB version.
* Improved scripts for testing speed.

# Version 0.1.4 [2012-08-21]

* BUG FIX: in optimizeLARS.R, the index in 'lambda' was off by 1.
* Added 'value' to the return values of 'optimizeLARS'.
* Renamed 'optimizeLARS' to 'gflars'.


# Version 0.1.3 [2012-08-20]

* BUG FIX: in optimizeLARS.R, 'Y' should be centered but not scaled.
* Added dpseg.R: finds the best set of change points from a(n over-)
  segmentation using joint dynamic programming.

# Version 0.1.2 [2012-08-17]

* SPEEDUP: replaced 'apply(*, 1, sum)' by 'rowSums(*)'.

# Version 0.1.1 [2012-08-15]

* Added test scripts for replicating MATLAB version.
* Added test scripts for speed trials.


# Version 0.1.0 [2012-08-01]

* Created.
