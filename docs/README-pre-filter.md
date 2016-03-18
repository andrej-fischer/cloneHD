# pre-filter command line arguments

The program `pre-filter` can be used to remove loci based on the observed read depth. It includes two heuristic filtering methods: loci are removed based on (i) their local variability and (ii) their being an outlier (see below). 

![pref](/images/prefilter-1.png "Pre-filtering of read depth via matched normal.")

The effect of `pre-filter` on read depth data: (A) Centromeric regions of real chromosomes often show huge large scale variability in their read depth. But there are also many small regions with very low read depth throughout. (B) After pre-filtering, the problematic regions are masked out. (C) Removing the same regions in the tumor data improves quality visibly while retaining biologically relevant features.


## Typical usage options

*  `--data [file]`  Input data to be pre-filtered. 
    
    The file format is the same as for cloneHD's `--cna` option. Only the first sample will be used for pre-filtering.

*  `--pre [string:"./out"]`  Set prefix for all output files. 

    The pre-filtered loci and data are print to a file named `pre.pref.txt`.

*  `--print-tracks [0/1:0]`  Print the window average and window variability. 
     
     The windowed tracks are used for pre-filtering. They are printed for all loci to a file named `pre.track.txt`. Use this to inspect and tune the pre-filter thresholds.

*  `--pick-from [file]`  Pre-filter data in this file by picking loci present in `match-to`. 

     Only loci are selected which fall into a bin also present in `match-to`. Bins in `match-to` are assumed to be of constant width with the given coordinate being the right bin end inclusive, e.g.

        1000  =  1-1000
        2000  =  1001-2000
        4000  =  3001-4000
        etc.

*  `--match-to [file]`  Use this file as reference to pick loci in `pick-from`. 

     Loci in this file are assumed to be equidistant (e.g. per 1 kb, not all bins need be present). The bin width is decided automatically by majority.

## Parameter options

*  `--window-size [int:100]` Set the window scale for smoothing (centered, +-size).

*  `--remove-outlier [double:3.0]`  Set the outlier threshold.

     All loci are removed, where the observed read depth is further than this value away from the local window-average (in units of sqrt(window-average), assuming Poisson distributed read depths). If set to `0.0`, filter is not applied. 

*  `--remove-variable [double:2.0]`  Set the variability threshold.

     All loci are removed, where the local window-variability exceeds this multiple of the global variability. Global (local) variability is defined as median (mean) of the absolute distance of observed read depths to the global median read depth. If set to `0.0`, filter is not applied. 
