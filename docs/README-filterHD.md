# filterHD command line arguments

## Typical usage options

*    `--data [file]`  Input data. 

     The file format is the same as below for `--cna`, `--baf` or
     `--snv`. Multiple amples are processed independently, one by one.

*    `--mode [1/2/3/4]`  Emission modes.

        1. Binomial (for SNV data and BAF data (use with `--reflect 1`))
        2. Beta-Binomial (over-dispersed Binomial)
        3: Poisson (for read depth data) 
        4: Negative-Binomial (over-dispersed Poisson)

    In modes 3/4, the range of the hidden emission rate is learned
    automatically. For modes 1/2, it is always in [0,1]. Reflective
    boundary conditions are used.

*    `--pre [string:"./out"]`  Prefix for all output files.

*    `--dist [0/1:0]`  Whether to print also the  posterior distribution. 
     
     The posterior mean, std-dev and jump probability are always printed  to files
     `pre.posterior-[int].txt`, one for each sample in the input. With 1, the
     whole posterior distribution is also printed, so files can be big. 

*    `--jumps [0/1:0]`  Whether to print posterior jump probability. 

     The posterior jump probability is compounded over all samples. It
     can be used with `--min-jump [double]` below, to consolidate jumps.

*    `--reflect [0/1:0]`  If 1, binomial observations `n in N` and
     `(N-n) in N` are assumed to be identical. Use this option for BAF data.

## Parameter options

The HMM underlying filterHD is determined by these four global
parameters. They can all be fixed, otherwise they are learned from the data.

*    `--jump [double]`   Fix the jump probability per length unit (bp).
*    `--sigma [double]`  Fix the diffusion constant. 
*    `--shape [double]`  Fix the shape parameter for modes 2/4. If >1000, use modes 1/3.
*    `--rnd [double]`    Fix the rate of random emissions.

For all of the above parameters, initial values for the numerical
optimization can be given. This might be useful if you suspect several
local optima and want to start in the neighbourhood of a particular one.

*    `--jumpi [double]`
*    `--sigmai [double]`
*    `--shapei [double]`
*    `--rndi [double]`

## Further advanced options

*    `--min-jump [double:0.0]`  Consolidate jumps down to `--min-jump`.

     The posterior jump probability track will be consolidated by merging neighboring jump events into
     unique jumps, down to the minimum value given here. Can only be used together with
     `--jumps 1`. 

*    `--filter-pVal [0/1:0]`  Use p-Value filter.

     Filter sites where the p-Value of the
     observation is below `10/nSites`, where `nSites` is the total number
     of sites in a sample.

*    `--filter-shortSeg [int:0]` Use short-segment filter.

     Filter sites within short segments between jumps. All filtered data will be in the file ending `pre.filtered.txt`, which will be in the same format as the input file.

*    `--grid [int:100]`  Set the grid size.

     The grid size for the internal representation of continuous distributions. For large ranges in
     mode 3/4, it can make sense to increase this resolution.

## filterHD output files

TBD
