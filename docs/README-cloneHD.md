# cloneHD command line arguments

## Typical usage options

Format of input files: the first two columns of all three input file
    types are always chromosome and coordinate of each observation. Both
    are expected to be integers (change X->23,Y->24). 

*   `--cna [file]` Read depth data file. 

    Format: For each sample, there are two additional columns with
    (i) the read depth and (ii) the number of independent observations
    this is the sum of. For human NGS data, use the mean read depth
    per 1 kb as the highest resolution. 

        1  1000  93   1  75   1  etc.
        1  2000  101  1  81   1
        1  3000  105  1  85   1
        1  5000  197  2  156  2
        etc.

*    `--baf [file]` B-allele read count data file. 

     Format: For each sample, there are two additional columns with
     (i) the number of reads of the minor allele and (ii) the total
     read depth at originally heterozygous loci.

        1  1036  43  90   28  72  etc.
        1  1287  47  99   32  80
        1  2877  30  100  36  82
        etc.

*    `--snv [file]` Somatic nucleotide variant read count data file.

     Format: For each sample, there are two additional columns with (i) the number 
     of reads of the somatic variant allele and (ii) the total read depth at that locus.

        1  1314  12  92   28  72  etc.
        1  1287  47  99   32  80
        1  2877  30  100  36  82
        etc.

*    `--pre [string:"./out"]`  Prefix for all output files.

*    `--bias [file]`  The bias field for the read depth data. 

     This must be a filterHD `pre.posterior-[int].txt` file, typically from a
     filterHD run on matched-normal read depth data. This is used to take into account
     the technical read depth modulation. Make sure that this is comparable between 
     matched normal and tumor (see below).

*    `--max-tcn [int]`  The maximum total copy number to be available.

     This is used as an upper limit for the total copy number genome wide (in all chr).
     If not specified, the normal copy number is used as limit for each chromosome.
     This number should be chosen conservatively, since it increases the
     HMM dimensionality and can open up the possibility for spurious solutions. 
     It is also possible to specify upper limits per chr and subclone (see below).

*    `--nmax [int:2]`  The maximum number of subclones to be tried.

     All subclone numbers `n` from `0..nmax` will be used and the one
     with maximum BIC chosen for output. BIC is affected by model complexity, 
     which is a combination of `n` and `max-tcn` (see below).

*    `--force [int]`  Fix the number `n` of subclones to be used.

*    `--trials [int:1]`  The number of independent optimizations.

     Global parameters, fractions `f` and mass `M`, are found numerically by local maximization of
     the total log-likelihood. The best result out of `trials` independent,
     newly seeded, runs will be used.

### snv-mode options

*    `--mean-tcn [file]`  Use a fixed mean total copy number for SNV data. 

     For a SNV data analysis, the cloneHD output file
     ending `pre.mean-tcn.txt` from a CNA(+BAF) run can be supplied here. Since
     the subclonal decomposition can be different for SNVs, this option
     ensures that a reasonable mean total copy number is used.

*    `--avail-cn [file]`  Use SNV copy number availablility constraint.

     For a SNV data analysis, the cloneHD output file
     ending `pre.avail-cn.txt` from a CNA(+BAF) analysis can be supplied here. Since
     the subclonal decomposition can be different for SNVs, this option
     ensures that the SNV genotype is consistent with the fraction of
     cells in which this number of copies is available, at all.
     This can only be used together with `--mean-tcn [file]` option. In
     combination, this is usually a much stronger constraint than using
     `mean-tcn` alone.

### Fuzzy segmentation options

For data with persistence along the genome, a fuzzy segmentation can
be used based on the filterHD posterior jump probability (must be a
`pre.jumps.txt` file). Data between potential jump sites, with a jump
probability of at least `min-jump`, is collapsed. The jump probability
is used in the HMM transition.

*    `--cna-jumps [file]`
*    `--baf-jumps [file]`
*    `--snv-jumps [file]`
*    `--min-jump [double:0.01]` 

## Parameter options

The shape parameter for the over-dispersed emission models
(Negative-Binomial or Beta-Binomial). If not specified, the normal
models are used (Poisson or Binomial).

*    `--cna-shape [double:inf]`
*    `--baf-shape [double:inf]`
*    `--snv-shape [double:inf]`

The rate for indiviual random emissions per data set. 

*    `--cna-rnd [double:1.0e-6]`
*    `--baf-rnd [double:1.0e-6]`
*    `--snv-rnd [double:1.0e-6]`

Both can be learned with filterHD for data with persistence.

## Advanced options

*    `--clones [file]` Use fixed mass(es) `M` and/or subclonal fractions `f`. 

     Either all mass parameters, or all subclonal fractions, or both 
     can be given (for each sample in the input data). The likelihoods
     and posteriors will be computed under these conditions. 
     Remaining parameters will be learned.
  
     Format: One line per sample. The first column, if greater than
     1.0, is interpreted as mass; the remaining as subclonal fractions.

        30.0 0.64 0.12
        28.0 0.31 0.23

    More than one parameter set can be given (as a continued list). Then,
    only the likelihoods are computed and printed to a file ending
    `pre.llh-values.txt`.  Useful for mapping the log-likelihood surface
    or comparing several given solutions.

*    `--purity [file]`  Use fixed purities, i.e. lower bounds for the sum of
     subclonal frequencies. One line per sample.

*    `--restarts [int:10]`  The number of perturbations in local random
     search mode.

     This simple random search routine is used: after finding a local
     maximum of LLH, the best solution is perturbed and a new optimum
     is sought. 

*    `--seed [int]`  A fixed seed to make inferences reproducible.

*    `--mass-gauging [0/1:1]`  Whether to use mass-gauging.

     The optimization in the space of masses (seq depths per haploid
     DNA) and subclonal frequencies can suffer from many local
     optima.  To fix the mass(es), one can, for a given solution,
     assume that an occupied state is actually all-normal. All
     occupied states will be proposed to fix the mass(es) 

*    `--min-occ [double:0.01]`  The minimum occupancy of levels to be
     used for the mass gauging.

*    `--print-all [0/1:0]`  If 1, the posterior for every observation
     is printed to files ending `*[cna/baf/snv].posterior.txt`. 
     If 0, only one line for each segment is printed.

*    `--max-tcn [file]`  The maximum total copy number per chr and subclone.

    This file should have the format: chr max1 max2 max3 etc., e.g.

        1  2
        2  2
        3  8  2
        4  2
        etc.

    The first column is the chromosome, the next columns are the limits to be used for subclone 1, 2 etc. For subclones not specified, the limit in the last column is used. In the example above, subclone 1 has an upper limit of 8 total copies in chr3, for all other subclones and in all other chromosomes, the upper limit is 2. If only SNV data is provided (and `--avail-cn [file]` is not given), this is used to fix the total number of copies. If `--max-tcn` is not given, cloneHD uses the normal copy number for each chr.

*    `--learn-priors [0/1:0]` For snv-mode only: if 1, then the parameters
     for the multiplicative SNV genotype priors are learned.

*    `--chr [file]`  Set normal copy numbers.

     The normal copy number for every single chromosome can be specified. This is needed only for non-human DNA. If not given, human DNA is assumed and the sex is inferred from the presence or absence of chr 24 (= chr Y) in the input data.

*    `--snv-fprate [double:1.0e-4]`  Set the false positive rate for SNVs.

*    `--snv-fpfreq [double:0.01]`  The typical frequency of false positive SNVs.

*    `--snv-pen-mult [double:0.01]`  Set the penalty against multiple SNVs.

*    `--snv-pen-high [doube:0.5]`  Set the penalty against higher genotypes.

*    `--cna-pen-zero [double:0.9]`  Set the penalty against zero total copies.

*    `--cna-pen-norm [double:1.0]`  Set the penalty against non-normal total c.n.

*    `--cna-pen-diff [double:1.0]`  Set the penalty against different total c.n.

*    `--baf-pen-comp [double:1.0]`  The penalty against complex minor allele status.

*    `--cna-jump [double:-1.0]`

*    `--baf-jump [double:-1.0]`

*    `--snv-jump [double:-1.0]`

    A constant jump probability per base pair can be set. If set to `-1.0`, then observations are uncorrellated along the genome (not the same as `1.0`). Can be learned with filterHD. No fuzzy data segmentation is performed.  Useful in combination with `--clones`, where very high-definition information available. Using this option will change the posterior output file format.

## Bulk options

These options are only needed if the sequenced cell population is a mixture of a
diverse bulk, with known allele frequency profile, and a number of
subclones with unknown genotypes and frequencies. Allele frequency
data is input with `--snv`. Data segmentation can be used with
`--snv-jumps`.  Read depth data can also be specified with `--cna`. 

*    `--bulk-mean [file]`  The bulk allele frequency profile. 

     Must be a filterHD `*posterior-[int].txt` file. Only the posterior mean is used.

*    `--bulk-prior [file]`  The bulk allele frequency profile. 

     Must be a filterHD `*posterior-[int].txt` file. The whole posterior
     distribution is used (run filterHD with `--dist 1` to obtain it).

*    `--bulk-updates [int:0]`  The number of (Bayesian) updates of the
     bulk allele frequency profile (if `--bulk-prior` was used).

*    `--bulk-fix [double:0.0]`  Use a flat and fixed bulk allele
     frequency profile.

## Technical options

The grid sizes for the pre-computed emission probabilities if fuzzy data segmentation is used.

*    `--cna-grid [int:300]`  
*    `--baf-grid [int:100]` 
*    `--snv-grid [int:100]`
