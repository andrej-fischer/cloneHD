# cloneHD and filterHD software

The current stable release, including pre-compiled executable binaries of filterHD and cloneHD for Mac OS X (64bit), can be found at <ftp://ftp.sanger.ac.uk/pub/teams/153/cloneHD/>. 

# Run test with simulated data

If you download cloneHD from the ftp site above, you can test both filterHD and cloneHD by running

`$ sh run-example.sh`

where you can see a typical workflow of analysing read depth and BAF
data with a matched normal. All command line arguments are explained below.

# Compilation  

To compile cloneHD yourself, you need the GNU scientific library (GSL) v1.15 or later. Change the paths in the Makefile to your GSL installation location (if non-standard). Then type 

`$ make -f Makefile_cloneHD`

in the source directory. The two executables, `filterHD` and `cloneHD`, will be in `./build`.

# filterHD command line arguments

## Typical usage options

*    `--data [file]`  Input data. File format: same as below for
     `--cna`, `--baf` or `--snv`. Samples are processed one by one.

*    `--mode [1/2/3/4]`  Emission model modes.

        1. Binomial (for SNV data and BAF data (with `--reflect 1`))
        2. Beta-Binomial (over-dispersed Binomial)
        3: Poisson (for read depth data) 
        4: Negative-Binomial (over-dispersed Poisson)

    In modes 3/4, the range of the hidden track is learned
    automatically. For modes 1/2, its always [0,1]. All with reflective
    boundary conditions.

*    `--pre [string:"./out"]`  Prefix for all output files.

*    `--dist [0/1:0]`  Whether to print the whole posterior
     distribution.

*    `--jumps [0/1:0]`  Whether to print a separate posterior jump
     probability track, compounded over all samples.

*    `--reflect [0/1:0]`  If 1, binomial observations `n in N` and
     `N-n in N` are assumed identical. Use for BAF data.

## Parameter options

*    `--jump [double]`
*    `--sigma [double]`
*    `--shape [double]`
*    `--rnd [double]`

*    `--jumpi [double]`
*    `--sigmai [double]`
*    `--shapei [double]`
*    `--rndi [double]`

## Advanced options

*    `--min-jump [double:0.0]`

*    `--filter-pVal [0/1:0]`

*    `--filter-shortSeg [int:0]`



# cloneHD command line arguments

## Typical usage options

Format of input files: the first two columns of all three input file
    types are always chromosome and coordinate of each observation. Chromosomes
    are expected to be integers (X->23,Y->24). 

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

     Format: For each sample, there are two additional columns with (i) the number of reads of the somatic variant allele and (ii) the total read depth at that locus.

        1  1314  12  92   28  72  etc.
        1  1287  47  99   32  80
        1  2877  30  100  36  82
        etc.

*    `--pre [string:"./out"]`  Prefix for all output files.

*    `--bias [file]`  The bias field for the read depth data. 

     This must be a filterHD `*posterior.*.txt` file, typically from a
     filterHD run on matched-normal read depth data, to estimate the
     technical read depth modulation.

*    `--maxcn [int:4]`  The maximum possible total copy number genome
     wide.

     This number should be chosen conservatively, since it increases the
     HMM dimensionality and can open the possibility for spurious solutions. 

*    `--nmax [int:2]`  The maximum number of subclones to be tried.

     All subclone numbers from 0 to `nmax` will be used and the one
     with maximum BIC chosen for output.

*    `--force [int]`  Force the number of subclones to be used.

*    `--trials [int:1]`  The number of independent optimization
     trials.

     Global parameters are found numerically by local maximization of
     the total log-likelihood. The best result out of `trials` independent,
     randomly seeded,  runs will be used.

*    `--copynumber [file]`

### Fuzzy segmentation options

For data with persistence along the genome, a fuzzy segmentation can
be used based on the filterHD posterior jump probability. Data between
potential jump sites, with a jump probability of at least `min-jump`, is
collapsed. The jump probability is used in the HMM transition.

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

The rate for indiviual random emissions per data set. Can be learned
with filterHD for data with persistence.

*    `--cna-rnd [double:0.0]`
*    `--baf-rnd [double:0.0]`
*    `--snv-rnd [double:0.0]`

A constant jump probability per base pair. If -1, then observations are
uncorrellated along the genome. Can be learned with filterHD. No fuzzy
data segmentation is performed. Very high-definition information
available. Useful in combination with `--clones`.

*    `--cna-jump [double:-1.0]`
*    `--baf-jump [double:-1.0]`
*    `--snv-jump [double:-1.0]`

## Advanced options

*    `--clones [file]`  Use fixed mass(es) and/or subclonal
     frequencies. Either all mass parameters, or all
     subclonal frequencies, or both can be given (for each sample in
     the data input). The likelihoods and posteriors will be computed
     under these conditions. Remaining parameters will be learned.
  
     Format: One line per sample. The first column, if greater than
     1.0, is interpreted as mass; the remaining as subclonal frequencies.

        30.0 0.64 0.12 0.03
        28.0 0.31 0.23 0.11 

    More than one parameter set can be given (appended below). Then,
    only likelihoods are computed and output to a file ending
    `*llh-values.txt`. 

*    `--purity [file]`  Fixed purities, lower bounds for the sum of
     subclonal frequencies. One line per sample.

*    `--restarts [int:10]`  The number of perturbations in local random
     search mode.

     This simple random search routine is used: after finding a local
     maximum of LLH, the best solution is perturbed and a new optimum
     is sought. 

*    `--seed [int]`  A fixed seed to make inferences reproducible.

*    `--mass-gauging [0/1:1]`  Whether to fix the mass(es) by assuming
     occupied states to be actually all-normal. 

*    `--min-occ [double:0.01]`  The minimum occupancy of levels to be
     used for the mass gauging.

*    `--print-all [0/1:0]`  If 1, the posterior for every observation
     is printed. If 0, only one line for each segment is printed.

*    `--learn-priors [0/1:0]` snv-mode only: if 1, then the parameters
     for the multiplicative genotype priors are learned.

*    `--maxcn-mask [file]`  A maximum total copy number for individual
     chromosomes can be specified. This is useful, if large copy
     numbers are expected in only in part of the genome.

*    `--chr [file]`  The normal copy number for every single
     chromosome can be specified. Needed for non-human DNA. If not
     given, human DNA is assumed and the sex inferred from the
     presence of chr 24 (= chr Y) in the input data.

*    `--snv-fpr [double:1.0e-4]`  The false positive rate for SNVs,
     i.e. rate of SNV data points of genotype all-0.

*    `--snv-err [double:0.01]`  The typical frequency of mis-called variant reads.

*    `--snv-pen [double:0.01]`  The penalty for higher than expected
     genotypes.

*    `--baf-pen [double:1.0]`  The penalty for complex minor allele status.

## Bulk options

These options are only needed if the population is a mixture of a
diverse bulk with known allele frequency profile and a number of
subclones. Allele frequency data is input with `--snv`.

*    `--bulk-mean [double]`  The bulk allele frequency profile. Must be a filterHD `*posterior.*.txt` file. Only the posterior mean is used.

*    `--bulk-prior [file]`  The bulk allele frequency profile. Must be
     a filterHD `*posterior.*.txt` file. The whole posterior
     distribution is used (run filterHD with `--dist 1`).

*    `--bulk-updates [int:0]`  The number of Bayesian updates of the
     bulk allele frequency profile.

*    `--bulk-fix [double:0.0]`  Using a flat fixed bulk allele
     frequency profile.

## Technical options

*    `--grid [int:300]`  The grid size for the pre-computed emission
     probabilities if fuzzy segmentation is used.
