# How to get cloneHD and filterHD?

The current stable release, including pre-compiled executable binaries
of filterHD and cloneHD for Mac OS X (64bit), can be found at
<ftp://ftp.sanger.ac.uk/pub/teams/153/cloneHD/>.  The source code can
also be downloaded here.

# Run a test with simulated data

If you download cloneHD from the ftp site above, you can test both filterHD and cloneHD by running

`$ sh run-example.sh`

where you can see a typical workflow of analysing read depth and BAF
data with a matched normal. All command line arguments are explained below.

# Compilation  

To compile cloneHD yourself, you need the GNU scientific library ([GSL](http://www.gnu.org/software/gsl/)) v1.15 or later. Change the paths in the Makefile to your GSL installation location (if non-standard). Then type 

`$ make`

in the source directory. The two executables, `filterHD` and
`cloneHD`, will be in `./build`.

# What are cloneHD and filterHD for?

cloneHD is a software for reconstructing the subclonal structure of a
population from next-generation short-read sequencing data. Read depth
data, B-allele count data and somatic nucleotide variant (SNV) data can be
used for the inference. cloneHD can find the number of subclonal
populations, their copy number profiles, their B-allele status and all
SNV genotypes with high resolution.

filterHD is a general purpose probabilistic filtering algorithm for one-dimensional
discrete data, similar in spirit to a Kalman filter. It is a continuous state
space Hidden Markov model with Poisson or Binomial emissions and a
jump-diffusion propagator. It can be used for scale-free smoothing, 
fuzzy data segmentation and data filtering. 

![cna gof](/images/cna.gof.png "CNA goodness of fit")
![cna gof](/images/cna.post.png "CNA posterior")
![baf gof](/images/baf.gof.png "BAF goodness of fit")
![cna gof](/images/baf.post.png "BAF posterior")

Visualization of the cloneHD output for the simulated data set. From
top to bottom: (i) the bias corrected read depth data and the cloneHD
posterior mean emission rate (ii) the total copy number posterior
distribution for subclone 1 with f1=0.52 and subclone 2 with f2=0.07
(iii) the BAF and (iv) the minor allele posterior. (Plots created with Wolfram
[Mathematica](http://www.wolfram.com/mathematica/).)

# filterHD command line arguments

## Typical usage options

*    `--data [file]`  Input data. 

     The file format is the same as below for `--cna`, `--baf` or
     `--snv`. Samples are processed one by one.

*    `--mode [1/2/3/4]`  Emission modes.

        1. Binomial (for SNV data and BAF data (use with `--reflect 1`))
        2. Beta-Binomial (over-dispersed Binomial)
        3: Poisson (for read depth data) 
        4: Negative-Binomial (over-dispersed Poisson)

    In modes 3/4, the range of the hidden emission rate is learned
    automatically. For modes 1/2, it is always in [0,1]. Reflective
    boundary conditions are used.

*    `--pre [string:"./out"]`  Prefix for all output files.

*    `--dist [0/1:0]`  Whether to print the  posterior distribution. 

     Files can be big. The posterior mean, std-dev and
     jump probability are printed in all cases to files
     `*posterior.*.txt`, one for each sample in the input.

*    `--jumps [0/1:0]`  Whether to print posterior jump probability. 

     The posterior jump probability is compounded over all samples. It
     can be used with `--min-jump [double]` below, to consolidate jumps..

*    `--reflect [0/1:0]`  If 1, binomial observations `n in N` and
     `(N-n) in N` are assumed to be identical. Use this option for BAF data.

## Parameter options

The HMM underlying filterHD is determined by these four
parameters. They can all be fixed, otherwise they are learned from the data.

*    `--jump [double]`  The jump probability per length unit (bp).
*    `--sigma [double]`  The diffusion constant. 
*    `--shape [double]`  The shape parameter for modes 2/4. If >1000, use modes 1/3.
*    `--rnd [double]`  The rate of random emissions.

For all of the above parameters, initial values for the numerical
optimization can be given. This might be useful if you suspect several
local optima.

*    `--jumpi [double]`
*    `--sigmai [double]`
*    `--shapei [double]`
*    `--rndi [double]`

## Advanced options

*    `--min-jump [double:0.0]`  Consolidate jumps down to `--min-jump`.

     The posterior jump probability track will be consolidated by merging neighboring jump events into
     unique jumps, down to the minimum value given here. Can only be used with
     `--jumps 1`. 

*    `--filter-pVal [0/1:0]`  Use p-Val filter.

     Filter sites where the p-Value of the
     observation is below 10/nSites, where nSites is the total number
     of sites in a sample.

*    `--filter-shortSeg [int:0]` Use short segment filter.

     Filter sites within short segments between jumps. All filtered data will be in the file ending `*filtered.txt`.

*    `--grid [int:100]`  Set grid size.

     The grid size for the internal representation of continuous distributions. For large ranges in
     mode 3/4, it can make sense to increase resolution.



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

*    `--max-tcn [file/int]`  The maximum total copy number.

     If a number is given, this is used as an upper limit for the total copy number genome wide (in all chr).
     If a file is given, it should have the format: chr max1 max2 max3 etc., e.g.

        1  2
        2  2
        3  8  2
        4  2
        etc.

     The first column is the chromosome, the next columns are the limits to be used for subclone 1, 2 etc.
     For subclones not specified, the limit in the last column is used. In the example above, subclone 1 has an upper limit of 8 total copies in chr3, for all other subclones and in all other chromosomes, the upper limit is 2. If only SNV data is provided (and `--avail-cn [file]` is not given), this is used to fix the total number of copies. If `--max-tcn` is not given, cloneHD uses the normal copy number for each chr.

     This number should be chosen conservatively, since it increases the
     HMM dimensionality and can open the possibility for spurious solutions. 

*    `--nmax [int:2]`  The maximum number of subclones to be tried.

     All subclone numbers from 0 to `nmax` will be used and the one
     with maximum BIC chosen for output.

*    `--force [int]`  Fix the number of subclones to be used.

*    `--trials [int:1]`  The number of independent optimizations.

     Global parameters are found numerically by local maximization of
     the total log-likelihood. The best result out of `trials` independent,
     randomly seeded,  runs will be used.

*    `--copynumber [file]`  Use copy number constraint for SNV data. 

     For a SNV data analysis, the cloneHD output file
     ending `*copynumber.txt` from a CNA(+BAF) run can be given here. Since
     the subclonal decomposition can be different for SNVs, this option
     ensures that the found solution is still consistent with the copy
     number profile.

### Fuzzy segmentation options

For data with persistence along the genome, a fuzzy segmentation can
be used based on the filterHD posterior jump probability (must be
`*jumps.txt` file). Data between potential jump sites, with a jump
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

*    `--clones [file]`  Use fixed mass(es) and/or subclonal frequencies. 

     Either all mass parameters, or all subclonal frequencies, or both 
     can be given (for each sample in the data input). The likelihoods
     and posteriors will be computed under these conditions. 
     Remaining parameters will be learned.
  
     Format: One line per sample. The first column, if greater than
     1.0, is interpreted as mass; the remaining as subclonal frequencies.

        30.0 0.64 0.12
        28.0 0.31 0.23

    More than one parameter set can be given (as a continued list). Then,
    only the likelihoods are computed and printed to a file ending
    `*llh-values.txt`.  Useful for mapping the log-likelihood surface
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

*    `--learn-priors [0/1:0]` For snv-mode only: if 1, then the parameters
     for the multiplicative genotype priors are learned.

*    `--maxcn-mask [file]`  Use max copy number mask.

     A maximum total copy number for individual
     chromosomes can be specified. This is useful, if large copy
     numbers are expected in only part of the genome. In
     chromosomes not specified in this file, the value of `--maxcn` is
     used as upper limit.

*    `--chr [file]`  Set normal copy numbers.

     The normal copy number for every single
     chromosome can be specified. This is needed only for non-human DNA. If not
     given, human DNA is assumed and the sex is inferred from the
     presence or absence of chr 24 (= chr Y) in the input data.

*    `--snv-fpr [double:1.0e-4]`  The false positive rate for SNVs,
     i.e. rate of SNV data points of genotype all-0.

*    `--snv-err [double:0.01]`  The typical frequency of mis-called variant reads.

*    `--snv-pen [double:0.01]`  The penalty for higher than expected
     genotypes.

*    `--baf-pen [double:1.0]`  The penalty for complex minor allele status.

## Bulk options

These options are only needed if the sequenced cell population is a mixture of a
diverse bulk, with known allele frequency profile, and a number of
subclones with unknown genotypes and frequencies. Allele frequency
data is input with `--snv`. Data segmentation can be used with
`--snv-jumps`.  Read depth data can also be specified with `--cna`. 

*    `--bulk-mean [double]`  The bulk allele frequency profile. 

     Must be a filterHD `*posterior.*.txt` file. Only the posterior mean is used.

*    `--bulk-prior [file]`  The bulk allele frequency profile. 

     Must be a filterHD `*posterior.*.txt` file. The whole posterior
     distribution is used (run filterHD with `--dist 1` to obtain it).

*    `--bulk-updates [int:0]`  The number of Bayesian updates of the
     bulk allele frequency profile (if `--bulk-prior` was used).

*    `--bulk-fix [double:0.0]`  Use a flat and fixed bulk allele
     frequency profile.

## Technical options

*    `--grid [int:300]`  The grid size for the pre-computed emission
     probabilities if fuzzy data segmentation is used.

# Tips and tricks

