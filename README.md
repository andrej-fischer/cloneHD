# cloneHD and filterHD software

The current stable release, including pre-compiled executable binaries of filterHD and cloneHD for Mac OS X (64bit), can be found at <ftp://ftp.sanger.ac.uk/pub/teams/153/cloneHD/>. 

# Run test with simulated data

If you download cloneHD from the ftp site above, you can test both filterHD and cloneHD by running

`$ sh run-example.sh`

where you can see a typical workflow of analysing read depth and BAF data with a matched normal. To view a complete list of the command line options, type e.g. 

`$ ./build/filterHD --print-options`  or  `$ ./build/cloneHD --print-options`.

# Compilation  

To compile cloneHD yourself, you need the GNU scientific library (GSL) v1.15 or later. Change the paths in the Makefile to your GSL installation location (if non-standard). Then type 

`$ make -f Makefile_cloneHD`

in the source directory. The two executables, `filterHD` and `cloneHD`, will be in `./build`.

# filterHD command line arguments

## Data input 

## Typical usage options

## Parameter options

# cloneHD command line arguments

## Typical usage options

Format of input files: the first two columns of all three input file
    types are always chromosome and coordinate of each observation. Chromosomes
    are expected to be integers (X->23,Y->24). 

*   `--cna [file]` Read depth data file. 

    Format: For each sample, there are two additional columns with  (i) the read depth and (ii) the number of independent observations this is the sum of. For NGS data, use the median read depth per 1 kb as the highest resolution.

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

*    `--cna-jumps [file]`
*    `--baf-jumps [file]`
*    `--snv-jumps [file]`

## Parameter options

*    `--cna-rnd [double:0.0]`
*    `--baf-rnd [double:0.0]`
*    `--snv-rnd [double:0.0]`

*    `--cna-jump [double:-1.0]`
*    `--baf-jump [double:-1.0]`
*    `--snv-jump [double:-1.0]`

*    `--cna-shape [double:inf]`
*    `--baf-shape [double:inf]`
*    `--snv-shape [double:inf]`

*    `--snv-err [double:-1.0]`
*    `--snv-fpr [double:-1.0]`

*    `--snv-pen [double:-1.0]`
*    `--baf-pen [double:-1.0]`

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
    only likelihoods are computed and output to file ending
    `*llh-values.txt`. 

*    `--purity [file]`  Fixed purities, lower bounds for the sum of
     subclonal frequencies. One line per sample.

*    `--restarts [int:10]`  The number of perturbations in local random
     search mode.

     This simple random search routine is used: after finding a local
     maximum of LLH, the best solution is perturbed and a new optimum
     is sought. 

*    `--seed [int]`  A fixed seed to make inferences reproducible.

*    `--min-occ [double:0.01]`
*    `--min-jump [double:0.01]`

*    `--print-all [0/1:0]`
*    `--mass-gauging [0/1:1]`
*    `--learn-priors [0/1:1]`

*    `--maxcn-mask [file]`
*    `--chr [file]`

## Bulk options

*    `--bulk-mean [double]`
*    `--bulk-prior [file]`
*    `--bulk-updates [int]`
*    `--bulk-fix [double:0.0]`

## Technical options

*    `--grid [int:300]`
