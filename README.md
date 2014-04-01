# How to get cloneHD and filterHD?

The current stable release, including pre-compiled executable binaries
of filterHD and cloneHD for Mac OS X (64bit), can be found at:

`ftp://ftp.sanger.ac.uk/pub/teams/153/cloneHD/`  

The source code can also be downloaded here.

# Run a test with simulated data

If you download cloneHD from the ftp site above, you can test both filterHD and cloneHD by running

`$ sh run-example.sh`

where you can see a typical workflow of analysing read depth and BAF
data with a matched normal. All command line arguments are explained below.

# Compilation  

For Mac OS X (64bit), there are pre-compiled binaries available on the ftp server above. To compile cloneHD yourself, you need the GNU scientific library ([GSL](http://www.gnu.org/software/gsl/)) v1.15 or later. Change the paths in the Makefile to your GSL installation location (if non-standard). Then type 

`$ make`

in the source directory. The two executables, `filterHD` and
`cloneHD`, will be in `./build`. For debugging with gdb, use `make -f Makefile.debug`.

# What are cloneHD and filterHD for?

cloneHD is a software for reconstructing the subclonal structure of a
population from short-read sequencing data. Read depth
data, B-allele count data and somatic nucleotide variant (SNV) data can be
used for the inference. cloneHD can estimate the number of subclonal
populations, their fractions in the sample, their individual total copy number profiles, 
their B-allele status and all the SNV genotypes with high resolution.

filterHD is a general purpose probabilistic filtering algorithm for one-dimensional
discrete data, similar in spirit to a Kalman filter. It is a continuous state
space Hidden Markov model with Poisson or Binomial emissions and a
jump-diffusion propagator. It can be used for scale-free smoothing, 
fuzzy data segmentation and data filtering. 

![cna gof](/images/cna.gof.png "CNA goodness of fit")
![baf gof](/images/baf.gof.png "BAF goodness of fit")
![cna post](/images/cna.post.png "CNA posterior")
![cna real](/images/cna.real.png "CNA real profile")
![baf post](/images/baf.post.png "BAF posterior")
![baf real](/images/baf.real.png "BAF real profile")
![snv gof](/images/snv.gof.png "SNV goodness of fit")

Visualization of the cloneHD output for the simulated data set. From
top to bottom: 
(i) The bias corrected read depth data and the cloneHD
prediction (red).
(ii) The BAF (B-allele frequency), reflected at 0.5 and the cloneHD prediction (red).
(iii) The total copy number posterior.
(iv) The real total copy number profile.
(v) The minor copy number posterior.
(vi) The real minor copy number profile.
(vii) The observed SNV frequencies, corrected for local ploidy, and per genotype (SNVs are assigned ramdomly according to the cloneHD SNV posterior).
(All plots are created with Wolfram [Mathematica](http://www.wolfram.com/mathematica/).)

# Full documentation

The full documentation can be found in the `/docs/` subfolder. Click below.

*  [pre-filter](/docs/README-pre-filter.md)
*  [filterHD](/docs/README-filterHD.md)
*  [cloneHD](/docs/README-cloneHD.md)

# Tips and tricks

*  Pre-filtering of data can be very important. If filterHD predicts
   many more jumps than you would expect, it might be necessary to
   filter the data, removing very short segments (with
   `--filter-shortSeg 10`).

*  Make sure that the bias field for the tumor CNA data is
   meaningful. If a matched normal sample was sequenced with the same
   pipeline, its read depth profile, as predicted by filterHD, can be used as a
   bias field for the tumor CNA data. Follow the logic of the example
   given here.

*  filterHD can sometimes run into local optima. It might be useful to
   fix initial values for the parameters via `--jumpi [double]` etc.

*  By default, cloneHD runs with mass gauging enabled. This seems like
   an overkill, but is actually quite useful because you can see some
   alternative explanations during the course of the analysis.

*  Don't put too much weight on the BIC criterion. It was calibrated
   using simulated data. For real data, it should be supplied with
   common sense and biological knowledge. Use `--force [int]` to use a
   fixed number of subclones.

*  For exome sequencing data, the read depth bias can be enormous. Use rather, if
   available, the jumps seen in the BAF data for both CNA and BAF.


