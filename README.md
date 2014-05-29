# How to get cloneHD and filterHD?

The current stable release, as well as pre-compiled executable binaries 
for Mac OS X and GNU Linux (64bit), can be found [here](https://github.com/andrej-fischer/cloneHD/releases). The cloneHD software is undergoing rapid development. Watch/Star this repo to receive updates.

# Run a test with simulated data

After downloading cloneHD from the release site, you can test both filterHD and cloneHD by running

`$ sh run-example.sh`

where you can see a typical workflow of analysing read depth and BAF
data with a matched normal. All command line arguments are explained below.

# Compilation  

For Mac OS X and GNU Linux (64bit), pre-compiled binaries are available [here](https://github.com/andrej-fischer/cloneHD/releases). To compile cloneHD yourself, you need the GNU scientific library ([GSL](http://www.gnu.org/software/gsl/)) v1.15 or later. Change the paths in the Makefile to point to your local GSL installation (if non-standard). Then type 

`$ make`

in the `src` directory. The executables will be in `build`. For debugging with gdb, use `make -f Makefile.debug`.

# Report bugs

To report bugs, use the [issue](https://github.com/andrej-fischer/cloneHD/issues) interface of github.

# Full documentation

The full documentation can be found in the `/docs/` subfolder. Click below.

*  [pre-filter](/docs/README-pre-filter.md)
*  [filterHD](/docs/README-filterHD.md)
*  [cloneHD](/docs/README-cloneHD.md)

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

# Tips and tricks

*  Note: all input files are assumed to be sorted by genomic coordinate. With Unix, this
   can be guaranteed with `sort -k1n,1 -k2n,2 file.txt > sorted-file.txt`.

*  Pre-filtering of data can be very important. If filterHD predicts
   many more jumps than you would expect, it might be necessary to
   pre-filter the data, removing variable regions, outliers or very short 
   segments (use programs the `pre-filter` and `filterHD`).

*  Make sure that the bias field for the tumor CNA data is
   meaningful. If a matched normal sample was sequenced with the same
   pipeline, its read depth profile, as predicted by filterHD, can be used as a
   bias field for the tumor CNA data. Follow the logic of the example
   given here.

*  If the matched-normal sample was sequenced at lower coverage than the tumor sample, 
   it might be necessary to run filterHD with a higher-than-optimal diffusion constant 
   (set with `--sigma [double]`) to obtain a more faithful bias field. Otherwise, the 
   filterHD solution is too stiff and you loose bias detail.

*  filterHD can sometimes run into local optima. In this case, it might be useful to
   set initial values for the parameters via `--jumpi [double]` etc.

*  By default, cloneHD runs with mass-gauging enabled. This seems wasteful,
   but is actually quite useful because you can see some alternative explanations
   during the course of the analysis.

*  Don't put too much weight on the BIC criterion. It was calibrated
   using simulated data. For real data, it should be supplemented with
   common sense and biological knowledge. Use `--force [int]` to use a
   fixed number of subclones and `--max-tcn [int]` to set the maximum possible total
   copy number.

*  If high copy numbers are expected only in a few chromosomes, you can increase performance
   by using the `--max-tcn [file]` option to specify per-chromosome upper limits.

*  For exome sequencing data, the read depth bias can be enormous. The filterHD estimate 
   of the bias field might not be very useful, especially in segmenting the tumor data.
   Use rather, if available, the jumps seen in the BAF data for both CNA and BAF data
   (give the BAF jumps file to both `--cna-jumps` and `--baf-jumps`).

#How to cite

The cloneHD and filterHD software is free under the GNU General Public License v3.
If you use this software in your work, please cite the accompanying publication:

Andrej Fischer, Ignacio Vazquez-Garcia, Christopher J.R. Illingworth and Ville Mustonen. High-definition reconstruction of subclonal composition in cancer. Cell Reports (2014), http://dx.doi.org/ 10.1016/j.celrep.2014.04.055
