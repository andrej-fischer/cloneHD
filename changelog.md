# changelog for cloneHD/filterHD

## v1.17.2 / 27.03.2014

*  new output: posterior per subclone, goodness of fit (GOF) per
   segment
*  changed file name `*clonal.txt` -> `*summary.txt`
*  filterHD STDOUT includes now GOF per sample
*  cloneHD `*summary.txt` includes now GOF per sample
*  changed `_` to `-` in all file names
*  fixed bug: BAF now symmetrized only in per-subclone-posterior
*  new CNA prior to penalize homozygous deletions `--cna-pen [double:0.9]`

## v1.17.1 / 01.03.2014

*  BAF posterior symmetrized for output
*  CNA transition matrix penalizes clones with zero copies of a segment
*  fixed bug in SNV prior computation
*  added pre-processor directives for conditional openMP compilation

## v1.17.0 / 25.02.2014 major release

### changed the way SNV priors are computed:

*  if CNA given: SNV prior informed by CNA posterior
*  if CNA+BAF given, SNV prior informed by BAF+CNA posterior
*  if SNV only and `--max-tcn` not given, assumes all chr to be
   all-normal, mean total c.n. to be normal; SNV prior parameters can
   be learned with `--learn-priors 1`.
*  if SNV only and `--max-tcn [int/file]` is given, this data is used
   to fix the total c.n. per chr and subclone; mean total c.n. is
   calculated on the fly; SNV prior parameters can be learned with
   `--learn-priors 1`.
*  if SNV only and `--max-tcn [int/file]` and `--avail-cn [file]` are
   given, SNV prior is calculated according to c.n. availability.

### more changes

*  changed option `--copynumber [file]` to  `--mean-tcn [file]`
*  new option  `--avail-cn [file]`
*  changed option `--maxcn [int:4]` to `--max-tcn [file/int]`
*  changed option `--snv-err [double]` to `--snv-fpfreq [double]`
*  changed option `--snv-fpr [double]` to `--snv-fprate [double]`
*  output file `*used-tcn.txt` to `*used_mean_tcn.txt`
*  output file `*copynumber.txt` to `*mean_tcn.txt`
*  new output file `*available_cn.txt`
*  changed `sample` to `chr` in cloneHD output files
*  slimmed down output of `--print-options`.
*  split clone.cpp into components clone-*.cpp
*  split off cloneHD-inference.cpp
*  new Makefile

## v1.16.7 / 19.02.2014

*  fixed bug in SNV w/ corr mode when --bulk-fix is used
*  introduced different grid sizes for CNA, BAF and SNV
*  fixed bug in Clone::get_interpolation(), at the boundaries
*  fixed bug in Clone::trapezoidal() (affected --bulk-prior vs --bulk-mean consistency)

## v1.16.6 / 12.02.2014

*  fixed major bug for SNV false positive emission rate and prior
*  introduced new functions:  Clone::update_snv_site_ncorr/fixed/nfixed()
*  fixed bug in SNV prior from CNA/BAF posterior computation (BAF normalization)
*  false positive SNV prior now includes P(c=all-zero)
*  fixed bug in used cn output
*  all-zero "observations" in SNV input (w/o corr) are ignored (and not printed!)
*  fixed bug in filterHD: all-zero observations are always retained.

## v1.16.5 / 07.04.2014

*  fixed major bug when CNA, BAF and SNV data used with males (X,Y with only one copy)
*  fixed bug in Clone::snv_prior_from_cna_baf_post()
*  fixed bug in posterior output for BAF and SNV
*  introduced prior masking for all update functions
*  introduced `--maxcn_mask [file]` option to limit total c.n. per chromosome
*  static linking of both libgcc and libstdc++ for increased portability

## v1.16.4 / 30.01.2014

*  filterHD: if `--reflect 1`, use only posterior in [0,0.5] for mean/std-dev
*  fixed bug with `--bulk-fix 0.0`

## v1.16.3 / 13.01.2014

*  introduced the option `--mass-gauging [0/1:1]` to switch off the mass gauging for cna data.

## v1.16.2 / 12.01.2104

*  snp -> snv and cnv -> cna in all code
*  introduced `--chr [file]`, candidate masses are computed via majority normal copy number

## v1.16.1 / 10.01.2014

*  cnv to cna  for all command line options
*  cnv to cna in all output file names and content
*  filterHD stdout modified

## v1.16 / 03.01.2014

*  first stable release of cloneHD
