# changelog for cloneHD/filterHD

v1.17.0 / to come

*  major release
*  changed the SNV prior computation from `--copynumber
	[file]` to available copy numbers from `--avail-cn [file]` and
	mean total copy numbers from `--mean-tcn [file]`.
*  slimmed down output of `--print-options`.
*  `--snv-err` to `--snv-fpfreq` and `--snv-fpr` to `--snv-fprate`
*  normal_copy: int* -> map<int,int>
*  introduced maxcn_per_clone

v1.16.7 / 19.02.2014

*  fixed bug in SNV+corr mode when --bulk-fix is used
*  introduced different grid sizes for CNA, BAF and SNV
*  fixed bug in Clone::get_interpolation(), at the boundaries
*  fixed bug in Clone::trapezoidal() (affected --bulk-prior vs --bulk-mean consistency)
	
v1.16.6 / 12.02.2014

*  fixed major bug for SNV false positive emission rate and prior
*  introduced new functions update_snv_site_ncorr/fixed/nfixed()
*  fixed bug in SNV prior from CNA/BAF posterior computation
*  false positive SNV prior now includes P(c=all-zero)
*  fixed bug in used cn output
*  all-zero "observations" in SNV input are ignored (and not printed!)
*  fixed bug in filterHD: all-zero obs are always retained.

v1.16.5 / 07.04.2014

*  fixed major bug when CNA, BAF and SNV data used with males (X,Y with only one copy)
*  fixed bug in Clone::snv_prior_from_cna_baf_post()
*  fixed bug in posterior output for BAF and SNV
*  introduced prior masking for all update functions
*  introduced --maxcn_mask [file] option to limit total c.n. per chromosome
*  static linking of both libgcc and libstdc++ for increased portability

v1.16.4 / 30.01.2014

*  filterHD: if --reflect 1, use only posterior in [0,0.5] for mean/std-dev
*  fixed bug with --bulk-fix 0.0

v1.16.3 / 13.01.2014

*  introduced the option --mass-gauging [0/1:1] to switch off the mass gauging for cna data.

v1.16.2 / 12.01.2104

*  snp -> snv and cnv -> cna in all code
*  introduced `--chr [file]`, candidate masses are computed via majority normal copy number

v1.16.1 / 10.01.2014

*  cnv to cna  for all command line options
*  cnv to cna in all output file names and content
*  filterHD stdout modified

v1.16 / 03.01.2014

*  first stable release of cloneHD
