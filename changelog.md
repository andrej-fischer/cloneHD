# changelog for cloneHD/filterHD

1.16.4 / ???
* filterHD: if --reflect 1, use only posterior in [0,0.5] for mean/std-dev 
* fixed bug with --bulk-fix 0.0

1.16.3 / 13.01.2014
* introduced the option --mass-gauging [0/1:1] to switch off the mass gauging for cna data.

1.16.2 / 12.01.2104
* snp -> snv and cnv -> cna in all code
* introduced `--chr [file]`, candidate masses are computed via majority normal copy number

1.16.1 / 10.01.2014
* --cnv* to --cna*  for command line options
* cnv to cna in all output file names and content
* filterHD stdout modified

1.16 / 03.01.2014
* first stable release of cloneHD
