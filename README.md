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

# Command line options

## input data
   
The first two columns of all three input file types are always chromosome and coordinate of each observation. Chromosomes are expected to be integers (X->23,Y->24).

* `--cna [file]` Read depth data file. 

  Format: For each sample, there are two additional columns with  (i) the read depth and (ii) the number of independent observations this is the sum of. For NGS data, use the median read depth per 1 kb as the highest resolution.

1 1000	  93  1	    75 	1
1 2000	  101 1	    81	1
1 3000	  105 1	    85  1
1 5000	  197 2	    156 2
1 6000	  91  1	    79  1

* `baf [file]` B-allele read count data file.

  Format: For each sample, there are two additional columns with (i) the number of reads of the minor allele and (ii) the total read depth at originally heterozygous loci.

1 1036	  43  90   28	   72
1 1287	  47  99   32	   80
1 2877	  30  100  36	   82
etc.

* `snv [file]` Somatic nucleotide variant read count data file.

  Format: For each sample, there are two additional columns with (i) the number of reads of the somatic variant allele and (ii) the total read depth at that locus.

1 1314	  12  92   28	   72
1 1287	  47  99   32	   80
1 2877	  30  100  36	   82
etc.

