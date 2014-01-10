# cloneHD and filterHD software

The current stable release, including pre-compiled executable binaries of filterHD and cloneHD for Mac OS X (64bit), can be found at <ftp://ftp.sanger.ac.uk/pub/teams/153/cloneHD/>. 

# Run test with simulated data

In the source directory, you can test filterHD and cloneHD by running

`$ sh run-example.sh`

where you can see a typical workflow of analysing read depth and BAF data with a matched normal. To view a complete list of the command line options, type e.g. 

`$ ./build/filterHD --print-options`  or  `$ ./build/cloneHD --print-options`.

# Compilation  

To compile cloneHD yourself, you need the GNU scientific library (GSL) v1.15 or later. Change the paths in the Makefile to your GSL installation location (if non-standard). Then type 

`$ make -f Makefile_cloneHD`

in the source directory. The two executables, `filterHD` and `cloneHD`, will be in `./build`.

