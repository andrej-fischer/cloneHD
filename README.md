# cloneHD and filterHD software

Pre-compiled executable binaries of filterHD and cloneHD for Mac OS X (64bit) can be found at <ftp://ftp.sanger.ac.uk/pub/teams/153/cloneHD/>. 

# Compilation  

To compile cloneHD, you need the GNU scientific library (GSL) v1.15 or later. Change the paths in the Makefile to your GSL installation location (if non-standard). Then type 

`$ make -f Makefile_cloneHD`

in the source directory. There will be two executables, `filterHD` and `cloneHD`, in `./build`. To view all the command line options, type e.g. `./build/cloneHD --print-options`.

