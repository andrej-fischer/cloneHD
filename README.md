To compile cloneHD, you need the GNU scientific library (GSL) v1.15 or later. Change the paths in the Makefile to your GSL installation location (if non-standard). Then type 

>make -f Makefile_cloneHD

in this directory. There will be two executables, filterHD and cloneHD, in ./build. To view all the command line options, type ./build/cloneHD --print-options.

