# RUN filterHD & cloneHD FOR A SIMULATED EXAMPLE DATA SET

# fix the number of threads
export OMP_NUM_THREADS=4;

# input data
data="./test/data/"
results="./test/results/"
filterHD="./build/filterHD"
cloneHD="./build/cloneHD"

### filterHD ###
echo "*** filterHD ***"
#emission modes: 
# 1: Binomial
# 2: Beta-Binomial
# 3: Poisson
# 4: Negative Binomial

# normal read depth
normalRD="${data}/normal.RD.txt"
cmd="$filterHD --data $normalRD --mode 3 --pre ${results}/normal.RD --grid 100"
echo $cmd
$cmd
echo

# tumour read depth without bias
tumorRD="${data}/tumor.RD.txt"
cmd="$filterHD --data $tumorRD --mode 3 --pre ${results}/tumor.RD --grid 100"
echo $cmd
$cmd
echo

# tumour read depth with bias from normal
bias="${results}/normal.RD.posterior.1.txt"
cmd="$filterHD --data $tumorRD --mode 3 --pre ${results}/tumor.RD.bias --bias $bias --grid 100 --sigma 0 --jumps 1"
echo $cmd
$cmd
echo
tumorRDjumps="${results}/tumor.RD.bias.jumps.txt"

# tumour BAF
tumorBAF="${data}/tumor.BAF.txt"
cmd="$exec --data $tumorBAF --mode 1 --pre ${results}/tumour.BAF --grid 100 --sigma 0 --jumps 1 --reflect 1 --dist 1"
echo $cmd
$cmd
echo
tumorBAFjumps="${results}/tumor.BAF.jumps.txt"


### cloneHD ###
echo "*** cloneHD ***"

cmd="$cloneHD --cnv $tumorRD --baf $tumorBAF --pre ${results}/tumor --bias $bias --seed 123 --trials 1 --nmax 3 --force --maxcn 4 --cnv-jumps $tumorRDjumps --cnv-rnd 0.001 --baf-rnd 0.001 --min-occ 0.01 --min-jump 0.01 --print-all 0 --restarts 10"
echo $cmd
$cmd
echo