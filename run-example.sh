# RUN filterHD & cloneHD FOR A SIMULATED EXAMPLE DATA SET

# fix the number of threads
export OMP_NUM_THREADS=1;

# input data
data="./test/data/"
results="./test/results/"
filterHD="./build/filterHD"
cloneHD="./build/cloneHD"

### filterHD ###
echo "*** filterHD ***"
echo

#emission modes: 
# 1: Binomial
# 2: Beta-Binomial
# 3: Poisson
# 4: Negative Binomial

# normal read depth
normalCNA="${data}/normal.cna.txt"
cmd="$filterHD --data $normalCNA --mode 3 --pre ${results}/normal.cna --grid 100 --rnd 0"
echo $cmd
$cmd
echo

# tumour read depth without bias
tumorCNA="${data}/tumor.cna.txt"
cmd="$filterHD --data $tumorCNA --mode 3 --pre ${results}/tumor.cna --grid 100 --rnd 0"
echo $cmd
$cmd
echo

# tumour read depth with bias from normal
bias="${results}/normal.cna.posterior.1.txt"
cmd="$filterHD --data $tumorCNA --mode 3 --pre ${results}/tumor.cna.bias --bias $bias --grid 100 --sigma 0 --rnd 0 --jumps 1"
echo $cmd
$cmd
echo
tumorCNAjumps="${results}/tumor.cna.bias.jumps.txt"

# tumour BAF
tumorBAF="${data}/tumor.baf.txt"
cmd="$filterHD --data $tumorBAF --mode 1 --pre ${results}/tumor.baf --grid 100 --sigma 0 --jumps 1 --reflect 1 --dist 1"
echo $cmd
$cmd
echo
tumorBAFjumps="${results}/tumor.baf.jumps.txt"


### cloneHD ###
echo "*** cloneHD ***"
echo

cmd="$cloneHD --cna $tumorCNA --baf $tumorBAF --pre ${results}/tumor --bias $bias --seed 123 --trials 1 --nmax 3 --force --maxcn 4 --cna-jumps $tumorCNAjumps --cna-rnd 0.0 --baf-rnd 0.0 --min-occ 0.01 --min-jump 0.01 --print-all 0 --restarts 20"
echo $cmd
$cmd
echo
cat ${results}/tumor.clonal.txt