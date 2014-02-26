# RUN filterHD & cloneHD FOR A SIMULATED EXAMPLE DATA SET

# fix the number of threads
export OMP_NUM_THREADS=4;

part=$1

# input data
data="./test/data/"
results="./test/results/"
filterHD="./build/filterHD"
cloneHD="./build/cloneHD"

normalCNA="${data}/normal.cna.txt"
tumorCNA="${data}/tumor.cna.txt"
tumorBAF="${data}/tumor.baf.txt"
bias="${results}/normal.cna.posterior.1.txt"
tumorCNAjumps="${results}/tumor.cna.bias.jumps.txt"
tumorBAFjumps="${results}/tumor.baf.jumps.txt"

### filterHD ###
if [ -z $part ] || [ $part -eq 1 ]
then
    echo "*** filterHD ***"
    echo    
#emission modes: 
# 1: Binomial
# 2: Beta-Binomial
# 3: Poisson
# 4: Negative Binomial
#    
# normal read depth
    cmd="$filterHD --data $normalCNA --mode 3 --pre ${results}/normal.cna --rnd 0"
    echo $cmd
    $cmd
    echo
# tumor read depth without bias
    cmd="$filterHD --data $tumorCNA --mode 3 --pre ${results}/tumor.cna --rnd 0"
    echo $cmd
    $cmd
    echo
# tumor read depth with bias from normal
    cmd="$filterHD --data $tumorCNA --mode 3 --pre ${results}/tumor.cna.bias --bias $bias --sigma 0 --rnd 0 --jumps 1"
    echo $cmd
    $cmd
    echo
# tumor BAF
    cmd="$filterHD --data $tumorBAF --mode 1 --pre ${results}/tumor.baf --sigma 0 --jumps 1 --reflect 1 --dist 1 --rnd 0"
    echo $cmd
    $cmd
    echo
fi

if [ -z $part ] || [ $part -eq 2 ] 
then
### cloneHD ###
    echo "*** cloneHD ***"
    echo
    cmd="$cloneHD --cna $tumorCNA --baf $tumorBAF --pre ${results}/tumor --bias $bias --seed 123 --trials 3\
 --nmax 3 --force --max-tcn 4 --cna-jumps $tumorCNAjumps --baf-jumps $tumorBAFjumps --min-jump 0.01 --restarts 20"    
    echo $cmd
    $cmd
    echo
    cat ${results}/tumor.clonal.txt
fi