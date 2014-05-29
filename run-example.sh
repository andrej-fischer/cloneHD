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
tumorSNV="${data}/tumor.snv.txt"
bias="${results}/normal.cna.posterior-1.txt"
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
 
# The normal read depth is analysed to estimate the technical read depth modulation. This will be later used to account
# for the bias field in cloneHD. In principal, jumps are not expected (so could set --jump 0). The simulations do not have 
# random emissions. 
    cmd="$filterHD --data $normalCNA --mode 3 --pre ${results}/normal.cna --rnd 0"
    echo $cmd
    $cmd
    echo

# The tumor read depth is first analysed without bias to get a benchmark for the LLH value. The result will not be used later. 
# In the tumor data, we do expect jumps, but we actually would like to learn the jumps only accounting for the bias field (below).
    cmd="$filterHD --data $tumorCNA --mode 3 --pre ${results}/tumor.cna --rnd 0"
    echo $cmd
    $cmd
    echo

# The tumor read depth is now analysed with the bias field from the matched normal. The diffusion constant is set to zero. 
# If left free, it should converge to a very small value. The jump rate could be slightly higher. The LLH should be higher than
# for the run above indicating the presence of the bias field. Now we are interested in the jumps.
    cmd="$filterHD --data $tumorCNA --mode 3 --pre ${results}/tumor.cna.bias --bias $bias --sigma 0 --rnd 0 --jumps 1"
    echo $cmd
    $cmd
    echo

# The tumor BAF data is analysed, mainly to get the emission parameters (shape, rnd) and jumps. In principle, there could be jumps
# visible in the BAF data, but not in the read depth (copy number neutral LOH within chromosomes). Diffusion should be switched off.
    cmd="$filterHD --data $tumorBAF --mode 1 --pre ${results}/tumor.baf --sigma 0 --jumps 1 --reflect 1 --dist 1 --rnd 0"
    echo $cmd
    $cmd
    echo
fi

if [ -z $part ] || [ $part -eq 2 ] 
then
### cloneHD ###
    echo "*** cloneHD ***"
    echo "True mass and cell fractions:" `cat test/data/clones.txt` 
    echo
    # The CNA and BAF data is analysed for subclonality.
    # Try varying the --min-jump, --force and --max-tcn values and try --mass-gauging 0. 
    # Try adding the SNV data to the mix.
    cmd="$cloneHD --cna $tumorCNA --baf $tumorBAF --pre ${results}/tumor --bias $bias --seed 123 --trials 2\
 --nmax 3 --force --max-tcn 4 --cna-jumps $tumorCNAjumps --baf-jumps $tumorBAFjumps --min-jump 0.01 --restarts 10 --mass-gauging 1" 
    echo $cmd
    $cmd
    echo
    cat ${results}/tumor.summary.txt
    echo

    # Using the information from above, the SNV data is analysed. Try what happens removing the --avail-cn and --mean-tcn options.
    cmd="$cloneHD --snv $tumorSNV --pre ${results}/tumorSNV --seed 123 --trials 2\
 --nmax 3 --force --max-tcn 4 --restarts 10 --mean-tcn ${results}/tumor.mean-tcn.txt --avail-cn ${results}/tumor.avail-cn.txt"    
    echo $cmd
    $cmd
    echo
    cat ${results}/tumorSNV.summary.txt
fi