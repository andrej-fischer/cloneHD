# To do list for filterHD/cloneHD

# Bugs/issues to be fixed

## filterHD

## cloneHD

*  with `--purity` and `--clones` (mass only), `--nmax` is ignored
*  if using BAF-jumps for CNA: mapping to be checked
*  print full if `--*jump` is (also) given, only events if `--*jumps`
*  check memory leaks

# Features to be added in future releases

## filterHD

*  filter loci incompatible with emission model via 2-state HMM
*  do Baum-Welch
*  `--filter-shortSeg [int]` via posterior jump-prob (not via pmean)!

## cloneHD

*  bias field as distribution
*  re-think update_snp_fixed (hashing)
*  different nLevels per chr via max-tcn (for memory only)

## both

*  string chromosome ids
