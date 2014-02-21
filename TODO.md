# To do list for filterHD/cloneHD

# Bugs/issues to be fixed

## filterHD

## cloneHD

*  check that CNA entry prior is used
*  with `--purity` and `--clones` (mass only), `--nmax` is ignored
*  check purity message ("lambda") 
*  check posterior output ("samples" vs. chr)
*  if use BAF-jumps for CNA: mapping to be checked
*  print full if `--*jump` is (also) given, only events if `--*jumps`
*  check memory leaks
*  check BIC complexity

# Features to be added in future releases

## filterHD

*  filter loci incompatible with emission model via 2-state HMM
*  do Baum-Welch
*  `--filter-shortSe`g via posterior jump-prob (not via pmean)!

## cloneHD

*  for output, sort clones by size in sample 1
*  change to --mean-tcn --avail-cn for the snv-mode
*  bias field as distribution
*  print also posterior per subclone
*  re-think update_snp_fixed (hashing)
