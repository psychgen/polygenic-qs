# Calculating and varying the SESOI
library(tidyverse)
library(effectsize)
library(parameters)

## PART 1: Setting and scaling an effect size of interest based on theory

# First, specify the effect size for a "perfect" PGS in an appropriate format
# e.g., if classification is the aim could use AUC

# Should we provide recommended benchmarks? E.g., for AUC, try to find a baseline model
# such as family history, and use the AUC from that as a SESOI? Of course, possible that
# PGS contributes additional meaningful prediction, even if not > on it's own so this
# would be quite strict...
# e.g., https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5809323/
# Parental history of ADHD symptoms

SESOI_perfect_PGS_OR <- 4.4   

# Or can simply use cohen's d benchmark effect sizes (0.2 is a small effect,
# so Krushke https://doi.org/10.1177/2515245918771304 assumes that half of that
# could be considered the lower bound for meaningfulness)
SESOI_perfect_PGS_d <- 0.1

#analysis defines the output format of the SESOI
#default values are for ADHD
scale_SESOI <- function(value,
                        stat,
                        target_trait_pgs_r2 = 0.055,
                        snp_h2 = 0.22,
                        prop_SNP_h2_exp_GWAS=0.1,
                        sample_size_GWAS=20000,
                        analysis="linear_regression",
                        scale=TRUE){
  if(analysis == "linear_regression"){
    
    if(stat=="OR"){
      d_eff <- effectsize::oddsratio_to_d(value, log= FALSE)
    }else if(stat=="r"){
      d_eff <- effectsize::r_to_d(value)
    }else if(stat == "d"){
      d_eff = value
    }else{
      stop(paste0("Effect size type: ", stat, " not recognised."))
    }
    sesoi <- effectsize::d_to_r(d_eff)
  }
  
  # Whatever form the sesoi is in, scale relative to some estimation of how far the 
  # current PGS is from some perfect PGS - I think the ideal for doing this is to weight
  # by some heuristically derivable version of the key elements of this paper:
  # https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008612
  # proportion of SNP heritability explained, discoverability, and polygenicity
  
  # So...
  # (sesoi * prop_SNP_h2_exp_GWAS)
  # ...would be the scaled SESOI if all PGS essentially maxxed out their prediction at the 
  # genomewide significant level - i.e., polygenicity was high and discoverability very low
  
  # What about using the r2 from the (relatively standard) out-of-sample prediction included
  # in the GWAS paper as a proportion of total SNP_h2 (or even twin h2?)
  if(scale==TRUE){
  scaled_sesoi = sesoi * (target_trait_pgs_r2/snp_h2)
  }else{
    scaled_sesoi=sesoi
  }
  return(scaled_sesoi)
}

scale_SESOI(SESOI_perfect_PGS_OR, stat="OR")
scale_SESOI(SESOI_perfect_PGS_d, stat="d")
## PART 2: Setting and scaling an effect size of interest based on "small telescopes"
# NB, the effect may first be calculated, then scaled up to the perfect PGS, then
# scaled back down, if different GWAS are used in the two studies 
 