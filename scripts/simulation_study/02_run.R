#02_run.R

#The goal of this script is to run the analyses in the simulation study
#part of the "polygenic questions" project
#This involves NHST and equivalence tests in a polygenic score analysis
#and an MR analysis, on data simulated in 01_simulate.R


library(tidyverse)
library(parameters)

# Source an external script with a procedure for scaling SESOIs for PGS
source("./scripts/sesoi_scale.R")

# Load simulated data
load("data/simstudy_dat.RData")

# Create an object to hold results from the loop overall
outer_res = data.frame()

# We loop through scenarios and replicates (simulated datasets)

for(scen in unique(alldat$scenario)){
  
  message(paste0("Current scenario is ", scen))
  
  #Create an object to hold results from the inner loop
  inner_res = data.frame()
  
  for(repl in 1:max(alldat$replicate)){
    if(round(repl,-2)>round(repl-1,-2)){
      message(paste0("Working on replicates in range ",round(repl,-2)-100, " to ", round(repl,-2)))
    }
    
    dat <- alldat %>%  
      filter(replicate==repl,
             scenario==scen)
    
    ### Polygenic score analyses
    # 1A) NHST: Linear regression (adhd_sx, bmi and neg1 on adhd_prs)
    ######
    
    lmods<-list(
      
      lm(scale(adhd_sx) ~ scale(adhd_prs), data =dat) ,
      lm(scale(bmi) ~ scale(adhd_prs), data =dat) ,
      lm(scale(neg1) ~ scale(adhd_prs), data =dat)
      
    )
    
    lmods_cis<- map(lmods, function(x){
      as.data.frame(confint(x))
    })
    lmods_summs <- map(lmods, summary)
    lmods_res <- map(lmods_summs, function(x){
      as.data.frame(x$coefficients)
    })
    
    lmods_res_ts<-lmods_res %>% 
      bind_rows() %>% 
      rownames_to_column("PGS") %>% 
      filter(str_detect(PGS, "prs" )) %>% 
      mutate(PGS= paste0(str_remove_all(PGS, "scale|[[:punct:]]|[[:digit:]]|prs"),"_PGS"),
             outcome = c("adhd_sx","bmi","neg1"),
             replicate = unique(dat$replicate),
             analysis = "pgs",
             test = "two-sided") %>% 
      select(PGS, outcome, everything()) %>% 
      rename(p=`Pr(>|t|)`) %>% 
      bind_cols(lmods_cis %>% 
                  bind_rows() %>% 
                  rownames_to_column("PGS") %>% 
                  filter(!str_detect(PGS,"Intercept")) %>% 
                  select(-PGS))
    
    ## Compute p-values for one-sided tests:
    
    lmods_res_os <- lmods_res_ts %>% 
      mutate(p = ifelse(Estimate>0,p/2,1-(p/2)),
             test ="one-sided")
    
    ## combine
    
    lmods_res <- lmods_res_ts %>% 
      bind_rows(lmods_res_os)
    
    # 1B) Equivalence test: Linear regression (adhd_sx, bmi and neg1 on adhd_prs)
    ######
    # For ADHD-ADHD, scaled version of family history benchmark
    # e.g., https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5809323/
    # Parental history of ADHD symptoms
    
    SESOI_perfect_PGS_OR <- 4.4 
    
    sesoi_adhd <- scale_SESOI(value=SESOI_perfect_PGS_OR, stat="OR")
    
    
    equiv_adhd <- parameters::equivalence_test(lmods[[1]], range= c(-sesoi_adhd,sesoi_adhd),p_values = T, rule="classic")  %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL) 
    
    # For ADHD-BMI, small telescopes from 10.1093/ije/dyaa214 (3320 DZ pairs in TEDS - 
    # we are ignoring the multilevel structure of their data for now)
    # unscaled (i.e., PGS is from the same GWAS and the presumption is that we
    # care about the size of this effect _now_, not in some hypothetical future
    # with a perfect PGS [cf. ADHD-ADHD above]):
    library(pwr)
    small_tel <-pwr.r.test(n = 6640, sig.level = 0.05, power = 0.33, alternative = "two.sided")
    sesoi_bmi <- small_tel$r
    
    
    equiv_bmi <- parameters::equivalence_test(lmods[[2]], range= c(-sesoi_bmi,sesoi_bmi),p_values = T, rule="classic")  %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL) 
    
    # For the unspecified negative control, we can use the generic version of SESOI
    # from Kruschke & Liddell 10.3758/s13423-016-1221-4; that is, half of a Cohen's
    # D "small effect" = 0.2/2 = 0.1
    # This time, again we do not scale (we care about the here and now):
    sesoi_neg <- effectsize::d_to_r(0.1)
    
    equiv_neg <- parameters::equivalence_test(lmods[[3]], range= c(-sesoi_neg,sesoi_neg),p_values = T, rule="classic") %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL) 
    
    equiv_res <- equiv_adhd %>% 
      bind_rows(equiv_bmi) %>% 
      bind_rows(equiv_neg) %>% 
      filter(Parameter!="(Intercept)")
    
    
    ### Mendelian randomization analyses
    # 2A) NHST: 2-stage least squares regression (ADHD on BMI, Neg1)
    ######
    
    ##Stage one
    stage1 <- lm(adhd_sx~scale(adhd_prs_gws), data=dat)
    #generate predicted values of X given Z
    dat$pred <- predict(stage1)
    
    
    mrmods <- list(
      ##Stage two - BMI
      stage2_bmi <-lm(scale(bmi)~scale(pred), data=dat),
      
      ##Stage two - Neg1
      stage2_neg1 <-lm(scale(neg1)~scale(pred), data=dat)
    )
    
    mrmods_cis<- map(mrmods, function(x){
      as.data.frame(confint(x))
    })
    mrmods_summs <- map(mrmods, summary)
    mrmods_res <- map(mrmods_summs, function(x){
      as.data.frame(x$coefficients)
    })
    
    mrmods_res<-mrmods_res %>% 
      bind_rows() %>% 
      rownames_to_column("exposure") %>% 
      filter(str_detect(exposure, "pred" )) %>% 
      mutate(exposure= ifelse(str_detect(exposure,"pred"), 
                              "PGS_predicted_adhd_sx",
                              str_remove_all(exposure, "[[:punct:]]|[[:digit:]]|prs")),
             outcome = c("bmi","neg1"),
             replicate = unique(dat$replicate),
             analysis = "mr",
             test = "two-sided") %>% 
      select(exposure, outcome, everything()) %>% 
      rename(p=`Pr(>|t|)`) %>% 
      bind_cols(mrmods_cis %>% 
                  bind_rows() %>% 
                  rownames_to_column("exposure") %>% 
                  filter(!str_detect(exposure,"Intercept")) %>% 
                  select(-exposure))
    
    # 2B) Equivalence test: MR 2SLS (bmi and neg1 on adhd)
    ######
    
    # For the ADHD-BMI analysis, we use a (hypothesized) threshold of clinical interest
    # i.e., on consultation with clinical experts, it is determined that a
    # 1 symptom change in ADHD predicting a 0.1kg/m2 change in BMI is the minimal effect
    # that would be clinically relevant. 
    # This is on the observed scale, and needs to be standardised for use in the equiv test:
    sesoi_mr_bmi <- (0.1/sd(dat$bmi, na.rm=T))*sd(dat$pred, na.rm=T)
    # NB: we use the sd of the _PGS-predicted_ ADHD_sx variable for this standardisation,
    # *not* the sd of the observed ADHD_sx. This is because the sd of the observed variable
    # is larger and we do not (in the context of MR) want to penalise the model (in terms
    # of making the SESOI a higher bar) on the basis of the limitations of PGS prediction of
    # the exposure. For equivalence testing in MR, the SESOI is set on the basis of the 
    # minimally relevant relationship between exposure and outcome (this is also why we don't
    # need to consider scaling the SESOI on the basis of how far along in the genomic discovery
    # process we are for this trait, as in the PGS analyses above); the fact that exposure
    # is instrumented by a polygenic score (and how well this is done) should not affect this.
    # As per MR guidelines, as long as the PGS-exposure association is robust enough, its
    # size does not matter
    
    mr_equiv_bmi<- parameters::equivalence_test(mrmods[[1]], range= c(-sesoi_mr_bmi,sesoi_mr_bmi),p_values = T, rule="classic")   %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL) 
    
    
    # For the unspecified negative control, we can use the generic version of SESOI
    # from Kruschke & Liddell 10.3758/s13423-016-1221-4; that is, half of a Cohen's
    # D "small effect" = 0.2/2 = 0.1
    
    sesoi_mr_neg <- effectsize::d_to_r(0.1)
    
    mr_equiv_neg<- parameters::equivalence_test(mrmods[[2]], range= c(-sesoi_mr_neg,sesoi_mr_neg),p_values = T, rule="classic")     %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL)
    
    mr_equiv_res <- mr_equiv_bmi %>% 
      bind_rows(mr_equiv_neg) %>% 
      filter(Parameter!="(Intercept)")
    
    ## Join together all results and prepare for output
    
    all_res <- lmods_res %>% 
      select(exposure=PGS, everything()) %>% 
      filter(!(outcome=="adhd_sx"&test=="two-sided"),
             !(outcome!="adhd_sx"&test=="one-sided")) %>% 
      arrange(outcome) %>% 
      bind_rows(mrmods_res) %>% 
      rename(ci95_low=`2.5 %`,
             ci95_high=`97.5 %`) %>% 
      bind_cols(equiv_res %>% 
                  bind_rows(mr_equiv_res) %>% 
                  select(ci90_low = CI_low, ci90_high = CI_high, ROPE_low, ROPE_high, ROPE_Percentage, ROPE_Equivalence) ) %>% 
      mutate(scenario = scen)
    
    inner_res <- rbind(inner_res,all_res)
    
    
    
  }
  
  outer_res <- rbind(outer_res, inner_res)
}


save(outer_res, file="./output/simulation_study_results.RData")
