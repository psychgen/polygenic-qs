#simulate_data.R

#The goal of this script is to simulate data (in different scenarios) on which to prepare and test the analysis
#scripts for real MoBa data in the Informative Polygenic Questions Framework

#Variables for analysis should be:

# ADHD prs (continuous, normal)
# MDD prs (continuous, normal)
# ADHD symptoms (continuous, skewed)
# ADHD diagnoses (binary)
# Education problems (categorical)
# BMI (continuous, normal)
# Neg control outcome (continuous, normal)
# Neg control outcome 2 (binary)

set.seed(219057)

library(tidyverse)
library(faux)

# Start by simulating a df based on a correlation matrix for 8 continuous variables 
# (these will be the basis of liabilities for categorical/binary variables)

# adhd_prs mdd_prs adhd_sx adhd_dx edu_probs bmi neg1 neg2
simulate_our_data <- function(n=7000,
                              adhd_prs_sx_corr = 0.25,
                              adhd_prs_bmi_corr = 0.17,
                              adhd_prs_edu_corr = 0.1,
                              prop_sample_adhd = 0.07,
                              corr_prs = 0.1,
                              prop_miss = 0.05,
                              replicates =1){
  
  
  df <- data.frame()
  
  for(i in 1:replicates){
    
    message(paste0("Simulating dataset ",i," of ",replicates))
    
    # mean SD values from ALSPC studies (We don't want to be too exact to MoBa)
    
    dat <- rnorm_multi(n = n, 
                       mu = c(0, 10, 17.7),
                       sd = c(1, 5, 2.8),
                       r = c(adhd_prs_sx_corr,adhd_prs_bmi_corr,0.20)  , 
                       varnames = c("adhd_prs", "adhd_sx","bmi"),
                       empirical = FALSE)
    
    # make a genomewide significant SNPs version of the PGS (for MR)
    # specify this to correlate at 0.47 with "main" PGS: this comes from 
    # my paper (https://doi.org/10.1017/S0033291721003330) in which we checked the correlation between
    # adhd PGS at 10-8 with a principal component PGS (PGS-PC) used in analyses
    
    dat <- dat %>% 
      mutate(adhd_prs_gws = rnorm_pre(adhd_prs, mu = 0, sd = 1, r = 0.47))
    
    # make adhd_dx liability variable as a function of symptoms plus noise
    # skew distribution of sx by adding floor effect (will slightly attenuate corr)
    
    dat <- dat %>% 
      mutate(adhd_dx = adhd_sx + rnorm(nrow(dat),0, 0.1),
             adhd_sx = round(ifelse(adhd_sx<0,0,adhd_sx),0))
    
    # assign adhd dx on basis of liability and prop sample adhd
    
    thresh_adhd_dx <- quantile(dat$adhd_dx, 1-prop_sample_adhd)
    
    dat <- dat %>% 
      mutate(adhd_dx = ifelse(adhd_dx>= thresh_adhd_dx, 1, 0))
    
    #add mdd_prs edu_probs liability and BMI based on desired correlation with adhd:prs
    
    dat <- dat %>% 
      mutate(mdd_prs = rnorm_pre(adhd_prs, mu = 0, sd = 1, r = corr_prs),
             bmi = rnorm_pre(adhd_prs, mu = 17.7, sd = 2.8, r = adhd_prs_bmi_corr),
             edu_probs = rnorm_pre(adhd_prs, mu = 0, sd = 1, r = adhd_prs_edu_corr))
    
    #Simulate the covariates (one binary and one continuous)
    
    dat <- dat %>% 
      mutate(sex = factor(rbinom(nrow(dat),1,0.51), levels = c(0,1), labels= c("Female","Male")),
             age = round(rnorm_pre(adhd_sx, mu = 8, sd = 0.5, r = 0),1))
    
    
    #Remove the impossible BMI values
    dat$bmi[dat$bmi<12]<- NA
    
    
    #convert edu liability to categorical
    
    prob = c('low' = 80, 
             'moderate' = 100, 
             'high' = 40, 
             'very high' = 20)
    
    dat <- dat %>% 
      mutate(edu_probs = norm2likert(edu_probs,prob))
    
    #add in negative control variables, specifying their correlation at around 0
    
    dat <- dat %>% 
      mutate(neg1 = rnorm_pre(adhd_prs, mu = 0, sd = 1, r = rnorm(1,0,0.03)),
             neg2 = rnorm_pre(adhd_prs, mu = 0, sd = 1, r = rnorm(1,0,0.05))) %>% 
      mutate(neg2 = ifelse(neg2>quantile(neg2,0.8), 1,0))
    
    
    #add in some missingness at a variable proportion
    
    #messy_dat <- messy(dat, prop_miss, names(dat)) # fails with "some columns not in data table" error - why?
    
    #check observed correlations between numeric vars
    
    dat %>% select_if(is.numeric) %>% cor(use = "pairwise.complete.obs") %>% round(2)
    
    dat <- dat %>% 
      mutate(iid= seq(1:nrow(dat)), 
             replicate=i) %>% 
      select(iid,sex, age, adhd_prs, adhd_prs_gws, mdd_prs, adhd_sx, adhd_dx, bmi, edu_probs, neg1, neg2, replicate)
    
    df <- rbind(df,dat)
  }
  return(df)
  
}

#Simstudy data:

# Make this programmatic - so that everything is defined relative to the assoc. between ADHD_prs and ADHD_sx

adhd_prs_sx_corr = 0.10
adhd_prs_bmi_corr = adhd_prs_sx_corr*(1/3)
adhd_prs_edu_corr = adhd_prs_sx_corr*(1/2)

alldat=data.frame()
for(n in c(1000,3000,5000,7000)){
  
  message(paste0("working on scenario N=",n))
  dat <- simulate_our_data(n=n,
                           adhd_prs_sx_corr=adhd_prs_sx_corr,
                           adhd_prs_bmi_corr=adhd_prs_bmi_corr,
                           adhd_prs_edu_corr=adhd_prs_edu_corr,
                           replicates=1000)
  dat <- dat %>% 
    mutate(scenario=paste0("N=",n))
  alldat<-rbind(alldat,dat)
}




save(alldat, file= "data/simstudy_dat.RData")
