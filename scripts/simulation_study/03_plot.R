#03_plot.R

#The purpose of this script is to visualise the results of the simulation study
#for which data are generated in 01_simulate and analysed in 02_run

library(tidyverse)
library(wesanderson)
names(wes_palettes)

load("./output/simulation_study_results.RData")

set.seed(213490)

# Not convinced that the automated decision-making in equivalence_test is doing what I think it should...
# I believe it is using 90%CIs for the NHST, not only for the equiv test - this is wrong (except in the
# case of ADHD Sx, which is a one-sided test. So will do my own... 


outer_res <- outer_res %>% 
  mutate(ci95_low=ifelse(outcome=="adhd_sx",ci90_low,ci95_low), #Adjust for the fact that ADHD is a one-sided test
         ci95_high=ifelse(outcome=="adhd_sx",ci90_high,ci95_high),
         `NHST+Equivalence`=factor(case_when(p<0.05 & (ci90_low<ROPE_low|ci90_high>ROPE_high) ~ "H0 rejected",
                                      p<0.05 & (ci90_low>ROPE_low&ci90_high<ROPE_high) ~ "H0 accepted",
                                      p>=0.05 & (ci90_low<ROPE_low|ci90_high>ROPE_high) ~ "Undecided",
                                      p>=0.05 & (ci90_low>ROPE_low&ci90_high<ROPE_high) ~ "H0 accepted")))
         

# Plot the estimates, 95 & 90% CIs, ROPE boundaries, and equivalence test conclusion for a random
# sample of replicates from the PGS analysis 

mini= outer_res %>% 
  filter(replicate%in%(1:50),
         analysis=="pgs") %>% 
  mutate(outcome=factor(outcome, 
                        levels = c("adhd_sx","bmi","neg1"),
                        labels = c("ADHD\nsymptoms", "BMI", "Negative\ncontrol")))

pd<- 0.2

ggplot(data=mini, aes(colour=`NHST+Equivalence`,fill=`NHST+Equivalence`), alpha = 0.8)+
  geom_rect(aes(xmin = ROPE_low, xmax = ROPE_high, ymin = -Inf, ymax = Inf),colour=NA,fill="grey80", alpha = 0.4)+
  geom_vline(aes(xintercept=0), linetype=4, colour="grey40" )+
  geom_vline(aes(xintercept=ROPE_low), linetype=2, colour="grey60")+
  geom_vline(aes(xintercept=ROPE_high), linetype=2, colour="grey60")+
  geom_errorbarh(aes(xmin=ci95_low, xmax=ci95_high, y=replicate),size=0.4,height=0)+
  geom_errorbarh(aes(xmin=ci90_low, xmax=ci90_high, y=replicate),size=0.8,height=0)+
  geom_point(aes(x=Estimate,y=replicate),shape=21,colour ="black", size=1)+
  facet_grid(scenario~outcome, scales = "free", switch = "y")+
  theme_minimal()+
  theme(text = element_text(size=18),
        axis.text.y = element_blank(),
        axis.text.x =element_text(angle=90),
        axis.title.y = element_blank(),
        strip.text.y = element_text(angle=0),
        legend.position = "bottom")+
  scale_color_discrete(type= wes_palette("Darjeeling1", 3, type = c("discrete")))+
  scale_fill_discrete(type= wes_palette("Darjeeling1", 3, type = c("discrete")))+
  scale_x_continuous("Effect estimate\n(standardized beta)")

# Plot the proportion of each conclusion in a NHST vs NHST+Equivalence situation for all
# scenarios, both analyses, and provide the average ROPE% (or SGPV)


summd= outer_res %>% 
  group_by(scenario, analysis, outcome, `NHST+Equivalence`) %>% 
  summarize(prop=n()/1000) %>%  
  mutate(test="NHST+Equivalence") %>% 
  rename(conclusion=`NHST+Equivalence`)

summd2= outer_res %>% 
  mutate(`conclusion` = factor(case_when(p<0.05~ "H0 rejected", 
                                        TRUE ~"Fail to\nreject H0"))) %>% 
  group_by(scenario, analysis, outcome, conclusion) %>% 
  summarize(prop=n()/1000) %>% 
  mutate(test="NHST only") 

summd3= outer_res %>% 
  group_by(scenario, analysis, outcome) %>% 
  summarize(avg_SGPV=mean(ROPE_Percentage),
            se_SGPV=sd(ROPE_Percentage)/sqrt(length(unique(outer_res$replicate))),
            `p:SGPV`=mean(p)/mean(ROPE_Percentage)) %>% 
  mutate(outcome=factor(outcome, 
                        levels = c("adhd_sx","bmi","neg1"),
                        labels = c("ADHD\nsymptoms", "BMI", "Negative\ncontrol"))) %>% 
  filter(analysis=="mr") 
  
summd_all = summd %>% 
  bind_rows(summd2)%>% 
  ungroup() %>% 
  mutate(outcome=factor(outcome, 
                        levels = c("adhd_sx","bmi","neg1"),
                        labels = c("ADHD\nsymptoms", "BMI", "Negative\ncontrol"))) %>% 
  filter(analysis=="mr") %>% 
  mutate(conclusion=factor(conclusion, levels= c("Fail to\nreject H0","H0 accepted","Undecided", "H0 rejected")),ordered=T)

dput(c(wes_palette("Darjeeling1", 3, type = c("discrete")), "grey70" ))

pal <-c("#FF0000", "#00A08A", "#F2AD00",  "grey70")
names(pal)<- c("H0 accepted", "H0 rejected", 
               "Undecided",  "Fail to\nreject H0")

ggplot()+
  geom_bar(data = summd_all %>% filter(test=="NHST only") %>% arrange(conclusion), 
           aes(x = 2, y = prop, fill = conclusion, group=test),alpha=0.8,stat = "identity")+
  geom_bar(data = summd_all %>% filter(test=="NHST+Equivalence")%>% arrange(conclusion), 
           aes(x = 2.8, y = prop, fill = conclusion, group=test),alpha=0.8,stat = "identity")+
  coord_polar("y",start=0) +
  geom_label(data= summd3, aes(x =4, y=0.4, label = round(avg_SGPV,2)), col = "grey20") +
  theme_void()+
  theme(text = element_text(size=18),
        panel.grid.major.y = element_line(colour="grey90")) +
  facet_grid(scenario~outcome, switch="y")+
  scale_fill_discrete(type= pal) +
  scale_x_continuous(limits = c(.2,4))+
  ylim(0,1)

