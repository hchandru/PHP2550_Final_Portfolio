knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "~/Downloads")

#Packages
install.packages("remotes", repos = "http://cran.us.r-project.org")
remotes::install_github("ddsjoberg/bstfun")

library(tidyverse)
library(naniar)
library(readr)
library(gtsummary)
library(kableExtra)
library(ggplot2)
library(stringr)
library(labelled)
library(bstfun)
library(forcats)
library(ggpubr)
library(corrplot)
library(gt)
#loading in the data
df <- read.csv("project1.csv")
#Data cleaning

lapply(df, class)

#One of the values of the mom_numcig variable is a range. We convert this to a mean of the low and high of the range. Another value of this column is "None". We change this to 0.

df[df$mom_numcig == "20-25", "mom_numcig"] <- mean(c(20,25))
df[df$mom_numcig == "None", "mom_numcig"] <- 0

#We will change the income variable and mom_numcig from character to numeric. 
char_cols <- c("income", "mom_numcig") 
df[char_cols] <- lapply(df[char_cols], parse_number)

#0's on SWAN inattentive and hyperactive are actually NA's
df$swan_hyperactive[df$swan_hyperactive == 0] <- NA
df$swan_inattentive[df$swan_inattentive == 0] <- NA

#Transforming the SWAN variables so that they indicate if the child is likely hyperactive/inattentive
df$swan_hyperactive[df$swan_hyperactive < 6] <- 0
df$swan_hyperactive[df$swan_hyperactive >= 6] <- 1
df$swan_inattentive[df$swan_inattentive < 6] <- 0
df$swan_inattentive[df$swan_inattentive >= 6] <- 1

#Transforming SDP and ETS variables so that they can be summarized easier
sdp_ets <- c("mom_smoke_16wk", "mom_smoke_22wk", "mom_smoke_32wk", "mom_smoke_pp1", "mom_smoke_pp2", "mom_smoke_pp12wk", "mom_smoke_pp6mo")

ind_func <- function(col){
  col_func <- df[,col]
  col_func <- case_when(col_func == "2=No" ~ 0,
                        col_func == "1=Yes" ~ 1)
  return(col_func)
}

df[sdp_ets] <- lapply(sdp_ets, ind_func)

#Turning any remaining blank values into NAs
df[df == ""] <- NA

#Change variables to factors and recode their levels
factors <- c("psex", "employ", "pedu", "tsex", "childasd")
df[,factors] <- lapply(factors, function(x){df[,x] <- as.factor(df[,x])})

df$psex <- fct_recode(df$psex, 
                      "Male" = "0", 
                      "Female" = "1")
df$employ <- fct_recode(df$employ, 
                        "No" = "0", 
                        "Part-Time" = "1", 
                        "Full-Time" = "2")
df$pedu <- fct_recode(df$pedu, 
                      "Some High School" = "0", 
                      "High School" = "1",
                      "GED" = "2",
                      "Some College" = "3",
                      "2 Year Degree" = "4",
                      "4 Year Degree" = "5",
                      "Postgraduate Degree" = "6")
df$tsex <- fct_recode(df$tsex, 
                      "Male" = "0", 
                      "Female" = "1")
df$childasd <- fct_recode(df$childasd,
                          "No" = "0",
                          "Diagnosed" = "1",
                          "Suspected" = "2")

#Convert "Prefer Not to Say" to NA
df$tethnic[df$tethnic == 2] <- NA

#Create variable labels
var_label(df) <- list(
  parent_id = "Parent ID",
  page = "Parent Age",
  psex = "Parent Sex",
  plang = "Another Language Spoken by Parent at Home",
  pethnic = "Parent Hispanic/Latino",
  paian = "Parent American Indian/Alaskan Native",
  pasian = "Parent Asian",
  pnhpi = "Parent Native Hawaiian or Pacific Islander",
  pblack = "Parent Black",
  pwhite = "Parent White",
  prace_other = "Parent Other Race",
  employ = "Parent Employed",
  pedu = "Highest Level of Education of Parent",
  income = "Estimated Annual Household Income",
  childasd = "ASD in Child",
  nidaalc = "Parent Alcohol Use in the Past 6 Months",
  nidatob = "Parent Tobacco Product Use in the Past 6 Months",
  nidapres = "Parent Prescription Drug Use for Non-Medical Reasons in the Past 6 Months",
  nidaill = "Parent Illegal Drug Use in the Past 6 Months",
  momcig = "Number of Days Parent Smoked in the Past 30 Days",
  mom_numcig = "Cigarettes Smoked per Day by Parent",
  mom_smoke_16wk = "Parent Smoker at 16 Weeks Pregnant",
  mom_smoke_22wk = "Parent Smoker at 22 Weeks Pregnant",
  mom_smoke_32wk = "Parent Smoker at 32 Weeks Pregnant",
  mom_smoke_pp1 = "Parent Smoker at First Postpartum Visit",
  mom_smoke_pp2 = "Parent Smoker at Second Postpartum Visit",
  mom_smoke_pp12wk = "Parent Smoker at 12 Weeks Postpartum", 
  mom_smoke_pp6mo = "Parent Smoker at 6 Months Postpartum", 
  cotimean_34wk = "Urine Cotinine at 34 Weeks Gestation",
  cotimean_pp6mo_baby = "Urine Cotinine at 6 Months Postpartum (Baby)", 
  cotimean_pp6mo = "Urine Cotinine at 6 Months Postpartum (Parent)",
  swan_inattentive = "Child Likely ADHD - Inattentive",
  swan_hyperactive = "Child Likely ADHD - Hyperactive/Impulsive",
  bpm_att_p = "Parent-Reported Attention Problems Score (BPM)",
  bpm_ext_p = "Parent-Reported Externalizing Problems Score (BPM)",
  bpm_int_p = "Parent-Reported Internalizing Problems Score (BPM)",
  smoke_exposure_6mo = "Smoke-Exposed Between 0 to 6 Months",
  smoke_exposure_12mo = "Smoke-Exposed Between 7 to 12 Months",
  smoke_exposure_2yr = "Smoke-Exposed in 2nd Year",
  smoke_exposure_3yr = "Smoke-Exposed in 3rd Year",
  smoke_exposure_4yr = "Smoke-Exposed in 4th Year",
  smoke_exposure_5yr = "Smoke-Exposed in 5th Year",
  ppmq_parental_knowledge = "Parent-Reported Average Response on Parental Knowledge (PKQ)",
  ppmq_child_disclosure = "Parent-Reported Average Response on Child Disclosure (PKQ)",
  ppmq_parental_solicitation = "Parent-Reported Average Response on Parental Solicitation (PKQ)",
  ppmq_parental_control = "Parent-Reported Average Response on Parental Control (PKQ)",
  bpm_att_a = "Self-Reported Attention Score for Parent (BPM)",
  bpm_ext_a = "Self-Reported Externalizing Problems Score for Parent (BPM)",
  bpm_int_a = "Self-Reported Internalizing Problems Score for Parent (BPM)",
  erq_cog_a = "Self-Reported Average Cognitive Reappraisal Score for Parent (ERQ)",
  erq_exp_a = "Self-Reported Average Expressive Suppression Score for Parent (ERQ)",
  tage = "Child Age",
  tsex = "Child Sex",
  language = "Another Language Spoken by Child at Home",
  tethnic = "Child Hispanic/Latino",
  taian = "Child American Indian/Alaskan Native",
  tasian = "Child Asian",
  tnhpi = "Child Native Hawaiian or Pacific Islander",
  tblack = "Child Black",
  twhite = "Child White",
  trace_other = "Child Other Race",
  cig_ever = "Cigarette Smoked Ever",
  num_cigs_30 = "Cigarette Use Over Past 30 Days (Child)",
  e_cig_ever = "E-Cigarette or Vape Used Ever",
  num_e_cigs_30 = "E-Cigarette or Vape Use Over Past 30 Days (Child)",
  mj_ever = "Marijuana Used Ever",
  num_mj_30 = "Marijuana Use Over Past 30 Days (Child)",
  alc_ever = "Alcohol Ever Consumed by Child (Not as Part of Religious Ceremonies)",
  num_alc_30 = "Alcohol Consumption Over Past 30 Days (Child)",
  bpm_att = "Self-Reported Attention Problems Score (BPM)",
  bpm_ext = "Self-Reported Externalizing Problems Score (BPM)",
  bpm_int = "Self-Reported Internalizing Problems Score (BPM)",
  erq_cog = "Average Cognitive Reappraisal Score of Child (ERQ)",
  erq_exp = "Average Expressive Suppression Score of Child (ERQ)",
  pmq_parental_knowledge = "Child-Reported Average Response on Parental Knowledge (PKQ)",
  pmq_child_disclosure = "Child-Reported Average Response on Child Disclosure (PKQ)",
  pmq_parental_solicitation = "Child-Reported Average Response on Parental Solicitation (PKQ)",
  pmq_parental_control = "Child-Reported Average Response on Parental Control (PKQ)"
)

#Remove variables not used in analysis
df <- df %>%
  select(-c(nidaalc, nidatob, nidapres, nidaill, momcig, mom_numcig, erq_cog_a, erq_exp_a, pmq_parental_knowledge, pmq_child_disclosure, pmq_parental_solicitation, pmq_parental_control, ppmq_parental_knowledge, ppmq_child_disclosure, ppmq_parental_solicitation, ppmq_parental_control, bpm_att_a, bpm_int_a, bpm_ext_a))

#Missing values
label_func <- function(var){
  return(attr(df[,var], "label"))
}

### Exploratory Data Analysis

#### A Summary of the Missing Data
df_miss <- df 
colnames(df_miss) <- lapply(colnames(df_miss), label_func)  

gg_miss_var(df_miss, show_pct = TRUE) +
  theme(axis.text.y = element_text(size = 5)) +
  ggtitle("Figure 1. Percent Missing Observations by Variable") +
  theme(plot.title = element_text(size = 10))
#### Demographics
var_label(df) <- list(
  pethnic = "Hispanic/Latino",
  paian = "American Indian/Alaskan Native",
  pasian = "Asian",
  pnhpi = "Native Hawaiian or Pacific Islander",
  pblack = "Black",
  pwhite = "White",
  prace_other = "Other Race",
  tethnic = "Hispanic/Latino",
  taian = "American Indian/Alaskan Native",
  tasian = "Asian",
  tnhpi = "Native Hawaiian or Pacific Islander",
  tblack = "Black",
  twhite = "White",
  trace_other = "Other Race"
)

#Tables summarizing parent demographics
pdems <- c("page", "psex", "plang", "employ", "pedu", "income")
prace_vars <- c("pethnic", "paian", "pasian", "pnhpi", "pblack", "pwhite", "prace_other")

theme_gtsummary_compact(set_theme=TRUE, font_size = 10)

df %>%
  select(all_of(prace_vars), all_of(pdems)) %>%
  tbl_summary(missing = "no") %>%
  add_variable_grouping("Parent Race" = prace_vars) %>%
  as_gt() %>%
  gt::tab_header(title = "Table 1. Parent Demographics") %>%
  gtsave(filename = "proj1_table1.png")

knitr::include_graphics('proj1_table1.png')
#Table summarizing child demographics
cdems <- c("tage", "tsex", "language")
crace_vars <- c("tethnic", "taian", "tasian", "tnhpi", "tblack", "twhite", "trace_other")

df %>%
  select(all_of(crace_vars), all_of(cdems)) %>%
  tbl_summary(missing = "no", type = list(tage ~ "continuous")) %>%
  add_variable_grouping("Child Race" = crace_vars) %>%
  as_gt() %>%
  gt::tab_header(title = "Table 2. Child Demographics") %>%
  gtsave(filename = "proj1_table2.png")

knitr::include_graphics('proj1_table2.png')
#### Substance Use in Child

su_ever_child <- c("cig_ever", "e_cig_ever", "mj_ever", "alc_ever")
su_dosage_child <- c("num_cigs_30", "num_e_cigs_30", "num_mj_30", "num_alc_30")

var_label(df) <- list(
  cig_ever = "Cigarette",
  e_cig_ever = "E-Cigarette or Vape",
  mj_ever = "Marijuana",
  alc_ever = "Alcohol (Not as Part of Religious Ceremonies)",
  num_cigs_30 = "Cigarette",
  num_e_cigs_30 = "E-Cigarette or Vape",
  num_mj_30 = "Marijuana",
  num_alc_30 = "Alcohol"
)

df %>%
  select(all_of(su_ever_child), all_of(su_dosage_child)) %>%
  tbl_summary(missing = "no",
              type = list(c("num_cigs_30", "num_e_cigs_30", "num_mj_30", "num_alc_30") ~ "continuous")) %>%
  add_variable_grouping("Substance Used (Ever)" = su_ever_child, "Number of Days Used in the Past 30 Days" = su_dosage_child) %>%
  as_gt() %>%
  gt::tab_header(title = "Table 3. Substance Use in Child") %>%
  gtsave(filename = "proj1_table3.png")

knitr::include_graphics('proj1_table3.png')
#### Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD in Child

sr_child <- c("erq_cog", "erq_exp")
ip_child <- c("bpm_int", "bpm_int_p")
ep_child <- c("bpm_ext", "bpm_ext_p")
ap_child <- c("bpm_att", "bpm_att_p")
adhd <- c("swan_hyperactive", "swan_inattentive")

var_label(df) <- list(
  erq_cog =  "Average Cognitive Reappraisal Score",
  erq_exp = "Average Expressive Suppression Score",
  bpm_int = "Self-Reported Score",
  bpm_int_p = "Parent-Reported Score",
  bpm_ext = "Self-Reported Score",
  bpm_ext_p = "Parent-Reported Score",
  bpm_att = "Self-Reported Score",
  bpm_att_p = "Parent-Reported Score",
  swan_hyperactive = "Likely ADHD - Hyperactive/Impulsive",
  swan_inattentive = "Likely ADHD - Inattentive"
)

df %>%
  select(all_of(sr_child), all_of(ip_child), all_of(ep_child), all_of(ap_child), childasd, all_of(adhd)) %>%
  tbl_summary(missing = "no",
              type = list(c(bpm_ext, bpm_ext_p, bpm_att, bpm_att_p) ~ "continuous")) %>%
  add_variable_grouping("Self-Regulation" = sr_child, 
                        "Internalizing Problems" = ip_child,
                        "Externalizing Problems" = ep_child,
                        "Attention Problems" = ap_child,
                        "ADHD in Child" = adhd) %>%
  as_gt() %>%
  gt::tab_header(title = "Table 4. Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD in Child") %>%
  gtsave(filename = "proj1_table4.png")

knitr::include_graphics('proj1_table4.png')
#Differing distributions of parent-reported and child-reported IP, EP and AP scores
scores_df <- df %>%
  select(all_of(ip_child), all_of(ep_child), all_of(ap_child)) %>%
  pivot_longer(cols = c(all_of(ip_child), all_of(ep_child), all_of(ap_child)), names_to = "score_type", values_to = "bpm_score") %>%
  mutate(reporter = case_when(score_type == "bpm_int_p" ~ "Parent",
                              score_type == "bpm_ext_p" ~ "Parent",
                              score_type == "bpm_att_p" ~ "Parent",
                              score_type == "bpm_int" ~ "Child (Self-Report)",
                              score_type == "bpm_ext" ~ "Child (Self-Report)",
                              score_type == "bpm_att" ~ "Child (Self-Report)"),
         score_type = case_when(score_type == "bpm_int_p" ~ "IP",
                                score_type == "bpm_ext_p" ~ "EP",
                                score_type == "bpm_att_p" ~ "AP",
                                score_type == "bpm_int" ~ "IP",
                                score_type == "bpm_ext" ~ "EP",
                                score_type == "bpm_att" ~ "AP"))

ggplot(data = scores_df, aes(x = reporter, y = bpm_score, fill = reporter)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = .1, fill = "black") +
  stat_summary(fun.y = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  ylab("Internalizing Problems Score") +
  labs(fill = "Reported By") +
  scale_fill_brewer() +
  facet_wrap(~score_type) +
  ggtitle("Figure 2. Difference in Distributions of Child-Reported and Parent-Reported Scores") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), title = element_text(size = 5))
#### Smoking During Pregnancy

sdp <- c("mom_smoke_16wk", "mom_smoke_22wk", "mom_smoke_32wk")

var_label(df) <- list(
  mom_smoke_16wk = "Smoker At 16 Weeks",
  mom_smoke_22wk = "Smoker At 22 Weeks",
  mom_smoke_32wk = "Smoker At 32 Weeks",
  cotimean_34wk = "UC At 34 Weeks Gestation"
)

df %>%
  select(all_of(sdp), cotimean_34wk) %>%
  tbl_summary(missing = "no") %>%
  as_gt() %>%
  gt::tab_header(title = "Table 5. Smoking During Pregnancy") %>%
  gtsave(filename = "proj1_table5.png")

knitr::include_graphics('proj1_table5.png')
#### Postpartum Smoking

pp_smoke <- c("mom_smoke_pp1", "mom_smoke_pp2", "mom_smoke_pp12wk", "mom_smoke_pp6mo")

var_label(df) <- list(
  mom_smoke_pp1 = "At First Postpartum Visit",
  mom_smoke_pp2 = "At Second Postpartum Visit",
  mom_smoke_pp12wk = "At 12 Weeks Postpartum", 
  mom_smoke_pp6mo = "At 6 Months Postpartum" 
)

df %>%
  select(all_of(pp_smoke)) %>%
  tbl_summary(missing = "no") %>%
  as_gt() %>%
  gt::tab_header(title = "Table 6. Postpartum Smoking") %>%
  gtsave(filename = "proj1_table6.png")

knitr::include_graphics('proj1_table6.png')
#### Environmental Tobacco Smoke

var_label(df) <- list(
  smoke_exposure_6mo = "0 to 6 Months",
  smoke_exposure_12mo = "7 to 12 Months", 
  smoke_exposure_2yr = "2nd Year",
  smoke_exposure_3yr = "3rd Year",
  smoke_exposure_4yr = "4th Year",
  smoke_exposure_5yr = "5th Year",
  cotimean_34wk = "At 34 Weeks Gestation",
  cotimean_pp6mo = "At 6 Months Postpartum (Parent)",
  cotimean_pp6mo_baby = "At 6 Months Postpartum (Baby)"
)

#Smoke exposure
ets_self <- c("smoke_exposure_6mo", "smoke_exposure_12mo", "smoke_exposure_2yr", "smoke_exposure_3yr", "smoke_exposure_4yr", "smoke_exposure_5yr")

#Urine cotinine
uc <- c("cotimean_pp6mo", "cotimean_pp6mo_baby")

df %>%
  select(all_of(ets_self), all_of(uc)) %>%
  tbl_summary(missing = "no") %>%
  add_variable_grouping("Smoke Exposure" = ets_self,
                        "UC levels" = uc) %>%
  as_gt() %>%
  gt::tab_header(title = "Table 7. Postpartum Environmental Tobacco Smoke Exposure") %>%
  gtsave(filename = "proj1_table7.png")

knitr::include_graphics('proj1_table7.png')
#### Association Between Timing of SDP and Outcomes

##### Substance Use by SDP Timing

#Outcomes - Substance use
df[,sdp] <- lapply(sdp, function(x){df[,x] <- as.factor(df[,x])})

df$mom_smoke_16wk <- fct_recode(df$mom_smoke_16wk,
                                "No" = "0",
                                "Yes" =  "1")
df$mom_smoke_22wk <- fct_recode(df$mom_smoke_22wk,
                                "No" = "0",
                                "Yes" =  "1")
df$mom_smoke_32wk <- fct_recode(df$mom_smoke_32wk,
                                "No" = "0",
                                "Yes" =  "1")

var_label(df) <- list(
  alc_ever = "Alcohol"
)

#Ever Used
df %>%
  select(all_of(sdp), all_of(su_ever_child)) %>%
  pivot_longer(cols = all_of(sdp), names_to = "sdp_timing", values_to = "exposed") %>%
  mutate(sdp_timing = case_when(sdp_timing == "mom_smoke_16wk" ~ "At 16 Weeks",
                                sdp_timing == "mom_smoke_22wk" ~ "At 22 Weeks",
                                sdp_timing == "mom_smoke_32wk" ~ "At 32 Weeks")) %>%
  tbl_strata(strata = sdp_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no")) %>%
  as_gt() %>%
  gt::tab_header(title = "Table 8. Substance Use (Ever) by Timing of SDP") %>%
  gtsave(filename = "proj1_table8.png")

knitr::include_graphics('proj1_table8.png')
#Frequency of Usage
df %>%
  select(all_of(sdp), all_of(su_dosage_child)) %>%
  pivot_longer(cols = all_of(sdp), names_to = "sdp_timing", values_to = "exposed") %>%
  mutate(sdp_timing = case_when(sdp_timing == "mom_smoke_16wk" ~ "At 16 Weeks",
                                sdp_timing == "mom_smoke_22wk" ~ "At 22 Weeks",
                                sdp_timing == "mom_smoke_32wk" ~ "At 32 Weeks")) %>%
  tbl_strata(strata = sdp_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no",
                           type = list(where(is.numeric) ~ "continuous"))) %>%
  as_gt() %>%
  gt::tab_header(title = "Table 9. Frequency of Substance Use by Timing of SDP") %>%
  gtsave(filename = "proj1_table9.png")

knitr::include_graphics('proj1_table9.png')
##### Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD by SDP Timing

var_label(df) <- list(
  erq_cog =  "Average CR Score",
  erq_exp = "Average ES Score",
  bpm_int = "Self-Reported Score",
  bpm_int_p = "Parent-Reported Score",
  bpm_ext = "Self-Reported Score",
  bpm_ext_p = "Parent-Reported Score",
  bpm_att = "Self-Reported Score",
  bpm_att_p = "Parent-Reported Score",
  swan_hyperactive = "Likely ADHD-H",
  swan_inattentive = "Likely ADHD-I"
)

df %>%
  select(all_of(sdp), all_of(sr_child), all_of(ip_child), all_of(ep_child), all_of(ap_child), childasd, all_of(adhd)) %>%
  pivot_longer(cols = all_of(sdp), names_to = "sdp_timing", values_to = "exposed") %>%
  mutate(sdp_timing = case_when(sdp_timing == "mom_smoke_16wk" ~ "At 16 Weeks",
                                sdp_timing == "mom_smoke_22wk" ~ "At 22 Weeks",
                                sdp_timing == "mom_smoke_32wk" ~ "At 32 Weeks")) %>%
  tbl_strata(strata = sdp_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no",
                           type = list(c(bpm_ext, bpm_ext_p, bpm_att, bpm_att_p, bpm_int, bpm_int_p) ~ "continuous"))) %>%
  add_variable_grouping("Self-Regulation" = sr_child, 
                        "Internalizing Problems" = ip_child,
                        "Externalizing Problems" = ep_child,
                        "Attention Problems" = ap_child,
                        "ADHD in Child" = adhd) %>%
  as_gt() %>%
  gt::tab_header("Table 10. Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD in Child by Timing of SDP") %>%
  gtsave(filename = "proj1_table10.png")

knitr::include_graphics('proj1_table10.png')
#### Associations Between Timing of Postpartum Smoking and Outcomes

##### Substance Use by Postpartum Smoking Timing
var_label(df) <- list(
  mom_smoke_pp1 = "At First PP Visit",
  mom_smoke_pp2 = "At Second PP Visit",
  mom_smoke_pp12wk = "At 12 Weeks PP", 
  mom_smoke_pp6mo = "At 6 Months PP" 
)

df[,pp_smoke] <- lapply(pp_smoke, function(x){df[,x] <- as.factor(df[,x])})

df$mom_smoke_pp1 <- fct_recode(df$mom_smoke_pp1,
                               "No" = "0",
                               "Yes" =  "1")
df$mom_smoke_pp2 <- fct_recode(df$mom_smoke_pp2,
                               "No" = "0",
                               "Yes" =  "1")
df$mom_smoke_pp12wk <- fct_recode(df$mom_smoke_pp12wk,
                                  "No" = "0",
                                  "Yes" =  "1")
df$mom_smoke_pp6mo <- fct_recode(df$mom_smoke_pp6mo,
                                 "No" = "0",
                                 "Yes" =  "1")

#Ever Used
df %>%
  select(all_of(pp_smoke), all_of(su_ever_child)) %>%
  pivot_longer(cols = all_of(pp_smoke), names_to = "pp_smoke_timing", values_to = "exposed") %>%
  mutate(pp_smoke_timing = case_when(pp_smoke_timing == "mom_smoke_pp1" ~ "1st PP Visit",
                                     pp_smoke_timing == "mom_smoke_pp2" ~ "2nd PP Visit",
                                     pp_smoke_timing == "mom_smoke_pp12wk" ~ "At 12 Weeks PP",
                                     pp_smoke_timing == "mom_smoke_pp6mo" ~ "At 6 Months PP")) %>%
  tbl_strata(strata = pp_smoke_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no")) %>%
  as_gt() %>%
  gt::tab_header("Table 11. Substance Use (Ever) by Timing of Postpartum Smoking") %>%
  gtsave(filename = "proj1_table11.png")

knitr::include_graphics('proj1_table11.png')


#Frequency of Usage
df %>%
  select(all_of(pp_smoke), all_of(su_dosage_child)) %>%
  pivot_longer(cols = all_of(pp_smoke), names_to = "pp_smoke_timing", values_to = "exposed") %>%
  mutate(pp_smoke_timing = case_when(pp_smoke_timing == "mom_smoke_pp1" ~ "1st PP Visit",
                                     pp_smoke_timing == "mom_smoke_pp2" ~ "2nd PP Visit",
                                     pp_smoke_timing == "mom_smoke_pp12wk" ~ "At 12 Weeks PP",
                                     pp_smoke_timing == "mom_smoke_pp6mo" ~ "At 6 Months PP")) %>%
  tbl_strata(strata = pp_smoke_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no",
                           type = list(where(is.numeric) ~ "continuous"))) %>%
  as_gt() %>%
  gt::tab_header("Table 12. Frequency of Substance Use by Timing of Postpartum Smoking") %>%
  gtsave(filename = "proj1_table12.png")

knitr::include_graphics('proj1_table12.png')
##### Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD

var_label(df) <- list(
  erq_cog =  "Average CR Score",
  erq_exp = "Average ES Score",
  bpm_int = "Self-Reported Score",
  bpm_int_p = "Parent-Reported Score",
  bpm_ext = "Self-Reported Score",
  bpm_ext_p = "Parent-Reported Score",
  bpm_att = "Self-Reported Score",
  bpm_att_p = "Parent-Reported Score",
  swan_hyperactive = "Likely ADHD-H",
  swan_inattentive = "Likely ADHD-I"
)

df %>%
  select(all_of(pp_smoke), all_of(sr_child), all_of(ip_child), all_of(ep_child), all_of(ap_child), childasd, all_of(adhd)) %>%
  pivot_longer(cols = all_of(pp_smoke), names_to = "pp_smoke_timing", values_to = "exposed") %>%
  mutate(pp_smoke_timing = case_when(pp_smoke_timing == "mom_smoke_pp1" ~ "1st PP Visit",
                                     pp_smoke_timing == "mom_smoke_pp2" ~ "2nd PP Visit",
                                     pp_smoke_timing == "mom_smoke_pp12wk" ~ "At 12 Weeks PP",
                                     pp_smoke_timing == "mom_smoke_pp6mo" ~ "At 6 Months PP")) %>%
  tbl_strata(strata = pp_smoke_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no",
                           type = list(c(erq_cog, erq_exp, bpm_ext, bpm_ext_p, bpm_att, bpm_att_p, bpm_int, bpm_int_p) ~ "continuous"))) %>%
  add_variable_grouping("Self-Regulation" = sr_child, 
                        "Internalizing Problems" = ip_child,
                        "Externalizing Problems" = ep_child,
                        "Attention Problems" = ap_child,
                        "ADHD in Child" = adhd) %>%
  as_gt() %>%
  gt::tab_header("Table 13. Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD in Child by Timing of Postpartum Smoking") %>%
  gtsave(filename = "proj1_table13.png")

knitr::include_graphics('proj1_table13.png')
#### Association Between Timing of ETS and Outcomes

##### Substance Use in Child by ETS Timing

var_label(df) <- list(
  smoke_exposure_6mo = "0 to 6 Months",
  smoke_exposure_12mo = "7 to 12 Months", 
  smoke_exposure_2yr = "2nd Year",
  smoke_exposure_3yr = "3rd Year",
  smoke_exposure_4yr = "4th Year",
  smoke_exposure_5yr = "5th Year",
  cotimean_34wk = "At 34 Weeks Gestation",
  cotimean_pp6mo = "At 6 Months Postpartum (Parent)",
  cotimean_pp6mo_baby = "At 6 Months Postpartum (Baby)"
)

df[,ets_self] <- lapply(ets_self, function(x){df[,x] <- as.factor(df[,x])})

df$smoke_exposure_6mo <- fct_recode(df$smoke_exposure_6mo,
                                    "No" = "0",
                                    "Yes" =  "1")
df$smoke_exposure_12mo <- fct_recode(df$smoke_exposure_12mo,
                                     "No" = "0",
                                     "Yes" =  "1")
df$smoke_exposure_2yr <- fct_recode(df$smoke_exposure_2yr,
                                    "No" = "0",
                                    "Yes" =  "1")
df$smoke_exposure_3yr <- fct_recode(df$smoke_exposure_3yr,
                                    "No" = "0",
                                    "Yes" =  "1")
df$smoke_exposure_4yr <- fct_recode(df$smoke_exposure_4yr,
                                    "No" = "0",
                                    "Yes" =  "1")
df$smoke_exposure_5yr <- fct_recode(df$smoke_exposure_5yr,
                                    "No" = "0",
                                    "Yes" =  "1")

#Ever Used
df %>%
  select(all_of(ets_self), all_of(su_ever_child)) %>%
  pivot_longer(cols = all_of(ets_self), names_to = "ets_timing", values_to = "exposed") %>%
  mutate(ets_timing = case_when(ets_timing == "smoke_exposure_6mo" ~ "Months 0 to 6",
                                ets_timing == "smoke_exposure_12mo" ~ "Months 7 to 12",
                                ets_timing == "smoke_exposure_2yr" ~ "Year 2",
                                ets_timing == "smoke_exposure_3yr" ~ "Year 3",
                                ets_timing == "smoke_exposure_4yr" ~ "Year 4",
                                ets_timing == "smoke_exposure_5yr" ~ "Year 5")) %>%
  tbl_strata(strata = ets_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no")) %>%
  modify_caption("**Substance Use (Ever) by Timing of ETS**")

#Frequency of Usage
df %>%
  select(all_of(ets_self), all_of(su_dosage_child)) %>%
  pivot_longer(cols = all_of(ets_self), names_to = "ets_timing", values_to = "exposed") %>%
  mutate(ets_timing = case_when(ets_timing == "smoke_exposure_6mo" ~ "Months 0 to 6",
                                ets_timing == "smoke_exposure_12mo" ~ "Months 7 to 12",
                                ets_timing == "smoke_exposure_2yr" ~ "Year 2",
                                ets_timing == "smoke_exposure_3yr" ~ "Year 3",
                                ets_timing == "smoke_exposure_4yr" ~ "Year 4",
                                ets_timing == "smoke_exposure_5yr" ~ "Year 5")) %>%
  tbl_strata(strata = ets_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no",
                           type = list(where(is.numeric) ~ "continuous"))) %>%
  modify_caption("**Frequency of Substance Use by Timing of ETS**") 
##### Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD in Child by ETS Timing

var_label(df) <- list(
  erq_cog =  "Average Cognitive Reappraisal Score",
  erq_exp = "Average Expressive Suppression Score",
  bpm_int = "Self-Reported Score",
  bpm_int_p = "Parent-Reported Score",
  bpm_ext = "Self-Reported Score",
  bpm_ext_p = "Parent-Reported Score",
  bpm_att = "Self-Reported Score",
  bpm_att_p = "Parent-Reported Score",
  swan_hyperactive = "Likely ADHD - Hyperactive/Impulsive",
  swan_inattentive = "Likely ADHD - Inattentive"
)

df %>%
  select(all_of(ets_self), all_of(sr_child), all_of(ip_child), all_of(ep_child), all_of(ap_child), childasd, all_of(adhd)) %>%
  pivot_longer(cols = all_of(ets_self), names_to = "ets_timing", values_to = "exposed") %>%
  mutate(ets_timing = case_when(ets_timing == "smoke_exposure_6mo" ~ "Months 0 to 6",
                                ets_timing == "smoke_exposure_12mo" ~ "Months 7 to 12",
                                ets_timing == "smoke_exposure_2yr" ~ "Year 2",
                                ets_timing == "smoke_exposure_3yr" ~ "Year 3",
                                ets_timing == "smoke_exposure_4yr" ~ "Year 4",
                                ets_timing == "smoke_exposure_5yr" ~ "Year 5")) %>%
  tbl_strata(strata = ets_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no",
                           type = list(c(erq_cog, erq_exp, bpm_ext, bpm_ext_p, bpm_att, bpm_att_p, bpm_int, bpm_int_p) ~ "continuous"))) %>%
  add_variable_grouping("Self-Regulation" = sr_child, 
                        "Internalizing Problems" = ip_child,
                        "Externalizing Problems" = ep_child,
                        "Attention Problems" = ap_child,
                        "ADHD in Child" = adhd) %>%
  modify_caption("**Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD in Child by Timing of ETS**")
#### Association Between SDP Dosage and Outcomes

##### SDP Dosage and Substance Use
uc_df_prenatal <- df %>%
  select(cotimean_34wk, all_of(su_ever_child), all_of(su_dosage_child), all_of(sr_child), all_of(ip_child), all_of(ep_child), all_of(ap_child), childasd, all_of(adhd))

#coverting to categorical
uc_df_prenatal$cotimean_34wk <- cut(uc_df_prenatal$cotimean_34wk, 
                                    breaks = c(0, 11, 30, 500, Inf), 
                                    labels = c("Less than 11", "11 to 30", "31 to 500", "Above 500"))

#Ever Used
uc_df_prenatal %>%
  select(cotimean_34wk, all_of(su_ever_child)) %>%
  tbl_summary(by = cotimean_34wk, missing = "no") %>%
  as_gt() %>%
  gt::tab_header("Table 14. Substance Use (Ever) by SDP Dosage") %>%
  gtsave(filename = "proj1_table14.png")

knitr::include_graphics('proj1_table14.png')

#Frequency of Usage
uc_df_prenatal %>%
  select(cotimean_34wk, all_of(su_dosage_child)) %>%
  tbl_summary(by = cotimean_34wk, missing = "no", type = list(where(is.numeric) ~ "continuous")) %>%
  as_gt() %>%
  gt::tab_header("Table 15. Frequency of Substance Use by SDP Dosage") %>%
  gtsave(filename = "proj1_table15.png")

knitr::include_graphics('proj1_table15.png')
##### Association between SDP Dosage and Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD
var_label(uc_df_prenatal) <- list(
  erq_cog =  "Average Cognitive Reappraisal Score",
  erq_exp = "Average Expressive Suppression Score",
  bpm_int = "Self-Reported Score",
  bpm_int_p = "Parent-Reported Score",
  bpm_ext = "Self-Reported Score",
  bpm_ext_p = "Parent-Reported Score",
  bpm_att = "Self-Reported Score",
  bpm_att_p = "Parent-Reported Score",
  swan_hyperactive = "Likely ADHD - Hyperactive/Impulsive",
  swan_inattentive = "Likely ADHD - Inattentive"
)

uc_df_prenatal %>%
  select(cotimean_34wk, all_of(sr_child), all_of(ip_child), all_of(ep_child), all_of(ap_child), childasd, all_of(adhd)) %>%
  tbl_summary(by = cotimean_34wk, missing = "no",
              type = list(c(erq_cog, erq_exp, bpm_ext, bpm_ext_p, bpm_att, bpm_att_p, bpm_int, bpm_int_p) ~ "continuous")) %>%
  add_variable_grouping("Self-Regulation" = sr_child, 
                        "Internalizing Problems" = ip_child,
                        "Externalizing Problems" = ep_child,
                        "Attention Problems" = ap_child,
                        "ADHD in Child" = adhd) %>%
  as_gt() %>%
  gt::tab_header("Table 16. Self-Regulation, Internalizing Problems, Externalizing Problems, Attention Problems, ASD and ADHD in Child by SDP Dosage") %>%
  gtsave(filename = "proj1_table16.png")

knitr::include_graphics('proj1_table16.png')
#### Interrelatedness of Prenatal and Postnatal Exposure

##### Self-Report Variables
#Exposure self-report
var_label(df) <- list(
  mom_smoke_pp1 = "Parent Smoker at First Postpartum Visit",
  mom_smoke_pp2 = "Parent Smoker at Second Postpartum Visit",
  mom_smoke_pp12wk = "Parent Smoker at 12 Weeks Postpartum", 
  mom_smoke_pp6mo = "Parent Smoker at 6 Months Postpartum", 
  smoke_exposure_6mo = "Smoke-Exposed Between 0 to 6 Months",
  smoke_exposure_12mo = "Smoke-Exposed Between 7 to 12 Months",
  smoke_exposure_2yr = "Smoke-Exposed in 2nd Year",
  smoke_exposure_3yr = "Smoke-Exposed in 3rd Year",
  smoke_exposure_4yr = "Smoke-Exposed in 4th Year",
  smoke_exposure_5yr = "Smoke-Exposed in 5th Year"
)

df %>%
  select(all_of(sdp), all_of(pp_smoke), all_of(ets_self)) %>%
  pivot_longer(cols = all_of(sdp), names_to = "sdp_timing", values_to = "exposed") %>%
  mutate(sdp_timing = case_when(sdp_timing == "mom_smoke_16wk" ~ "At 16 Weeks",
                                sdp_timing == "mom_smoke_22wk" ~ "At 22 Weeks",
                                sdp_timing == "mom_smoke_32wk" ~ "At 32 Weeks")) %>%
  tbl_strata(strata = sdp_timing, 
             .tbl_fun = 
               ~ .x %>%
               tbl_summary(by = exposed, missing = "no")) %>%
  as_gt() %>%
  gt::tab_header("Table 17. Postnatal Smoke Exposure by Exposure to SDP") %>%
  gtsave(filename = "proj1_table17.png")

knitr::include_graphics('proj1_table17.png')
##### Urine Cotinine
uc_df_prepost <- df %>%
  select(cotimean_34wk, cotimean_pp6mo_baby)

ggplot(data = uc_df_prepost, aes(x = cotimean_34wk, y = cotimean_pp6mo_baby)) +
  geom_point() +
  stat_cor(method = "pearson", label.x = -5, label.y = 40) 
#### Interrelations Among Self-Regulation Data
df1 <- read.csv("project1.csv")
df$swan_hyperactive <- df1$swan_hyperactive
df$swan_inattentive <- df1$swan_inattentive

scores_data <- df %>%
  select(all_of(sr_child), all_of(ip_child), all_of(ep_child), all_of(ap_child), all_of(adhd)) 

scores_data <- scores_data[complete.cases(scores_data), ]

var_label(df) <- list(
  erq_cog = "Average CR Score",
  erq_exp = "Average ES Score",
  bpm_att_p = "Parent-Reported AP Score",
  bpm_ext_p = "Parent-Reported EP Score",
  bpm_int_p = "Parent-Reported IP Score",
  bpm_att = "Self-Reported AP Score",
  bpm_ext = "Self-Reported EP Score",
  bpm_int = "Self-Reported IP Score",
  swan_inattentive = "SWAN Score - Inattentive",
  swan_hyperactive = "SWAN Score - Hyperactive"
)
colnames(scores_data) <- lapply(colnames(scores_data), label_func)

correlation_matrix <- cor(scores_data)
corrplot(correlation_matrix, method = "number", number.cex = 0.5, title = "Figure 3. Correlations Between Self-Regulation Variables", tl.cex = 0.5, mar=c(0,0,2,0))