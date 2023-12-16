knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "~/Downloads")

#Packages
install.packages("riskCommunicator", repos = "http://cran.us.r-project.org")
install.packages("nhanesA", repos = "http://cran.us.r-project.org")

library(riskCommunicator)
library(tidyverse)
library(tableone)
library(nhanesA)
library(naniar)
library(mice)
library(survey)
library(corrplot)
library(rsimsum)
library(knitr)
library(MASS)
data("framingham")

# The Framingham data has been used to create models for cardiovascular risk.
# The variable selection and model below are designed to mimic the models used
# in the paper General Cardiovascular Risk Profile for Use in Primary Care 
# This paper is available (cvd_risk_profile.pdf) on Canvas.

framingham_df <- framingham %>% dplyr::select(c(CVD, TIMECVD, SEX, TOTCHOL, AGE,
                                                SYSBP, DIABP, CURSMOKE, DIABETES, BPMEDS,
                                                HDLC, BMI))
framingham_df <- na.omit(framingham_df)

framingham_df$CURSMOKE <- as.factor(framingham_df$CURSMOKE)
framingham_df$DIABETES <- as.factor(framingham_df$DIABETES)
framingham_df$BPMEDS <- as.factor(framingham_df$BPMEDS)
framingham_df$CVD <- as.factor(framingham_df$CVD)

framingham_dem <- print(CreateTableOne(data=framingham_df, strata = c("SEX"))$ContTable)

# Get blood pressure based on whether or not on BPMEDS
framingham_df$SYSBP_UT <- ifelse(framingham_df$BPMEDS == 0, 
                                 framingham_df$SYSBP, 0)
framingham_df$SYSBP_T <- ifelse(framingham_df$BPMEDS == 1, 
                                framingham_df$SYSBP, 0)

# Looking at risk within 15 years - remove censored data
# dim(framingham_df)
framingham_df <- framingham_df %>%
  filter(!(CVD == 0 & TIMECVD <= 365*15)) %>%
  dplyr::select(-c(TIMECVD)) %>%
  mutate(S = 1)
# dim(framingham_df)

# Filter to each sex
framingham_df_men <- framingham_df %>% filter(SEX == 1) 
framingham_df_women <- framingham_df %>% filter(SEX == 2) 

# The NHANES data here finds the same covariates among this national survey data

# blood pressure, demographic, bmi, smoking, and hypertension info
bpx_2017 <- nhanes("BPX_J") %>% 
  dplyr::select(SEQN, BPXSY1 ) %>% 
  rename(SYSBP = BPXSY1)
demo_2017 <- nhanes("DEMO_J") %>% 
  dplyr::select(SEQN, RIAGENDR, RIDAGEYR) %>% 
  rename(SEX = RIAGENDR, AGE = RIDAGEYR)
bmx_2017 <- nhanes("BMX_J") %>% 
  dplyr::select(SEQN, BMXBMI) %>% 
  rename(BMI = BMXBMI)
smq_2017 <- nhanes("SMQ_J") %>%
  mutate(CURSMOKE = case_when(SMQ040 %in% c(1,2) ~ 1,
                              SMQ040 == 3 ~ 0, 
                              SMQ020 == 2 ~ 0)) %>%
  dplyr::select(SEQN, CURSMOKE)
bpq_2017 <- nhanes("BPQ_J") %>% 
  mutate(BPMEDS = case_when(
    BPQ020 == 2 ~ 0,
    BPQ040A == 2 ~ 0,
    BPQ050A == 1 ~ 1,
    TRUE ~ NA )) %>%
  dplyr::select(SEQN, BPMEDS) 
tchol_2017 <- nhanes("TCHOL_J") %>% 
  dplyr::select(SEQN, LBXTC) %>% 
  rename(TOTCHOL = LBXTC)
hdl_2017 <- nhanes("HDL_J") %>% 
  dplyr::select(SEQN, LBDHDD) %>% 
  rename(HDLC = LBDHDD)
diq_2017 <- nhanes("DIQ_J") %>% 
  mutate(DIABETES = case_when(DIQ010 == 1 ~ 1, 
                              DIQ010 %in% c(2,3) ~ 0, 
                              TRUE ~ NA)) %>%
  dplyr::select(SEQN, DIABETES) 

# Join data from different tables
df_2017 <- bpx_2017 %>%
  full_join(demo_2017, by = "SEQN") %>%
  full_join(bmx_2017, by = "SEQN") %>%
  full_join(hdl_2017, by = "SEQN") %>%
  full_join(smq_2017, by = "SEQN") %>%
  full_join(bpq_2017, by = "SEQN") %>%
  full_join(tchol_2017, by = "SEQN") %>%
  full_join(diq_2017, by = "SEQN") %>%
  filter(AGE < 63 & AGE > 29) #subset of data that meets the eligibility criteria for the Framingham study

df_2017$CURSMOKE <- as.factor(df_2017$CURSMOKE)
df_2017$DIABETES <- as.factor(df_2017$DIABETES)
df_2017$BPMEDS <- as.factor(df_2017$BPMEDS)

nhanes_dem <- print(CreateTableOne(data = df_2017[-1], strata = c("SEX"))$ContTable)

df_2017$SYSBP_UT <- ifelse(df_2017$BPMEDS == 0, 
                           df_2017$SYSBP, 0)
df_2017$SYSBP_T <- ifelse(df_2017$BPMEDS == 1, 
                          df_2017$SYSBP, 0)
framingham_dem %>%
  kable(caption = "Summary Statistics for Framingham Data")

nhanes_dem %>%
  kable(caption = "Summary Statistics for NHANES Data")
miss_var_summary(df_2017) %>%
  kable(caption = "Summary of Missing Data in NHANES")
set.seed(1234)

#data sampled from the source population
#to create test and train data sets
ignore_source_men <- sample(c(TRUE, FALSE), nrow(framingham_df_men), replace = TRUE, prob = c(0.3,0.7))
ignore_source_women <- sample(c(TRUE, FALSE), nrow(framingham_df_women), replace = TRUE, prob = c(0.3,0.7))

framingham_df_men_train <- framingham_df_men[!ignore_source_men,]
framingham_df_men_test <- framingham_df_men[ignore_source_men,]

framingham_df_women_train <- framingham_df_women[!ignore_source_women,]
framingham_df_women_test <- framingham_df_women[ignore_source_women,]

# Fit risk score models with log transforms for all continuous variables on training set from source population data
mod_men <- glm(CVD~log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                 log(SYSBP_T+1)+CURSMOKE+DIABETES,
               data= framingham_df_men_train, family= "binomial")

mod_women <- glm(CVD~log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                   log(SYSBP_T+1)+CURSMOKE+DIABETES,
                 data= framingham_df_women_train, family= "binomial")

#Brier score on test set from source population data
#men
framingham_df_men_test$pred_probs <- predict(mod_men, newdata = framingham_df_men_test, type = "response")
brier_men_source_test <- (sum((framingham_df_men_test$pred_probs - (as.numeric(framingham_df_men_test$CVD)-1))^2))/(nrow(framingham_df_men_test))

#women
framingham_df_women_test$pred_probs <- predict(mod_women, newdata = framingham_df_women_test, type = "response")
brier_women_source_test <- (sum((framingham_df_women_test$pred_probs - (as.numeric(framingham_df_women_test$CVD)-1))^2))/(nrow(framingham_df_women_test))

#table of brier scores
brier_df_source <- data.frame(
  Model = c("Men", "Women"),
  Brier_Score = c(brier_men_source_test, brier_women_source_test)
)

kable(brier_df_source, col.names = c("Model", "Brier Score"), caption = "Brier Score Estimates in the Test Set of the Source Population")
set.seed(1234)

#data sampled from the target population
complete_df_2017 <- df_2017[complete.cases(df_2017),] %>%
  mutate(S = 0)

#filtering by sex
complete_df_2017_men <- complete_df_2017 %>%
  filter(SEX == 1)
complete_df_2017_women <- complete_df_2017 %>%
  filter(SEX == 2)

#combining the testing data from both the source and target populations to get the inverse odds weights
men_test_full_df <- rbind(framingham_df_men_test[-c(1,6,15)], complete_df_2017_men[-1])
women_test_full_df <- rbind(framingham_df_women_test[-c(1,6,15)], complete_df_2017_women[-1])

#models for odds of belonging to the source population

#men
men_source_lo_test <- glm(S ~ log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                            log(SYSBP_T+1)+CURSMOKE+DIABETES,
                          data = men_test_full_df, family = "binomial")
men_test_full_df$odds <- predict(men_source_lo_test, newdata = men_test_full_df, type = "response")/(1- predict(men_source_lo_test, newdata = men_test_full_df, type = "response"))#obtaining odds
men_test_full_df$weights <- 1/men_test_full_df$odds #inverse odds weights
men_test_source <- men_test_full_df %>% filter(S == 1) #only source population data to fit model
men_test_source <- cbind(framingham_df_men_test$CVD, men_test_source) %>%
  rename("CVD" = "framingham_df_men_test$CVD")
men_test_source$CVD <- as.numeric(men_test_source$CVD)-1

#brier score
men_test_source$pred_probs <- predict(mod_men, newdata = men_test_source, type = "response") 
men_test_source$brier_num <- (men_test_source$weights)*((men_test_source$CVD - men_test_source$pred_probs)^2)

men_brier_score <- sum(men_test_source$brier_num)/nrow(complete_df_2017_men) 

#women
women_source_lo_test <- glm(S ~ log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                              log(SYSBP_T+1)+CURSMOKE+DIABETES,
                            data = women_test_full_df, family = "binomial")
women_test_full_df$odds <- (predict(women_source_lo_test, newdata = women_test_full_df, type = "response"))/(1 - predict(women_source_lo_test, newdata = women_test_full_df, type = "response")) #obtaining odds
women_test_full_df$weights <- 1/women_test_full_df$odds
women_test_source <- women_test_full_df %>% filter(S == 1)
women_test_source <- cbind(framingham_df_women_test$CVD, women_test_source) %>%
  rename("CVD" = "framingham_df_women_test$CVD")
women_test_source$CVD <- as.numeric(women_test_source$CVD)-1

#brier score
women_test_source$pred_probs <- predict(mod_women, newdata = women_test_source, type = "response") 
women_test_source$brier_num <- (women_test_source$weights)*((women_test_source$CVD - women_test_source$pred_probs)^2)

women_brier_score <- sum(women_test_source$brier_num)/nrow(complete_df_2017_women)

#table of brier scores
brier_df <- data.frame(
  Model = c("Men", "Women"),
  Brier_Score = c(men_brier_score, women_brier_score)
)

kable(brier_df, col.names = c("Model", "Brier Score"), caption = "Weighted Brier Score Estimates in the Target Population")
set.seed(1234)
sim_brier_function_m <- function(sample_size){
  sim_target_male_df <- data.frame(
    SYSBP = rnorm(n = sample_size, mean = 126.07, sd = 16.63),
    AGE = rnorm(n = sample_size, mean = 47.16, sd = 9.97),
    HDLC = rnorm(n = sample_size, mean = 47.45, sd = 14.54),
    CURSMOKE = rbinom(n = sample_size, size = 1, prob = 0.259),
    BPMEDS = rbinom(n = sample_size, size = 1, prob = 0.755),
    TOTCHOL = rnorm(n = sample_size, mean = 192.86, sd = 40.71),
    DIABETES = rbinom(n = sample_size, size = 1, prob = 0.131),
    S = 0
  )
  
  sim_target_male_df$CURSMOKE <- as.factor(sim_target_male_df$CURSMOKE)
  sim_target_male_df$DIABETES <- as.factor(sim_target_male_df$DIABETES)
  sim_target_male_df$BPMEDS <- as.factor(sim_target_male_df$BPMEDS)
  
  # Get blood pressure based on whether or not on BPMEDS
  sim_target_male_df$SYSBP_UT <- ifelse(sim_target_male_df$BPMEDS == 0, 
                                        sim_target_male_df$SYSBP, 0)
  sim_target_male_df$SYSBP_T <- ifelse(sim_target_male_df$BPMEDS == 1, 
                                       sim_target_male_df$SYSBP, 0)
  
  #combining the testing data from both the source and target populations to get the inverse odds weights
  men_test_full_df_sim <- rbind(framingham_df_men_test[,-c(1,2,6,11,15)], sim_target_male_df)
  
  #models for odds of belonging to the source population
  
  #men
  men_source_lo_test_sim <- glm(S ~ log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                                  log(SYSBP_T+1)+CURSMOKE+DIABETES,
                                data = men_test_full_df_sim, family = "binomial")
  men_test_full_df_sim$odds <- (predict(men_source_lo_test_sim, newdata = men_test_full_df_sim, type = "response"))/(1-(predict(men_source_lo_test_sim, newdata = men_test_full_df_sim, type = "response"))) #obtaining odds
  men_test_full_df_sim$weights <- 1/men_test_full_df_sim$odds #inverse odds weights
  men_test_source_sim <- men_test_full_df_sim %>% filter(S == 1) #only source population data to fit model
  men_test_source_sim <- cbind(framingham_df_men_test$CVD, men_test_source_sim) %>%
    rename("CVD" = "framingham_df_men_test$CVD")
  men_test_source_sim$CVD <- as.numeric(men_test_source_sim$CVD)-1
  
  #brier score
  men_test_source_sim$pred_probs <- predict(mod_men, newdata = men_test_source_sim, type = "response") 
  men_test_source_sim$brier_num <- (men_test_source_sim$weights)*((men_test_source_sim$CVD - men_test_source_sim$pred_probs)^2)
  
  men_brier_score_sim <- sum(men_test_source_sim$brier_num)/nrow(sim_target_male_df)
  
  return(men_brier_score_sim)
}

men_brier_sim_res <- replicate(1600, sim_brier_function_m(1500))
output_men_df <- data.frame(dataset = 1:1600, bs = men_brier_sim_res)
simsum_men <- simsum(data = output_men_df, estvarname = "bs", true = men_brier_score)
rsimsum::kable(simsum_men, stats = c("thetamean", "thetamedian", "bias", "empse", "mse"), caption = "Performance Measures for the Simulations Estimating Brier Scores for Men using DGM 1", col.names = c("Measure", "Estimate", "MCSE"))
set.seed(1234)
sim_brier_function_f <- function(sample_size){
  sim_target_female_df <- data.frame(
    SYSBP = rnorm(n = sample_size, mean = 122.46, sd = 18.76),
    AGE = rnorm(n = sample_size, mean = 46.75, sd = 9.85),
    HDLC = rnorm(n = sample_size, mean = 57.59, sd = 16.25),
    CURSMOKE = rbinom(n = sample_size, size = 1, prob = 0.171),
    BPMEDS = rbinom(n = sample_size, size = 1, prob = 0.8),
    TOTCHOL = rnorm(n = sample_size, mean = 195.40, sd = 38.95),
    DIABETES = rbinom(n = sample_size, size = 1, prob = 0.11),
    S = 0
  )
  
  sim_target_female_df$CURSMOKE <- as.factor(sim_target_female_df$CURSMOKE)
  sim_target_female_df$DIABETES <- as.factor(sim_target_female_df$DIABETES)
  sim_target_female_df$BPMEDS <- as.factor(sim_target_female_df$BPMEDS)
  
  # Get blood pressure based on whether or not on BPMEDS
  sim_target_female_df$SYSBP_UT <- ifelse(sim_target_female_df$BPMEDS == 0, 
                                          sim_target_female_df$SYSBP, 0)
  sim_target_female_df$SYSBP_T <- ifelse(sim_target_female_df$BPMEDS == 1, 
                                         sim_target_female_df$SYSBP, 0)
  
  #combining the testing data from both the source and target populations to get the inverse odds weights
  women_test_full_df_sim <- rbind(framingham_df_women_test[,-c(1,2,6,11,15)], sim_target_female_df)
  
  #models for odds of belonging to the source population
  
  #women
  women_source_lo_test_sim <- glm(S ~ log(HDLC)+log(TOTCHOL)+log(AGE)+log(SYSBP_UT+1)+
                                    log(SYSBP_T+1)+CURSMOKE+DIABETES,
                                  data = women_test_full_df_sim, family = "binomial")
  women_test_full_df_sim$odds <- (predict(women_source_lo_test_sim, newdata = women_test_full_df_sim, type = "response"))/(1-(predict(women_source_lo_test_sim, newdata = women_test_full_df_sim, type = "response"))) #obtaining odds
  women_test_full_df_sim$weights <- 1/women_test_full_df_sim$odds #inverse odds weights
  women_test_source_sim <- women_test_full_df_sim %>% filter(S == 1) #only source population data to fit model
  women_test_source_sim <- cbind(framingham_df_women_test$CVD, women_test_source_sim) %>%
    rename("CVD" = "framingham_df_women_test$CVD")
  women_test_source_sim$CVD <- as.numeric(women_test_source_sim$CVD)-1
  
  #brier score
  women_test_source_sim$pred_probs <- predict(mod_women, newdata = women_test_source_sim, type = "response") 
  women_test_source_sim$brier_num <- (women_test_source_sim$weights)*((women_test_source_sim$CVD - women_test_source_sim$pred_probs)^2)
  
  women_brier_score_sim <- sum(women_test_source_sim$brier_num)/nrow(sim_target_female_df)
  
  return(women_brier_score_sim)
}

women_brier_sim_res <- replicate(1600, sim_brier_function_f(1500))
output_women_df <- data.frame(dataset = 1:1600, bs = women_brier_sim_res)
simsum_women <- simsum(data = output_women_df, estvarname = "bs", true = women_brier_score)
rsimsum::kable(simsum_women, stats = c("thetamean", "thetamedian", "bias", "empse", "mse"), caption = "Performance Measures for the Simulations Estimating Brier Scores for Women using DGM 1", col.names = c("Measure", "Estimate", "MCSE"))
#sigma matrices
cor_mat_men <- cor(complete_df_2017_men[c("SYSBP", "AGE", "HDLC", "TOTCHOL")])
cor_mat_women <- cor(complete_df_2017_women[c("SYSBP", "AGE", "HDLC", "TOTCHOL")])

set.seed(1234)
sim_brier_function_m_2 <- function(sample_size){
  mv_res <- mvrnorm(n = sample_size, mu = c(126.07, 47.16, 47.45, 192.86),
                    Sigma = cor_mat_men)
  sim_target_male_df <- data.frame(
    SYSBP = mv_res[,1],
    AGE = mv_res[,2],
    HDLC = mv_res[,3],
    TOTCHOL = mv_res[,4],
    CURSMOKE = rbinom(n = sample_size, size = 1, prob = 0.259),
    BPMEDS = rbinom(n = sample_size, size = 1, prob = 0.755),
    DIABETES = rbinom(n = sample_size, size = 1, prob = 0.131),
    S = 0
  )
  
  sim_target_male_df$CURSMOKE <- as.factor(sim_target_male_df$CURSMOKE)
  sim_target_male_df$DIABETES <- as.factor(sim_target_male_df$DIABETES)
  sim_target_male_df$BPMEDS <- as.factor(sim_target_male_df$BPMEDS)
  
  # Get blood pressure based on whether or not on BPMEDS
  sim_target_male_df$SYSBP_UT <- ifelse(sim_target_male_df$BPMEDS == 0, 
                                        sim_target_male_df$SYSBP, 0)
  sim_target_male_df$SYSBP_T <- ifelse(sim_target_male_df$BPMEDS == 1, 
                                       sim_target_male_df$SYSBP, 0)
  
  #combining the testing data from both the source and target populations to get the inverse odds weights
  men_test_full_df_sim <- rbind(framingham_df_men_test[,-c(1,2,6,11,15)], sim_target_male_df)
  
  #models for odds of belonging to the source population
  
  #men
  men_source_lo_test_sim <- glm(S ~ HDLC + TOTCHOL + AGE + SYSBP_UT + SYSBP_T + CURSMOKE + DIABETES,
                                data = men_test_full_df_sim, family = "binomial")
  men_test_full_df_sim$odds <- (predict(men_source_lo_test_sim, newdata = men_test_full_df_sim, type = "response"))/(1-(predict(men_source_lo_test_sim, newdata = men_test_full_df_sim, type = "response"))) #obtaining odds
  men_test_full_df_sim$weights <- 1/men_test_full_df_sim$odds #inverse odds weights
  men_test_source_sim <- men_test_full_df_sim %>% filter(S == 1) #only source population data to fit model
  men_test_source_sim <- cbind(framingham_df_men_test$CVD, men_test_source_sim) %>%
    rename("CVD" = "framingham_df_men_test$CVD")
  men_test_source_sim$CVD <- as.numeric(men_test_source_sim$CVD)-1
  
  #brier score
  men_test_source_sim$pred_probs <- predict(mod_men, newdata = men_test_source_sim, type = "response") 
  men_test_source_sim$brier_num <- (men_test_source_sim$weights)*((men_test_source_sim$CVD - men_test_source_sim$pred_probs)^2)
  
  men_brier_score_sim <- sum(men_test_source_sim$brier_num)/nrow(sim_target_male_df)
  
  return(men_brier_score_sim)
}

men_brier_sim_res_2 <- replicate(1600, sim_brier_function_m_2(1500))
output_men_df_2 <- data.frame(dataset = 1:1600, bs = men_brier_sim_res_2)
simsum_men_2 <- simsum(data = output_men_df_2, estvarname = "bs", true = men_brier_score)
rsimsum::kable(simsum_men_2, stats = c("thetamean", "thetamedian", "bias", "empse", "mse"), caption = "Performance Measures for the Simulations Estimating Brier Scores for Men using DGM 2", col.names = c("Measure", "Estimate", "MCSE"))
set.seed(1234)
sim_brier_function_f_2 <- function(sample_size){
  mv_res <- mvrnorm(n = sample_size, mu = c(122.46, 46.75, 57.59, 195.40),
                    Sigma = cor_mat_women)
  sim_target_female_df <- data.frame(
    SYSBP = mv_res[,1],
    AGE = mv_res[,2],
    HDLC = mv_res[,3],
    TOTCHOL = mv_res[,4],
    CURSMOKE = rbinom(n = sample_size, size = 1, prob = 0.171),
    BPMEDS = rbinom(n = sample_size, size = 1, prob = 0.8),
    DIABETES = rbinom(n = sample_size, size = 1, prob = 0.11),
    S = 0
  )
  
  sim_target_female_df$CURSMOKE <- as.factor(sim_target_female_df$CURSMOKE)
  sim_target_female_df$DIABETES <- as.factor(sim_target_female_df$DIABETES)
  sim_target_female_df$BPMEDS <- as.factor(sim_target_female_df$BPMEDS)
  
  # Get blood pressure based on whether or not on BPMEDS
  sim_target_female_df$SYSBP_UT <- ifelse(sim_target_female_df$BPMEDS == 0, 
                                          sim_target_female_df$SYSBP, 0)
  sim_target_female_df$SYSBP_T <- ifelse(sim_target_female_df$BPMEDS == 1, 
                                         sim_target_female_df$SYSBP, 0)
  
  #combining the testing data from both the source and target populations to get the inverse odds weights
  women_test_full_df_sim <- rbind(framingham_df_women_test[,-c(1,2,6,11,15)], sim_target_female_df)
  
  #models for odds of belonging to the source population
  
  #women
  women_source_lo_test_sim <- glm(S ~ HDLC + TOTCHOL + AGE + SYSBP_UT + SYSBP_T + CURSMOKE + DIABETES,
                                  data = women_test_full_df_sim, family = "binomial")
  women_test_full_df_sim$odds <- (predict(women_source_lo_test_sim, newdata = women_test_full_df_sim, type = "response"))/(1-(predict(women_source_lo_test_sim, newdata = women_test_full_df_sim, type = "response"))) #obtaining odds
  women_test_full_df_sim$weights <- 1/women_test_full_df_sim$odds #inverse odds weights
  women_test_source_sim <- women_test_full_df_sim %>% filter(S == 1) #only source population data to fit model
  women_test_source_sim <- cbind(framingham_df_women_test$CVD, women_test_source_sim) %>%
    rename("CVD" = "framingham_df_women_test$CVD")
  women_test_source_sim$CVD <- as.numeric(women_test_source_sim$CVD)-1
  
  #brier score
  women_test_source_sim$pred_probs <- predict(mod_women, newdata = women_test_source_sim, type = "response") 
  women_test_source_sim$brier_num <- (women_test_source_sim$weights)*((women_test_source_sim$CVD - women_test_source_sim$pred_probs)^2)
  
  women_brier_score_sim <- sum(women_test_source_sim$brier_num)/nrow(sim_target_female_df)
  
  return(women_brier_score_sim)
}

women_brier_sim_res_2 <- replicate(1600, sim_brier_function_f_2(1500))
output_women_df_2 <- data.frame(dataset = 1:1600, bs = women_brier_sim_res_2)
simsum_women_2 <- simsum(data = output_women_df_2, estvarname = "bs", true = women_brier_score)
rsimsum::kable(simsum_women_2, stats = c("thetamean", "thetamedian", "bias", "empse", "mse"), caption = "Performance Measures for the Simulations Estimating Brier Scores for Women using DGM 2", col.names = c("Measure", "Estimate", "MCSE"))