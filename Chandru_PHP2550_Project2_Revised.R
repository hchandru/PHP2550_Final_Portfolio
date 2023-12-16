knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "~/Downloads")

#packages
install.packages("Hmisc", repos = "http://cran.us.r-project.org")
install.packages("performance", repos = "http://cran.us.r-project.org")
install.packages("webshot2", repos = "http://cran.us.r-project.org")
install.packages("rcompanion", repos = "http://cran.us.r-project.org")

library(tidyverse)
library(naniar)
library(mice)
library(gtsummary)
library(gt)
library(webshot2)
library(vcd)
library(corrplot)
library(readr)
library(labelled)
library(Hmisc)
library(glmnet)
library(pROC)
library(lme4)
library(knitr)
library(kableExtra)
#laoding the data
df <- read_csv("project2.csv")
#changing the values of the sga variable and removing the mat_race variable
df <- df %>%
  mutate(sga = case_when(sga == "Not SGA" ~ "No",
                         sga == "SGA" ~ "Yes")) %>%
  select(-mat_race)

length(unique(df$record_id)) #only 996 unique ids but we have 999 observations
df <- df[!duplicated(df$record_id),] #removing the duplicates

#getting the center number from the record id
df$center <- ifelse(nchar(as.character(df$record_id)) == 7, as.numeric(substr(df$record_id, 1, 1)), as.numeric(substr(df$record_id, 1, 2)))

#Change variables to factors
df$center <- as.factor(df$center)
df$mat_ethn <- as.factor(df$mat_ethn)
df$del_method <- as.factor(df$del_method)
df$prenat_ster <- as.factor(df$prenat_ster)
df$com_prenat_ster <- as.factor(df$com_prenat_ster)
df$mat_chorio <- as.factor(df$mat_chorio)
df$gender <- as.factor(df$gender)
df$sga <- as.factor(df$sga)
df$any_surf <- as.factor(df$any_surf)
df$ventilation_support_level.36 <- as.factor(df$ventilation_support_level.36)
df$med_ph.36 <- as.factor(df$med_ph.36)
df$ventilation_support_level_modified.44 <- as.factor(df$ventilation_support_level_modified.44)
df$med_ph.44 <- as.factor(df$med_ph.44)
df$Trach <- as.factor(df$Trach)
df$Death <- as.factor(df$Death)

#create variable labels
var_label(df) <- list(
  record_id = "Patient ID",
  center = "Medical Center",
  mat_ethn = "Maternal Ethnicity",
  bw = "Birth Weight (g)",
  ga = "Obstetrical Gestational Age",
  blength = "Birth Length (cm)",
  birth_hc = "Birth Head Circumference (cm)",
  del_method = "Delivery Method",
  prenat_ster = "Prenatal Corticosteroids",
  com_prenat_ster = "Complete Prenatal Steroids",
  mat_chorio = "Maternal Chorioamnionitis",
  gender = "Gender",
  sga = "Small for Gestational Age",
  any_surf = "Received Surfactant in the First 72 Hours",
  weight_today.36 = "Weight at 36 Weeks",
  ventilation_support_level.36 = "Ventilation Support Level at 36 Weeks",
  inspired_oxygen.36 = "Fraction of Inspired Oxygen at 36 Weeks",
  p_delta.36 = "Peak Inspiratory Pressure (cm H2O) at 36 Weeks",
  peep_cm_h2o_modified.36 = "Positive End Exploratory Pressure (cm H2O) at 36 Weeks",
  med_ph.36 = "Medication for Pulmonary Hypertension at 36 Weeks",
  weight_today.44 = "Weight at 44 Weeks",
  ventilation_support_level_modified.44 = "Ventilation Support Level at 44 Weeks",
  inspired_oxygen.44 = "Fraction of Inspired Oxygen at 44 Weeks",
  p_delta.44 = "Peak Inspiratory Pressure (cm H2O) at 44 Weeks",
  peep_cm_h2o_modified.44 = "Positive End Exploratory Pressure (cm H2O) at 44 Weeks",
  med_ph.44 = "Medication for Pulmonary Hypertension at 44 Weeks",
  hosp_dc_ga = "Hospital Discharge Gestational Age",
  Trach = "Tracheostomy",
  Death = "Death"
)
#missing data summary
df_miss <- df
colnames(df_miss) <- label(df)

gg_miss_var(df_miss, show_pct = TRUE) +
  theme(axis.text.y = element_text(size = 5)) +
  ggtitle("Figure 1. Percent Missing Observations by Variable") +
  theme(plot.title = element_text(size = 10))
theme_gtsummary_compact(set_theme=TRUE, font_size = 10)

#table of descriptive statistics - overall
df_miss <- df_miss[,-1]

overall_desc_stats <- df_miss %>%
  tbl_summary(missing = "no") %>%
  as_gt() %>%
  gt::tab_header(title = "Table 1. Overall Descriptive Statistics") %>%
  gtsave(filename = "proj2_table1.png")

knitr::include_graphics('proj2_table1.png')
#table of descriptive statistics by center
desc_stats_center <- df_miss %>%
  tbl_summary(by = "Medical Center", missing_text = "Missing") %>%
  as_gt() %>%
  gt::tab_header(title = "Table 2. Descriptive Statistics by Center") %>%
  gtsave(filename = "proj2_table2.png")

knitr::include_graphics('proj2_table2.png')

#Correlation matrix for the continuous variables
numeric_df <- df %>%
  select(c("bw", "ga", "blength", "birth_hc", "weight_today.36",
           "inspired_oxygen.36", "p_delta.36", "peep_cm_h2o_modified.36",
           "weight_today.44", "inspired_oxygen.44", "p_delta.44", "peep_cm_h2o_modified.44", "hosp_dc_ga"))

correlation_matrix <- cor(numeric_df[complete.cases(numeric_df),])
corrplot(correlation_matrix, method = "number", number.cex = 0.5, title = "Figure 2. Correlations Between Continuous Variables", tl.cex = 0.5, mar=c(0,0,2,0))
#ventilation support level at 44 weeks vs 36 weeks
assocstats(table(df$ventilation_support_level_modified.44, df$ventilation_support_level.36))$chisq_tests[2,3]

#medication for pulmonary hypertension at 44 weeks vs 36 weeks
assocstats(table(df$med_ph.44, df$med_ph.36))$chisq_tests[2,3]

#complete prenatal steroids vs prenatal steroids
assocstats(table(df$com_prenat_ster, df$prenat_ster))$chisq_tests[2,3]
set.seed(155)
#to create test and train data sets
ignore <- sample(c(TRUE, FALSE), 996, replace = TRUE, prob = c(0.3,0.7))

#data using tracheostomy as the outcome (removing death), and removing any_surf since it has a lot of missing data and cannot be imputed with information from any other variables
df <- df %>%
  select(-c(Death, any_surf, record_id))

#train and test data sets
train_df <- df[!ignore, ]
test_df <- df[ignore, ]

#imputing 5 test and train datasets
imp.train <- mice(train_df[,-c(1,25)], m = 5, print = FALSE, seed = 155)
imp.test <- mice.mids(imp.train, newdata = test_df[,-c(1,25)])

df_train_imp_36 <- vector("list", length = 5)
df_test_imp_36 <- vector("list", length = 5)
df_train_imp_44 <- vector("list", length = 5)
df_test_imp_44 <- vector("list", length = 5)

for(i in 1:5){
  df_train_imp_36[[i]] <- mice::complete(imp.train, i) %>%
    select(-c(weight_today.44, ventilation_support_level_modified.44, inspired_oxygen.44, p_delta.44, peep_cm_h2o_modified.44, med_ph.44))
  df_train_imp_36[[i]]$center <- train_df$center
  df_train_imp_36[[i]]$hosp_dc_ga <- train_df$hosp_dc_ga
  
  df_test_imp_36[[i]] <- mice::complete(imp.test, i) %>%
    select(-c(weight_today.44, ventilation_support_level_modified.44, inspired_oxygen.44, p_delta.44, peep_cm_h2o_modified.44, med_ph.44))
  df_test_imp_36[[i]]$center <- test_df$center
  df_test_imp_36[[i]]$hosp_dc_ga <- test_df$hosp_dc_ga
  
  df_train_imp_44[[i]] <- mice::complete(imp.train, i) %>%
    select(-c(weight_today.36, ventilation_support_level.36, inspired_oxygen.36, p_delta.36, peep_cm_h2o_modified.36, med_ph.36))
  df_train_imp_44[[i]]$center <- train_df$center
  df_train_imp_44[[i]]$hosp_dc_ga <- train_df$hosp_dc_ga
  
  df_test_imp_44[[i]] <- mice::complete(imp.test, i) %>%
    select(-c(weight_today.36, ventilation_support_level.36, inspired_oxygen.36, p_delta.36, peep_cm_h2o_modified.36, med_ph.36))
  df_test_imp_44[[i]]$center <- test_df$center
  df_test_imp_44[[i]]$hosp_dc_ga <- test_df$hosp_dc_ga
  
}

#storing all the imputed data 
#only those discharged after 36 weeks
df_train_imp_36_long <- mice::complete(imp.train, action = "long") %>%
  select(-c(weight_today.44, ventilation_support_level_modified.44, inspired_oxygen.44, p_delta.44, peep_cm_h2o_modified.44, med_ph.44))
df_train_imp_36_long$center <- rep(train_df$center, 5)
df_train_imp_36_long$hosp_dc_ga <- rep(train_df$hosp_dc_ga, 5)
df_train_imp_36_long <- df_train_imp_36_long %>%
  filter(hosp_dc_ga > 36)

#only those discharged after 36 weeks
df_test_imp_36_long <- mice::complete(imp.test, action = "long") %>%
  select(-c(weight_today.44, ventilation_support_level_modified.44, inspired_oxygen.44, p_delta.44, peep_cm_h2o_modified.44, med_ph.44))
df_test_imp_36_long$center <- rep(test_df$center, 5)
df_test_imp_36_long$hosp_dc_ga <- rep(test_df$hosp_dc_ga, 5)
df_test_imp_36_long <- df_test_imp_36_long %>%
  filter(hosp_dc_ga > 36)

#only those discharged after 44 weeks
df_train_imp_44_long <- mice::complete(imp.train, action = "long") %>%
  select(-c(weight_today.36, ventilation_support_level.36, inspired_oxygen.36, p_delta.36, peep_cm_h2o_modified.36, med_ph.36))
df_train_imp_44_long$center <- rep(train_df$center, 5)
df_train_imp_44_long$hosp_dc_ga <- rep(train_df$hosp_dc_ga, 5)
df_train_imp_44_long <- df_train_imp_44_long %>%
  filter(hosp_dc_ga > 44)

#only those discharged after 44 weeks
df_test_imp_44_long <- mice::complete(imp.test, action = "long") %>%
  select(-c(weight_today.36, ventilation_support_level.36, inspired_oxygen.36, p_delta.36, peep_cm_h2o_modified.36, med_ph.36))
df_test_imp_44_long$center <- rep(test_df$center, 5)
df_test_imp_44_long$hosp_dc_ga <- rep(test_df$hosp_dc_ga, 5)
df_test_imp_44_long <- df_test_imp_44_long %>%
  filter(hosp_dc_ga > 44)
lasso <- function(df) { 
  #' Runs 10-fold CV for lasso and returns corresponding coefficients 
  #' @param df, data set
  #' @return coef, coefficients for minimum cv error
  
  # Matrix form for ordered variables 
  x.ord <- model.matrix(Trach~., data = df[,-c(19,20)])[,-1] 
  y.ord <- df$Trach 
  
  # Generate folds
  k <- 10 
  set.seed(1) # consistent seeds between imputed data sets
  folds <- sample(1:k, nrow(df), replace=TRUE)
  
  # Lasso model
  lasso_mod_cv <- cv.glmnet(x.ord, y.ord, nfolds = 10, foldid = folds, 
                            alpha = 1, family = "binomial") 
  lasso_mod <- glmnet(x.ord, y.ord, nfolds = 10, alpha = 1, family = "binomial",
                      lambda = lasso_mod_cv$lambda.min)
  
  # Get coefficients 
  coef <- coef(lasso_mod) 
  return(coef) 
} 

# Find average lasso coefficients over imputed datasets - 36 week data
lasso_coef1_36 <- lasso(df_train_imp_36[[1]]) 
lasso_coef2_36 <- lasso(df_train_imp_36[[2]]) 
lasso_coef3_36 <- lasso(df_train_imp_36[[3]]) 
lasso_coef4_36 <- lasso(df_train_imp_36[[4]]) 
lasso_coef5_36 <- lasso(df_train_imp_36[[5]]) 
lasso_coef_36 <- cbind(lasso_coef1_36, lasso_coef2_36, lasso_coef3_36, 
                       lasso_coef4_36, lasso_coef5_36) 
avg_coefs_lasso_36 <- apply(lasso_coef_36, 1, mean) 

# Find average lasso coefficients over imputed datasets - 44 week data
lasso_coef1_44 <- lasso(df_train_imp_44[[1]]) 
lasso_coef2_44 <- lasso(df_train_imp_44[[2]]) 
lasso_coef3_44 <- lasso(df_train_imp_44[[3]]) 
lasso_coef4_44 <- lasso(df_train_imp_44[[4]]) 
lasso_coef5_44 <- lasso(df_train_imp_44[[5]]) 
lasso_coef_44 <- cbind(lasso_coef1_44, lasso_coef2_44, lasso_coef3_44, 
                       lasso_coef4_44, lasso_coef5_44) 
avg_coefs_lasso_44 <- apply(lasso_coef_44, 1, mean) 
wk36_vars <- names(avg_coefs_lasso_36[avg_coefs_lasso_36!=0])[-1]
wk36_or <- data.frame(Variable = c("Mother is not Hispanic/Latino", "Birth weight in grams", "Birth length in cm", "Birth head circumference", "Cesarean section delivery", "Prenatal steroids were administered", "Complete prenatal steroids were administered", "Maternal chorioamnionitis was present", "Infant was male", "Infant weight at 36 weeks", "Ventilation support through invasive positive pressure at 36 weeks", "Fraction of inspired oxygen at 36 weeks", "Peak inspiratory pressure at 36 weeks", "Positive end exploratory pressure at 36 weeks", "Medication for pulmonary hypertension was administered at 36 weeks"), OR = exp(as.numeric(avg_coefs_lasso_36[avg_coefs_lasso_36!=0]))[-1]) %>%
  kbl(caption = "Odds Ratios from the 36-Week Model", booktabs = T) 
wk36_or
wk44_vars <- names(avg_coefs_lasso_44[avg_coefs_lasso_44!=0])[-1]
wk44_or <- data.frame(Variable = c("Mother is not Hispanic/Latino", "Birth weight in grams", "Obstetrical gestational age", "Birth length in cm", "Birth head circumference", "Cesarean section delivery", "Prenatal steroids were administered", "Complete prenatal steroids were administered", "Maternal chorioamnionitis was present", "Infant was male", "Infant was small for gestational age", "Infant weight at 44 weeks", "Ventilation support through non-invasive positive pressure at 44 weeks", "Ventilation support through invasive positive pressure at 44 weeks", "Fraction of inspired oxygen at 44 weeks", "Peak inspiratory pressure at 44 weeks", "Positive end exploratory pressure at 44 weeks", "Medication for pulmonary hypertension was administered at 44 weeks"), OR = exp(as.numeric(avg_coefs_lasso_44[avg_coefs_lasso_44!=0]))[-1]) %>%
  kbl(caption = "Odds Ratios from the 44-Week Model", booktabs =T)
wk44_or
# Find predicted probabilities and predicted classes on long imputed test data 

#36 week model
x_vars_36 <- model.matrix(Trach~. , df_test_imp_36_long[,-c(1,2,21,22)])
df_test_imp_36_long$pred_probs <- as.vector(plogis(x_vars_36 %*% avg_coefs_lasso_36)) #predicted probabilities
df_test_imp_36_long$pred <- ifelse(df_test_imp_36_long$pred_probs > 0.1, 1, 0) #predicted classes

#44 week model
x_vars_44 <- model.matrix(Trach~. , df_test_imp_44_long[,-c(1,2,21,22)])
df_test_imp_44_long$pred_probs <- as.vector(plogis(x_vars_44 %*% avg_coefs_lasso_44)) #predicted probabilities
df_test_imp_44_long$pred <- ifelse(df_test_imp_44_long$pred_probs > 0.103, 1, 0) #predicted classes
#36 week model

#confusion matrix
confusion_matrix_36 <- table(Observed = df_test_imp_36_long$Trach, Predicted = df_test_imp_36_long$pred)

#accuracy
accuracy_36 <- sum(diag(confusion_matrix_36)) / sum(confusion_matrix_36)

#precision (positive predictive value)
precision_36 <- confusion_matrix_36[2, 2] / sum(confusion_matrix_36[, 2])

#recall (sensitivity / true positive rate)
recall_36 <- confusion_matrix_36[2, 2] / sum(confusion_matrix_36[2, ])

#specificity (true negative rate)
specificity_36 <- confusion_matrix_36[1, 1] / sum(confusion_matrix_36[1, ])

#f score
f_score_36 <- 2 * (precision_36 * recall_36) / (precision_36 + recall_36)

#roc and auc
roc_obj_36 <- roc(df_test_imp_36_long$Trach, df_test_imp_36_long$pred_probs, auc=TRUE)
# plot(roc_obj_36, print.thres = TRUE)
auc_roc_36 <- auc(roc_obj_36)

#brier score
brier_score_36 <- mean((df_test_imp_36_long$pred_probs - (as.numeric(df_test_imp_36_long$Trach)-1))^2)

#44 week model

#confusion matrix
confusion_matrix_44 <- table(Observed = df_test_imp_44_long$Trach, Predicted = df_test_imp_44_long$pred)

#accuracy
accuracy_44 <- sum(diag(confusion_matrix_44)) / sum(confusion_matrix_44)

#precision (positive predictive value)
precision_44 <- confusion_matrix_44[2, 2] / sum(confusion_matrix_44[, 2])

#recall (sensitivity / true positive rate)
recall_44 <- confusion_matrix_44[2, 2] / sum(confusion_matrix_44[2, ])

#specificity (true negative rate)
specificity_44 <- confusion_matrix_44[1, 1] / sum(confusion_matrix_44[1, ])

#f score
f_score_44 <- 2 * (precision_44 * recall_44) / (precision_44 + recall_44)

#roc and auc
roc_obj_44 <- roc(df_test_imp_44_long$Trach, df_test_imp_44_long$pred_probs)
# plot(roc_obj_44, print.thres = TRUE)
auc_roc_44 <- auc(roc_obj_44)

#brier score
brier_score_44 <- mean((df_test_imp_44_long$pred_probs - (as.numeric(df_test_imp_44_long$Trach)-1))^2)

perf_tab <- data.frame(
  Measure = c("Accuracy", "Precision", "Recall", "Specificity", "F Score", "AUC", "Brier Score"),
  Week36_Model = c(
    sprintf("%.4f", accuracy_36),
    sprintf("%.4f", precision_36),
    sprintf("%.4f", recall_36),
    sprintf("%.4f", specificity_36),
    sprintf("%.4f", f_score_36),
    sprintf("%.4f", auc_roc_36),
    sprintf("%.4f", brier_score_36)
  ),
  Week44_Model = c(
    sprintf("%.4f", accuracy_44),
    sprintf("%.4f", precision_44),
    sprintf("%.4f", recall_44),
    sprintf("%.4f", specificity_44),
    sprintf("%.4f", f_score_44),
    sprintf("%.4f", auc_roc_44),
    sprintf("%.4f", brier_score_44)
  )
)

perf_tab %>%
  kbl(col.names = c("Measure", "36 Week Model", "44 Week Model"), caption = "Performance Measures for the Two Models", booktabs = T) %>%
  kable_styling(latex_options = "hold_position")