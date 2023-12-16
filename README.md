# PHP2550_Final_Portfolio

This Repository contains the code, R Markdown files and PDF files of my final, revised reports of the three projects completed for PHP 2550 taught by Dr. Alice Paul in Fall 2023. 

## Project 1

### Objective

The main objective of this exploratory data analysis is to determine the presence of any associations between exposure to smoking during pregnancy (SDP) and environmental tobacco smoke (ETS), and child outcomes related to substance use, self-regulation, internalizing problems (IP), externalizing problems (EP), attention problems (AP), Attention Deficit Hyperactivity Disorder (ADHD), and Autism Spectrum Disorder (ASD).

### Data and Preprocessing

The data used in this project originates from a collaboration with Dr. Lauren Micalizzi, a researcher in the Department of Behavioral and Social Sciences at Brown University. Dr. Micalizzi's work focuses on investigating the impact of prenatal tobacco exposure on child outcomes, with a particular emphasis on exposure to smoking during pregnancy (SDP) and environmental tobacco smoke (ETS).

Exposure to SDP and ETS represents two of the most prevalent and detrimental environmental factors affecting children (Micalizzi). Studies indicate that between seven to fifteen percent of infants born each year are exposed to smoking during pregnancy, and more than twenty-five percent of children are exposed to household environmental tobacco smoke. The economic burden of SDP alone is substantial, imposing a $4 billion annual cost in healthcare expenses.

Previous research by Dr. Micalizzi involved a randomized controlled trial aimed at reducing SDP and ETS through a tailored video intervention. This trial, which included 738 low-income women, targeted smoke avoidance during pregnancy and sought to diminish children's exposure to ETS in the immediate postpartum period.

For the current project, a subset of adolescents (N=100) and their mothers was randomly selected from the larger cohort recruited in the previous study. The purpose of this investigation is to further understand the association between smoking during pregnancy (SDP) and environmental tobacco smoke (ETS) exposure and its impact on self-regulation, externalizing behavior, and substance use in children.

There are 49 observations of 78 variables in this data set. Each observation corresponds to a particular parent-child pair and the variables contain information on demographics, on exposure, and on outcomes in children collected from both parents and the children. The data also contains information on outcomes in the parents, but these are not of interest in our current study.

#### Variables

The variables included in this analysis include those on:

-   Parent demographics (race, age, sex, number of languages spoken, employment, highest level of education and annual household income),
-   Child demographics (race, age, sex, number of languages spoken),
-   Exposure to SDP (self-reported information from the parent about whether they smoked when 16 weeks, 22 weeks and 32 weeks pregnant and urine cotinine (UC) levels from parent at 34 weeks pregnant)
-   Exposure to ETS (self-reported information on whether the parent or partner smoked during the first 5 years of the child's life, on whether the parent smoked at the first and second postpartum visits, and at 12 weeks and 6 months postpartum, and urine cotinine levels from both mom and baby at 6 months postpartum)
-   Substance use in child (if cigarettes, e-cigarettes/vape, marijuana, or alcohol were ever used and if so, how frequently in the last 30 days),
-   Internalizing, attention and externalizing problems in child (self-reported and parent-reported scores on the Brief Problem Monitor),
-   Self-regulation in child (self-reported average cognitive reappraisal (CR) and expressive suppression (ES) scores on the Emotion Regulation Questionnaire),
-   ADHD in child (SWAN scores on the inattentive and hyperactive items),
-   ASD in child (absence, presence or suspicion)

### Files

`Chandru_PHP2550_Project1_Revised.Rmd` is an R Markdown file containing the code and report pertaining to this analysis.

`Chandru_PHP2550_Project1_Revised.pdf` is a pdf file containing the code and report pertaining to this analysis.

`Chandru_PHP2550_Project1_Revised.R` is an R script file containing the code used in this analysis.

## Project 2

### Introduction
This project, undertaken in partnership with Dr. Chris Schmid from the Biostatistics Department, addresses the ongoing challenges surrounding severe bronchopulmonary dysplasia (sBPD) in neonatal care. The uncertainty surrounding precise indication criteria and optimal timing for tracheostomy placement in neonates with sBPD motivates this research initiative. Existing studies suggest potential benefits associated with earlier tracheostomy placement, particularly concerning neonatal growth (Schmid). Notably, previous analyses of extensive databases have successfully predicted the likelihood of tracheostomy placement or death based on baseline demographics and clinical diagnoses. However, these analyses lacked the inclusion of detailed respiratory parameters and did not offer predictions at various postmenstrual ages (PMA). This project aims to fill this gap by incorporating a comprehensive set of respiratory variables measured at critical time points, specifically at 36 and 44 weeks PMA. 

### Data Overview
The dataset originates from the BPD Collaborative Registry, a multi-center consortium of interdisciplinary BPD programs across the United States and Sweden. It includes infants with gestational age below 32 weeks and diagnosed with sBPD according to the 2001 NHLBI criteria. Standard demographic and clinical data are collected at four critical time points: birth, 36 weeks postmenstrual age (PMA), 44 weeks PMA, and discharge. The dataset, covering the period between January 1 and July 19, 2021, was queried for patients with sBPD and complete growth data. Ten BPD Collaborative centers contributed data meeting the study's inclusion criteria.

The data analyzed here contained different types of variables, and consisted of 996 observations of 30 variables. Patient ID and center number were included. The demographic variables included were maternal race and ethnicity (1 = Hispanic/Latino, 2 = Not Hispanic or Latino). The birth variables included were birth weight, obstetrical gestational age, whether the infant was small for gestational age, birth length, birth head circumference, delivery method (1 = Vaginal Delivery, 2 = Cesarean Section), whether prenatal corticosteroids were administered, whether complete prenatal steroids were administered, whether maternal chorioamnionitis was present, gender, and whether the infant received any surfactant at any point in the first 72 hours. The weight of the infant, and respiratory support variables such as ventilation support (0 = No respiratory support or supplemental oxygen, 1 = Non-invasive positive pressure, 2 = Invasive positive pressure), fraction of inspired oxygen, peak inspiratory pressure, positive end exploratory pressure, and whether medication for pulmonary hypertension was administered, were recorded at both 36 and 44 weeks. In addition, data on the infant’s gestational age at the time of discharge and whether tracheostomy (0 = No, 1 = Yes) or death had occurred at this point were recorded.

### Files

`Chandru_PHP2550_Project2_Revised.Rmd` is an R Markdown file containing the code and report pertaining to this analysis.

`Chandru_PHP2550_Project2_Revised.pdf` is a pdf file containing the code and report pertaining to this analysis.

`Chandru_PHP2550_Project2_Revised.R` is an R script file containing the code used in this analysis.

## Project 3

### Introduction

Cardiovascular disease (CVD) remains a leading cause of morbidity and mortality worldwide, necessitating accurate risk prediction models for timely intervention and prevention (World Health Organization, 2021). In the realm of predictive modeling, the application of risk prediction models across diverse populations is of paramount importance. This collaborative project, conducted in partnership with Dr. Jon Steingrimsson from the Biostatistics Department, focuses on evaluating the transportability of a cardiovascular risk prediction model originally developed on a source population when applied to a target population both when data from the target population is available, and when it is simulated.

The predictive model under scrutiny originates from the Framingham Heart Study, a seminal investigation comprising participants aged 30-62 (National Heart, Lung, and Blood Institute). Historically, such models, while valuable in their original context, have demonstrated challenges when extrapolated to populations with different demographic compositions, as seen in the example of the Framingham ATP-III model's suboptimal generalization to multi-ethnic populations (Steingrimsson).

The objectives of this project are two-fold:

-   Evaluate performance of a cardiovascular risk prediction model (D'Agostino et al., 2008) in a target population underlying NHANES: The evaluation is conducted using the weighting estimator for the Brier score in the target population (Steingrimsson et al., 2022).
-   Conduct the same analysis when the target population is simulated: The simulation process involves utilizing summary statistics from the NHANES data to simulate the target population. The same estimator for weighted Brier score in the target population is used to evaluate the performance of the models in the simulated populations.

### Data and Preprocessing

The Framingham data is obtained from the `riskCommunicator` package, and it comes from the Framingham Heart Study, a study that began in 1948 by collecting data from individuals between the ages of 30-62 from Framingham, Massachusetts, with the aim of identifying risk factors for cardiovascular disease (National Heart, Lung, and Blood Institute).

NHANES data comes from annual surveys conducted by the NHANES program, whose participants are selected to be representative of the U.S. population belonging to all ages (Centers for Disease Control and Prevention). The NHANES data is obtained from the `nhanesA` package.

The preprocessing involves using data from the Framingham study to obtain a sample of source population data. From the Framingham data, the variables CVD (indicating whether myocardial infarction, fatal coronary heart disease, atherothrombotic infarction, cerebral embolism, intracerebral hemorrhage, or subarachnoid hemorrhage or fatal cerebrovascular disease occured during followup), TIMECVD (number of days from baseline exam to first CVD event during followup or number of days from baseline to censor date), SEX (participant sex), TOTCHOL (serum total cholesterol in mg/dL), AGE (age at exam), SYSBP (systolic blood pressure in mmHg; mean of last two of three measurements), DIABP (diastolic blood pressure in mmHg; mean of last two of three measurements), CURSMOKE (indicating whether participants are current smokers), DIABETES (indicating whether the participant is diabetic), BPMEDS (whether anti-hypertensive medication was being used at the time of exam), HDLC (high density lipoprotein cholesterol in mg/dL) and BMI (body mass index) are extracted, and any NA's are removed.

Different variables (`SYSBP_T` and `SYSBP_UT`) are created to store the blood pressure values depending on whether or not the individuals were on anti-hypertensive medication. In addition, any observations corresponding to individuals without cardiovascular events within 15 years were removed (removal of censored data). Two different datasets corresponding to males and females was then created.

Then, from the NHANES data for the year 2017-18, variables corresponding to those selected from the Framingham study (above) are similarly selected to create a separate dataset corresponding to a sample of the target population data. The observations are filtered to only include data corresponding to participants belonging to ages 30-62 so that the eligibility criteria for the source population (Framingham study participants) is met.

### Files

`Chandru_PHP2550_Project3_Revised.Rmd` is an R Markdown file containing the code and report pertaining to this analysis.

`Chandru_PHP2550_Project3_Revised.pdf` is a pdf file containing the code and report pertaining to this analysis.

`Chandru_PHP2550_Project3_Revised.R` is an R script file containing the code used in this analysis.
