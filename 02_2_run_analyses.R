# Example script: GROUP CFA, Chi-square trend test, SEM and multigroup SEM
library(tidyverse)
library(mice)
library(lavaan)
library(semTools)
library(purrr)


# GROUP CFA, Chi-square trend test, SEM and multigroup SEM

eatdis$cohort <- factor(eatdis$cohort)
levels(eatdis$cohort)
eatdis$cohort <- factor(eatdis$cohort, levels = c("non_pandemic", "pandemic"))

# Group CFA 
twofactor_model <-'f1 =~ Q1 + Q2 + Q3 + Q4 + Q5
                   f2 =~ Q6 + Q7 + Q8 + Q9
                  '

fit1 <- cfa(twofactor_model, data = eatdis, estimator = "WLSMV", group = "cohort", group.label =c("non_pandemic", "pandemic"), 
            ordered = TRUE, meanstructure = TRUE)


fit2 <- cfa(twofactor_model, data = eatdis, estimator = "WLSMV", group = "cohort", group.label =c("non_pandemic", "pandemic"), 
            ordered = TRUE, meanstructure = TRUE,
            group.equal = c("thresholds","loadings"))

summary(fit2, fit.measures = TRUE, standardized = TRUE)
anova(fit1,fit2)


# Multigroup SEM with Year answered

multigroup_SEM <-"f1 =~ Q1 + Q2 + Q3 + Q4 + Q5
                        f2 =~ Q6 + Q7 + Q8 + Q9

                        f1 ~  YA_2
                        f2 ~  YA_2"


fit <- sem(multigroup_SEM, data = eatdis, estimator = "WLSMV", group = "cohort", group.label =c("non_pandemic", "pandemic"), 
                 ordered = TRUE, meanstructure = TRUE,
                 group.equal = c("thresholds","loadings"))

summary(fit, fit.measures = TRUE, standardized = TRUE)


# create long file for this test

eatdis_1 <- eatdis %>%
  select(cohort, Q1:Q9)

eatdis_2 <- eatdis_1 %>%
  pivot_longer(matches("Q")) %>%
  count(cohort,name,value) %>% 
  group_by(cohort, name) %>%
  mutate(prop=prop.table(n))%>%
  na.omit()

# 
perform_test <- function(data, question, pandemic_col = "pandemic", non_pandemic_col = "non_pandemic") {
  unexp <- data %>%
    filter(cohort == non_pandemic_col, name == question) %>%
    .$n
  exp <- data %>%
    filter(cohort == pandemic_col, name == question) %>%
    .$n
  denom <- exp + unexp
  
  # Perform prop.trend.test
  prop.trend.test(exp, denom)
}

questions <- paste0("Q", 1:9)
results <- map(questions, ~ perform_test(eatdis_2, .x))


# Multiple imputation


# For continuous variables: "norm", for categorical: "logreg" or "polyreg" etc.
methodList <- c(
  
  Q1 = "polr", Q2 = "polr", Q3 = "polr", Q4 = "polr", Q5 = "polr",
  Q6 = "polr", Q7 = "polr", Q8 = "polr", Q9 = "polr",
  
  cohort = "",
  
  exercise = "pmm",
  
  
  socialmedia = "pmm", screen = "pmm", stress = "pmm", 
  
  
  scl_score = "norm", conflict_score = "pmm",anx.pgs.pc = "norm", ocd.pgs.pc = "norm", bmi.pgs.pc = "norm",
  mdd.pgs.pc = "norm", neurot2018.pgs.pc = "norm", an2019.pgs.pc = "norm", asd.pgs.pc = "norm",
  
  ED_M = "", BMI_overvekt = "logreg", BMI_undervekt = "logreg",
  
  YA = "", YA_2 = "",AGE_YRS_UB = "")

# Here the datasets are imputed. This can take a while
imputedData <- mice(dataset, method = methodList, m = 50)

# save dataset

saveRDS(imputedData ,file="./data/imputedData.RDS")

# load dataset
imputedData <- readRDS("N:/durable/projects/cme-eatdis/data/imputedData.RDS")

# Apply your analysis model to each imputed dataset

twofactor_SEM <-"f1 =~ Q1 + Q2 + Q3 + Q4 + Q5
                 f2 =~ Q6 + Q7 + Q8 + Q9

                 f1 + f2 ~ exercise + socialmedia + screen +stress + scl_score + conflict_score + anx.pgs.pc + ocd.pgs.pc + bmi.pgs.pc + 
                 mdd.pgs.pc + ED_M + neurot2018.pgs.pc + an2019.pgs.pc + asd.pgs.pc + BMI_overvekt + BMI_undervekt + YA + AGE_YRS_UB"


# pool the results using runMI from the semTools package
pooled_results_G_6 <- runMI(model = twofactor_SEM, data = imputedData, miPackage = "mice", 
                            fun = "sem", estimator = "WLSMV", ordered = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9"))

summary(pooled_results,fit.measures = TRUE, test = "D2",pool.robust = TRUE)

fitMeasures(pooled_results)

lavTestLRT.mi(pooled_results, fit.measures = TRUE, test = "D2", pool.robust = TRUE)

parameterEstimates(pooled_results)

# Multigroup SEM: with pandemic groups 

multigroup_SEM <-"f1 =~ Q1 + Q2 + Q3 + Q4 + Q5
                 f2 =~ Q6 + Q7 + Q8 + Q9

                 f1 ~ exercise + socialmedia + screen + stress + scl_score + conflict_score + anx.pgs.pc + ocd.pgs.pc + bmi.pgs.pc + mdd.pgs.pc +
                 neurot2018.pgs.pc + an2019.pgs.pc + asd.pgs.pc + ED_M + BMI_overvekt + BMI_undervekt + YA + AGE_YRS_UB

                 f2 ~ exercise + socialmedia + screen + stress + scl_score + conflict_score + anx.pgs.pc + ocd.pgs.pc + bmi.pgs.pc + mdd.pgs.pc +
                 neurot2018.pgs.pc + an2019.pgs.pc + asd.pgs.pc + ED_M + BMI_overvekt + BMI_undervekt + YA + AGE_YRS_UB"


pooled_multigroup <- runMI(model = multigroup_SEM, data = imputedData, miPackage = "mice", group = "cohort",
                               fun = "sem", estimator = "WLSMV", ordered = c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9"))


summary(pooled_multigroup,fit.measures = TRUE, test = "D2",pool.robust = TRUE)


