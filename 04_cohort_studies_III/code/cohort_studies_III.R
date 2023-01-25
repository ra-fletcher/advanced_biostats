#*******************************************************************************
#
# Course:  MPhil Population Health Sciences
# Module:  Advanced Biostatistics for Epidemiology
# Date:    06-Jan-2023
# Authors: Angela Wood & Rob Fletcher
# Contact: amw79@medschl.cam.ac.uk | raf69@medschl.cam.ac.uk
# Purpose: Analysis of Cohort Studies III
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# This is the code for Practical Session 4: Analysis of Cohort Studies III

# `|>` is the base R pipe operator. It will be used throughout this code instead
# of the pipe operator imported with the magrittr library `%>%`


# Install dependencies (if required) --------------------------------------

# Library names
libs <- c("epiR", "jtools", "multcomp", "tidymodels", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

lapply(libs, library, character.only = TRUE)


# Exercise 1; Task A ------------------------------------------------------

# Objective: Read data into R and explore data

# Define file path to practical directory (EDIT THIS)
dir <- "/Users/robertfletcher/Documents/phd/projects"

# Define practical directory (DO NOT EDIT THIS)
prac <- "advanced_biostats/04_cohort_studies_III"

# Read data
epic <- readr::read_csv(glue::glue("{dir}/{prac}/data/cohort_03_epic.csv"))

# Inspect data
print(epic, n = 10)

# Get number of cardiovascular disease deaths in each age group
epic |> 
  group_by(age) |> 
  summarise(events = sum(cvd))


# Exercise 1; Task B ------------------------------------------------------

# Objective: Explore association between hostility and cardiovascular mortality

# Compute unadjusted/crude rate ratio using Poisson regression
fit_crude <- 
  glm(cvd ~ hostility, data = epic, family = poisson, offset = log(pmonths))

summary(fit_crude)

# Get model fit metrics
jtools::summ(fit_crude, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_crude <- fit_crude |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)


# Exercise 1; Task C ------------------------------------------------------

# Objective: Explore associations of hostility and age with cardiovascular 
#            mortality and whether effect modification between these exposures
#            exists

# Create age factor variable (`age_fct`)
epic <- epic |> 
  mutate(
    age_fct = factor(
      age, levels = 1:3, labels = c("40-59 yrs", "60-69 yrs", "70-80 yrs")
    )
  )

# Fit Poisson model adjusted for age groups
fit_age <- 
  glm(
    cvd ~ hostility + age_fct, data = epic, family = poisson, 
    offset = log(pmonths)
  )

summary(fit_age)

# Get model fit metrics
jtools::summ(fit_age, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_age <- fit_age |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)

# Fit Poisson model with interaction term between `age_fct` and `hostility`
fit_ageint <- 
  glm(
    cvd ~ hostility * age_fct, data = epic, family = poisson, 
    offset = log(pmonths)
  )

summary(fit_ageint)

# Get model fit metrics
jtools::summ(fit_ageint, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_ageint <- fit_ageint |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)

# Compute LRT
anova(fit_age, fit_ageint, test = "LRT")


# Exercise 1; Task D ------------------------------------------------------

# Objective: Explore associations of hostility and smoking with cardiovascular 
#            mortality and whether effect modification between these exposures
#            exists

# Create smoking factor variable
epic <- epic |> 
  mutate(
    smoking_fct = factor(
      smoking, levels = 1:3, 
      labels = c("Current smoker", "Former smoker", "Never smoker")
    )
  )

# Fit Poisson model adjusted for smoking groups
fit_smok <-
  glm(
    cvd ~ hostility + smoking_fct, data = epic, family = poisson, 
    offset = log(pmonths)
  )

summary(fit_age)

# Get model fit metrics
jtools::summ(fit_smok, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_smok <- fit_smok |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)

# Fit Poisson model with interaction term between `smoking_fct` and `hostility`
fit_smokint <- 
  glm(
    cvd ~ hostility * smoking_fct, data = epic, family = poisson, 
    offset = log(pmonths)
  )

summary(fit_smokint)

# Get model fit metrics
jtools::summ(fit_smokint, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_smokint <- fit_smokint |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)

# Compute LRT
anova(fit_smok, fit_smokint, test = "LRT")

# Check further if smoking is a confounder for the relationship between
# hostility and cardiovascular mortality

# Check univariable association of smoking with cardiovascular mortality
fit_smok2 <- 
  glm(
    cvd ~ smoking_fct, data = epic, family = poisson, offset = log(pmonths)
  ) |> 
  broom::tidy(confint = TRUE, exponentiate = TRUE)

# Check univariable association of smoking with hostility
table(epic$hostility, epic$smoking)
chisq.test(table(epic$hostility, epic$smoking_fct))

glm(hostility ~ smoking_fct, data = epic, family = binomial) |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)


# Exercise 1; Task E ------------------------------------------------------

# Objective: Explore associations of hostility, age, and smoking

# Fit Poisson model adjusted for age and smoking
fit_agesmok <-
  glm(
    cvd ~ hostility + age_fct + smoking_fct, data = epic, family = poisson, 
    offset = log(pmonths) 
  )
summ(fit_agesmok, confint = TRUE, digits = 3, exp = TRUE)

# Fit Poisson model with interaction term between `hostility` and `age`
fit_ageint_smok <- 
  glm(
    cvd ~ hostility * age_fct + smoking_fct, data = epic, family = poisson, 
    offset = log(pmonths)
  )
summ(fit_ageint_smok, confint = TRUE, digits = 3, exp = TRUE)

# Compute LRT
anova(fit_agesmok, fit_ageint_smok, test = "LRT")
