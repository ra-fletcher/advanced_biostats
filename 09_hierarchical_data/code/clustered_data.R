#*******************************************************************************
#
# Course:  MPhil Population Health Sciences
# Module:  Advanced Biostatistics for Epidemiology
# Date:    11-Jan-2023
# Authors: Angela Wood & Rob Fletcher
# Contact: amw79@medschl.cam.ac.uk | raf69@medschl.cam.ac.uk
# Purpose: Analysis of Hierarchical Data; Clustered Data
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# This is the code for Practical Session 9: Analysis of Hierarchical Data; 
# Clustered Data

# `|>` is the base R pipe operator. It will be used throughout this code instead 
# of the pipe operator imported with the magrittr library `%>%`


# Install dependencies (if required) --------------------------------------

# Library names
libs <- c("jtools", "lmtest", "sandwich", "skimr", "survival", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

library(jtools)
library(lmtest)
library(sandwich)
library(skimr)
library(survival)
library(tidyverse)


# Exercise 1; Task A ------------------------------------------------------

# Objective: Read data into R

# Define file path to practical directory (EDIT THIS)
dir <- "/Users/robertfletcher/Documents/phd/projects"

# Define practical directory (DO NOT EDIT THIS)
prac <- "advanced_biostats/09_hierarchical_data"

# Read data
sga <- readr::read_csv(glue::glue("{dir}/{prac}/data/sga.csv"))

# Inspect data
print(sga, n = 10)
skimr::skim(sga)
dplyr::glimpse(sga)

# Define categorical variables as factors
sga <- sga |> 
  dplyr::mutate(
    depcat7 = as.factor(depcat7),
    hospital = as.factor(hospital),
    parcat3 = factor(
      parcat3, levels = 1:3, 
      labels = c("No previous", "Parity 1 to 2", "Parity >=3")
    ),
    sga_5th = factor(
      sga_5th, levels = 0:1, labels = c("Normal weight", "Low weight")
    ),
    smokcat3 = factor(
      smokcat3, levels = 1:3, labels = c("Never", "Smoker", "Ex-smoker")
    )
  )


# Exercise 1; Task B ------------------------------------------------------

# Objective: Explore data

# Number of unique hospital clusters
sga |> 
  dplyr::summarise(`unique clusters` = length(unique(hospital)))

# Number of deliveries per hospital
sga |> 
  dplyr::count(hospital)

# Number of births and mothers
sga |> 
  dplyr::summarise(births = n(), mothers = length(unique(id_mum)))


# Exercise 1; Task C ------------------------------------------------------

# Objective: Use logistic regression to explore association between smoking 
#.           status and SGA

# Fit simple logistic regression model
fit1 <-
  glm(sga_5th ~ smokcat3, data = sga, family = binomial) 

fit1 |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)

# Fit multiple logistic regression model to assess confounding
fit2 <- 
  glm(
    sga_5th ~ smokcat3 + height + weight + age + parcat3 + depcat7,
    data = sga, family = binomial
  ) 

fit2 |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)


# Exercise 1; Task D ------------------------------------------------------

# Objective: Obtain robust standard errors after clustering on id_mum

# Robust standard errors for model 1
lmtest::coeftest(fit1, vcov = sandwich::vcovCL, cluster = sga$id_mum)
summary(fit1)

# Robust standard errors for model 2
lmtest::coeftest(fit2, vcov = sandwich::vcovCL, cluster = sga$id_mum)
summary(fit2)


# Exercise 1; Task E ------------------------------------------------------

# Objective: Obtain robust standard errors after clustering on hospital

# Robust standard errors for model 1
lmtest::coeftest(fit1, vcov = sandwich::vcovCL, cluster = sga$hospital)

# Robust standard errors for model 2
lmtest::coeftest(fit2, vcov = sandwich::vcovCL, cluster = sga$hospital)

fit3 <- 
  glm(
    sga_5th ~ smokcat3 + height + weight + age + parcat3 + depcat7 + hospital, 
    data = sga, family = binomial
  ) |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)


# Exercise 2; Task A ------------------------------------------------------

# Objective: Read diabetic retinopathy data and convert to long format

# Read data
drs <- readr::read_csv(glue::glue("{dir}/{prac}/data/drs.csv"))

# Inspect data
skimr::skim(drs)
dplyr::glimpse(drs)

# Convert data to long format (one row per eye)
drs_l <- drs |> 
  tidyr::pivot_longer(
    -c(id, age_dx), names_sep = "\\.", names_to = c("param", "eye")
  ) |> 
  tidyr::pivot_wider(
    id_cols = c(id, age_dx, eye), names_from = "param", values_from = "value"
  )


# Exercise 2; Task B ------------------------------------------------------

# Objective: Produce Kaplan-Meier curves for each treatment arm

# Create survival object
eye_surv <- survival::Surv(drs_l$time, drs_l$blind)

# Create Kaplan-Meier object
km <- survival::survfit(eye_surv ~ trt, data = drs_l)

# Plot Kaplan-Meier curve
plot(
  km, xlab = "months", ylab = "Probability of survival (%)",
  main = "Kaplan-Meier survivor function", mark.time = TRUE, yscale = 100
)


# Exercise 2; Task C ------------------------------------------------------

# Objective: Fit unadjusted and adjusted Cox models assuming all the data are
#            independent

# Unadjusted
fit_un <- survival::coxph(eye_surv ~ trt, data = drs_l) |> 
  summary(exp = TRUE)

# Adjusted
fit_ad <- survival::coxph(eye_surv ~ trt + eye + age_dx, data = drs_l) |> 
  summary(exp=TRUE)


# Exercise 3; Task D ------------------------------------------------------

# Objective: Fit unadjusted and adjusted Cox models with robust standard errors 
#            after clustering on `id`

# Unadjusted
fit_un2 <- 
  survival::coxph(eye_surv ~ trt, cluster = id, data = drs_l) |> 
  summary(exp = TRUE)

# Adjusted
fit_ad2 <- 
  survival::coxph(eye_surv ~ trt + eye + age_dx, cluster = id, data = drs_l) |> 
  summary(exp = TRUE)
