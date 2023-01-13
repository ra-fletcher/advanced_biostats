#*******************************************************************************
#
# Course:  MPhil Population Health Sciences
# Module:  Advanced Biostatistics for Epidemiology
# Date:    06-Jan-2023
# Author:  Rob Fletcher
# Contact: raf69@medschl.cam.ac.uk
# Purpose: Analysis of Cohort Studies II
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# This is the code for Practical Session 3: Analysis of Cohort Studies II

# `|>` is the base R pipe operator. It will be used throughout this code instead
# of the pipe operator imported with the magrittr library `%>%`


# Install dependencies (if required) --------------------------------------

# Library names
libs <- c("epiR", "jtools", "multcomp", "tidymodels", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

library(epiR)
library(jtools)
library(multcomp)
library(tidymodels)
library(tidyverse)


# Exercise 1; Task A ------------------------------------------------------

# Objective: Create data

# Create new data frame
df <- 
  tibble::tibble(
    soc = seq(0, 6, by = 1),
    deaths = c(77, 272, 371, 192, 77, 16, 3),
    pyears = c(10911, 30296, 30416, 18463, 6155, 1321, 224)
  )

# Inspect data
print(df)


# Exercise 1; Task B ------------------------------------------------------

# Objective: Fit Poisson model with binary exposure variable

# Create binary exposure variable
df <- df |> 
  dplyr::mutate(soc2 = dplyr::if_else(soc <= 1, 1, 0))

# Aggregate data
aggr <- df |> 
  dplyr::group_by(soc2) |> 
  dplyr::summarise(dplyr::across(c(deaths, pyears), sum)) |> 
  dplyr::arrange(dplyr::desc(soc2))

# Compute rate ratio by classical method
results <- aggr |> 
  dplyr::select(deaths, pyears) |> 
  as.matrix() |> 
  epiR::epi.2by2(method = "cohort.time")

irr <- 
  tibble::tibble(summary(results)$massoc.summary) |>
  dplyr::filter(var == "Inc rate ratio")
p_val <- 
  summary(results)$massoc.detail$chi2.strata.uncor$p.value.2s

# Compute rate ratio using Poisson regression
fit <-
  glm(deaths ~ soc2, data = df, family = poisson, offset = log(pyears)) |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)


# Exercise 1; Task C ------------------------------------------------------

# Objective: Consider exposure variable with more than 2 levels

# Create new categorical exposure variable
df <- df |> 
  dplyr::mutate(soc4 = as.factor(c(0, 1, 2, 3, 3, 3, 3)))

# Aggregate data
aggr2 <- df |> 
  dplyr::group_by(soc4) |> 
  dplyr::summarise(dplyr::across(c(deaths, pyears), sum))

# Fit Poisson regression with categorical `soc4` variable
fit_cat <- 
  glm(deaths ~ soc4, data = aggr2, family = poisson, offset = log(pyears)) 

summary(fit_cat)

# Get model fit metrics
jtools::summ(fit_cat, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_cat <- fit_cat |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)

# Use intercept estimate to get incidence rate (per 1000 person years) for 
# reference group (`soc4` == 0) 
ref_rate <- coef_cat |> 
  dplyr::filter(term == "(Intercept)") |> 
  dplyr::pull(estimate) |> 
  purrr::map_dbl( ~ . * 1000)

# Use intercept estimate + `soc4` == 3 estimate to get incidence rate for 
# `soc4` == 3 
lvl3_rate <- coef_cat |> 
  dplyr::filter(term == "soc43") |> 
  dplyr::pull(estimate) |> 
  purrr::map_dbl( ~ . * ref_rate)

# To test for linear trend, use soc variable as continuous in a new model
# and compare the two models using likelihood ratio test (LRT)

# Create continuous `soc4` variable
aggr2 <- aggr2 |> 
  dplyr::mutate(soc4_c = as.numeric(as.character(soc4)))

# Fit Poisson with continuous `soc4` variable
fit_lin <- 
  glm(deaths ~ soc4_c, data = aggr2, family = poisson, offset = log(pyears))

summary(fit_lin)

# Get model fit metrics
jtools::summ(fit_lin, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_lin <- fit_lin |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)

# Compute LRT
anova(fit_lin, fit_cat, test = "LRT")


# Exercise 1; Task D ------------------------------------------------------

# Fit a Poisson model with a continuous exposure variable

# Use original `soc` variable with 7 levels as continuous exposure variable
fit_cont <- 
  glm(deaths ~ soc, data = df, family = poisson, offset = log(pyears)) 

summary(fit_cont)

# Get model fit metrics
jtools::summ(fit_cont, confint = TRUE, digits = 5, exp = TRUE)

# Get model results as tibble
coef_cont <- fit_cont |> 
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)

# Interpretation of risk per 3 point increase in `soc`
inc_3pt <- 
  confint(multcomp::glht(fit_cont, linfct = c("3 * soc = 0")))

# Convert results to tibble
tibble::as_tibble(inc_3pt$confint)

# Alternative
df <- df |> 
  mutate(soc_div3 = soc / 3)

# Fit Poisson regression model
glm(deaths ~ soc_div3, data = df, family = poisson, offset = log(pyears)) |> 
  # Get model results as tibble
  broom::tidy(conf.int = TRUE, exponentiate = TRUE)
