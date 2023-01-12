#*******************************************************************************
#
# Course:  MPhil Population Health Sciences
# Module:  Advanced Biostatistics for Epidemiology
# Date:    11-Jan-2023
# Author:  Rob Fletcher
# Contact: raf69@medschl.cam.ac.uk
# Purpose: Analysis of Hierarchical Data; Repeat Measures
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# This is the code for Practical Session X: Analysis of Hierarchical Data

# `|>` is the base R pipe operator. It will be used throughout this code in 
# favour of the pipe operator imported with the magrittr library `%>%`


# Install dependencies (if required) --------------------------------------

# Library names
libs <- c("jtools", "lme4", "skimr", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

library(jtools)
library(lme4)
library(skimr)
library(tidyverse)


# Exercise 1; Task A ------------------------------------------------------

# Objective: Read data into R and recode variables

# Define file path to practical directory (EDIT THIS)
dir <- "/Users/robertfletcher/Documents/phd/training"

# Define practical directory (DO NOT EDIT THIS)
prac <- "advanced_biostats/09_hierarchical_data"

# Read data
skin <- readr::read_csv(glue::glue("{dir}/{prac}/data/skin.csv"))

# Inspect data
print(skin, n = 10)

# Define categoricsal variables as factors
skin <- skin |> 
  dplyr::mutate(
    treat = factor(treat, levels = 0:1, labels = c("Placebo", "Treatment")),
    resp = factor(
      resp, levels = 0:1, labels = c("Excellent/good", "Fair/poor")
    )
  )


# Exercise 1; Task B ------------------------------------------------------

# Objective: Tabulate responses by time and treatment

# Define function to tabulate
tabulate <- function(.data, .var1, .var2) {
  tab <- .data |> 
  dplyr::group_by({{ .var1 }}, {{ .var2 }}) |> 
    dplyr::summarise(n = n(), .groups = "drop") |> 
    tidyr::pivot_wider(
      id_cols = {{ .var1 }}, names_from = {{ .var2 }}, values_from = n
    ) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(total = sum(dplyr::c_across(`3`:`12`))) |> 
    dplyr::ungroup() |> 
    dplyr::add_row(resp = "total") |> 
    dplyr::mutate(
      dplyr::across(
        tidyselect:::where(is.numeric), 
        ~ dplyr::case_when(resp == "total" ~ sum(., na.rm = TRUE), TRUE ~ .)
      ),
      dplyr::across(
        tidyselect:::where(is.numeric), 
        ~ paste0(., " (", round(. / (sum(.) / 2) * 100, digits = 1), "%)")
      )
    )
  return(tab)
}

# Tabulate for the entire data
skin |> 
  tabulate(resp, time)

# Split by treatment
skin |> 
  dplyr::group_split(treat) |> 
  purrr::map(\(x) tabulate(x, resp, time))


# Exercise 1; Task C ------------------------------------------------------

# Objective: Fit mixed-effects logistic regression models

# Fit logistic model to responses, with fixed terms for treatment and time, and
# random intercept term to account for the correlation within individuals
fit <- 
  lme4::glmer(resp ~ treat + time + (1 | id), data = skin, family = "binomial")

jtools::summ(fit, confint = TRUE, digits = 3, exp = TRUE)

# Explore non-linear association with time, and interaction between time and 
# treatment
skin <- skin |> 
  dplyr::mutate(time = as.factor(time))

fit_nl <- 
  lme4::glmer(
    resp ~ treat * time + (1 | id), data = skin, family = "binomial"
  ) |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)

# Add interaction term
fit_int <- 
  lme4::glmer(
    resp ~ treat * time + (1 | id), data = skin, family = "binomial"
  ) |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)


# Extension Exercise 2; Task A --------------------------------------------

# Objective: Read data into R and check for missing data

# Read data
bp <- readr::read_csv(glue::glue("{dir}/{prac}/data/mmbp.csv"))

# Check for missing data
bp |> 
  purrr::map_df(\(x) sum(is.na(x)))

# Define categorical variables as factors
bp <- bp |> 
  dplyr::mutate(
    sex = factor(sex, levels = 1:2, labels = c("male", "female")),
    agegrp = factor(
      agegrp, levels = 1:3, labels = c("75-79 yrs", "80-84 yrs", "85+ yrs")
    ),
    smoking = factor(
      smoking, levels = 0:2, 
      labels = c("never smoker", "ex-smoker", "current smoker")
    )
  )


# Extension Exercise 2; Task B --------------------------------------------

# Objective:	Plot observed responses in systolic blood pressure over time for 
#             the three different age groups

# Set seed to set the random number generator state
set.seed(69)

# Select 10 random individuals from each age group
bp_slice <- bp |> 
  tidyr::drop_na() |> 
  dplyr::group_by(id) |> 
  dplyr::filter(n() > 1) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(agegrp) |> 
  dplyr::select(id) |> 
  dplyr::slice_sample(n = 10) |> 
  dplyr::ungroup() |> 
  dplyr::left_join(bp, by = c("id", "agegrp"))

# Create spaghetti plot
ggplot2::ggplot(bp_slice, ggplot2::aes(x = time, y = sbp, group = id)) +
  ggplot2::geom_point(aes(colour = factor(id))) +
  ggplot2::geom_line(aes(colour = factor(id))) +
  ggplot2::facet_wrap(ggplot2::vars(agegrp)) +
  ggplot2::labs(
    x = "Time from first health check (months)",
    y = "Systolic blood pressure (mm Hg)"
  )


# Extension Exercise 2; Task C --------------------------------------------

# Objective: Centre continuous variables

# Centre `weight`
bp <- bp |> 
  dplyr::mutate(weight_c = weight - mean(weight, na.rm = TRUE))

# Plot histogram of `weight`
ggplot2::ggplot(bp, aes(x = weight, y = ..density..)) + 
  ggplot2::geom_histogram()

# Plot histogram of `weight_c`
ggplot2::ggplot(bp, aes(x = weight_c, y = ..density..)) + 
  ggplot2::geom_histogram()

# Scale `time` to years
bp <- bp |> 
  dplyr::mutate(time_years = time / 12)


# Extension Exercise 2; Task D --------------------------------------------

# Objective: Fit a linear mixed model to `sbp`, with fixed terms for `sex`, 
#            `agegrp`, `smoking`, `weight_c` and `time_years`, and random 
#            intercept and slope terms to account for the correlation within 
#            individuals

# Fit linear mixed model
fit <- 
  lme4::lmer(
    sbp ~ sex + agegrp + smoking + weight_c + time_years + 
    (1 + time_years | id), 
    data = bp
  )
summary(fit)


# Extension Exercise 3; Task E --------------------------------------------

# Objective: Investigate use of quadratic term for time

# Create new time variable
bp <- bp |> 
 dplyr::mutate(time_years2 = time_years * time_years)

# Fit linear mixed model
fit2 <- 
  lme4::lmer(
    sbp ~ sex + agegrp + smoking + weight_c + time_years + time_years2 + 
    ( 1 + time_years | id),
    data = bp
  )
summary(fit)


# Extension Exercise 3; Task F --------------------------------------------

# Objective: Plot the predicted curves for sbp over time 

# Get model predictions
bp <- bp |> 
  drop_na(-dbp) |> 
  mutate(pred = predict(fit2))

# Plot model predictions by age group
ggplot2::ggplot(bp, ggplot2::aes(x = time_years, y = pred, group = agegrp)) + 
  ggplot2::geom_smooth(
    ggplot2::aes(colour = agegrp), method = "lm", formula = y ~ x
  ) +
  ggplot2::labs(
    x = "Time from first health check (years)",
    y = "Predicted systolic blood pressure (mm Hg)"
  )
