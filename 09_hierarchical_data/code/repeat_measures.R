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
  group_split(treat) |> 
  purrr::map( ~ tabulate(., resp, time))


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
