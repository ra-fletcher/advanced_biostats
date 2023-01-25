#*******************************************************************************
#
# Course:  MPhil Population Health Sciences
# Module:  Advanced Biostatistics for Epidemiology
# Date:    12-Jan-2023
# Authors: Angela Wood & Rob Fletcher
# Contact: amw79@medschl.cam.ac.uk | raf69@medschl.cam.ac.uk
# Purpose: Missing Data and Multiple Imputation
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# This is the code for Practical Session 10: Missing Data and Multiple 
# Imputation

# `|>` is the base R pipe operator. It will be used throughout this code instead
# of the pipe operator imported with the magrittr library `%>%`


# Install dependencies (if required) --------------------------------------

# Library names
libs <- c("jtools", "mice", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

library(jtools)
library(mice)
library(tidyverse)


# Exercise 1; Task A ------------------------------------------------------

# Objective: Read data into R and explore pattern of missingness

# Define file path to practical directory (EDIT THIS)
dir <- "/Users/robertfletcher/Documents/phd/training"

# Define practical directory (DO NOT EDIT THIS)
prac <- "advanced_biostats/10_missing_data_multiple_imputation"

# Read data
kar <- readr::read_csv(glue::glue("{dir}/{prac}/data/karonga_2022.csv"))

# Recode variables and define categorical variables as factors
kar2 <- kar %>% 
  dplyr::mutate(
    d = case_when(
      d == "Control" ~ 0,
      d == "Case" ~ 1,
      TRUE ~ NA_real_
    ),
    age = factor(
      age, 
      levels = c("5-9", "10-14", "15-19", "20-24", "25-29", "30-44", "45+"), 
      labels = c("5-9", "10-14", "15-19", "20-24", "25-29", "30-44", "45+")
    ),
    sex = factor(
      sex, levels = c("Male", "Female"), labels = c("Male", "Female")
    ),
    bcg = factor(
      bcg, levels = c("Absent", "Present"), labels = c("Absent", "Present")
    ),
    house = factor(
      house, 
      levels = c(
        "wattle and daub", "temporary shelter", 
        "sun-dried bricks or pounded mud", "burnt brick"
      ),
      labels = c(
        "wattle and daub", "temporary shelter", 
        "sun-dried bricks or pounded mud", "burnt brick"
      )
    ),
    school = factor(
      school,
      levels = c(
       "none", "1-5yr primary", "6-8yr primary", "secondary/tertiary"
      ),
      labels = c(
        "none", "1-5yr primary", "6-8yr primary", "secondary/tertiary"
      )
    )
  )

# Inspect data
print(kar2, n = 10)

# Explore extent of missingness in each variable (simple)
kar |> 
  purrr::map_df(\(x) sum(is.na(x)))

# Explore extent of missingness in each variable (a little more complex)
kar |> 
  purrr::map_df(\(x) sum(is.na(x))) |> 
  tidyr::pivot_longer(
    dplyr::everything(), names_to = "variable", values_to = "missing"
  ) |> 
  dplyr::mutate(percentage = round(missing / nrow(kar) * 100, digits = 1))

# Create indicator variables for variables with missing data
kar2 <- kar2 |> 
  dplyr::mutate(
    dplyr::across(
      c(school, house), \(x) dplyr::if_else(is.na(x), 1, 0),
      .names = "{.col}_na"
    )
  )

# Explore association of missingness with other variables

# Univariable approach
kar2 |>
  dplyr::select(d, age, sex, bcg, house) |> 
  purrr::map(
    \(x) glm(school_na ~ x, data = kar2) |> 
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)
  )

kar2 |>
  dplyr::select(d, age, sex, bcg, school) |> 
  purrr::map(
    \(x) glm(house_na ~ x, data = kar2) |> 
      broom::tidy(exponentiate = TRUE, conf.int = TRUE)
  )

# Multivariable approach
school_assoc <- 
  glm(
    school_na ~ d + age + sex + bcg + house, data = kar2, family = "binomial"
  ) |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)

house_assoc <-   
  glm(
    house_na ~ d + age + sex + bcg + school, data = kar2, family = "binomial"
  ) |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)


# Exercise 1; Task B ------------------------------------------------------

# Objective: Perform complete case analysis

# Fit logistic model of `bcg` against `d`
glm(d ~ bcg, data = kar2, family = "binomial") |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)

# Fit logistic model of `bcg` against `d` adjusting for other variables
glm(d ~ bcg + age + sex + house + school, data = kar2, family = "binomial") |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)

# Fit logistic model of `bcg` against `d` removing individuals with any missing 
# data
glm(
  d ~ bcg + age + sex + house + school, data = drop_na(kar2), 
  family = "binomial"
) |> 
  jtools::summ(confint = TRUE, digits = 3, exp = TRUE)


# Exercise 1; Task C ------------------------------------------------------

# Objective: Conduct basic multiple imputation

# Multiple imputation using the `mice` package
imputed <- 
  mice::mice(
    kar2[,2:7], method = c("", "", "", "", "polr", "polyreg"), 
    m = 5, seed = 2022
  )
summary(imputed)
attributes(imputed)


# Exercise 1; Task D ------------------------------------------------------

# Objective: Analyse multiply imputed data 

# Run logistic regression across five datasets
fit <- 
  mice::with(
    data = imputed, 
    exp = glm(d ~ age + sex + bcg + house + school, family = "binomial")
  )

# Print output
model_fit

# Get coefficients for `bcg`
coef <- 
  model_fit[["analyses"]] |> 
  purrr::map(\(x) broom::tidy(x) |> dplyr::filter(str_detect(term, "bcg"))) |> 
  purrr::map_dbl(\(x) dplyr::pull(x, estimate))

sum(coef)/5
exp(sum(coef)/5)

# Pool model results using Rubin’s rules 
mice::pool(model_fit)
