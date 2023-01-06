#*******************************************************************************
#
# Course:  MPhil Population Health Sciences
# Module:  Advanced Biostatistics for Epidemiology
# Date:    06-Jan-2023
# Author:  Robert A Fletcher
# Purpose: Analysis of cohort studies I
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# `|>` is the base R pipe operator


# Install dependencies (if required) --------------------------------------

# Library names
libs <- c("epi", "epitools", "popEpi", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

library(Epi)
library(epiR)
library(epitools)
library(popEpi)
library(tidyverse)


# Exercise 1; Task A ------------------------------------------------------

# Objective: Read data into R

# Define working directory (EDIT THIS)
dir <- 
  "/Users/robertfletcher/Documents/phd/training"

# Define practical directory (DO NOT EDIT THIS)
prac <- 
  "advanced_biostats/02_cohort_studies_I/"

# Read data
bc <- readr::read_csv(glue::glue("{dir}/{prac}/data/cohort_01_bc_rates.csv"))
at <- readr::read_csv(glue::glue("{dir}/{prac}/data/cohort_01_at_relatives.csv"))

# Inspect data
print(bc, n = 20)
print(at, n = 20)


# Exercise 1; Task B ------------------------------------------------------

# Objective: Explore breast cancer data among relatives of AT patients

# Rate of breast cancer between 1970-75 in women aged 60-64 years
bc |> 
  dplyr::filter(yearband == 1970, ageband == 60) |> 
  dplyr::pull(ratesper100000)

# Rate of breast cancer between 1990-95 in women aged 65-69 years
bc |>
  dplyr::filter(yearband == 1990, ageband == 65) |>
  dplyr::pull(ratesper100000)

# Numbers of each relative
at |>
  dplyr::count(relative)

# Numbers of breast cancers diagnosed
at |>
  dplyr::count(brcancer)


# Exercise 1; Task C ------------------------------------------------------

# Objective: Summarise cohort data

# Code date of entry
at <- at |> 
  dplyr::mutate(
    dentry = dplyr::case_when(
      dob > 1970 ~ dob,
      TRUE ~ 1970
    )
  )

# Aggregate data
at_aggr <- 
  popEpi::lexpand(
    at, birth = dob, entry = dentry, exit = dexit, status = brcancer == 1, 
    breaks = list(age = seq(0, 90, 5), per = seq(1970, 2000, 5)),
    aggre = list(age, year = per)
  ) |> 
  tibble::tibble()

# Get number of events
at_aggr |>
  dplyr::summarise(events = sum(from0to1))


# Exercise 1; Task D ------------------------------------------------------

# Objective: Merge data with population rates and calculate standardised 
#            mortality ratio (SMR)

# Rename variables in the breast cancer data
bc2 <- bc |>
  dplyr::rename(year = yearband, age = ageband)

# Calculate SMR
dplyr::left_join(at_aggr, bc2, by = c("year", "age")) |> 
  dplyr::mutate(expected = pyrs * ratesper100000/100000) |> 
  dplyr::summarise(
    smr = sum(from0to1, na.rm = TRUE) / sum(expected, na.rm = TRUE),
    lci = exp(log(smr) - 1.96 * (1 / sqrt(20))),
    uci = exp(log(smr) + 1.96 * (1 / sqrt(20)))
  ) |> 
  dplyr::mutate(
    dplyr::across(tidyr::everything(), ~ round(., digits = 2)),
    res = paste0(smr, " (95% CI: ", lci, " to ", uci, ")")
  )


# Exercise 1; Task E ------------------------------------------------------

# Objective: Calculate the SMR with the `sir()` function

# Recode rates
bc2 <- bc2 |>
  dplyr::mutate(rate = ratesper100000/100000)

# Calculate SMR using `sir()` function
popEpi::sir(
  coh.data = at_aggr, coh.obs = from0to1, coh.pyrs = pyrs, ref.data = bc2, 
  ref.rate = rate, adjust = c("age", "year"), conf.type= "univariate"
)


# Exercise 1; Task F ------------------------------------------------------

# Objective: What are the results (observed, expected, SMR and 95% CI) for aunts 
#            and for grandmothers? 

# Re-aggregate data but now including `relative`
at_aggr_rel <- 
  popEpi::lexpand(
    at, birth = dob, entry = dentry, exit = dexit, status = brcancer == 1, 
    breaks = list(age = seq(0, 90, 5), per = seq(1970, 2000, 5)),
    aggre = list(age, year = per, relative)
  ) |> 
  tibble::tibble()

# Aunts
popEpi::sir(
  coh.data = at_aggr_rel[at_aggr_rel$relative == "Aunt", ], coh.obs = from0to1, 
  coh.pyrs = pyrs, ref.data = bc2 , ref.rate = rate,
  adjust = c("age", "year"), conf.type = "univariate"
)

# Grandmothers
popEpi::sir(
  coh.data = at_aggr_rel[at_aggr_rel$relative == "Grandmother", ], 
  coh.obs = from0to1, coh.pyrs = pyrs, ref.data = bc2 , ref.rate = rate,
  adjust = c("age", "year"), conf.type = "univariate"
)

# Power calculation to detect SMR = 2
alpha = 0.05
-qnorm(alpha)
power = 0.80
beta = 1 - power
-qnorm(beta)
e <- round(((-qnorm(alpha) + -qnorm(beta))^2) / (4 * (sqrt(2) - 1)^2), 2)

# E should be greater than
paste0("E should be greater than ", e)


# Exercise 2; Task A ------------------------------------------------------

# Objective: Create data

# Create new data frame
df <- 
  tibble::tibble(
    deaths = c(47, 185, 117, 117, 350, 192),
    py = c(3896, 17801, 19510, 6559, 22671, 27349),
    smoking = c(1, 2, 3, 1, 2, 3),
    soc = c(1, 1, 1, 0, 0, 0)
  )


# Exercise 2; Task B ------------------------------------------------------

# Objective: Use classical methods to compute the incidence rate ratio (IRR)

# Order data by smoking and exposed/unexposed SOC
df <- df |>
  dplyr::arrange(smoking, dplyr::desc(soc))

# Get number of deaths and total person years for exposed and unexposed 
df2 <- df |> 
  dplyr::group_by(soc) |> 
  dplyr::summarise(dplyr::across(tidyr::everything(), sum)) |> 
  dplyr::arrange(dplyr::desc(soc))

# Compute the rate ratio and 95% CI
results <- 
  epiR::epi.2by2(cbind(df2$deaths, df2$py), method = "cohort.time")
irr <- 
  tibble::tibble(summary(results)$massoc.summary) |>
  dplyr::filter(var == "Inc rate ratio")
p_val <- 
  summary(results)$massoc.detail$chi2.strata.uncor$p.value.2s

# Create a 2x2 contingency table for each smoking stratum
tab <- df |>
  dplyr::group_split(smoking) |>
  purrr::map( ~ cbind(.$deaths, .$py)) |>
  unlist() |> 
  array(dim = c(2, 2, 3)) 

# Compute the rate ratio and 95% CI for each smoking stratum and 
# Mantel–Haenszel-adjusted rate ratio
mh <- tab |> 
  epiR::epi.2by2(method = "cohort.time") |> 
  summary()

names <- tibble::tibble(stratum = c("1", "2", "3", "MH"))

stratum_irr <- tibble::tibble(mh$massoc.detail$IRR.strata.wald)
mh_irr <- tibble::tibble(mh$massoc.detail$IRR.mh.wald)
comb_irr <- bind_rows(stratum_irr, mh_irr)

stratum_pval <- tibble::tibble(p = mh$massoc.detail$chi2.strata.uncor$p.value.2s)
mh_pval <- tibble::tibble(p = mh$massoc.detail$chi2.mh$p.value.2s)
comb_pval <- bind_rows(stratum_pval, mh_pval)

bind_cols(names, comb_irr, comb_pval)


# Exercise 3; Task A ------------------------------------------------------

# Load data
ht <- readr::read_csv(glue::glue("{wd}/data/cohort_01_heart.csv"))

# Format the data into five 2 by 2 tables for each age group
ht_tab <- ht |>
  dplyr::group_split(age) |>
  purrr::map( ~ as.matrix(dplyr::select(., deaths, pyears))) |> 
  unlist() |> 
  array(dim = c(2, 2, 5)) 

# Compute the rate ratio and 95% CI for each age stratum and 
# Mantel–Haenszel-adjusted rate ratio

mh_ht <- ht_tab |> 
  epiR::epi.2by2(method = "cohort.time") |> 
  summary()

names_ht <- tibble::tibble(stratum = c("1", "2", "3", "4", "5", "MH"))

stratum_irr_ht <- tibble::tibble(mh_ht$massoc.detail$IRR.strata.wald)
mh_irr_ht <- tibble::tibble(mh_ht$massoc.detail$IRR.mh.wald)
comb_irr_ht <- bind_rows(stratum_irr_ht, mh_irr_ht)

stratum_pval_ht <- tibble::tibble(p = mh_ht$massoc.detail$chi2.strata.uncor$p.value.2s)
mh_pval_ht <- tibble::tibble(p = mh_ht$massoc.detail$chi2.mh$p.value.2s)
comb_pval_ht <- bind_rows(stratum_pval_ht, mh_pval_ht)

bind_cols(names_ht, comb_irr_ht, comb_pval_ht)
