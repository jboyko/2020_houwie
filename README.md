# 2020_houwie
joint modeling of discrete and continuous traits

hOUwieNode.R - all of the main functions for hOUwie necessary to fit the likelihoods and run ancestral state reconstructions. these functions were mainly moved to the OUwie R package and are no longer maintained. 

Utils.R - many utility functions for simulating datasets. most of these functions are obsolete and were implemented formally into hOUwieNode.R

01_sim-hOUwie-data.R - simulates the data as described in the manuscript. 

02_fit-hOUwie-model.R - fits the hOUwie model to all of the data sets constructed in 01.

03_analyze-hOUwie-model.R - summarizes the model fits from 02 and constructs some of the figures in the manuscript.

04_empirical-seed-dispersal.R - fits the set of hOUwie models to our seed dispersal dataset as described in Vasoncelos et al. 2021.

05_forward-simulation-figure.R - generates forward simulations to show the different behaviours of the various houwie models. 

06_compare-simmap-generation.R - for testing our different conditional sampling approaches and demonstrating how adaptive sampling improves things. 

07_empirical-analysis.R - summarizes the results of our empirical test case in 04.

biggest_table.Rsave - the complete datafiles are too large to upload to GitHub so this summary table summarizes the resulting parameter estimates. necessary for script 03.

checks/ - a folder containing scripts which test some of the properties of the hOUwie model and ensure it is working correctly.

doc/ - a folder containing various stages of the writing of the manuscript.

empirical_fit/ - a folder containing all of the fit houwie models to the empirical dataset. unfortunately, the bootstrap result is too large to upload to GitHub, so contact me directly for that (jboyko [at] uark.edu). 

old_code/ - various old code that is no longer needed for the manuscript.

tables/ - tables that appear in the manuscript.

lit-review-houwie.xlsx - an excel worksheet containing brief summaries of similar models 

doc/Appendix_note-on-path-probabilities.pdf - a note describing the calculation differences between what is used in hOUwie and what is described by May and Moore (2020).

