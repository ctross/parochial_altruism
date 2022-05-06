############################################################# Set your directory
setwd("C:\\Users\\Mind Is Moving\\Dropbox\\Open Papers\\Parochial altruism\\BBS\\Workflow\\")

################################################################# Load Libraries
library(igraph)
library(reshape2)
library(plyr)
library(kinship2)
library(geosphere)
library(GGally)
library(network)
library(ggplot2)
library(rethinking)
library(colorspace)
library(parallel)
library(Cairo)
library(qgraph)
options(mc.cores = parallel::detectCores())

################################################################# Build Database
load("Data\\ColombianDataWithImputations_Site1.RData")       # Loads anonymized and rescaled
load("Data\\ColombianDataWithImputations_Site2.RData")       # data with hard-coded
                                                             # imputations of missings

################################################################# Fit models
source("Code\\fit_models.R")
source("Code\\plot_results.R")

################################################################# Other estimates
source("Code\\print_descriptives.R")
source("Code\\plot_scatters.R")
source("Code\\plot_networks.R")





