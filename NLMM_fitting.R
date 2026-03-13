#!/usr/bin/env RScript

#command-line executable:  chmod +x NLMM_fitting.R
# run as: ./NLMM_fitting.R

# Script for fitting the growth curves for knockout and wildtype to a Non-Linear Mixed effects Model (NLMM) and extracting the mean log2 fold-change model estimates in growth rate and maximum capacity.

TodaysDate <- ""
AnalysisName <- "Systematic"

################################################## Loading Necessary Libraries and Source Code #########################################################################
source("./R/analysis_group_functions.R")
source("./R/nlmm_functions.R")

library(nlme)
library(tidyr)
library(emmeans)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(magrittr)
library("RMySQL")
library(plyr)
library(dplyr)
library(cowplot)
library(Hmisc)
library(gdata)
library(patchwork)
library(gsubfn)
library(gridExtra)
library(tidyverse)
library(qpcR)
library(zoo)
library(data.table)
library(rlang)
library('dbx')
library(openxlsx)
library(mgsub)

############################################################# Connecting to Database and Loading Desired Tables #############################################################

#upload Incucyte confluency data from MySQL database:
mysqlconnection = dbConnect(RMySQL::MySQL(),
                            dbname='',
                            host='',
                            user='',
                            password='')

#list out existing tables in my database
dbListTables(mysqlconnection)

#write query to access confluency data from table
IncucyteConfluencyTableContents = dbSendQuery(mysqlconnection, "select * from Incucyte_Confluency_Dataset")

#save entire dataset
IncucyteDataDF=fetch(IncucyteConfluencyTableContents, n = -1)

#write query to access meta data from table
IncucyteMetaDataTableContents = dbSendQuery(mysqlconnection, "select * from Incucyte_Meta_Data")

#save entire metadata
IncucyteMetaDataDF=fetch(IncucyteMetaDataTableContents, n = -1)

#write query to access condition IDs from table
IncucyteConditionTableContents = dbSendQuery(mysqlconnection, "select * from Incucyte_Condition_Table")

#save entire condition table
IncucyteConditionDF=fetch(IncucyteConditionTableContents, n = -1)

############################################################# Data Wrangling #############################################################

# labelling of bio and technical replicates by combining meta data with raw confluency data
IncucyteDataAndMetaData <- IncucyteDataDF %>% left_join(IncucyteMetaDataDF, by = c("project_ID","plate_ID","well_ID"))

# join condition df to previous dataframe in order to include cell line, treatments, etc info
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaData %>% left_join(IncucyteConditionDF, by = c("condition_ID"))

# remove outer wells/outliers
IncucyteDataAndMetaDataAndConditionsDF= subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$exclude==0)

# generate bio_rep column:join condition_ID, project_ID and plate_ID
IncucyteDataAndMetaDataAndConditionsDF$bio_rep_id<- paste(IncucyteDataAndMetaDataAndConditionsDF$condition_ID, IncucyteDataAndMetaDataAndConditionsDF$project_ID, IncucyteDataAndMetaDataAndConditionsDF$plate_ID, sep = '|')

# average across technical replicates per biological replicate
IncucyteDataAndMetaDataAndConditionsDF<- IncucyteDataAndMetaDataAndConditionsDF %>%
  group_by(bio_rep_id,timepoint,condition_ID,plate_ID, cell_line, cell_line_modifications, cell_number, treatment_0, concentration_0,treatment_1, concentration_1) %>%
  dplyr::summarise(mean_conf=mean(confluency))

# rename cell lines
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>%
  mutate(cell_line_modifications = case_when(
    cell_line_modifications=="C631"  ~ "WT",
    TRUE ~cell_line_modifications))

# replace NA in treatment_0 and concentration_0 with Untreated and none, respectively
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF%>%
  mutate(treatment_0 = case_when(is.na(treatment_0) & is.na(treatment_1)~ "Untreated",
                                 is.na(treatment_0) & treatment_1=="BRM014"~ "BRM014",
                                 TRUE ~ treatment_0)) %>%
  mutate(concentration_0 = case_when(is.na(concentration_0) & is.na(concentration_1) ~ "none",
                                     is.na(concentration_0) & !(is.na(concentration_1)) ~ concentration_1,
                                     TRUE ~ concentration_0))


# CL_treat_conc_id: generate concatenated string of cell line modification, treatment, concentration and cell number????
IncucyteDataAndMetaDataAndConditionsDF$CL_treat_conc_id<- paste(IncucyteDataAndMetaDataAndConditionsDF$cell_line_modifications, IncucyteDataAndMetaDataAndConditionsDF$treatment_0,
                                                                IncucyteDataAndMetaDataAndConditionsDF$concentration_0, IncucyteDataAndMetaDataAndConditionsDF$cell_number, sep = '|')

# filter for conditions with more than 1 biorep
IncucyteDataAndMetaDataAndConditionsDF1 <- IncucyteDataAndMetaDataAndConditionsDF %>%
  group_by(condition_ID) %>%
  dplyr::distinct(bio_rep_id) %>%
  dplyr::filter(n()>1)

IncucyteDataAndMetaDataAndConditionsDF <- subset(IncucyteDataAndMetaDataAndConditionsDF, bio_rep_id %in% IncucyteDataAndMetaDataAndConditionsDF1$bio_rep_id)


############################################################# Assign Analysis ID to each Analysis Group #############################################################

# rename cell line
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>%
  mutate(cell_line_modifications = case_when(
    cell_line_modifications=="WT" & treatment_0=="BRM014"  ~ "WT-",
    TRUE ~cell_line_modifications)) %>%
  mutate(treatment_0 = case_when(treatment_0=="BRM014"~ "Untreated",
                                 TRUE ~ treatment_0)) %>%
  mutate(treatment_1 = case_when(treatment_1=="BRM014"~ NA,
                                 TRUE ~ treatment_1)) %>%
  mutate(concentration_0 = case_when(!(is.na(concentration_0)) ~ "none",
                                     TRUE ~ concentration_0)) %>%
  mutate(concentration_1 = case_when(!(is.na(concentration_1)) ~ NA,
                                     TRUE ~ concentration_1))


# subset for untreated conditions with starting cell number of 4500
IncucyteDataAndMetaDataAndConditionsDFFClones <- subset(IncucyteDataAndMetaDataAndConditionsDF,treatment_0 %in% c("BRM014","Untreated") & cell_number==4500 & !(cell_line_modifications !="WT-"&treatment_0=="BRM014"&cell_number==4500))

# create analysis ID per analysis comparison group (KO & WT)
IncucyteDataAndMetaDataAndConditionsLabeledDF <- getSamePlatesPerUntreatedComparisonGroup(IncucyteDataAndMetaDataAndConditionsDFFClones,wildtype = "WT",AnalysisName)

# remove columns pertaining to second treatment information (as it's not applicable) - treatment_1 & concentration_1
IncucyteDataAndMetaDataAndConditionsLabeledDF <- IncucyteDataAndMetaDataAndConditionsLabeledDF[,-c(9:10)]

############################################################# Running NLMM & extracting Model Estimates & p-values #############################################################

log_growth_pars <- function(pars) {

  # Transforms the parameters Asym, xmid, scal of the logistic function SSlogis to parameters log2K, y0, log2r of the logistic growth function SSlogisGrowth

  setNames(c(log2(pars[1]), pars[1]/(1+exp(pars[2]/pars[3])), log2(1/pars[3])),

           nm = c("l2K", "y0", "l2r"))
}

# function defining initial parameters - in log scale
initlogisGrowth <- function(mCall, data, LHS, ...) {

  xy <- sortedXyData(mCall[["x"]], LHS, data)

  pars <- getInitial(y ~ SSlogis(x, Asym, xmid, scal), data = xy)

  return(log_growth_pars(pars))

}

# self-start function for model - in log scale
logSSlogisGrowth <- selfStart(model = ~ (2^l2K)*y0*exp((2^l2r)*x) / ((2^l2K) + y0*(exp((2^l2r)*x) - 1)),

                              initial = initlogisGrowth,

                              parameter = c("l2K", "y0", "l2r"))


ListOfNLMMFitFailed <- list()

ListOfFailedParameterExtract <- list()

ListOfFailedPlots <- list()

ListOfNLMMAnalysisID = unique(IncucyteDataAndMetaDataAndConditionsLabeledDF$NLMM_Analysis_ID)

for(k in 1:length(ListOfNLMMAnalysisID)){

  IncucyteDatasetSubset <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, IncucyteDataAndMetaDataAndConditionsLabeledDF$NLMM_Analysis_ID==ListOfNLMMAnalysisID[k])#k

  # fit conditions (KO & WT) for one given analysis group to NLMM
  if(isTRUE(k==1)==TRUE){
    ListOfNLMMtTestResults <- runNLMM(IncucyteDatasetSubset,treatmentendtime = 120,ListOfControlCL = "WT",ListOfNLMMFitFailed,ListOfFailedParameterExtract,ListOfFailedPlots,TodaysDate, AnalysisName)
  } else{ListOfNLMMtTestResults <- runNLMM(IncucyteDatasetSubset,treatmentendtime = 120,ListOfControlCL = "WT",ListOfNLMMtTestResults$ListOfNLMMFitFailed,ListOfNLMMtTestResults$ListOfFailedParameterExtract,
                                                ListOfNLMMtTestResults$ListOfFailedPlots,TodaysDate, AnalysisName)
  }

  #export NLMM significance statistics (p-values) to MySQL table
  dbWriteTable(conn = mysqlconnection, name = 'Incucyte_NLMM_Untreat_Interaction_Table', value = ListOfNLMMtTestResults$NLMMEstimatesTable, append = TRUE, header = FALSE, row.names = FALSE)

  SubDir <- sprintf("./Modeling-dependent_Phenotypes/NLMM_Fitted_Curves/%s/",AnalysisName)

  # check if directory exists - create if it doesn't
  if(isTRUE(file.exists(SubDir))==TRUE){
    print(paste("Plotting an NLMM model fit for", unique(IncucyteDatasetSubset$NLMM_Analysis_ID)))
  } else {
    dir.create(SubDir, recursive = TRUE)
    print(paste("Plotting an NLMM model fit for", unique(IncucyteDatasetSubset$NLMM_Analysis_ID)))
  }

  IncucyteDatasetSubsetNoWT <- subset(IncucyteDatasetSubset, !(cell_line_modifications%in%ListOfControlCL))

  cell_line <- unique(IncucyteDatasetSubsetNoWT$cell_line_modifications)

  treat_type <- "Untreated"

  # save model fit plots
  new_filename = sprintf("./Modeling-dependent_Phenotypes/NLMM_Fitted_Curves/%s/%s_%s_%s_%s_%s.pdf", AnalysisName, TodaysDate,unique(IncucyteDatasetSubset$NLMM_Analysis_ID),
                         treat_type,unique(IncucyteDatasetSubset$cell_number),cell_line)
  ggsave(filename=new_filename, plot=ListOfNLMMtTestResults$ModelFitPlot,width = 25, height = 15, dpi = 150, units = "in", device='pdf')
}

print("ListOfNLMMFitFailed t-test: \n")
print(ListOfNLMMtTestResults$ListOfNLMMFitFailed)

print("ListOfFailedParameterExtract t-test: \n")
print(ListOfNLMMtTestResults$ListOfFailedParameterExtract)

print("ListOfFailedPlots t-test: \n")
print(ListOfNLMMtTestResults$ListOfFailedPlots)

