#!/usr/bin/env RScript

#command-line executable:  chmod +x assigning_phenotypes.R
# run as: ./assigning_phenotypes.R

# This script processes the confluency data and generates the following phenotype readouts: confluency after 24h, 48h, 72h, 96h and 96h/48h of treatment.

################################################## Loading Necessary Libraries and Source Code #########################################################################
source("./R/splitpoint.R")

library("RMySQL")
library(nlme)
library(emmeans)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(Hmisc)
library(gdata)
library(plyr)
library(dplyr)
library(ggpubr)
library(patchwork)
library(gsubfn)
library(gridExtra)
library(tidyverse)
library(qpcR)
library(zoo)
library(data.table)

############################################################# Connecting to Database and Loading Desired Tables #############################################################

#upload Incucyte confluency data from MySQL database:
mysqlconnection = dbConnect(RMySQL::MySQL(),
                            dbname='',
                            host='',
                            user='',
                            password='')

#list out existing tables in my database
dbListTables(mysqlconnection)

#load confluency dataset
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
#subset dataset via metadata
#labeling of bio and technical replicates by combining meta data with df

IncucyteDataAndMetaData <- IncucyteDataDF %>% left_join(IncucyteMetaDataDF, by = c("project_ID","plate_ID","well_ID"))

#join condition df to previous dataframe in order to include cell line, treatments, etc info
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaData %>% left_join(IncucyteConditionDF, by = c("condition_ID"))

#remove outer wells/outliers
IncucyteDataAndMetaDataAndConditionsDF= subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$exclude==0)

#generate bio_rep column
IncucyteDataAndMetaDataAndConditionsDF$bio_rep<- paste(IncucyteDataAndMetaDataAndConditionsDF$condition_ID, IncucyteDataAndMetaDataAndConditionsDF$project_ID, IncucyteDataAndMetaDataAndConditionsDF$plate_ID, sep = '|')


############################################################# Calculating First Derivative #############################################################
IncucyteDataFirstDerivativeDF <- getFirstDerivative(IncucyteDataAndMetaDataAndConditionsDF)

############################################################# Calculating Second Derivative  #############################################################

IncucyteFirstAndSecondDerivativeWithSplitPointsDF <- getSecondDerivativeAndSplitPoint(IncucyteDataFirstDerivativeDF,treatmentstarttime = 24, treatmentendtime = 120)

############################################################# Create Splitpoint Table and Save to Database #############################################################
IncucyteSplitPointTable <- createSplitPointTable(IncucyteFirstAndSecondDerivativeWithSplitPointsDF,treatmentstarttime, treatmentendtime)

# project_ID
# plate_ID
# condition_ID
# logistic_growth?
# logistic_growth_duration
# second_phase?
# second_phase_duration
# split_timepoint
# split_first_conf_difference
# after_24h_confluency
# after_48h_confluency
# after_72h_confluency
# last_confluency
# conf_96hr_48hr_ratio
# phenotype

#upload IncucyteSplitPointTable to MySQL DB
dbWriteTable(conn = mysqlconnection, name = 'Incucyte_Split_Point_Table', value = IncucyteSplitPointTable, append = TRUE, header = FALSE, row.names = FALSE)
