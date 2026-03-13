#!/usr/bin/env RScript

#command-line executable:  chmod +x growth_curves.R
# run as: ./growth_curves.R

# This script plots the growth curves for each analysis group (knockout treated, knockout untreated, wildtype treated, widltype untreated) across all biological replicates.

TodaysDate <- ""

################################################## Loading Necessary Libraries and Source Code #########################################################################
source("./R/analysis_group_functions.R")
source("./R/plot_functions.R")

library("RMySQL")
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(Hmisc)
library(gdata)
library(dplyr)
library(plyr)
library(ggpubr)
library(patchwork)
library(gsubfn)
library(gridExtra)
library(tidyverse)
library(qpcR)
library(zoo)
library(data.table)
library(rlang)

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

#write query to access phenotype readouts
IncucyteSplitPointTableContents = dbSendQuery(mysqlconnection, "select * from Incucyte_Split_Point_Table")

#save entire phenotype table
IncucyteSplitPointDF=fetch(IncucyteSplitPointTableContents, n = -1)

#get DMSO dilution table
IncucyteDMSODilutionTableContents = dbSendQuery(mysqlconnection, "select * from Incucyte_Compound_DMSO_Dilution_Table")

IncucyteDMSODilutionDF=fetch(IncucyteDMSODilutionTableContents, n = -1)


############################################################# Data Wrangling #############################################################
# #subset dataset via metadata
#labelling of bio and technical replicates by combining meta data with df
#join has to include proj_ID
IncucyteDataAndMetaData <- IncucyteDataDF %>% left_join(IncucyteMetaDataDF, by = c("project_ID","plate_ID","well_ID"))

#join condition df to previous dataframe in order to include cell line, treatments, etc info
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaData %>% left_join(IncucyteConditionDF, by = c("condition_ID"))

#remove outer wells/outliers
IncucyteDataAndMetaDataAndConditionsDF= subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$exclude==0)

#generate bio_rep column:join condition_ID, project_ID and plate_ID
IncucyteDataAndMetaDataAndConditionsDF$bio_rep_id<- paste(IncucyteDataAndMetaDataAndConditionsDF$condition_ID, IncucyteDataAndMetaDataAndConditionsDF$project_ID, IncucyteDataAndMetaDataAndConditionsDF$plate_ID, sep = '|')

#bio_rep_id
#change composite ID to: join condition_ID, project_ID and plate_ID = get bio rep for each condition
IncucyteSplitPointDF$bio_rep_id<- paste(IncucyteSplitPointDF$condition_ID, IncucyteSplitPointDF$project_ID, IncucyteSplitPointDF$plate_ID, sep = '|')

#get first phase of growth curve by merging IncucyteSplitPointDF
IncucyteDataAndMetaDataAndConditionsDF <- merge(IncucyteDataAndMetaDataAndConditionsDF, IncucyteSplitPointDF, by = c("bio_rep_id","project_ID","plate_ID","condition_ID"), all.x = TRUE, all.y = TRUE)

#remove outer wells/outliers
IncucyteDataAndMetaDataAndConditionsDF= subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$exclude==0)

IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(treatment = case_when(is.na(treatment_0) ~ "Untreated",
                                                                                                                  is.na(treatment_1) & is.na(treatment_2) ~ treatment_0,
                                                                                                                  !(is.na(treatment_1)) & is.na(treatment_2) ~ paste(treatment_0, treatment_1, sep='_'),
                                                                                                                  !(is.na(treatment_1)) & !(is.na(treatment_2)) ~ paste(treatment_0, treatment_1,treatment_2, sep='_')))



IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(treat_conc_id = case_when(is.na(treatment_0) &  is.na(treatment_1)~ paste("Untreated", cell_number, sep='_'),
                                                                                                                      is.na(treatment_0) & !(is.na(treatment_1)) ~ paste(treatment_1, concentration_1, sep='_'),
                                                                                                                      is.na(treatment_1) & is.na(treatment_2) ~ paste(treatment_0, concentration_0, sep='_'),
                                                                                                                      !(is.na(treatment_1)) & is.na(treatment_2) ~ paste(treatment_0, concentration_0,treatment_1, concentration_1, sep='_'),
                                                                                                                      !(is.na(treatment_1)) & !(is.na(treatment_2)) ~ paste(treatment_0, concentration_0,treatment_1, concentration_1,treatment_2, concentration_2, sep='_')))

#generate cell line modification and compound concatenated label column
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(CL_treat_conc_id = paste(gsub("_"," ",cell_line_modifications), gsub("_"," ",treat_conc_id), sep=' | '))


############################################################# Assign Analysis ID to each Analysis Group #############################################################
#ensuring same plates are fitted per comparison group by assigning them an analysis ID
IncucyteDataAndMetaDataAndConditionsLabeledDF <- getSamePlatesPerComparisonGroupQC(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDMSODilutionDF,ListOfSpecialCompounds = "HydroxyUrea" , wildtype = "C631")

#left_join IncucyteDataAndMetaDataAndConditionsLabeledDF to get mean_conf values back
IncucyteDataAndMetaDataAndConditionsDF2 <- subset(IncucyteDataAndMetaDataAndConditionsDF, bio_rep_id %in% IncucyteDataAndMetaDataAndConditionsLabeledDF$bio_rep_id)


IncucyteDataAndMetaDataAndConditionsLabeledDF <- merge(IncucyteDataAndMetaDataAndConditionsLabeledDF, IncucyteDataAndMetaDataAndConditionsDF2, by = c("bio_rep_id","condition_ID","cell_line", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0","treatment_timepoint_0","treatment_duration_0",
                                                                                                                                                      "treatment_1", "concentration_1","treatment_timepoint_1","treatment_duration_1","treatment_2", "concentration_2","treatment_timepoint_2","treatment_duration_2",
                                                                                                                                                      "treatment_3", "concentration_3","treatment_timepoint_3","treatment_duration_3", "project_ID", "plate_ID","well_ID", "condition_ID", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
                                                                                                                                                      "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
                                                                                                                                                      "split_first_conf_difference","phenotype","treatment","CL_treat_conc_id","treat_conc_id"))

#only want to keep conditions labeled with NLMM Analysis ID...
IncucyteDataAndMetaDataAndConditionsLabeledDF <- IncucyteDataAndMetaDataAndConditionsLabeledDF %>% drop_na(NLMM_Analysis_ID)

#average across tech reps
AvgIncucyteDataAndMetaDataAndConditionsLabeledDF<- IncucyteDataAndMetaDataAndConditionsLabeledDF %>%
  group_by(condition_ID,cell_line, cell_line_modifications, cell_number, treatment_0, concentration_0, treatment_1, concentration_1, treatment_2, concentration_2, treatment_3, concentration_3, NLMM_Analysis_ID,bio_rep_id,timepoint, split_timepoint,plate_ID,second_phase,treatment,CL_treat_conc_id,treat_conc_id) %>%
  dplyr::summarise_at(c("confluency"),list(mean = mean, sd = sd, se = ~ sd(.)/ sqrt(length(.))), na.rm=TRUE)

#label treatment & cell line
IncucyteDataAndMetaDataAndConditionsLabeledDF <- IncucyteDataAndMetaDataAndConditionsLabeledDF%>%
  mutate(CL_TRT = paste(cell_line_modifications, treatment_0,sep=" | "))

AvgIncucyteDataAndMetaDataAndConditionsLabeledDF <- AvgIncucyteDataAndMetaDataAndConditionsLabeledDF%>%
  mutate(CL_TRT = paste(cell_line_modifications, treatment_0,sep=" | "))

############################################################# Visualize each Analysis Group per condition with  and without splitpoint #############################################################
# create necessary colors and list of conditions that are to be plotted

CLTRTList <- c("ARID1A KO","ARID2 KO","ARID1B KO","BRD7 KO","BRD9 KO","PBRM1 KO","PHF10 KO","SMARCA4 KO","SMARCC1 KO","SMARCD1 KO","C631")

ColorListTreat <- c("#C52030", "#176533", "#8A181A","#566B30","#D2AE2A","#93C13E", "#4C491D", "#80519B","#3569A1", "#9ABCDC", "darkgrey")

ColorListControl <- c("#C52030", "#176533", "#8A181A","#566B30","#D2AE2A","#93C13E", "#4C491D", "#80519B","#3569A1", "#9ABCDC", "darkgrey")

# list mapping cell line to line type for plot
LinetypeListTreat <- c("solid","solid","solid","solid","solid", "solid","solid","solid", "solid","solid","solid")

LinetypeListControl <- c("dotted","dotted","dotted","dotted","dotted", "dotted","dotted","dotted", "dotted","dotted", "dotted")

treatmentendtime <- 120

ListOfNLMMAnalysisIDs = unique(AvgIncucyteDataAndMetaDataAndConditionsLabeledDF$NLMM_Analysis_ID)

for(j in 1:length(ListOfNLMMAnalysisIDs)){
  IncucyteDataAndMetaDataAndConditionsSubset= subset(AvgIncucyteDataAndMetaDataAndConditionsLabeledDF, AvgIncucyteDataAndMetaDataAndConditionsLabeledDF$NLMM_Analysis_ID==ListOfNLMMAnalysisIDs[j])

  IncucyteDataAndMetaDataAndConditionsSubset <- IncucyteDataAndMetaDataAndConditionsSubset %>%
    mutate(cell_line_modifications = gsub("_", " ", cell_line_modifications))%>%
    mutate(CL_TRT = gsub("_", " ", CL_TRT))%>%
    mutate(plate_ID = gsub("_", " ", plate_ID))%>%
    mutate(treat_conc_id = gsub("_", " ", treat_conc_id))

  TypeOfTreatmentData <- "Treated_Data"

  AnalysisValue <- as.numeric(gsub("nlmm_analysis_", "", unique(IncucyteDataAndMetaDataAndConditionsSubset$NLMM_Analysis_ID)))

  IncucyteDataAndMetaDataAndConditionsSubsetValue <- subset(IncucyteDataAndMetaDataAndConditionsSubset, IncucyteDataAndMetaDataAndConditionsSubset$condition_ID == paste("condition_",AnalysisValue,sep=""))

  CompoundName <- unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$treatment)

  CellLineModification <- gsub(" ", "_", unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$cell_line_modifications))

  # plot analysis group for each condition
  PlotObject <- plotGrowthCurves(IncucyteDataAndMetaDataAndConditionsSubset,CLTRTList,ColorListTreat,ColorListControl,wildtype,treatmentendtime, AnalysisValue,LinetypeListTreat=LinetypeListTreat,LinetypeListControl=LinetypeListControl)


  #make sure directory exists
  SubDir <- sprintf("./Confluency_Curves/%s/%s/Individual_QC_Plots/analysis_group_plots_per_condition/%s/",TypeOfTreatmentData, CompoundName,gsub(" ", "_",CellLineModification))

  if(isTRUE(file.exists(SubDir))==TRUE){
    print(paste("Plotting analysis group growth curves of",CellLineModification, "treated with", CompoundName))
  } else {
    dir.create(SubDir, recursive = TRUE)
    print(paste("Plotting analysis group growth curves of",CellLineModification, "treated with", CompoundName))
  }

  FilenameForPlotObject = sprintf("./Confluency_Curves/%s/%s/Individual_QC_Plots/analysis_group_plots_per_condition/%s/%s_%s_all_plates.pdf",TypeOfTreatmentData, CompoundName,gsub(" ","_",CellLineModification), TodaysDate, unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$treat_conc_id))

  ggsave(filename=FilenameForPlotObject, plot=PlotObject,width = 30, height = 20, dpi = 300, units = "in", device='pdf')


}

