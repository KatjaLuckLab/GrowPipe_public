#!/usr/bin/env RScript

#command-line executable:  chmod +x unclustered_heatmap.R
# run as: ./unclustered_heatmap.R

# Script is part of the modeling-independent analysis and calculates the area under the curve (AUC) values from the dose-response curves for each cell line, treatment and phenotype followed by a paired t-test comparing KO to WT.
# These values are visualized in a heatmap for all tested knockout cell lines against a panel of compounds.Here the computed AUC ratio compound/DMSO for a given KO cell line is normalized to WT in order to understand the effect of knocking out a BAF subunit has on the treated state.

TodaysDate <- ""
AnalysisName<-"Systematic"

################################################## Loading Necessary Libraries and Source Code #########################################################################
source("./R/analysis_group_functions.R")
source("./R/phenotype_quantification_functions.R")
source("./R/plot_functions.R")
source("./R/clustering_functions.R")

library(tidyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(magrittr)
library("RMySQL")
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
library("ComplexHeatmap")
library("colorRamp2")
library(rstatix)
library(openxlsx)
library(factoextra)
library(fpc)
library(NbClust)
library(clValid)
library(GGally)
library(ggbiplot)
library(ppclust)
library(limma)
library(Silhouette)
library(clusterSim)

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

#save entire condition table
IncucyteSplitPointDF=fetch(IncucyteSplitPointTableContents, n = -1)


#get DMSO dilution table
IncucyteDMSODilutionTableContents = dbSendQuery(mysqlconnection, "select * from Incucyte_Compound_DMSO_Dilution_Table")

IncucyteDMSODilutionDF=fetch(IncucyteDMSODilutionTableContents, n = -1)


############################################################# Data Wrangling ############################################################
# labeling of biological and technical replicates by combining meta data with raw confluency data
IncucyteDataAndMetaData <- IncucyteDataDF %>% left_join(IncucyteMetaDataDF, by = c("project_ID","plate_ID","well_ID"))

# join condition df to previous dataframe in order to include cell line, treatments, etc info
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaData %>% left_join(IncucyteConditionDF, by = c("condition_ID"))

# remove outer wells/outliers
IncucyteDataAndMetaDataAndConditionsDF= subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$exclude==0)

# generate bio_rep column:join condition_ID, project_ID and plate_ID
IncucyteDataAndMetaDataAndConditionsDF$bio_rep_id<- paste(IncucyteDataAndMetaDataAndConditionsDF$condition_ID, IncucyteDataAndMetaDataAndConditionsDF$project_ID, IncucyteDataAndMetaDataAndConditionsDF$plate_ID, sep = '|')

# generate bio_rep_id: join condition_ID, project_ID and plate_ID = get bio rep for each condition
IncucyteSplitPointDF$bio_rep_id<- paste(IncucyteSplitPointDF$condition_ID, IncucyteSplitPointDF$project_ID, IncucyteSplitPointDF$plate_ID, sep = '|')

#get phenotypes of growth curve by merging split point data
IncucyteDataAndMetaDataAndConditionsDF <- merge(IncucyteDataAndMetaDataAndConditionsDF, IncucyteSplitPointDF, by = c("bio_rep_id","condition_ID","plate_ID"), all.x = TRUE, all.y = TRUE)

# create treatment label
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(treatment = case_when(is.na(treatment_0) ~ "Untreated",
                                                                                                                  is.na(treatment_1) & is.na(treatment_2) ~ treatment_0,
                                                                                                                  !(is.na(treatment_1)) & is.na(treatment_2) ~ paste(treatment_0, treatment_1, sep='_'),
                                                                                                                  !(is.na(treatment_1)) & !(is.na(treatment_2)) ~ paste(treatment_0, treatment_1,treatment_2, sep='_')))

# remove unwanted treatments
IncucyteDataAndMetaDataAndConditionsDF <- subset(IncucyteDataAndMetaDataAndConditionsDF, !(treatment %in% c("MNNG","MNNG_O6BG")))

# create treat_conc_id label
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(treat_conc_id = case_when(is.na(treatment_0) & is.na(treatment_1) ~ paste("Untreated", cell_number, sep='_'),
                                                                                                                      is.na(treatment_0) & !(is.na(treatment_1)) ~ paste(treatment_1, concentration_1, sep='_'),
                                                                                                                      is.na(treatment_1) & is.na(treatment_2) ~ paste(treatment_0, concentration_0, sep='_'),
                                                                                                                      !(is.na(treatment_1)) & is.na(treatment_2) ~ paste(treatment_0, concentration_0,treatment_1, concentration_1, sep='_'),
                                                                                                                      !(is.na(treatment_1)) & !(is.na(treatment_2)) ~ paste(treatment_0, concentration_0,treatment_1, concentration_1,treatment_2, concentration_2, sep='_')))

#generate cell line modification and compound concatenated label column
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(CL_treat_conc_id = paste(gsub("_"," ",cell_line_modifications), gsub("_"," ",treat_conc_id), sep=' | '))

# filter for conditions with more than one biological replicate
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>%
  group_by(condition_ID) %>%
  dplyr::distinct(bio_rep_id,plate_ID, cell_line_modifications, cell_number, treatment_0, concentration_0,treatment_timepoint_0,treatment_duration_0,
                  treatment_1, concentration_1,treatment_timepoint_1,treatment_duration_1,treatment_2, concentration_2,treatment_timepoint_2,treatment_duration_2,
                  treatment_3, concentration_3,treatment_timepoint_3,treatment_duration_3, split_timepoint, logistic_growth_phase, logistic_growth_duration, second_phase, second_phase_duration,
                  after_24h_confluency,  after_48h_confluency,after_72h_confluency,last_confluency,conf_96hr_48hr_ratio,
                  split_first_conf_difference, phenotype,treatment,CL_treat_conc_id,treat_conc_id) %>%
  dplyr::filter(n()>1)


################################ Assigning Identifier per Analysis Group ################################
ListOfSpecialCompounds <- c("HydroxyUrea")
wildtype <- "C631"

# generate analysis ID for KO treated, WT treated, KO untreated, and WT untreated for each compound and concentration
IncucyteDataAndMetaDataAndConditionsLabeledDF <- getSamePlatesPerComparisonGroup(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDMSODilutionDF,ListOfSpecialCompounds = ListOfSpecialCompounds,wildtype = wildtype)


################################ Calculate Treatment:DMSO ratio for each phenotype ################################
ListOfSpecialCompounds <- c("HydroxyUrea")

controlcompound <- "DMSO"

# normalize phenotypes to control treatment
PhenotypesNormToDMSODF <- normalizeToDMSO(IncucyteDataAndMetaDataAndConditionsLabeledDF,controlcompound,ListOfSpecialCompounds)

# rename cell lines
PhenotypesNormToDMSODF <- PhenotypesNormToDMSODF %>%
  mutate(cell_line_modification = replace(cell_line_modification, cell_line_modification == "C631", "WT"))

# update wildtype variable
wildtype<- "WT"

# filter for DRC of bioreps for a compound with more than two concentrations
PhenotypesNormToDMSODF <- PhenotypesNormToDMSODF %>% dplyr::group_by(treatment_0,cell_line_modification, plate_ID) %>% dplyr::filter(n()>2)

################################ Calculating AUC & performing paired t-test ################################

# calculating area under the curve (AUC) of the dose-response curves normalized to DMSO per condition & phenotype
AUCsPerConditionDF <- getAUCsPerConditionNormalizedToDMSO(PhenotypesNormToDMSODF,wildtype)

# Apply pairwise t-Test for each condition & phenotype across all biological replicates - using limma package
CompleteAUCDF <- getAUCSignificance(AUCsPerConditionDF,wildtype)

################################ Quality Control: Plotting AUC values for KO and paired WT per Phenotype ################################

# Create necessary metadata for AUC Quality Control Plots
################################################################################

# remove underscores from cell line names - will use these for plots
AUCsPerConditionDF<- AUCsPerConditionDF %>% mutate(cell_line_modification = gsub("_"," ",cell_line_modification))

# create BAFSubtypeDF - will map subtype label to each cell line
CL <- unique(AUCsPerConditionDF$cell_line_modification)

BAFSubtypeDF <- data.frame(CL)

BAFSubtypeDF <- BAFSubtypeDF %>% mutate(Subtype=case_when(
  CL %in% c("SMARCA4 KO","SMARCC1 KO", "SMARCD1 KO") ~ "general BAF",
  CL %in% c("ARID1A KO", "ARID1B KO") ~ "cBAF",
  CL %in% c("ARID2 KO", "BRD7 KO", "PBRM1 KO", "PHF10 KO") ~ "PBAF",
  CL =="BRD9 KO" ~ "ncBAF",
  CL == "WT" ~ "control"))

# create list of growth curve phenotypes to parse through for plotting
ListOfPhenotypes <- c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency","First_Phase_Duration","First_Phase_Change_in_Confluency")

################################ Plotting Heatmap for each Phenotype ################################

#Create necessary metadata for heatmaps
################################################################################

# create meta data linking compounds to their mechanism of action

CompleteAUCDF<- CompleteAUCDF %>%mutate(treatment = gsub("_"," + ",treatment))

Compound <- unique(CompleteAUCDF$treatment)

CompoundMOADF <- data.frame(Compound)

CompoundMOADF <- CompoundMOADF %>% mutate(MOA=case_when(Compound %in% c("Palbociclib") ~ "G1/S phase",
                                                        Compound %in% c("Aphidicolin","Gemcitabine","HydroxyUrea") ~ "S phase",
                                                        Compound %in% c("Adavosertib","BI-2536","Paclitaxel") ~ "G2/M phase",
                                                        Compound %in% c("Camptothecin","Topotecan") ~ "TopI inhibitor",
                                                        Compound %in% c("Doxorubicin","Etoposide") ~ "TopII inhibitor",
                                                        Compound %in% c("IlludinS","Methyl-Methanesulfonate") ~ "DNA Alkylating",
                                                        Compound == "Niraparib" ~ "PARP inhibitor",
                                                        Compound %in% c("MLN4924","NMS-873") ~ "Protein homeostasis",
                                                        Compound %in% c("Lapatinib","Sapitinib") ~ "EGFR inhibitor",
                                                        Compound %in% c("Cobimetinib","Mirdametinib") ~ "MEK inhibitor",
                                                        Compound =="BIBR1532" ~ "Telomerase inhibitor"))


#create list of compound MOA linked to their respective color - used in heatmap
MOAColors = list("MOA" = c("G1/S phase"= "lightgreen","S phase"= "green","S phase"= "green","S phase"= "green",
                           "G2/M phase"="yellowgreen","G2/M phase"="yellowgreen","G2/M phase"="yellowgreen",
                           "DNA Alkylating"="pink","DNA Alkylating"="pink","TopI inhibitor"="violet",
                           "TopI inhibitor"="violet", "TopII inhibitor"="purple","TopII inhibitor"="purple",
                           "PARP inhibitor" = "orange","Protein homeostasis"="red","Protein homeostasis"='red',
                           "EGFR inhibitor"="blue","EGFR inhibitor"="blue",
                           "MEK inhibitor"='lightblue', "MEK inhibitor"="lightblue", "Telomerase inhibitor"="black"))


# update BAF subtype metadata to exclude WT - as mean log2 fold change of AUC per cell line is KO vs WT
BAFSubtypeDF <- subset(BAFSubtypeDF, CL != "WT")

# create list of BAF subtype linked to their respective color - used in heatmap
SubtypeColors <- list ("Subtype"=c("cBAF"="darkred", "PBAF"="darkgreen","ncBAF"="gold", "general BAF"="blue"))

# list providing desired ordering of compounds in heatmap
ListOfCompoundOrder = c("Palbociclib", "Aphidicolin","Gemcitabine","HydroxyUrea","Adavosertib","BI-2536", "Paclitaxel",
                        "IlludinS","Methyl-Methanesulfonate", "Camptothecin", "Topotecan","Doxorubicin","Etoposide","Niraparib","MLN4924",
                        "NMS-873","Lapatinib", "Sapitinib","Cobimetinib","Mirdametinib","BIBR1532")

# list providing desired ordering of cell lines in heatmap
ListOfCLOrder = c("ARID1A KO","ARID1B KO","ARID2 KO","BRD7 KO","PHF10 KO","PBRM1 KO", "BRD9 KO","SMARCD1 KO","SMARCC1 KO","SMARCA4 KO")


# desired min and max AUC values for color range
scaling <- "scaled"
minColorRange <- -0.5
maxColorRange <- 0.5

# how KO and WT should be labeled in heatmap
kolabel<-"KO"
wtlabel <- "WT"

################################################################################

###### Plot unclustered & clustered heatmaps for each phenotype ######

for(n in 1:length(ListOfPhenotypes)){

  phenotype <- ListOfPhenotypes[n]

  PhenotypeName <- gsub("_", " ", phenotype)

  # subset dataframe containing mean log2 fold change of AUC KO vs WT per compound for one given growth curve phenotype
  CompleteAUCDFSubset <- subset(CompleteAUCDF, Phenotype==phenotype)

  # check if directory exists - create if it doesn't
  SubDir <- sprintf("./Modeling-independent_Phenotypes/Growth_Curve_Phenotypes/%s_%s/Heatmaps/Unclustered/", n, phenotype)

  if(isTRUE(file.exists(SubDir))==TRUE){
    print(paste("Plotting an unclustered heatmap for", PhenotypeName))
  } else {
    dir.create(SubDir, recursive = TRUE)
    print(paste("Plotting an unclustered heatmap for", PhenotypeName))
  }

  NewFilename = sprintf("./Modeling-independent_Phenotypes/Growth_Curve_Phenotypes/%s_%s/Heatmaps/Unclustered/%s_%s_unclustered_%s_%s.pdf", n, phenotype, TodaysDate, phenotype, AnalysisName,scaling)

  pdf(NewFilename, width = 7, height = 7)#8

  #Plot unclustered heatmap and save
  comprehensive_heatmap_unclust <-plotAUCHeatmapUnclustered(CompleteAUCDFSubset,CompoundMOADF,MOAColors,BAFSubtypeDF,SubtypeColors,ListOfCompoundOrder,ListOfCLOrder,kolabel,wtlabel,scaling,minColorRange,maxColorRange)

  draw(comprehensive_heatmap_unclust,merge_legend = TRUE, column_dend_side = "top")
  dev.off()

}
