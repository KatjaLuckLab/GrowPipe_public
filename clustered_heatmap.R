#!/usr/bin/env RScript

#command-line executable:  chmod +x clustered_heatmap.R
# run as: ./clustered_heatmap.R

#Script for clustering across the 5 modeling-independent phenotypes (Confluency after 24hr, 48hr, 72hr, & 96hr along with the relative change in confluency in the second phase
# (confluency after 96hr/48hr)) testing 4 methods: PAM, hierarchical, k-means and c-means. Clustered summary heatmaps using the various clustering methods are produced, along with plots demonstrating the clustering accuracy for each.
#  Here the computed AUC ratio compound/DMSO for a given KO cell line is normalized to WT in order to understand the effect of knocking out a BAF subunit has on the treated state.

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
#labeling of bio and technical replicates by combining meta data with df
IncucyteDataAndMetaData <- IncucyteDataDF %>% left_join(IncucyteMetaDataDF, by = c("project_ID","plate_ID","well_ID"))

#join condition df to previous dataframe in order to include cell line, treatments, etc info
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaData %>% left_join(IncucyteConditionDF, by = c("condition_ID"))

#remove outer wells/outliers
IncucyteDataAndMetaDataAndConditionsDF= subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$exclude==0)#0

#generate bio_rep column:join condition_ID, project_ID and plate_ID
IncucyteDataAndMetaDataAndConditionsDF$bio_rep_id<- paste(IncucyteDataAndMetaDataAndConditionsDF$condition_ID, IncucyteDataAndMetaDataAndConditionsDF$project_ID, IncucyteDataAndMetaDataAndConditionsDF$plate_ID, sep = '|')

#bio_rep_id
#change composite ID to: join condition_ID, project_ID and plate_ID = get bio rep for each condition
IncucyteSplitPointDF$bio_rep_id<- paste(IncucyteSplitPointDF$condition_ID, IncucyteSplitPointDF$project_ID, IncucyteSplitPointDF$plate_ID, sep = '|')

#get first phase of growth curve by merging IncucyteSplitPointDF
IncucyteDataAndMetaDataAndConditionsDF <- merge(IncucyteDataAndMetaDataAndConditionsDF, IncucyteSplitPointDF, by = c("bio_rep_id","condition_ID","plate_ID"), all.x = TRUE, all.y = TRUE)

# create treatment label
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(treatment = case_when(is.na(treatment_0) ~ "Untreated",
                                                                                                                  is.na(treatment_1) & is.na(treatment_2) ~ treatment_0,
                                                                                                                  !(is.na(treatment_1)) & is.na(treatment_2) ~ paste(treatment_0, treatment_1, sep='_'),
                                                                                                                  !(is.na(treatment_1)) & !(is.na(treatment_2)) ~ paste(treatment_0, treatment_1,treatment_2, sep='_')))

# remove unwanted treatments
IncucyteDataAndMetaDataAndConditionsDF <- subset(IncucyteDataAndMetaDataAndConditionsDF, !(treatment %in% c("MNNG","MNNG_O6BG")))

#generate compound and concentration concatenated label column
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

# normalize phenotypes to DMSO
PhenotypesNormToDMSODF <- normalizeToDMSO(IncucyteDataAndMetaDataAndConditionsLabeledDF,controlcompound,ListOfSpecialCompounds)

# rename cell lines
PhenotypesNormToDMSODF <- PhenotypesNormToDMSODF %>%
  mutate(cell_line_modification = replace(cell_line_modification, cell_line_modification == "C631", "WT"))

# update wildtype variable
wildtype<- "WT"

#make sure only create DRC of bioreps for a compound with more than two concentrations
PhenotypesNormToDMSODF <- PhenotypesNormToDMSODF %>% dplyr::group_by(treatment_0,cell_line_modification, plate_ID) %>% dplyr::filter(n()>2)


################################ Calculating AUC & performing paired t-test ################################

# calculating area under the curve (AUC) of the dose-response curves normalized to DMSO per condition
AUCsPerConditionDF <- getAUCsPerConditionNormalizedToDMSO(PhenotypesNormToDMSODF,wildtype)

# Apply pairwise t-Test for each condition for each phenotype across all biological replicates
CompleteAUCDF <- getAUCSignificance(AUCsPerConditionDF,wildtype)

#Create necessary metadata for heatmaps
################################################################################

# create metadata mapping compound to mechanism of action
Compound <- unique(CompleteAUCDF$treatment)

CompoundMOADF <- data.frame(Compound)


CompoundMOADF <- CompoundMOADF %>% mutate(MOA=case_when(Compound %in% c("Palbociclib","Aphidicolin","Gemcitabine","HydroxyUrea","Adavosertib","BI-2536","Paclitaxel") ~ "Cell cycle regulation",
                                                        Compound %in% c("Camptothecin","Topotecan","Doxorubicin","Etoposide","IlludinS","Methyl-Methanesulfonate","Niraparib") ~ "DNA damage",
                                                        Compound %in% c("MLN4924","NMS-873") ~ "Protein homeostasis",
                                                        Compound =="BIBR1532" ~ "Telomerase inhibitor",
                                                        Compound %in% c("Cobimetinib","Mirdametinib","Lapatinib","Sapitinib") ~ "Growth signaling"))


MOAColors = list("MOA" = c("Cell cycle regulation"= "skyblue","DNA damage"= "darkorange","Protein homeostasis"= "gold","Telomerase inhibitor"= "black","Growth signaling"="green"))

# create metadata mapping cell line to BAF subtype
CL <- gsub("_"," ",unique(AUCsPerConditionDF$cell_line_modification))

BAFSubtypeDF <- data.frame(CL)

BAFSubtypeDF <- BAFSubtypeDF %>% mutate(Subtype=case_when(
  CL %in% c("SMARCA4 KO","SMARCC1 KO", "SMARCD1 KO") ~ "general BAF",
  CL %in% c("ARID1A KO", "ARID1B KO") ~ "cBAF",
  CL %in% c("ARID2 KO", "BRD7 KO", "PBRM1 KO", "PHF10 KO") ~ "PBAF",
  CL =="BRD9 KO" ~ "ncBAF"))

BAFSubtypeDF <- subset(BAFSubtypeDF, CL != "WT")

# list mapping BAF subtype to color for heatmap
SubtypeColors <- list ("Subtype"=c("cBAF"="darkred", "PBAF"="darkgreen","ncBAF"="gold","general BAF"="blue"))

# list of desired compound order for unclustered heatmap
ListOfCompoundOrder = c("Gemcitabine","Lapatinib", "Palbociclib", "Aphidicolin","HydroxyUrea","Adavosertib","BI-2536", "Paclitaxel","Camptothecin", "Topotecan","Doxorubicin","Etoposide", "IlludinS","Methyl-Methanesulfonate","MLN4924","BIBR1532",
                        "Cobimetinib","Mirdametinib","NMS-873","Niraparib","Sapitinib")

# list of desired cell line order for unclustered heatmap
ListOfCLOrder = c("SMARCA4 KO","SMARCC1 KO","SMARCD1 KO","ARID1A KO","ARID1B KO","ARID2 KO","BRD7 KO","PHF10 KO","PBRM1 KO", "BRD9 KO")

################################################################################

############# Preparing Data for Clustering across 5 Phenotypes (Conf after 24h, 48h, 72h, 96 and relative change in conf in second phase (96h/48h)) ################

# create necessary input data -- transform AUCDF into matrix where column names are phenotypes and rownames are conditions
AUCDF5PhenotypesDF <- getSummaryAUCMatrix(CompleteAUCDF)

#do the same for AUCSignificanceDF
PValue5PhenotypesDF <- getSummaryPValueMatrix(CompleteAUCDF)


#Visualize manner transformed data groups based on already existing labels
#creating various labels that will be used in subsequent external cluster validation steps
################################################################################

AUCDF5PhenotypesDF <- AUCDF5PhenotypesDF %>%
  mutate(symbol = case_when(treatment=="IlludinS" ~ "I",
                            treatment=="MLN4924" ~ "M",
                            treatment=="Etoposide" ~ "E",
                            treatment=="Methyl-Methanesulfonate" ~ "MM",
                            treatment=="Camptothecin" ~ "C",
                            treatment=="BIBR1532" ~ "B",
                            treatment=="BI-2536" ~ "BI",
                            treatment=="Topotecan" ~ "T",
                            treatment=="Lapatinib" ~ "L",
                            treatment=="Palbociclib" ~ "P",
                            treatment=="Paclitaxel" ~ "Px",
                            treatment=="Aphidicolin" ~ "Ap",
                            treatment=="Gemcitabine" ~ "G",
                            treatment=="HydroxyUrea" ~ "H",
                            treatment=="Doxorubicin" ~ "D",
                            treatment=="Adavosertib" ~ "A",
                            treatment=="Cobimetinib" ~ "CB",
                            treatment=="Mirdametinib" ~"MI",
                            treatment=="NMS-873" ~ "NM",
                            treatment=="Niraparib" ~ "N",
                            treatment=="Sapitinib" ~ "S"

  ))%>%
  mutate(colorkey = case_when(cell_line_modification =="BRD9_KO" ~ "gold",
                              cell_line_modification =="BRD7_KO" ~ "darkgreen",
                              cell_line_modification =="ARID2_KO" ~ "darkgreen",
                              cell_line_modification =="PHF10_KO" ~ "darkgreen",
                              cell_line_modification =="PBRM1_KO" ~ "darkgreen",
                              cell_line_modification =="SMARCA4_KO" ~ "blue",
                              cell_line_modification =="SMARCC1_KO" ~ "blue",
                              cell_line_modification =="SMARCD1_KO" ~ "blue",
                              cell_line_modification =="ARID1A_KO" ~ "darkred",
                              cell_line_modification =="ARID1B_KO" ~ "darkred"))%>%
  mutate(BAF_Subtype = case_when(cell_line_modification =="BRD9_KO" ~ "ncBAF",
                                 cell_line_modification =="BRD7_KO" ~ "PBAF",
                                 cell_line_modification =="ARID2_KO" ~ "PBAF",
                                 cell_line_modification =="PHF10_KO" ~ "PBAF",
                                 cell_line_modification =="PBRM1_KO" ~ "PBAF",
                                 cell_line_modification =="SMARCA4_KO" ~ "general BAF",
                                 cell_line_modification =="SMARCC1_KO" ~ "general BAF",
                                 cell_line_modification =="SMARCD1_KO" ~ "general BAF",
                                 cell_line_modification =="ARID1A_KO" ~ "cBAF",
                                 cell_line_modification =="ARID1B_KO" ~ "cBAF")) %>%
  mutate(Compound_MOA = case_when(treatment %in% c("Palbociclib","Aphidicolin","Gemcitabine","HydroxyUrea","Adavosertib","BI-2536","Paclitaxel") ~ "Cell cycle regulation",
                                  treatment %in% c("Camptothecin","Topotecan","Doxorubicin","Etoposide","IlludinS","Methyl-Methanesulfonate","Niraparib") ~ "DNA damage",
                                  treatment%in% c("MLN4924","NMS-873") ~ "Protein homeostasis",
                                  treatment%in% c("Cobimetinib","Mirdametinib","Lapatinib","Sapitinib") ~ "Growth signaling",
                                  treatment =="BIBR1532" ~ "Telomerase inhibitor"))%>%
  mutate(Domain = case_when(
    cell_line_modification ==  "ARID1A_KO"  ~ "ARID",
    cell_line_modification ==  "ARID1B_KO"  ~ "ARID",
    cell_line_modification ==  "ARID2_KO"  ~ "ARID",
    cell_line_modification ==  "BRD7_KO"  ~ "Bromo",
    cell_line_modification ==  "BRD9_KO"  ~ "Bromo",
    cell_line_modification ==  "PBRM1_KO"  ~ "Bromo",
    cell_line_modification ==  "PHF10_KO"  ~ "PHD_Zinc",
    cell_line_modification ==  "SMARCC1_KO"  ~ "SWIRM",
    cell_line_modification ==  "SMARCD1_KO"  ~ "SWIB",
    cell_line_modification ==  "SMARCA4_KO"  ~ "HSA")) %>%
  mutate(Module = case_when(
    cell_line_modification ==  "ARID1A_KO"  ~ "accessory",
    cell_line_modification ==  "ARID1B_KO"  ~ "accessory",
    cell_line_modification ==  "ARID2_KO"  ~ "accessory",
    cell_line_modification ==  "BRD7_KO"  ~ "accessory",
    cell_line_modification ==  "BRD9_KO"  ~ "accessory",
    cell_line_modification ==  "PBRM1_KO"  ~ "accessory",
    cell_line_modification ==  "PHF10_KO"  ~ "accessory",
    cell_line_modification ==  "SMARCC1_KO"  ~ "core",
    cell_line_modification ==  "SMARCD1_KO"  ~ "core",
    cell_line_modification ==  "SMARCA4_KO"  ~ "atpase"))%>%
  mutate(Subtype_MOA = paste(BAF_Subtype, Compound_MOA, sep="|")) %>%
  mutate(Compound_CL = paste(symbol, cell_line_modification, sep="|")) %>%
  mutate(CL_Compound = paste(cell_line_modification,treatment, sep="|")) %>%
  mutate(Subtype_Compound = paste(BAF_Subtype, treatment, sep="|")) %>%
  mutate(CL_MOA = paste(cell_line_modification, Compound_MOA, sep="|")) %>%
  mutate(Module_MOA = paste(Module, Compound_MOA, sep="|")) %>%
  mutate(Module_Compound = paste(Module, treatment, sep="|")) %>%
  mutate(Domain_MOA = paste(Domain, Compound_MOA, sep="|")) %>%
  mutate(Domain_Compound = paste(Domain, treatment, sep="|"))

# converting labels into numerical values so they may be used in subsequent external cluster validation

AUCDF5PhenotypesDF <- AUCDF5PhenotypesDF %<>%
  mutate(Subtype_MOA_Class = case_when(
    Subtype_MOA ==  "cBAF|Cell cycle regulation"  ~ 1,
    Subtype_MOA ==  "PBAF|Cell cycle regulation" ~ 2,
    Subtype_MOA == "ncBAF|Cell cycle regulation"~ 3,
    Subtype_MOA == "general BAF|Cell cycle regulation"  ~ 4,
    Subtype_MOA == "cBAF|DNA damage"  ~ 5,
    Subtype_MOA == "PBAF|DNA damage" ~ 6,
    Subtype_MOA == "ncBAF|DNA damage" ~ 7,
    Subtype_MOA == "general BAF|DNA damage" ~ 8,
    Subtype_MOA == "cBAF|Protein homeostasis" ~ 9,
    Subtype_MOA == "PBAF|Protein homeostasis" ~ 10,
    Subtype_MOA == "ncBAF|Protein homeostasis"  ~ 11,
    Subtype_MOA == "general BAF|Protein homeostasis" ~ 12,
    Subtype_MOA == "cBAF|Growth signaling" ~ 13,
    Subtype_MOA == "PBAF|Growth signaling" ~ 14,
    Subtype_MOA == "ncBAF|Growth signaling"    ~ 15,
    Subtype_MOA == "general BAF|Growth signaling" ~ 16,
    Subtype_MOA == "cBAF|Telomerase inhibitor"  ~ 17,
    Subtype_MOA == "PBAF|Telomerase inhibitor" ~ 18,
    Subtype_MOA == "ncBAF|Telomerase inhibitor" ~ 19,
    Subtype_MOA == "general BAF|Telomerase inhibitor" ~ 20

  )) %<>%
  mutate(Subtype_Compound_Class = case_when(
    Subtype_Compound ==  "cBAF|Adavosertib"  ~ 1,
    Subtype_Compound ==  "PBAF|Adavosertib" ~ 2,
    Subtype_Compound == "ncBAF|Adavosertib" ~ 3,
    Subtype_Compound == "general BAF|Adavosertib"  ~ 4,
    Subtype_Compound == "cBAF|Aphidicolin"  ~ 5,
    Subtype_Compound == "PBAF|Aphidicolin" ~ 6,
    Subtype_Compound == "ncBAF|Aphidicolin"  ~ 7,
    Subtype_Compound == "general BAF|Aphidicolin" ~ 8,
    Subtype_Compound == "cBAF|BIBR1532" ~ 9,
    Subtype_Compound == "PBAF|BIBR1532"  ~ 10,
    Subtype_Compound == "ncBAF|BIBR1532"  ~ 11,
    Subtype_Compound == "general BAF|BIBR1532"  ~ 12,
    Subtype_Compound == "cBAF|BI-2536" ~ 13,
    Subtype_Compound == "PBAF|BI-2536" ~ 14,
    Subtype_Compound == "ncBAF|BI-2536"    ~ 15,
    Subtype_Compound == "general BAF|BI-2536"   ~ 16,
    Subtype_Compound == "cBAF|Camptothecin" ~ 17,
    Subtype_Compound == "PBAF|Camptothecin"  ~ 18,
    Subtype_Compound == "ncBAF|Camptothecin"  ~ 19,
    Subtype_Compound == "general BAF|Camptothecin"  ~ 20,
    Subtype_Compound ==  "cBAF|Doxorubicin"  ~ 21,
    Subtype_Compound ==  "PBAF|Doxorubicin" ~ 22,
    Subtype_Compound ==  "ncBAF|Doxorubicin"~ 23,
    Subtype_Compound ==  "general BAF|Doxorubicin"~ 24,
    Subtype_Compound ==  "cBAF|Etoposide" ~ 25,
    Subtype_Compound == "PBAF|Etoposide" ~ 26,
    Subtype_Compound == "ncBAF|Etoposide" ~ 27,
    Subtype_Compound == "general BAF|Etoposide"   ~ 28,
    Subtype_Compound =="cBAF|Gemcitabine"  ~ 29,
    Subtype_Compound ==   "PBAF|Gemcitabine" ~ 30,
    Subtype_Compound ==  "ncBAF|Gemcitabine"   ~ 31,
    Subtype_Compound == "general BAF|Gemcitabine"  ~ 32,
    Subtype_Compound == "cBAF|HydroxyUrea" ~ 33,
    Subtype_Compound == "PBAF|HydroxyUrea"  ~ 34,
    Subtype_Compound == "ncBAF|HydroxyUrea"   ~ 35,
    Subtype_Compound == "general BAF|HydroxyUrea"  ~ 36,
    Subtype_Compound == "cBAF|IlludinS" ~ 37,
    Subtype_Compound == "PBAF|IlludinS" ~ 38,
    Subtype_Compound == "ncBAF|IlludinS" ~ 39,
    Subtype_Compound == "general BAF|IlludinS" ~ 40,
    Subtype_Compound ==  "cBAF|Lapatinib"   ~ 41,
    Subtype_Compound ==  "PBAF|Lapatinib"  ~ 42,
    Subtype_Compound == "ncBAF|Lapatinib" ~ 43,
    Subtype_Compound == "general BAF|Lapatinib"  ~ 44,
    Subtype_Compound == "cBAF|MLN4924"   ~ 45,
    Subtype_Compound == "PBAF|MLN4924"  ~ 46,
    Subtype_Compound == "ncBAF|MLN4924"~ 47,
    Subtype_Compound == "general BAF|MLN4924"  ~ 48,
    Subtype_Compound == "cBAF|Methyl-Methanesulfonate"  ~ 49,
    Subtype_Compound == "PBAF|Methyl-Methanesulfonate"  ~ 50,
    Subtype_Compound ==  "ncBAF|Methyl-Methanesulfonate"  ~ 51,
    Subtype_Compound ==  "general BAF|Methyl-Methanesulfonate" ~ 52,
    Subtype_Compound == "cBAF|Palbociclib" ~ 53,
    Subtype_Compound == "PBAF|Palbociclib"  ~ 54,
    Subtype_Compound == "ncBAF|Palbociclib"   ~ 55,
    Subtype_Compound == "general BAF|Palbociclib" ~ 56,
    Subtype_Compound == "cBAF|Paclitaxel"  ~ 57,
    Subtype_Compound == "PBAF|Paclitaxel" ~ 58,
    Subtype_Compound == "ncBAF|Paclitaxel" ~ 59,
    Subtype_Compound == "general BAF|Paclitaxel" ~ 60,
    Subtype_Compound ==  "cBAF|Topotecan"    ~ 61,
    Subtype_Compound ==  "PBAF|Topotecan"  ~ 62,
    Subtype_Compound == "ncBAF|Topotecan"   ~ 63,
    Subtype_Compound == "general BAF|Topotecan"  ~ 64,

    Subtype_Compound ==  "cBAF|Cobimetinib"    ~ 65,
    Subtype_Compound ==  "PBAF|Cobimetinib"  ~ 66,
    Subtype_Compound == "ncBAF|Cobimetinib"   ~ 67,
    Subtype_Compound == "general BAF|Cobimetinib"  ~ 68,


    Subtype_Compound ==  "cBAF|Mirdametinib"    ~ 69,
    Subtype_Compound ==  "PBAF|Mirdametinib"  ~ 70,
    Subtype_Compound == "ncBAF|Mirdametinib"   ~ 71,
    Subtype_Compound == "general BAF|Mirdametinib"  ~ 72,


    Subtype_Compound ==  "cBAF|NMS-873"    ~ 73,
    Subtype_Compound ==  "PBAF|NMS-873"  ~ 74,
    Subtype_Compound == "ncBAF|NMS-873"   ~ 75,
    Subtype_Compound == "general BAF|NMS-873"  ~ 76,

    Subtype_Compound ==  "cBAF|Niraparib"    ~ 77,
    Subtype_Compound ==  "PBAF|Niraparib"  ~ 78,
    Subtype_Compound == "ncBAF|Niraparib"   ~ 79,
    Subtype_Compound == "general BAF|Niraparib"  ~ 80,

    Subtype_Compound ==  "cBAF|Sapitinib"    ~ 81,
    Subtype_Compound ==  "PBAF|Sapitinib"  ~ 82,
    Subtype_Compound == "ncBAF|Sapitinib"   ~ 83,
    Subtype_Compound == "general BAF|Sapitinib"  ~ 84

  ))%<>%
  mutate(Subtype_Class = case_when(
    BAF_Subtype ==  "cBAF"  ~ 1,
    BAF_Subtype ==  "PBAF" ~ 2,
    BAF_Subtype == "ncBAF" ~ 3,
    BAF_Subtype == "general BAF"  ~ 4))%<>%
  mutate(Compound_Class = case_when(
    treatment ==  "Adavosertib"  ~ 1,
    treatment ==  "Aphidicolin" ~ 2,
    treatment == "BIBR1532" ~ 3,
    treatment == "BI-2536"  ~ 4,
    treatment ==  "Camptothecin"  ~ 5,
    treatment ==  "Doxorubicin" ~ 6,
    treatment == "Etoposide" ~ 7,
    treatment == "Gemcitabine"  ~ 8,
    treatment == "HydroxyUrea"  ~ 9,
    treatment ==  "IlludinS"  ~ 10,
    treatment ==  "Lapatinib" ~ 11,
    treatment == "MLN4924" ~ 12,
    treatment == "Methyl-Methanesulfonate"  ~ 13,
    treatment ==  "Palbociclib" ~ 14,
    treatment == "Paclitaxel" ~ 15,
    treatment == "Topotecan"  ~ 16,

    treatment ==  "Cobimetinib" ~ 17,
    treatment == "Mirdametinib" ~ 18,
    treatment == "NMS-873"  ~ 19,
    treatment ==  "Niraparib" ~ 20,
    treatment == "Sapitinib"  ~ 21

  ))%<>%
  mutate(CL_MOA_Class = case_when(
    CL_MOA ==  "ARID1A_KO|Cell cycle regulation"   ~ 1,
    CL_MOA ==  "ARID1B_KO|Cell cycle regulation" ~ 2,
    CL_MOA == "ARID2_KO|Cell cycle regulation"  ~ 3,
    CL_MOA == "BRD7_KO|Cell cycle regulation"  ~ 4,
    CL_MOA ==  "BRD9_KO|Cell cycle regulation"  ~ 5,
    CL_MOA ==  "PBRM1_KO|Cell cycle regulation" ~ 6,
    CL_MOA == "PHF10_KO|Cell cycle regulation"  ~ 7,
    CL_MOA == "SMARCA4_KO|Cell cycle regulation"   ~ 8,
    CL_MOA == "SMARCC1_KO|Cell cycle regulation"  ~ 9,
    CL_MOA ==  "SMARCD1_KO|Cell cycle regulation"  ~ 10,
    CL_MOA ==  "DPF2_KO|Cell cycle regulation"  ~ 11,
    CL_MOA ==  "ARID1A_KO|DNA damage"  ~ 12,
    CL_MOA == "ARID1B_KO|DNA damage" ~ 13,
    CL_MOA == "ARID2_KO|DNA damage"   ~ 14,
    CL_MOA ==  "BRD7_KO|DNA damage"  ~ 15,
    CL_MOA == "BRD9_KO|DNA damage"~ 16,
    CL_MOA == "PBRM1_KO|DNA damage"  ~ 17,
    CL_MOA ==  "PHF10_KO|DNA damage" ~ 18,
    CL_MOA == "SMARCA4_KO|DNA damage"    ~ 19,
    CL_MOA == "SMARCC1_KO|DNA damage" ~ 20,
    CL_MOA ==  "SMARCD1_KO|DNA damage"   ~ 21,
    CL_MOA ==  "DPF2_KO|DNA damage"   ~ 22,

    CL_MOA ==  "ARID1A_KO|Telomerase inhibitor"    ~ 23,
    CL_MOA ==  "ARID1B_KO|Telomerase inhibitor" ~ 24,
    CL_MOA == "ARID2_KO|Telomerase inhibitor"  ~ 25,
    CL_MOA == "BRD7_KO|Telomerase inhibitor"  ~ 26,
    CL_MOA ==  "BRD9_KO|Telomerase inhibitor"  ~ 27,
    CL_MOA ==  "PBRM1_KO|Telomerase inhibitor" ~ 28,
    CL_MOA == "PHF10_KO|Telomerase inhibitor"  ~ 29,
    CL_MOA == "SMARCA4_KO|Telomerase inhibitor"   ~ 30,
    CL_MOA == "SMARCC1_KO|Telomerase inhibitor"  ~ 31,
    CL_MOA ==  "SMARCD1_KO|Telomerase inhibitor"  ~ 32,

    CL_MOA ==  "ARID1A_KO|Protein homeostasis"   ~ 33,
    CL_MOA ==  "ARID1B_KO|Protein homeostasis" ~ 34,
    CL_MOA == "ARID2_KO|Protein homeostasis"  ~ 35,
    CL_MOA == "BRD7_KO|Protein homeostasis"  ~ 36,
    CL_MOA ==  "BRD9_KO|Protein homeostasis"  ~ 37,
    CL_MOA ==  "PBRM1_KO|Protein homeostasis" ~ 38,
    CL_MOA == "PHF10_KO|Protein homeostasis"  ~ 39,
    CL_MOA == "SMARCA4_KO|Protein homeostasis"   ~ 40,
    CL_MOA == "SMARCC1_KO|Protein homeostasis"  ~ 41,
    CL_MOA ==  "SMARCD1_KO|Protein homeostasis"  ~ 42,
    CL_MOA ==  "DPF2_KO|Protein homeostasis"  ~ 43,
    CL_MOA ==  "ARID1A_KO|Growth signaling"   ~ 44,
    CL_MOA ==  "ARID1B_KO|Growth signaling" ~ 45,
    CL_MOA == "ARID2_KO|Growth signaling"  ~ 46,
    CL_MOA == "BRD7_KO|Growth signaling"  ~ 47,
    CL_MOA ==  "BRD9_KO|Growth signaling"  ~ 48,
    CL_MOA ==  "PBRM1_KO|Growth signaling" ~ 49,
    CL_MOA == "PHF10_KO|Growth signaling"  ~ 50,
    CL_MOA == "SMARCA4_KO|Growth signaling"   ~ 51,
    CL_MOA == "SMARCC1_KO|Growth signaling"  ~ 52,
    CL_MOA ==  "SMARCD1_KO|Growth signaling"  ~ 53,
    CL_MOA ==  "DPF2_KO|Growth signaling"  ~ 54
  )) %<>%
  mutate(Compound_MOA_Class = case_when(
    Compound_MOA ==  "Cell cycle regulation"   ~ 1,
    Compound_MOA ==  "DNA damage" ~ 2,
    Compound_MOA == "Telomerase inhibitor"  ~ 3,
    Compound_MOA == "Protein homeostasis"  ~ 4,
    Compound_MOA ==  "Growth signaling"  ~ 5


  )) %<>%
  mutate(CL_Class = case_when(
    cell_line_modification ==  "ARID1A_KO"   ~ 1,
    cell_line_modification ==  "ARID1B_KO" ~ 2,
    cell_line_modification == "ARID2_KO"  ~ 3,
    cell_line_modification == "BRD7_KO"  ~ 4,
    cell_line_modification ==  "BRD9_KO"  ~ 5,
    cell_line_modification ==  "PBRM1_KO" ~ 6,
    cell_line_modification == "PHF10_KO"  ~ 7,
    cell_line_modification == "SMARCA4_KO"   ~ 8,
    cell_line_modification == "SMARCC1_KO"  ~ 9,
    cell_line_modification == "SMARCD1_KO"   ~ 10,
    cell_line_modification == "DPF2_KO"   ~ 11

  )) %<>%
  mutate(Domain_Class = case_when(
    Domain == "ARID" ~ 1,
    Domain == "Bromo" ~ 2,
    Domain ==  "PHD_Zinc" ~ 3,
    Domain == "SWIRM" ~ 4,
    Domain == "SWIB" ~5,
    Domain == "HSA" ~6))%<>%
  mutate(Domain_MOA_Class = case_when(
    Domain_MOA == "ARID|Cell cycle regulation" ~ 1,
    Domain_MOA == "Bromo|Cell cycle regulation" ~ 2,
    Domain_MOA ==  "PHD_Zinc|Cell cycle regulation" ~ 3,
    Domain_MOA == "SWIRM|Cell cycle regulation" ~ 4,
    Domain_MOA == "SWIB|Cell cycle regulation" ~5,
    Domain_MOA == "HSA|Cell cycle regulation" ~6,
    Domain_MOA == "ARID|DNA damage" ~ 7,
    Domain_MOA == "Bromo|DNA damage" ~ 8,
    Domain_MOA ==  "PHD_Zinc|DNA damage" ~ 9,
    Domain_MOA == "SWIRM|DNA damage" ~ 10,
    Domain_MOA == "SWIB|DNA damage" ~11,
    Domain_MOA == "HSA|DNA damage" ~12,
    Domain_MOA == "ARID|Telomerase inhibitor" ~ 13,
    Domain_MOA == "Bromo|Telomerase inhibitor" ~ 14,
    Domain_MOA ==  "PHD_Zinc|Telomerase inhibitor" ~ 15,
    Domain_MOA == "SWIRM|Telomerase inhibitor" ~ 16,
    Domain_MOA == "SWIB|Telomerase inhibitor" ~17,
    Domain_MOA == "HSA|Telomerase inhibitor" ~18,
    Domain_MOA == "ARID|Protein homeostasis" ~ 19,
    Domain_MOA == "Bromo|Protein homeostasis" ~ 20,
    Domain_MOA ==  "PHD_Zinc|Protein homeostasis" ~ 21,
    Domain_MOA == "SWIRM|Protein homeostasis" ~ 22,
    Domain_MOA == "SWIB|Protein homeostasis" ~23,
    Domain_MOA == "HSA|Protein homeostasis" ~24,
    Domain_MOA == "ARID|Growth signaling" ~ 25,
    Domain_MOA == "Bromo|Growth signaling" ~ 26,
    Domain_MOA ==  "PHD_Zinc|Growth signaling" ~ 27,
    Domain_MOA == "SWIRM|Growth signaling" ~ 28,
    Domain_MOA == "SWIB|Growth signaling" ~29,
    Domain_MOA == "HSA|Growth signaling" ~30
  )) %<>%
  mutate(Module_Class = case_when(
    Module ==  "accessory"  ~ 1,
    Module ==  "core"  ~ 2,
    Module ==  "atpase"  ~ 3))%<>%
  mutate(Module_MOA_Class = case_when(
    Module_MOA == "accessory|Cell cycle regulation" ~ 1,
    Module_MOA == "atpase|Cell cycle regulation" ~ 2,
    Module_MOA ==  "core|Cell cycle regulation" ~ 3,
    Module_MOA == "accessory|DNA damage" ~ 4,
    Module_MOA == "atpase|DNA damage" ~5,
    Module_MOA == "core|DNA damage" ~6,
    Module_MOA == "accessory|Telomerase inhibitor" ~ 7,
    Module_MOA == "atpase|Telomerase inhibitor" ~ 8,
    Module_MOA ==  "core|Telomerase inhibitor" ~ 9,
    Module_MOA == "accessory|Protein homeostasis" ~ 10,
    Module_MOA == "atpase|Protein homeostasis" ~11,
    Module_MOA == "core|Protein homeostasis" ~12,
    Module_MOA == "accessory|Growth signaling" ~ 13,
    Module_MOA == "atpase|Growth signaling" ~ 14,
    Module_MOA ==  "core|Growth signaling" ~ 15))%>%
  mutate(Module_Compound_Class = case_when(
    Module_Compound ==  "accessory|Adavosertib"  ~ 1,
    Module_Compound ==  "atpase|Adavosertib" ~ 2,
    Module_Compound == "core|Adavosertib" ~ 3,
    Module_Compound == "accessory|Aphidicolin"  ~ 4,
    Module_Compound == "atpase|Aphidicolin" ~ 5,
    Module_Compound == "core|Aphidicolin"  ~ 6,
    Module_Compound == "accessory|BIBR1532" ~ 7,
    Module_Compound == "atpase|BIBR1532"  ~ 8,
    Module_Compound == "core|BIBR1532"  ~ 9,
    Module_Compound == "accessory|BI-2536" ~ 10,
    Module_Compound == "atpase|BI-2536" ~ 11,
    Module_Compound == "core|BI-2536"    ~ 12,
    Module_Compound == "accessory|Camptothecin" ~ 13,
    Module_Compound == "atpase|Camptothecin"  ~ 14,
    Module_Compound == "core|Camptothecin"  ~ 15,
    Module_Compound ==  "accessory|Doxorubicin"  ~ 16,
    Module_Compound ==  "atpase|Doxorubicin" ~ 17,
    Module_Compound ==  "core|Doxorubicin"~ 18,
    Module_Compound ==  "accessory|Etoposide" ~ 19,
    Module_Compound == "atpase|Etoposide" ~ 20,
    Module_Compound == "core|Etoposide" ~ 21,
    Module_Compound =="accessory|Gemcitabine"  ~ 22,
    Module_Compound ==   "atpase|Gemcitabine" ~ 23,
    Module_Compound ==  "core|Gemcitabine"   ~ 24,
    Module_Compound == "accessory|HydroxyUrea" ~ 25,
    Module_Compound == "atpase|HydroxyUrea"  ~ 26,
    Module_Compound == "core|HydroxyUrea"   ~ 27,
    Module_Compound == "accessory|IlludinS" ~ 28,
    Module_Compound == "atpase|IlludinS" ~ 29,
    Module_Compound == "core|IlludinS" ~ 30,
    Module_Compound ==  "accessory|Lapatinib"   ~ 31,
    Module_Compound ==  "atpase|Lapatinib"  ~ 32,
    Module_Compound == "core|Lapatinib" ~ 33,
    Module_Compound == "accessory|MLN4924"   ~ 34,
    Module_Compound == "atpase|MLN4924"  ~ 35,
    Module_Compound == "core|MLN4924"~ 36,
    Module_Compound == "accessory|Methyl-Methanesulfonate"  ~ 37,
    Module_Compound == "atpase|Methyl-Methanesulfonate"  ~ 38,
    Module_Compound ==  "core|Methyl-Methanesulfonate"  ~ 39,
    Module_Compound == "accessory|Palbociclib" ~ 40,
    Module_Compound == "atpase|Palbociclib"  ~ 41,
    Module_Compound == "core|Palbociclib"   ~ 42,
    Module_Compound == "accessory|Paclitaxel"  ~ 43,
    Module_Compound == "atpase|Paclitaxel" ~ 44,
    Module_Compound == "core|Paclitaxel" ~ 45,
    Module_Compound ==  "accessory|Topotecan"    ~ 46,
    Module_Compound ==  "atpase|Topotecan"  ~ 47,
    Module_Compound == "core|Topotecan"   ~ 48,

    Module_Compound ==  "accessory|Cobimetinib"    ~ 49,
    Module_Compound ==  "atpase|Cobimetinib"  ~ 50,
    Module_Compound == "core|Cobimetinib"   ~ 51,

    Module_Compound ==  "accessory|Mirdametinib"    ~ 52,
    Module_Compound ==  "atpase|Mirdametinib"  ~ 53,
    Module_Compound == "core|Mirdametinib"   ~ 54,

    Module_Compound ==  "accessory|NMS-873"    ~ 55,
    Module_Compound ==  "atpase|NMS-873"  ~ 56,
    Module_Compound == "core|NMS-873"   ~ 57,

    Module_Compound ==  "accessory|Niraparib"    ~ 58,
    Module_Compound ==  "atpase|Niraparib"  ~ 59,
    Module_Compound == "core|Niraparib"   ~ 60,

    Module_Compound ==  "accessory|Sapitinib"    ~ 61,
    Module_Compound ==  "atpase|Sapitinib"  ~ 62,
    Module_Compound == "core|Sapitinib"   ~ 63
  ))%<>%
  mutate(Domain_Compound_Class = case_when(
    Domain_Compound ==  "ARID|Adavosertib"  ~ 1,
    Domain_Compound ==  "Bromo|Adavosertib" ~ 2,
    Domain_Compound == "PHD_Zinc|Adavosertib" ~ 3,
    Domain_Compound == "SWIRM|Adavosertib" ~ 4,
    Domain_Compound == "SWIB|Adavosertib" ~ 5,
    Domain_Compound == "HSA|Adavosertib" ~ 6,

    Domain_Compound == "ARID|Aphidicolin"  ~ 7,
    Domain_Compound == "Bromo|Aphidicolin" ~ 8,
    Domain_Compound == "PHD_Zinc|Aphidicolin"  ~ 9,
    Domain_Compound == "SWIRM|Aphidicolin"  ~ 10,
    Domain_Compound == "SWIB|Aphidicolin"  ~ 11,
    Domain_Compound == "HSA|Aphidicolin"  ~ 12,

    Domain_Compound == "ARID|BIBR1532" ~ 13,
    Domain_Compound == "Bromo|BIBR1532"  ~ 14,
    Domain_Compound == "PHD_Zinc|BIBR1532"  ~ 15,
    Domain_Compound == "SWIRM|BIBR1532"  ~ 16,
    Domain_Compound == "SWIB|BIBR1532"  ~ 17,
    Domain_Compound == "HSA|BIBR1532"  ~ 18,

    Domain_Compound == "ARID|BI-2536" ~ 19,
    Domain_Compound == "Bromo|BI-2536" ~ 20,
    Domain_Compound == "PHD_Zinc|BI-2536"    ~ 21,
    Domain_Compound == "SWIRM|BI-2536"    ~ 22,
    Domain_Compound == "SWIB|BI-2536"    ~ 23,
    Domain_Compound == "HSA|BI-2536"    ~ 24,

    Domain_Compound == "ARID|Camptothecin" ~ 25,
    Domain_Compound == "Bromo|Camptothecin"  ~ 26,
    Domain_Compound == "PHD_Zinc|Camptothecin"  ~ 27,
    Domain_Compound == "SWIRM|Camptothecin"  ~ 28,
    Domain_Compound == "SWIB|Camptothecin"  ~ 29,
    Domain_Compound == "HSA|Camptothecin"  ~ 30,

    Domain_Compound ==  "ARID|Doxorubicin"  ~ 31,
    Domain_Compound ==  "Bromo|Doxorubicin" ~ 32,
    Domain_Compound ==  "PHD_Zinc|Doxorubicin"~ 33,
    Domain_Compound ==  "SWIRM|Doxorubicin"~ 34,
    Domain_Compound ==  "SWIB|Doxorubicin"~ 35,
    Domain_Compound ==  "HSA|Doxorubicin"~ 36,

    Domain_Compound ==  "ARID|Etoposide" ~ 37,
    Domain_Compound == "Bromo|Etoposide" ~ 38,
    Domain_Compound == "PHD_Zinc|Etoposide" ~ 39,
    Domain_Compound == "SWIRM|Etoposide" ~ 40,
    Domain_Compound == "SWIB|Etoposide" ~ 41,
    Domain_Compound == "HSA|Etoposide" ~ 42,

    Domain_Compound =="ARID|Gemcitabine"  ~ 43,
    Domain_Compound ==   "Bromo|Gemcitabine" ~ 44,
    Domain_Compound ==  "PHD_Zinc|Gemcitabine"   ~ 45,
    Domain_Compound ==  "SWIRM|Gemcitabine"   ~ 46,
    Domain_Compound ==  "SWIB|Gemcitabine"   ~ 47,
    Domain_Compound ==  "HSA|Gemcitabine"   ~ 48,

    Domain_Compound == "ARID|HydroxyUrea" ~ 49,
    Domain_Compound == "Bromo|HydroxyUrea"  ~ 50,
    Domain_Compound == "PHD_Zinc|HydroxyUrea"   ~ 51,
    Domain_Compound == "SWIRM|HydroxyUrea"   ~ 52,
    Domain_Compound == "SWIB|HydroxyUrea"   ~ 53,
    Domain_Compound == "HSA|HydroxyUrea"   ~ 54,

    Domain_Compound == "ARID|IlludinS" ~ 55,
    Domain_Compound == "Bromo|IlludinS" ~ 56,
    Domain_Compound == "PHD_Zinc|IlludinS" ~ 57,
    Domain_Compound == "SWIRM|IlludinS" ~ 58,
    Domain_Compound == "SWIB|IlludinS" ~ 59,
    Domain_Compound == "HSA|IlludinS" ~ 60,

    Domain_Compound ==  "ARID|Lapatinib"   ~ 61,
    Domain_Compound ==  "Bromo|Lapatinib"  ~ 62,
    Domain_Compound == "PHD_Zinc|Lapatinib" ~ 63,
    Domain_Compound == "SWIRM|Lapatinib" ~ 64,
    Domain_Compound == "SWIB|Lapatinib" ~ 65,
    Domain_Compound == "HSA|Lapatinib" ~ 66,

    Domain_Compound == "ARID|MLN4924"   ~ 67,
    Domain_Compound == "Bromo|MLN4924"  ~ 68,
    Domain_Compound == "PHD_Zinc|MLN4924"~ 69,
    Domain_Compound == "SWIRM|MLN4924"~ 70,
    Domain_Compound == "SWIB|MLN4924"~ 71,
    Domain_Compound == "HSA|MLN4924"~ 72,

    Domain_Compound == "ARID|Methyl-Methanesulfonate"  ~ 73,
    Domain_Compound == "Bromo|Methyl-Methanesulfonate"  ~ 74,
    Domain_Compound ==  "PHD_Zinc|Methyl-Methanesulfonate"  ~ 75,
    Domain_Compound ==  "SWIRM|Methyl-Methanesulfonate"  ~ 76,
    Domain_Compound ==  "SWIB|Methyl-Methanesulfonate"  ~ 77,
    Domain_Compound ==  "HSA|Methyl-Methanesulfonate"  ~ 78,

    Domain_Compound == "ARID|Palbociclib" ~ 79,
    Domain_Compound == "Bromo|Palbociclib"  ~ 80,
    Domain_Compound == "PHD_Zinc|Palbociclib"   ~ 81,
    Domain_Compound == "SWIRM|Palbociclib"   ~ 82,
    Domain_Compound == "SWIB|Palbociclib"   ~ 83,
    Domain_Compound == "HSA|Palbociclib"   ~ 84,

    Domain_Compound == "ARID|Paclitaxel"  ~ 85,
    Domain_Compound == "Bromo|Paclitaxel" ~ 86,
    Domain_Compound == "PHD_Zinc|Paclitaxel" ~ 87,
    Domain_Compound == "SWIRM|Paclitaxel" ~ 88,
    Domain_Compound == "SWIB|Paclitaxel" ~ 89,
    Domain_Compound == "HSA|Paclitaxel" ~ 90,

    Domain_Compound ==  "ARID|Topotecan"    ~ 91,
    Domain_Compound ==  "Bromo|Topotecan"  ~ 92,
    Domain_Compound == "PHD_Zinc|Topotecan"   ~ 93,
    Domain_Compound == "SWIRM|Topotecan"   ~ 94,
    Domain_Compound == "SWIB|Topotecan"   ~ 95,
    Domain_Compound == "HSA|Topotecan"   ~ 96,

    Domain_Compound ==  "ARID|Cobimetinib"    ~ 97,
    Domain_Compound ==  "Bromo|Cobimetinib"  ~ 98,
    Domain_Compound == "PHD_Zinc|Cobimetinib"   ~ 99,
    Domain_Compound == "SWIRM|Cobimetinib"   ~ 100,
    Domain_Compound == "SWIB|Cobimetinib"   ~ 101,
    Domain_Compound == "HSA|Cobimetinib"   ~ 102,

    Domain_Compound ==  "ARID|Mirdametinib"    ~ 103,
    Domain_Compound ==  "Bromo|Mirdametinib"  ~ 104,
    Domain_Compound == "PHD_Zinc|Mirdametinib"   ~ 105,
    Domain_Compound == "SWIRM|Mirdametinib"   ~ 106,
    Domain_Compound == "SWIB|Mirdametinib"   ~ 107,
    Domain_Compound == "HSA|Mirdametinib"   ~ 108,

    Domain_Compound ==  "ARID|NMS-873"    ~ 109,
    Domain_Compound ==  "Bromo|NMS-873"  ~ 110,
    Domain_Compound == "PHD_Zinc|NMS-873"   ~ 111,
    Domain_Compound == "SWIRM|NMS-873"   ~ 112,
    Domain_Compound == "SWIB|NMS-873"   ~ 113,
    Domain_Compound == "HSA|NMS-873"   ~ 114,

    Domain_Compound ==  "ARID|Niraparib"    ~ 115,
    Domain_Compound ==  "Bromo|Niraparib"  ~ 116,
    Domain_Compound == "PHD_Zinc|Niraparib"   ~ 117,
    Domain_Compound == "SWIRM|Niraparib"   ~ 118,
    Domain_Compound == "SWIB|Niraparib"   ~ 119,
    Domain_Compound == "HSA|Niraparib"   ~ 120,

    Domain_Compound ==  "ARID|Sapitinib"    ~ 121,
    Domain_Compound ==  "Bromo|Sapitinib"  ~ 122,
    Domain_Compound == "PHD_Zinc|Sapitinib"   ~ 123,
    Domain_Compound == "SWIRM|Sapitinib"   ~ 124,
    Domain_Compound == "SWIB|Sapitinib"   ~ 125,
    Domain_Compound == "HSA|Sapitinib"   ~ 126
  ))%<>%
  dplyr::rename(Compound = treatment)

################################################################################

################################ Running Clustering Methods across the 5 Phenotypes ################################

#set range of K clusters to test
KRange <- 2:10

# specify data labels to use for external cluster validation
LabelTypesList <- c("Subtype_MOA_Class", "Subtype_Compound_Class", "Module_MOA_Class","Module_Compound_Class", "Domain_MOA_Class","Domain_Compound_Class","CL_MOA_Class","Subtype_Class","Domain_Class",
                    "Module_Class","Compound_MOA_Class","Compound_Class","CL_Class")

# KO & WT labels for heatmap
kolabel<-"KO"
wtlabel <- "WT"

# text size for heatmap
rowtextsize<-2

set.seed(123)

UnsupervisedResultsDF <- data.frame(matrix(ncol=10,nrow=0))
colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Label_Type","CRI","MVI")

for(K in KRange){

  # run K-means clustering
  RowInput <- runKMeans(AUCDF5PhenotypesDF[,1:5],K,AUCDF5PhenotypesDF,LabelTypesList)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)
  # }

  # run PAM clustering
  RowInput <- runPAM(AUCDF5PhenotypesDF[,1:5], K, AUCDF5PhenotypesDF,LabelTypesList)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)

  #run Hierarchical clustering
  RowInput <- runHierarchicalClustering(AUCDF5PhenotypesDF[,1:5], K, AUCDF5PhenotypesDF,LabelTypesList)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)

  #run C-means clustering
  RowInput <- runCMeans(AUCDF5PhenotypesDF[,1:5],K,AUCDF5PhenotypesDF,LabelTypesList)

  UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)


}

################################ Comparing Performance of Tested Clustering Methods ################################

#creating necessary directories for storing cluster method comparisons plots
SubDir <- sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/",DirectoryName,AnalysisName)

if (file.exists(SubDir)){
  print("Comparing clustering algorithms")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("Comparing clustering algorithms")
}

SubDir <- sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/C-means/",DirectoryName,AnalysisName)

if(isTRUE(file.exists(SubDir))==TRUE){
  print("...")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("...")
}

SubDir <- sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/",DirectoryName,AnalysisName)

if(isTRUE(file.exists(SubDir))==TRUE){
  print("...")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("...")
}

SubDir <- sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/",DirectoryName,AnalysisName)

if(isTRUE(file.exists(SubDir))==TRUE){
  print("...")
} else {
  dir.create(SubDir, recursive = TRUE)
  print("...")
}


# run performance comparison across methods tested on specified label types

UnsupervisedResultsDF <- as.data.frame(UnsupervisedResultsDF)

for(Label in unique(UnsupervisedResultsDF$Label_Type)){

  #Compare performance of the tested clustering algorithms w/ specified number of clusters for one given label type
  ListOfUnsupervisedMethodComps <- runUnsupervisedAlgorithmComparison(AUCDF5PhenotypesDF[,1:5], UnsupervisedResultsDF, Label)

  # save plots for comparing all methods
  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_silhouette2_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_dunn_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_connectivity_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_connectivity_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_ch_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_wss_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$WSS,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


  # save plots for comparing k clusters using hierarchical clustering

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/%s_hc_silhouette_%s.pdf",DirectoryName,AnalysisName, TodaysDate, AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$HC$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/%s_hc_dunn_%s.pdf",DirectoryName,AnalysisName, TodaysDate, AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$HC$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/%s_hc_connectivity_%s.pdf",DirectoryName,AnalysisName,TodaysDate, AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$HC$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/%s_hc_ch_%s.pdf",DirectoryName, AnalysisName,TodaysDate, AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$HC$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/%s_hc_wss_%s.pdf",DirectoryName,AnalysisName,TodaysDate, AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$HC$WSS,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


  # save plots for comparing k clusters using k-means clustering
  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/%s_km_silhouette_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$KM$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/%s_km_dunn_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$KM$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/%s_km_connectivity_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$KM$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/%s_km_ch_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$KM$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/%s_km_wss_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$KM$WSS,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


  #save plots for comparing k clusters using PAM clustering
  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/%s_pam_silhouette_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$PAM$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/%s_PAM_dunn_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$PAM$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/%s_pam_connectivity_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$PAM$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/%s_pam_ch_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$PAM$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/%s_pam_wss_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$PAM$WSS,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  # save plots for comparing k clusters using c-means clustering

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/C-means/%s_cm_silhouette_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$CM$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/C-means/%s_cm_dunn_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$CM$Dunn,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/C-means/%s_cm_connectivity_%s.pdf",DirectoryName,AnalysisName, TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$CM$Connectivity,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/C-means/%s_cm_ch_%s.pdf",DirectoryName, AnalysisName,TodaysDate,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$CM$CH,width = 15, height = 10, dpi = 150, units = "in", device='pdf')


  #saving all plots demonstrating accuracy of label types

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_%s_corr_rand_%s.pdf",DirectoryName,AnalysisName, TodaysDate, Label, AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$LabelAccuracy$CRI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/%s_%s_mvi_%s.pdf",DirectoryName,AnalysisName,TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$AllMethods$LabelAccuracy$MVI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/%s_%s_hc_corr_rand_%s.pdf",DirectoryName,AnalysisName, TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$HC$LabelAccuracy$CRI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/Hierarchical/%s_%s_hc_mvi_%s.pdf",DirectoryName, AnalysisName,TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$HC$LabelAccuracy$MVI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/%s_%s_km_corr_rand_%s.pdf",DirectoryName,AnalysisName, TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$KM$LabelAccuracy$CRI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/K-means/%s_%s_km_mvi_%s.pdf",DirectoryName,AnalysisName, TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$KM$LabelAccuracy$MVI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/%s_%s_pam_corr_rand_%s.pdf",DirectoryName,AnalysisName, TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$PAM$LabelAccuracy$CRI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/PAM/%s_%s_pam_mvi_%s.pdf",DirectoryName, AnalysisName,TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$PAM$LabelAccuracy$MVI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/C-means/%s_%s_cm_corr_rand_%s.pdf",DirectoryName,AnalysisName, TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$CM$LabelAccuracy$CRI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

  NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Algorithm_Comparison/%s/C-means/%s_%s_cm_mvi_%s.pdf",DirectoryName,AnalysisName, TodaysDate, Label,AnalysisName)

  ggsave(filename=NewFilename, plot=ListOfUnsupervisedMethodComps$CM$LabelAccuracy$MVI,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

}


################################ Plotting Heatmap clustered across 5 Phenotypes ################################

# setting specific color range for plots
scaling <- "scaled"
minColorRange <- -0.5
maxColorRange <- 0.5

columnk <- 0

plotlist<- getAUCSummaryHeatmapClustered(AUCDF5PhenotypesDF[,1:5],AUCDF5PhenotypesDF,PValue5PhenotypesDF,ListOfUnsupervisedMethodComps$bestClustMethod,BAFSubtypeDF,CompoundMOADF,SubtypeColors,MOAColors,kolabel,wtlabel,
                                         columnk,rowtextsize,scaling,minColorRange,maxColorRange)

#creating necessary directories to store Cmeans-related plots
SubDir <- sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Heatmaps/",DirectoryName,ListOfUnsupervisedMethodComps$bestClustMethod$Cluster_Num)

if(isTRUE(file.exists(SubDir))==TRUE){
  print(paste("Plotting",ListOfUnsupervisedMethodComps$bestClustMethod$Method,"clustered heatmap with" ,ListOfUnsupervisedMethodComps$bestClustMethod$Cluster_Num, "clusters"))
} else {
  dir.create(SubDir, recursive = TRUE)
  print(paste("Plotting",ListOfUnsupervisedMethodComps$bestClustMethod$Method,"clustered heatmap with" ,ListOfUnsupervisedMethodComps$bestClustMethod$Cluster_Num, "clusters"))
}

# save clustered heatmap
NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Heatmaps/%s_%s_%s_heatmap_%s.pdf",DirectoryName, TodaysDate, ListOfUnsupervisedMethodComps$bestClustMethod$Method,
                      ListOfUnsupervisedMethodComps$bestClustMethod$Cluster_Num,AnalysisName)

pdf(NewFilename, width = 8, height = 8)

draw(plotlist$Heatmap,merge_legend = TRUE, column_dend_side = "top") #
dev.off()

#save average silhouette plot
NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Heatmaps/%s_%s_%s_avg_sil_%s.pdf",DirectoryName, TodaysDate, ListOfUnsupervisedMethodComps$bestClustMethod$Method,
                      ListOfUnsupervisedMethodComps$bestClustMethod$Cluster_Num,AnalysisName)

ggsave(filename=NewFilename, plot=plotlist$AvgSil,width = 15, height = 10, dpi = 150, units = "in", device='pdf')

# save scatter plot
NewFilename = sprintf("./Modeling-independent_Phenotypes/%s/Unsupervised_Clustering/Heatmaps/%s_%s_%s_scatter_plot_%s.pdf",DirectoryName, TodaysDate, ListOfUnsupervisedMethodComps$bestClustMethod$Method,
                      ListOfUnsupervisedMethodComps$bestClustMethod$Cluster_Num,AnalysisName)

ggsave(filename=NewFilename, plot=plotlist$ClusterViz,width = 15, height = 10, dpi = 150, units = "in", device='pdf')





