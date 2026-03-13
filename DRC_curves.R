#!/usr/bin/env RScript

#command-line executable:  chmod +x DRC_curves.R
# run as: ./DRC_curves.R

# Script for generating dose-response curves (DRC) normalized to a control treatment (i.e. DMSO) for all modeling-independent phenotypes averaged across all biological replicates per condition.
# Here the computed ratio is compound/DMSO for each cell line in order to understand the effect of knocking out a BAF subunit has on the treated state.

TodaysDate <- ""
AnalysisName<-"Systematic"

################################################## Loading Necessary Libraries and Source Code #########################################################################
source("./R/analysis_group_functions.R")
source("./R/plot_functions.R")
source("./R/phenotype_quantification_functions.R")

library(tidyr)
library(emmeans)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(magrittr)
library("RMySQL")
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(Hmisc)
library(gdata)
library(ggpubr)
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
library(lme4)
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
# labeling of bio and technical replicates by combining meta data with raw confluency data
IncucyteDataAndMetaData <- IncucyteDataDF %>% left_join(IncucyteMetaDataDF, by = c("project_ID","plate_ID","well_ID"))

# join condition df to previous dataframe in order to include cell line, treatments, etc info
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaData %>% left_join(IncucyteConditionDF, by = c("condition_ID"))

# remove outer wells/outliers
IncucyteDataAndMetaDataAndConditionsDF= subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$exclude==0)

# generate bio_rep column:join condition_ID, project_ID and plate_ID
IncucyteDataAndMetaDataAndConditionsDF$bio_rep_id<- paste(IncucyteDataAndMetaDataAndConditionsDF$condition_ID, IncucyteDataAndMetaDataAndConditionsDF$project_ID, IncucyteDataAndMetaDataAndConditionsDF$plate_ID, sep = '|')

# create bio_rep_id: join condition_ID, project_ID and plate_ID = get bio rep for each condition
IncucyteSplitPointDF$bio_rep_id<- paste(IncucyteSplitPointDF$condition_ID, IncucyteSplitPointDF$project_ID, IncucyteSplitPointDF$plate_ID, sep = '|')

# get phenotypes of growth curve by merging IncucyteSplitPointDF
IncucyteDataAndMetaDataAndConditionsDF <- merge(IncucyteDataAndMetaDataAndConditionsDF, IncucyteSplitPointDF, by = c("bio_rep_id","condition_ID","plate_ID"), all.x = TRUE, all.y = TRUE)


# filter out certain cell lines
IncucyteDataAndMetaDataAndConditionsDF <- subset(IncucyteDataAndMetaDataAndConditionsDF, !(cell_line_modifications %in% c("SMARCD1_KO","DPF2_KO","BRD7_KO","PBRM1_KO", "PHF10_KO")))

# create treatment label
IncucyteDataAndMetaDataAndConditionsDF <- IncucyteDataAndMetaDataAndConditionsDF %>% mutate(treatment = case_when(is.na(treatment_0) ~ "Untreated",
                                                                                                                  is.na(treatment_1) & is.na(treatment_2) ~ treatment_0,
                                                                                                                  !(is.na(treatment_1)) & is.na(treatment_2) ~ paste(treatment_0, treatment_1, sep='_'),
                                                                                                                  !(is.na(treatment_1)) & !(is.na(treatment_2)) ~ paste(treatment_0, treatment_1,treatment_2, sep='_')))


# create treat_conc_id summarizing compound and concentration
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

# generate analysis group ID for KO treated, WT treated, KO untreated, WT untreated for one given compound and concentration
IncucyteDataAndMetaDataAndConditionsLabeledDF <- getSamePlatesPerComparisonGroup(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDMSODilutionDF,ListOfSpecialCompounds = ListOfSpecialCompounds,wildtype = wildtype)


################################ Calculate Treatment:DMSO ratio for each phenotype ################################
ListOfSpecialCompounds <- c("HydroxyUrea")
controlcompound <- "DMSO"

# normalize phenotypes to control treatment
PhenotypesNormToDMSODF <- normalizeToDMSO(IncucyteDataAndMetaDataAndConditionsLabeledDF,controlcompound,ListOfSpecialCompounds)

# rename cell line
PhenotypesNormToDMSODF <- PhenotypesNormToDMSODF %>%
  mutate(cell_line_modification = replace(cell_line_modification, cell_line_modification == "C631", "WT"))

# update wildtype variable
wildtype<- "WT"

# filter for DRC of bioreps for a compound with more than two concentrations
PhenotypesNormToDMSODF <- PhenotypesNormToDMSODF %>% dplyr::group_by(treatment_0,cell_line_modification, plate_ID) %>% dplyr::filter(n()>2)

################################ Plot Dose-response curves normalized to DMSO ################################
# Prepare Treatment:DMSO ratio df to be plotted as dose-response curve

# average ratio per condition
PhenotypesNormToDMSODF<- PhenotypesNormToDMSODF %>%
  dplyr::group_by(condition_ID, cell_line_modification, treatment_0, concentration_0,treatment,CL_treat_conc_id,treat_conc_id,NLMM_Analysis_ID) %>%
  dplyr::summarise_at(c("Confluency_after_24h_ratio","Confluency_after_48h_ratio","Confluency_after_72h_ratio","Confluency_after_96h_ratio","Second_Phase_Relative_Change_in_Confluency_ratio",
                        "First_Phase_Duration_ratio","First_Phase_Change_in_Confluency_ratio"), list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(length(.))), na.rm=TRUE)

#find KO cell line each WT is paired with based on NLMM_Analysis_ID
PhenotypesNormToDMSODFKO <- subset(PhenotypesNormToDMSODF, cell_line_modification!="WT")

namevector <- "KO_partner"
PhenotypesNormToDMSODFKO[ , namevector] <- NA

PhenotypesNormToDMSODFWT <- subset(PhenotypesNormToDMSODF, cell_line_modification=="WT")

PhenotypesNormToDMSODFWT[ , namevector] <- NA

SharedNLMMAnalysisIDDF <- merge(PhenotypesNormToDMSODFKO,PhenotypesNormToDMSODFWT, by = "NLMM_Analysis_ID")

SharedNLMMAnalysisIDList <-SharedNLMMAnalysisIDDF$NLMM_Analysis_ID

#subset df by NLMM analysis IDs that have both WT and KO
PhenotypesNormToDMSODF2 <- subset(PhenotypesNormToDMSODF, NLMM_Analysis_ID %in% SharedNLMMAnalysisIDList)

# create KO_partner label to link KO with its respective WT
NLMM_groups <- PhenotypesNormToDMSODF2 %>%
  dplyr::arrange(factor(cell_line_modification, levels= c("ARID1A_KO","ARID1B_KO","ARID2_KO", "BRD9_KO","SMARCA4_KO","SMARCC1_KO", "WT"))) %>%
  dplyr::group_by(NLMM_Analysis_ID) %>%
  dplyr::summarise(KO_partner = paste0(cell_line_modification, collapse="|"))

PhenotypesNormToDMSODF <- merge(PhenotypesNormToDMSODF2, NLMM_groups , by=c("NLMM_Analysis_ID"), all.x = TRUE, all.y = TRUE)


#Create Metadata for Dose-response curves normalized to DMSO
################################################################################

# create metadata mapping compound to mechanism of action
Compound <- unique(PhenotypesNormToDMSODF$treatment)

CompoundMOADF <- data.frame(Compound)

CompoundMOADF <- CompoundMOADF %>% mutate(MOA=case_when(Compound %in% c("Gemcitabine","Lapatinib","Palbociclib") ~ "G1/S phase",
                                                        Compound %in% c("Aphidicolin","HydroxyUrea") ~ "S phase",
                                                        Compound %in% c("Adavosertib","BI-2536","Paclitaxel") ~ "G2/M phase",
                                                        Compound %in% c("Camptothecin","Topotecan") ~ "TopI inhibitor",
                                                        Compound %in% c("Doxorubicin","Etoposide") ~ "TopII inhibitor",
                                                        Compound %in% c("IlludinS","Methyl-Methanesulfonate") ~ "DNA Alkylating",
                                                        Compound == "MLN4924" ~ "NAE inhibitor",
                                                        Compound =="BIBR1532" ~ "Telomerase inhibitor",
                                                        Compound %in% c("Cobimetinib","Mirdametinib") ~ "MEK inhibitor",
                                                        Compound =="NMS-873" ~ "VCP/p97 inhibitor",
                                                        Compound == "Niraparib" ~ "PARP inhibitor",
                                                        Compound == "Sapitinib" ~ "EGFR inhibitor"))

# create metadata mapping cell lines to BAF subtype
CL <- unique(PhenotypesNormToDMSODF$cell_line_modification)

BAFSubtypeDF <- data.frame(CL)

BAFSubtypeDF <- BAFSubtypeDF %>% mutate(Subtype=case_when(
  CL %in% c("SMARCA4_KO","SMARCC1_KO") ~ "general BAF",
  CL %in% c("ARID1A_KO", "ARID1B_KO") ~ "cBAF",
  CL %in% c("ARID2_KO") ~ "PBAF",
  CL =="BRD9_KO" ~ "ncBAF")) %>%
  mutate(CL_name=gsub("_", " ", CL))

BAFSubtypeDF <- subset(BAFSubtypeDF, CL != wildtype)

# update cell line names to match all.group.lines
BAFSubtypeDF <- BAFSubtypeDF %>% mutate(CL_name=CL) %>%
  mutate(CL=gsub(" ", "_", CL_name)) %>%
  mutate(CL_name=gsub("_", " ", CL_name))

BAFSubtypeDF$CL_name <- factor(BAFSubtypeDF$CL_name,
                               levels=c("SMARCA4 KO","SMARCC1 KO","ARID1A KO","ARID1B KO","ARID2 KO", "BRD9 KO"))

# create KO partner label
BAFSubtypeDF <- BAFSubtypeDF %>%mutate(KO_partner_label=paste(CL,"|WT")) %>% mutate(KO_partner_label=gsub(" ", "", KO_partner_label))%>%mutate(KO_partner_name=gsub("_", " ", KO_partner_label))


BAFSubtypeDF <- BAFSubtypeDF  %>% mutate(KO_partner_label = fct_relevel(KO_partner_label,
                                                                        c("SMARCA4_KO|WT","SMARCC1_KO|WT","ARID1A_KO|WT","ARID1B_KO|WT","ARID2_KO|WT","BRD9_KO|WT")))%>% dplyr::arrange(KO_partner_label)

rownames(BAFSubtypeDF) <- 1:nrow(BAFSubtypeDF)
PhenotypesNormToDMSODF$KO_partner <- factor(PhenotypesNormToDMSODF$KO_partner,
                                            levels=c("SMARCA4_KO|WT","SMARCC1_KO|WT","ARID1A_KO|WT", "ARID1B_KO|WT","ARID2_KO|WT", "BRD9_KO|WT"))

# list mapping cell line to color for plot
all.group.colors <- c("SMARCA4 KO"= "#80519B","SMARCC1 KO" = "#3569A1","ARID1A KO"= "#C52030", "ARID1B KO" = "#8A181A","ARID2 KO"= "#176533","BRD9 KO"= "#D2AE2A","WT"="darkgrey")


# list mapping cell line to line type for plot
all.group.lines <- c("WT" = "solid","SMARCA4 KO"= "solid", "SMARCC1 KO" = "solid","ARID1A KO"= "solid", "ARID1B KO" = "solid","ARID2 KO"= "solid","BRD9 KO"= "solid")


# list of compounds to parse through
ListOfCompounds <- unique(PhenotypesNormToDMSODF$treatment)

# list of phenotypes to parse through
ListOfPhenotypes <- c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency","First_Phase_Duration","First_Phase_Change_in_Confluency")

# list of BAF subtypes to parse through
BAFSubtypeDFNoWT2 <- subset(BAFSubtypeDF, CL!=wildtype)

BAFSubtypeDFNoControl <- subset(BAFSubtypeDF, Subtype!="control")

SubtypeList <- unique(BAFSubtypeDFNoControl$Subtype)

# not setting specific y-axis range for plots
scaling = "unscaled"

################################################################################

for(k in 1:length(ListOfCompounds)){

  # create treatment label for plot
  treatmentlabel<- ListOfCompounds[k]

  # create control treatment label for plot
  if(isTRUE(treatmentlabel=="HydroxyUrea")==TRUE){
    controllabel<-"Untreated"
  }else{controllabel<-"DMSO"}

  for(n in 1:length(ListOfPhenotypes)){

    for(l in 1:length(SubtypeList)){

      #make sure directory exists - create if it doesn't
      SubDir <- sprintf("./Modeling-independent_Phenotypes/Growth_Curve_Phenotypes/%s_%s/Dose_response_curves/%s/BAF_Subtype_Specific/%s/",n,ListOfPhenotypes[n],
                        ListOfCompounds[k],gsub(" ", "_",SubtypeList[l]))

      if(isTRUE(file.exists(SubDir))==TRUE){
        print(paste("Plotting dose-response curves of",ListOfPhenotypes[n],"for",SubtypeList[l], "treated with", unique(ListOfCompounds[k]),"normalized to DMSO"))
      } else {
        dir.create(SubDir, recursive = TRUE)
        print(paste("Plotting dose-response curves of",ListOfPhenotypes[n],"for",SubtypeList[l], "treated with", unique(ListOfCompounds[k]),"normalized to DMSO"))
      }

      # plot averaged dose-response curve per KO cell line of a given BAF subtype along with their paired WT for a given phenotype
      DRCPlotObject<- plotDoseResponseCurves(PhenotypesNormToDMSODF,BAFSubtypeDF,all.group.colors,all.group.lines,wildtype,compound=ListOfCompounds[k],treatmentlabel,controllabel,phenotype=ListOfPhenotypes[n],subtype=SubtypeList[l])

      #save plot
      NewFilename = sprintf("./Modeling-independent_Phenotypes/Growth_Curve_Phenotypes/%s_%s/Dose_response_curves/%s/BAF_Subtype_Specific/%s/%s_%s_%s_%s_%s.pdf",n,ListOfPhenotypes[n],ListOfCompounds[k],gsub(" ", "_",SubtypeList[l]), TodaysDate,
                            ListOfCompounds[k],ListOfPhenotypes[n],AnalysisName, scaling)

      tryCatch(

        expr = {

          ggsave(filename=NewFilename, plot=DRCPlotObject,width=3, height=3, dpi=300, units = "in", device='pdf')

        },error=function(e){
          return(print(paste0(unique(ListOfCompounds[k],ListOfPhenotypes[n], "Plot for",SubtypeList[l],"is empty"))))
        },warning = function(w) {
          ggsave(filename=NewFilename, plot=DRCPlotObject,width = 20, height = 15, dpi = 300, units = "in", device='pdf')
        }
      )
    }
  }
}
