#' Assigning Identifiers for Analysis Groups - Quality Control

#'
#' TThis function assigns an analysis ID (i.e. nlmm_analysis_1) to every analysis group: KO + treatment, KO + control, WT + treatment, and WT + control, where the goal is to compare compound treatment to control treatment and
#' knockout (KO) cell line to the wildtype (WT) cell line. This version is used before plotting the growth curves of all members of one analysis group on the same graph during the quality control step of the pipeline.
#' @param IncucyteDataAndMetaDataAndConditionsDF Dataframe that combines Incucyte Confluency Dataset, Meta Data and Condition ID table.
#' @param IncucyteDMSODilutionDF Dataframe that contains the contents of the Incucyte DMSO Dilution Table. This is used for finding the respective control treatment and concentration for a given compound.
#' @param ListOfSpecialCompounds List of compounds that are normalized to untreated cells rather than DMSO-treated cells.
#' @param ListofNewClones List of cell line modifications that use a different DMSO concentration.
#' @param wildtype Variable representing name of wildtype cell line.
#' @return
#'\item{IncucyteDataAndMetaDataAndConditionsLabeledDF}{ Dataframe that combines Incucyte Confluency Dataset, Meta Data, Condition IDs and analysis IDs.}
#' @author Caroline Barry
#' @import dplyr
#' @export

getSamePlatesPerComparisonGroupQC <- function(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDMSODilutionDF,ListOfSpecialCompounds = NULL,ListofNewClones=NULL, wildtype){

  IncucyteDataAndMetaDataAndConditionsLabeledDF <- data.frame()

  i = 0

  IncucyteDataAndMetaDataAndConditionsDFNoWTNoDMSONoUntreated <- subset(IncucyteDataAndMetaDataAndConditionsDF, !(cell_line_modifications == wildtype) & treatment_0!="DMSO" & treatment!="Untreated")

  ListOfConditionID = unique(IncucyteDataAndMetaDataAndConditionsDFNoWTNoDMSONoUntreated$condition_ID)

  ListOfDataFrames <- vector(mode = "list", length = length(ListOfConditionID))

  for(k in 1:length(ListOfConditionID)){

    # get that KO TreatDF
    IncucyteKOTreatDatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$condition_ID==ListOfConditionID[k])

    # get that WT TreatDF
    IncucyteWTTreatDatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, cell_line_modifications==wildtype & treatment_0!="DMSO" & treatment!="Untreated")

    IncucyteWTTreatDatasetSubsetOriginal <- subset(IncucyteWTTreatDatasetSubsetOriginal, treatment_0 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_0 & concentration_0 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_0 &treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1 & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & plate_ID %in% IncucyteKOTreatDatasetSubsetOriginal$plate_ID)

    # get that KO DMSO DF
    CompoundName <-unique(IncucyteKOTreatDatasetSubsetOriginal$treatment_0)

    if(isTRUE(CompoundName %in% ListOfSpecialCompounds)==TRUE){

      # need to get untreated cell line and normalize treated to that
      # now subset by all_plates column
      IncucyteKODMSODatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, cell_line_modifications %in% IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications & cell_number %in% IncucyteKOTreatDatasetSubsetOriginal$cell_number & treatment =="Untreated"& treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1 & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1)
      IncucyteWTDMSODatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, cell_line_modifications ==wildtype & cell_number %in% IncucyteKOTreatDatasetSubsetOriginal$cell_number & treatment =="Untreated" & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1 & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1)

    }
    if(isTRUE(CompoundName %in% ListOfSpecialCompounds)==FALSE){
      DMSOConcentration <- subset(IncucyteDMSODilutionDF, IncucyteDMSODilutionDF$Compound_name %in%IncucyteKOTreatDatasetSubsetOriginal$treatment_0)
      if(isTRUE(grepl("KO_[0-9]",unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications))|unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications) %in% ListofNewClones)==TRUE){

        DMSOTreatedDFSubset <- subset(IncucyteDataAndMetaDataAndConditionsDF, concentration_0%in%DMSOConcentration$DMSO_concentration_new_clones_0 | concentration_0%in%DMSOConcentration$DMSO_concentration_new_clones_1 | concentration_0%in%DMSOConcentration$DMSO_concentration_new_clones_2)

        DMSOTreatedDFSubset <- subset(DMSOTreatedDFSubset, treatment_0=="DMSO" & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1)
      }
      if(isTRUE(grepl("KO_[0-9]",unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications))|unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications) %in% ListofNewClones)==FALSE){
        DMSOTreatedDFSubset <- subset(IncucyteDataAndMetaDataAndConditionsDF, concentration_0%in%DMSOConcentration$DMSO_concentration_0 & treatment_0=="DMSO"& concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1)
      }
      IncucyteKODMSODatasetSubsetOriginal <- subset(DMSOTreatedDFSubset, cell_line_modifications %in% IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications & plate_ID %in% IncucyteKOTreatDatasetSubsetOriginal$plate_ID)
      IncucyteWTDMSODatasetSubsetOriginal <- subset(DMSOTreatedDFSubset, cell_line_modifications==wildtype& plate_ID %in% IncucyteKOTreatDatasetSubsetOriginal$plate_ID)}

    IncucyteDataAndMetaDataAndConditionsKOWTTreatDF <- merge(IncucyteKOTreatDatasetSubsetOriginal, IncucyteWTTreatDatasetSubsetOriginal, by=c("bio_rep_id", "timepoint","condition_ID","cell_line", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0","treatment_timepoint_0","treatment_duration_0",
                                                                                                                                              "treatment_1", "concentration_1","treatment_timepoint_1","treatment_duration_1","treatment_2", "concentration_2","treatment_timepoint_2","treatment_duration_2",
                                                                                                                                              "treatment_3", "concentration_3","treatment_timepoint_3","treatment_duration_3","confluency","exclude", "project_ID", "plate_ID","well_ID", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
                                                                                                                                              "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
                                                                                                                                              "split_first_conf_difference", "phenotype","treatment","CL_treat_conc_id","treat_conc_id"), all.x = TRUE, all.y = TRUE)

    IncucyteDataAndMetaDataAndConditionsKOWTTreatDMSODF <- merge(IncucyteDataAndMetaDataAndConditionsKOWTTreatDF, IncucyteKODMSODatasetSubsetOriginal, by=c("bio_rep_id", "timepoint","condition_ID","cell_line", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0","treatment_timepoint_0","treatment_duration_0",
                                                                                                                                                            "treatment_1", "concentration_1","treatment_timepoint_1","treatment_duration_1","treatment_2", "concentration_2","treatment_timepoint_2","treatment_duration_2",
                                                                                                                                                            "treatment_3", "concentration_3","treatment_timepoint_3","treatment_duration_3","confluency","exclude", "project_ID", "plate_ID","well_ID", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
                                                                                                                                                            "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
                                                                                                                                                            "split_first_conf_difference", "phenotype","treatment","CL_treat_conc_id","treat_conc_id"), all.x = TRUE, all.y = TRUE)

    IncucyteDataAndMetaDataAndConditionsCompGroupDF <- merge(IncucyteDataAndMetaDataAndConditionsKOWTTreatDMSODF, IncucyteWTDMSODatasetSubsetOriginal, by=c("bio_rep_id", "timepoint","condition_ID","cell_line", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0","treatment_timepoint_0","treatment_duration_0",
                                                                                                                                                            "treatment_1", "concentration_1","treatment_timepoint_1","treatment_duration_1","treatment_2", "concentration_2","treatment_timepoint_2","treatment_duration_2",
                                                                                                                                                            "treatment_3", "concentration_3","treatment_timepoint_3","treatment_duration_3","confluency","exclude", "project_ID", "plate_ID","well_ID", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
                                                                                                                                                            "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
                                                                                                                                                            "split_first_conf_difference","phenotype","treatment","CL_treat_conc_id","treat_conc_id"), all.x = TRUE, all.y = TRUE)


    IncucyteDataAndMetaDataAndConditionsCompGroupDF <- IncucyteDataAndMetaDataAndConditionsCompGroupDF[,-c(2,23,24)]


    ConditionsCheckDF <- IncucyteDataAndMetaDataAndConditionsCompGroupDF%>%
      group_by(plate_ID) %>%
      distinct(condition_ID) %>%
      dplyr::filter(n()>3)

    if(isTRUE(nrow(ConditionsCheckDF)>0)==TRUE){

      IncucyteDataAndMetaDataAndConditionsCompGroupDF <- subset(IncucyteDataAndMetaDataAndConditionsCompGroupDF, IncucyteDataAndMetaDataAndConditionsCompGroupDF$plate_ID %in% ConditionsCheckDF$plate_ID)

      # assign NLMM Analysis ID

      ConditionValue <- as.numeric(gsub("condition_", "", unique(IncucyteKOTreatDatasetSubsetOriginal$condition_ID)))

      IncucyteDataAndMetaDataAndConditionsCompGroupDF <- IncucyteDataAndMetaDataAndConditionsCompGroupDF %>% mutate(NLMM_Analysis_ID = paste0("nlmm_analysis_", ConditionValue,sep=""))

      print(unique(IncucyteDataAndMetaDataAndConditionsCompGroupDF))

      ListOfDataFrames[[k]] <- data.frame(unique(IncucyteDataAndMetaDataAndConditionsCompGroupDF))

    }
  }

  IncucyteDataAndMetaDataAndConditionsLabeledDF <- do.call("rbind", ListOfDataFrames)

  return(IncucyteDataAndMetaDataAndConditionsLabeledDF)

}

#########################################################################################################################################################################################################################################################################

#' Assigning Identifiers for Analysis Groups

#' This function assigns an analysis ID (i.e. nlmm_analysis_1) to every analysis group: KO + treatment, KO + control, WT + treatment, and WT + control, where the goal is to compare compound treatment to control treatment and
#' knockout (KO) cell line to the wildtype (WT) cell line. This version is used before running the modeling-independent (phenotypes) analysis, as it accounts for the additional columns from the Split Point Table that previously weren't there when getSamePlatesPerComparisonGroupQC() was ran.
#' @param IncucyteDataAndMetaDataAndConditionsDF Dataframe that combines Incucyte Confluency Dataset, Meta Data and Condition ID table.
#' @param IncucyteDMSODilutionDF Dataframe that contains the contents of the Incucyte DMSO Dilution Table. This is used for finding the respective control treatment and concentration for a given compound.
#' @param ListOfSpecialCompounds List of compounds that are normalized to untreated cells rather than DMSO-treated cells.
#' @param ListofNewClones List of cell line modifications that use a different DMSO concentration.
#' @param wildtype Variable representing name of wildtype cell line.
#' @return
#'\item{IncucyteDataAndMetaDataAndConditionsLabeledDF}{  Dataframe that combines Incucyte Confluency Dataset, Meta Data, Condition IDs and analysis IDs.}
#' @author Caroline Barry
#' @import dplyr
#' @export

getSamePlatesPerComparisonGroup <- function(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDMSODilutionDF, ListOfSpecialCompounds=NULL,ListofNewClones=NULL, wildtype){

  IncucyteDataAndMetaDataAndConditionsLabeledDF <- data.frame()

  IncucyteDataAndMetaDataAndConditionsDFNoWTNoDMSONoUntreated <- subset(IncucyteDataAndMetaDataAndConditionsDF, !(cell_line_modifications == wildtype) & treatment_0!="DMSO" & treatment!="Untreated")

  ListOfConditionID = unique(IncucyteDataAndMetaDataAndConditionsDFNoWTNoDMSONoUntreated$condition_ID)

  ListOfDataFrames <- vector(mode = "list", length = length(ListOfConditionID))

  for(k in 1:length(ListOfConditionID)){

    #get that KO TreatDF
    IncucyteKOTreatDatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$condition_ID==ListOfConditionID[k])

    #get that WT TreatDF
    IncucyteWTTreatDatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, cell_line_modifications==wildtype& treatment_0!="DMSO" & treatment!="Untreated")

    IncucyteWTTreatDatasetSubsetOriginal <- subset(IncucyteWTTreatDatasetSubsetOriginal, treatment_0 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_0 & concentration_0 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_0 & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1 & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & plate_ID %in% IncucyteKOTreatDatasetSubsetOriginal$plate_ID)

    #get that KO DMSO DF
    CompoundName <-unique(IncucyteKOTreatDatasetSubsetOriginal$treatment_0)

    if(isTRUE(CompoundName %in% ListOfSpecialCompounds)==TRUE){

      #need to get untreated cell line and normalize treated to that
      #wanna now subset by all_plates column
      IncucyteKODMSODatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, cell_line_modifications %in% IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications & cell_number %in% IncucyteKOTreatDatasetSubsetOriginal$cell_number & treatment =="Untreated" &concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1)
      IncucyteWTDMSODatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, cell_line_modifications ==wildtype & cell_number %in% IncucyteKOTreatDatasetSubsetOriginal$cell_number & treatment =="Untreated" &concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1)

    }else{
      DMSOConcentration <- subset(IncucyteDMSODilutionDF, IncucyteDMSODilutionDF$Compound_name %in%IncucyteKOTreatDatasetSubsetOriginal$treatment_0)
      if(isTRUE(grepl("KO_[0-9]",unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications))|unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications) %in% ListofNewClones)==TRUE){

        DMSOTreatedDFSubset <- subset(IncucyteDataAndMetaDataAndConditionsDF, concentration_0%in%DMSOConcentration$DMSO_concentration_new_clones_0 | concentration_0%in%DMSOConcentration$DMSO_concentration_new_clones_1 | concentration_0%in%DMSOConcentration$DMSO_concentration_new_clones_2)

        DMSOTreatedDFSubset <- subset(DMSOTreatedDFSubset, treatment_0=="DMSO" & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1)

      }
      if(isTRUE(grepl("KO_[0-9]",unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications))|unique(IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications) %in% ListofNewClones)==FALSE){
        DMSOTreatedDFSubset <- subset(IncucyteDataAndMetaDataAndConditionsDF, concentration_0%in%DMSOConcentration$DMSO_concentration_0 & treatment_0=="DMSO" & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 & treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1)

      }

      IncucyteKODMSODatasetSubsetOriginal <- subset(DMSOTreatedDFSubset, cell_line_modifications %in% IncucyteKOTreatDatasetSubsetOriginal$cell_line_modifications & plate_ID %in% IncucyteKOTreatDatasetSubsetOriginal$plate_ID)
      IncucyteWTDMSODatasetSubsetOriginal <- subset(DMSOTreatedDFSubset, cell_line_modifications==wildtype & plate_ID %in% IncucyteKOTreatDatasetSubsetOriginal$plate_ID)}

    IncucyteDataAndMetaDataAndConditionsKOWTTreatDF <- merge(IncucyteKOTreatDatasetSubsetOriginal, IncucyteWTTreatDatasetSubsetOriginal, by=c("bio_rep_id", "condition_ID","plate_ID", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0","treatment_timepoint_0","treatment_duration_0",
                                                                                                                                              "treatment_1", "concentration_1","treatment_timepoint_1","treatment_duration_1","treatment_2", "concentration_2","treatment_timepoint_2","treatment_duration_2",
                                                                                                                                              "treatment_3", "concentration_3","treatment_timepoint_3","treatment_duration_3", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
                                                                                                                                              "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
                                                                                                                                              "split_first_conf_difference", "phenotype","treatment","CL_treat_conc_id","treat_conc_id"), all.x = TRUE, all.y = TRUE)

    IncucyteDataAndMetaDataAndConditionsKOWTTreatDMSODF <- merge(IncucyteDataAndMetaDataAndConditionsKOWTTreatDF, IncucyteKODMSODatasetSubsetOriginal, by=c("bio_rep_id", "condition_ID","plate_ID", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0","treatment_timepoint_0","treatment_duration_0",
                                                                                                                                                            "treatment_1", "concentration_1","treatment_timepoint_1","treatment_duration_1","treatment_2", "concentration_2","treatment_timepoint_2","treatment_duration_2",
                                                                                                                                                            "treatment_3", "concentration_3","treatment_timepoint_3","treatment_duration_3", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
                                                                                                                                                            "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
                                                                                                                                                            "split_first_conf_difference", "phenotype","treatment","CL_treat_conc_id","treat_conc_id"), all.x = TRUE, all.y = TRUE)


    IncucyteDataAndMetaDataAndConditionsCompGroupDF <- merge(IncucyteDataAndMetaDataAndConditionsKOWTTreatDMSODF, IncucyteWTDMSODatasetSubsetOriginal, by=c("bio_rep_id", "condition_ID","plate_ID", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0","treatment_timepoint_0","treatment_duration_0",
                                                                                                                                                            "treatment_1", "concentration_1","treatment_timepoint_1","treatment_duration_1","treatment_2", "concentration_2","treatment_timepoint_2","treatment_duration_2",
                                                                                                                                                            "treatment_3", "concentration_3","treatment_timepoint_3","treatment_duration_3", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
                                                                                                                                                            "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
                                                                                                                                                            "split_first_conf_difference", "phenotype","treatment","CL_treat_conc_id","treat_conc_id"), all.x = TRUE, all.y = TRUE)


    IncucyteDataAndMetaDataAndConditionsCompGroupDF <- IncucyteDataAndMetaDataAndConditionsCompGroupDF[,-c(8,9,12,13,14:17)]

    #make sure WT untreat, WT treat, KO untreat, KO treat are included for each plate
    ConditionsCheckDF <- IncucyteDataAndMetaDataAndConditionsCompGroupDF%>%
      group_by(plate_ID) %>%
      distinct(condition_ID) %>%
      dplyr::filter(n()>3)

    if(isTRUE(nrow(ConditionsCheckDF)>0)==TRUE){

      IncucyteDataAndMetaDataAndConditionsCompGroupDF <- subset(IncucyteDataAndMetaDataAndConditionsCompGroupDF, IncucyteDataAndMetaDataAndConditionsCompGroupDF$plate_ID %in% ConditionsCheckDF$plate_ID)

      #assign NLMM Analysis ID

      ConditionValue <- as.numeric(gsub("condition_", "", unique(IncucyteKOTreatDatasetSubsetOriginal$condition_ID)))

      IncucyteDataAndMetaDataAndConditionsCompGroupDF <- IncucyteDataAndMetaDataAndConditionsCompGroupDF %>% mutate(NLMM_Analysis_ID = paste0("nlmm_analysis_", ConditionValue,sep=""))


      print(unique(IncucyteDataAndMetaDataAndConditionsCompGroupDF))

      ListOfDataFrames[[k]] <- data.frame(unique(IncucyteDataAndMetaDataAndConditionsCompGroupDF))

    }
  }

  IncucyteDataAndMetaDataAndConditionsLabeledDF <- do.call("rbind", ListOfDataFrames)


  return(IncucyteDataAndMetaDataAndConditionsLabeledDF)

}

#########################################################################################################################################################################################################################################################################

#' Assigning Identifiers for Analysis Groups - NLMM Analysis on KO vs WT

#'
#' Ths function assigns an untreated analysis ID (i.e. untreat_nlmm_analysis_1) to every untreated analysis group: KO + untreated and WT + untreated, where the goal is to compare the knockout (KO) cell line to the wildtype (WT)
#'  cell line in the untreated state. This function is used before running NLMM on untreated data for the KO and WT cell lines.
#' @param IncucyteDataAndMetaDataAndConditionsDF Dataframe that combines Incucyte Confluency Dataset, Meta Data and Condition ID table.
#' @param wildtype A variable representing the wildtype cell line.
#' @param AnalysisName User-defined unique name of the analysis being conducted.
#' @return
#'\item{IncucyteDataAndMetaDataAndConditionsLabeledDF}{  Dataframe that combines Incucyte Confluency Dataset, Meta Data, Condition IDs and analysis IDs.}
#' @author Caroline Barry
#' @import dplyr
#' @export

getSamePlatesPerUntreatedComparisonGroup <- function(IncucyteDataAndMetaDataAndConditionsDF, wildtype, AnalysisName){

  IncucyteDataAndMetaDataAndConditionsLabeledDF <- data.frame()

  IncucyteDataAndMetaDataAndConditionsDFNoWTNoDMSONoUntreated <- subset(IncucyteDataAndMetaDataAndConditionsDF, !(cell_line_modifications == wildtype))

  ListOfConditionID = unique(IncucyteDataAndMetaDataAndConditionsDFNoWTNoDMSONoUntreated$condition_ID)

  for(k in 1:length(ListOfConditionID)){

    #get Untreated KO DF
    IncucyteKOTreatDatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, IncucyteDataAndMetaDataAndConditionsDF$condition_ID==ListOfConditionID[k])

    #get that Untreated WT DF

    IncucyteWTTreatDatasetSubsetOriginal <- subset(IncucyteDataAndMetaDataAndConditionsDF, cell_line_modifications==wildtype)

    IncucyteWTTreatDatasetSubsetOriginal <- subset(IncucyteWTTreatDatasetSubsetOriginal, treatment_0 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_0 & concentration_0 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_0 &
                                                     treatment_1 %in% IncucyteKOTreatDatasetSubsetOriginal$treatment_1 & concentration_1 %in% IncucyteKOTreatDatasetSubsetOriginal$concentration_1 &
                                                     plate_ID %in% IncucyteKOTreatDatasetSubsetOriginal$plate_ID & cell_number %in% IncucyteKOTreatDatasetSubsetOriginal$cell_number)

    IncucyteDataAndMetaDataAndConditionsKOWTTreatDF <- merge(IncucyteKOTreatDatasetSubsetOriginal, IncucyteWTTreatDatasetSubsetOriginal, by=c("bio_rep_id", "timepoint","condition_ID","cell_line", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0",
                                                                                                                                              "treatment_1", "concentration_1","mean_conf", "plate_ID","CL_treat_conc_id"), all.x = TRUE, all.y = TRUE)


    IncucyteDataAndMetaDataAndConditionsCompGroupDFBiorepList <- IncucyteDataAndMetaDataAndConditionsKOWTTreatDF[,-c(2,11)]


    IncucyteDataAndMetaDataAndConditionsCompGroupDFBiorepList <- IncucyteDataAndMetaDataAndConditionsCompGroupDFBiorepList%>%
      dplyr::distinct() %>%
      group_by(plate_ID) %>%
      dplyr::filter(n()>1)

    if(isTRUE(nrow(IncucyteDataAndMetaDataAndConditionsCompGroupDFBiorepList)>0)==TRUE){

      IncucyteDataAndMetaDataAndConditionsKOWTTreatDF <- subset(IncucyteDataAndMetaDataAndConditionsKOWTTreatDF, bio_rep_id %in% IncucyteDataAndMetaDataAndConditionsCompGroupDFBiorepList$bio_rep_id)

      #assign NLMM Analysis ID

      ConditionValue <- as.numeric(gsub("condition_", "", unique(IncucyteKOTreatDatasetSubsetOriginal$condition_ID)))

      IncucyteDataAndMetaDataAndConditionsKOWTTreatDF <- IncucyteDataAndMetaDataAndConditionsKOWTTreatDF %>% mutate(NLMM_Analysis_ID = paste0(AnalysisName,"_nlmm_analysis_", ConditionValue,sep=""))

      IncucyteDataAndMetaDataAndConditionsLabeledDF <- rbind(IncucyteDataAndMetaDataAndConditionsLabeledDF, IncucyteDataAndMetaDataAndConditionsKOWTTreatDF)}
  }

  return(IncucyteDataAndMetaDataAndConditionsLabeledDF)

}



