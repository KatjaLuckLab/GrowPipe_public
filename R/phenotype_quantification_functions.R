#' Normalize Phenotypes to Control Treatment

#'
#' This function normalizes all non-NLMM phenotypes (i.e. logistic phase duration) to the control treatment (i.e. DMSO).
#' @param IncucyteDataAndMetaDataAndConditionsLabeledDF Dataframe containing condition IDs, analysis group IDs and non-NLMM phenotypes.
#' @param controlcompound Variable representing the compound all treatments should be normalized to (i.e. DMSO).
#' @param ListOfSpecialCompounds List of compounds that are normalized to untreated cells rather than DMSO-treated cells.
#' @param normalizationtype Variable representing how conditions should be normalized: 'global' (default), all conditions normalized to same compound; 'firsttreatment', user-specified conditions containing more than 1 compound that should be normalized to the first compound.
#' @param respectconditions Variable detailing whether treated conditions are normalized to the control compound that also has the same second treatment at the same concentration, i.e. treatment_1 and concentration_1, ('TRUE', default)
#' or whether the second treatment is ignored when normalizing to the control compound ('FALSE').
#' @param ListOfNormalizationPairs List of compound pairs where the first is to be normalized to the second. Default NULL.

#' @return
#'\item{PhenotypesNormToDMSODF}{ Dataframe containing condition IDs, analysis group IDs and phenotypes normalized to the control treatment.}
#' @author Caroline Barry
#' @export

normalizeToDMSO <- function(IncucyteDataAndMetaDataAndConditionsLabeledDF,controlcompound, ListOfSpecialCompounds = NULL, normalizationtype = "global",respectconditions =TRUE, ListOfNormalizationPairs = NULL){

  PhenotypesNormToDMSODF <- data.frame(stringsAsFactors = FALSE)


  # if all conditions get normalized to the same compound (i.e. DMSO):
  if(isTRUE(normalizationtype=="global")==TRUE){

    # IncucyteSplitPointPlateIDsDFDMSO <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment_0 ==controlcompound)

    # IncucyteSplitPointPlateIDsDFUntreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment =="Untreated")

    ListOfSpecialCompounds2 <- NULL


    #if second treatment is taken into account when normalizing:
    if(isTRUE(respectconditions == TRUE)==TRUE){

      IncucyteSplitPointPlateIDsDFTreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment_0 !=controlcompound & treatment !="Untreated")

      IncucyteSplitPointPlateIDsDFDMSO <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment_0 ==controlcompound)

      # if there are user-specified compounds that get normalized to the untreated condition:
      if(isTRUE(!(is.null(ListOfSpecialCompounds)))==TRUE){

        IncucyteSplitPointPlateIDsDFUntreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment =="Untreated")}

    }else {
      #if second treatment is ignored when normalizing:
      IncucyteSplitPointPlateIDsDFTreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment !=controlcompound & treatment !="Untreated")

      IncucyteSplitPointPlateIDsDFDMSO <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment ==controlcompound)

      # if there are user-specified compounds that get normalized to the untreated condition:
      if(isTRUE(!(is.null(ListOfSpecialCompounds)))==TRUE){

        IncucyteSplitPointPlateIDsDFUntreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment =="Untreated")
      }

    }
    # else{
    #   if(isTRUE(respectconditions == TRUE)==TRUE){IncucyteSplitPointPlateIDsDFTreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment_0 !=controlcompound)
    #   IncucyteSplitPointPlateIDsDFDMSO <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment_0 ==controlcompound)
    #
    #   }else{
    #     IncucyteSplitPointPlateIDsDFTreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment !=controlcompound)
    #
    #     IncucyteSplitPointPlateIDsDFDMSO <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment ==controlcompound)
  }

  # }
  # }

  # conditions where more than one compound was used (treatment_0 & treatment_1), get normalized to the first compound (treatment_0) as user-specified in ListOfNormalizationPairs
  if(isTRUE(normalizationtype=="firsttreatment")==TRUE){

    #find conditions that were only treated with treatment_0 that is same as controlcompound
    IncucyteSplitPointPlateIDsDFFirstCmpnd <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment %in% controlcompound)

    ListOfSpecialCompounds2 <- sapply(ListOfNormalizationPairs, `[`, 1)

    NormalizationPairsDF <- as.data.frame(ListOfNormalizationPairs)

    NormalizationPairsDF <- t(NormalizationPairsDF)


    # IncucyteSplitPointPlateIDsDFUntreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment =="Untreated")

    #if there are certain compounds that need to be normalized to the untreated condition
    if(isTRUE(!(is.null(ListOfSpecialCompounds)))==TRUE){

      IncucyteSplitPointPlateIDsDFUntreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, treatment =="Untreated")

      IncucyteSplitPointPlateIDsDFTreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, !(treatment %in%controlcompound) & treatment !="Untreated")
    } else{
      #if there are no compounds that need to be normalized to the untreated condition
      IncucyteSplitPointPlateIDsDFTreat <- subset(IncucyteDataAndMetaDataAndConditionsLabeledDF, !(treatment %in%controlcompound))}

  } # & is.na(treatment_1)


  # loop row-wise through conditions that need to be normalized
  for(i in 1:nrow(IncucyteSplitPointPlateIDsDFTreat)){

    #get plate_ID
    TreatPartner <- as.character(IncucyteSplitPointPlateIDsDFTreat[i,3])#9

    #if compound matches list for those that need to be normalized to untreated condition
    if(isTRUE(as.character(IncucyteSplitPointPlateIDsDFTreat[i,26]) %in% ListOfSpecialCompounds)==TRUE){

      IncucyteSplitPointPlateIDsDFDMSOSubset <- subset(IncucyteSplitPointPlateIDsDFUntreat, cell_line_modifications==as.character(IncucyteSplitPointPlateIDsDFTreat[i,4]) & cell_number==as.character(IncucyteSplitPointPlateIDsDFTreat[i,5])&
                                                         plate_ID==TreatPartner & NLMM_Analysis_ID==as.character(IncucyteSplitPointPlateIDsDFTreat[i,29]))
    }

    #if compound matches list for those that need to be normalized to first condition
    if(isTRUE(as.character(IncucyteSplitPointPlateIDsDFTreat[i,26]) %in% ListOfSpecialCompounds2)==TRUE){

      cntrlcmpnd <- subset(NormalizationPairsDF, NormalizationPairsDF[,1]==IncucyteSplitPointPlateIDsDFTreat[i,26])
      # if(isTRUE(normalizationtype=="firsttreatment")==TRUE){
      IncucyteSplitPointPlateIDsDFDMSOSubset <- subset(IncucyteSplitPointPlateIDsDFFirstCmpnd, cell_line_modifications==as.character(IncucyteSplitPointPlateIDsDFTreat[i,4]) & cell_number==as.character(IncucyteSplitPointPlateIDsDFTreat[i,5])&
                                                         plate_ID==TreatPartner& treatment ==cntrlcmpnd[,2]  & concentration_0 == as.character(IncucyteSplitPointPlateIDsDFTreat[i,7]))
    }#treatment_0 ==as.character(IncucyteSplitPointPlateIDsDFTreat[i,6])  & concentration_0 == as.character(IncucyteSplitPointPlateIDsDFTreat[i,7]) & NLMM_Analysis_ID==as.character(IncucyteSplitPointPlateIDsDFTreat[i,29])

    #if compound doesn't list for those that need to be normalized to untreated OR first condition AND all conditions get normalized to same control compound:
    if(isTRUE(!(as.character(IncucyteSplitPointPlateIDsDFTreat[i,26]) %in% ListOfSpecialCompounds) & normalizationtype=="global")==TRUE){
      # DMSOConcentration <- subset(IncucyteDMSODilutionDF, IncucyteDMSODilutionDF$Compound_name==as.character(IncucyteSplitPointPlateIDsDFTreat[i,6]))

      # if(isTRUE(grepl("KO_[0-9]",as.character(IncucyteSplitPointPlateIDsDFTreat[i,4])))==TRUE){
      #   IncucyteSplitPointPlateIDsDFDMSOSubset <- subset(IncucyteSplitPointPlateIDsDFDMSO, concentration_0==DMSOConcentration$DMSO_concentration_1)
      #   # DMSOConcentration <- DMSOConcentration$DMSO_concentration_1
      #   # subset(IncucyteDMSODilutionDF, IncucyteDMSODilutionDF$DMSO_concentration_1%in% IncucyteKOTreatDatasetSubsetOriginal$treatment_0 )
      # }
      # if(isTRUE(grepl("KO_[0-9]",as.character(IncucyteSplitPointPlateIDsDFTreat[i,4])))==FALSE){
      #   IncucyteSplitPointPlateIDsDFDMSOSubset <- subset(IncucyteSplitPointPlateIDsDFDMSO, concentration_0==DMSOConcentration$DMSO_concentration_0)
      #   # DMSOConcentration <- DMSOConcentration$DMSO_concentration_0
      #   # subset(IncucyteDMSODilutionDF, IncucyteDMSODilutionDF$Compound_name%in% IncucyteKOTreatDatasetSubsetOriginal$treatment_0 )
      # }
      IncucyteSplitPointPlateIDsDFDMSOSubset <- subset(IncucyteSplitPointPlateIDsDFDMSO,cell_line_modifications ==as.character(IncucyteSplitPointPlateIDsDFTreat[i,4])& plate_ID==TreatPartner & NLMM_Analysis_ID==as.character(IncucyteSplitPointPlateIDsDFTreat[i,29]))

      # # IncucyteSplitPointPlateIDsDFDMSOSubset <- subset(IncucyteSplitPointPlateIDsDFDMSO, concentration_0%in%DMSOConcentration$DMSO_concentration)
      # IncucyteSplitPointPlateIDsDFDMSOSubset <- subset(IncucyteSplitPointPlateIDsDFDMSOSubset, cell_line_modifications ==as.character(IncucyteSplitPointPlateIDsDFTreat[i,4]))
      # IncucyteSplitPointPlateIDsDFDMSOSubset <-subset(IncucyteSplitPointPlateIDsDFDMSOSubset,plate_ID==TreatPartner)
    }

    if(isTRUE(nrow(IncucyteSplitPointPlateIDsDFDMSOSubset)>1)==TRUE){

      IncucyteSplitPointPlateIDsDFDMSOSubset <- IncucyteSplitPointPlateIDsDFDMSOSubset[1,]#subset(IncucyteSplitPointPlateIDsDFDMSOSubset, NLMM_Analysis_ID ==as.character(IncucyteSplitPointPlateIDsDFTreat[i,29]))

    }

    #populate empty dataframe
    if(nrow(IncucyteSplitPointPlateIDsDFDMSOSubset)>0) {

      PhenotypesNormToDMSODF[i,1]<- IncucyteSplitPointPlateIDsDFTreat[i,3]#9
      PhenotypesNormToDMSODF[i,2]<- IncucyteSplitPointPlateIDsDFTreat[i,2]#2
      PhenotypesNormToDMSODF[i,3]<- IncucyteSplitPointPlateIDsDFTreat[i,4]#4
      PhenotypesNormToDMSODF[i,4]<- IncucyteSplitPointPlateIDsDFTreat[i,6]#6
      PhenotypesNormToDMSODF[i,5]<- IncucyteSplitPointPlateIDsDFTreat[i,7]
      PhenotypesNormToDMSODF[i,6]<- IncucyteSplitPointPlateIDsDFTreat[i,8]#6
      PhenotypesNormToDMSODF[i,7]<- IncucyteSplitPointPlateIDsDFTreat[i,9]

      PhenotypesNormToDMSODF[i,8]<- IncucyteSplitPointPlateIDsDFTreat[i,16]
      PhenotypesNormToDMSODF[i,9]<- IncucyteSplitPointPlateIDsDFTreat[i,16]/IncucyteSplitPointPlateIDsDFDMSOSubset$logistic_growth_duration
      PhenotypesNormToDMSODF[i,10]<- IncucyteSplitPointPlateIDsDFTreat[i,24]
      PhenotypesNormToDMSODF[i,11]<- IncucyteSplitPointPlateIDsDFTreat[i,24]/IncucyteSplitPointPlateIDsDFDMSOSubset$split_first_conf_difference
      PhenotypesNormToDMSODF[i,12]<- IncucyteSplitPointPlateIDsDFTreat[i,19]
      PhenotypesNormToDMSODF[i,13]<- IncucyteSplitPointPlateIDsDFTreat[i,19]/IncucyteSplitPointPlateIDsDFDMSOSubset$after_24h_confluency
      PhenotypesNormToDMSODF[i,14]<- IncucyteSplitPointPlateIDsDFTreat[i,20]
      PhenotypesNormToDMSODF[i,15]<- IncucyteSplitPointPlateIDsDFTreat[i,20]/IncucyteSplitPointPlateIDsDFDMSOSubset$after_48h_confluency
      PhenotypesNormToDMSODF[i,16]<- IncucyteSplitPointPlateIDsDFTreat[i,21]
      PhenotypesNormToDMSODF[i,17]<- IncucyteSplitPointPlateIDsDFTreat[i,21]/IncucyteSplitPointPlateIDsDFDMSOSubset$after_72h_confluency
      PhenotypesNormToDMSODF[i,18]<- IncucyteSplitPointPlateIDsDFTreat[i,22]
      PhenotypesNormToDMSODF[i,19]<- IncucyteSplitPointPlateIDsDFTreat[i,22]/IncucyteSplitPointPlateIDsDFDMSOSubset$last_confluency
      PhenotypesNormToDMSODF[i,20]<- IncucyteSplitPointPlateIDsDFTreat[i,23]
      PhenotypesNormToDMSODF[i,21]<- IncucyteSplitPointPlateIDsDFTreat[i,23]/IncucyteSplitPointPlateIDsDFDMSOSubset$conf_96hr_48hr_ratio
      PhenotypesNormToDMSODF[i,22]<- IncucyteSplitPointPlateIDsDFTreat[i,26]
      PhenotypesNormToDMSODF[i,23]<- IncucyteSplitPointPlateIDsDFTreat[i,27]
      PhenotypesNormToDMSODF[i,24]<- IncucyteSplitPointPlateIDsDFTreat[i,28]
      PhenotypesNormToDMSODF[i,25]<- IncucyteSplitPointPlateIDsDFTreat[i,29]
    }
  }


  colnames(PhenotypesNormToDMSODF) <- c("plate_ID", "condition_ID", "cell_line_modification", "treatment_0", "concentration_0","treatment_1", "concentration_1",
                                        "First_Phase_Duration",
                                        "First_Phase_Duration_ratio", "First_Phase_Change_in_Confluency", "First_Phase_Change_in_Confluency_ratio","Confluency_after_24h","Confluency_after_24h_ratio",
                                        "Confluency_after_48h","Confluency_after_48h_ratio","Confluency_after_72h","Confluency_after_72h_ratio",
                                        "Confluency_after_96h","Confluency_after_96h_ratio","Second_Phase_Relative_Change_in_Confluency","Second_Phase_Relative_Change_in_Confluency_ratio",
                                        "treatment","CL_treat_conc_id","treat_conc_id","NLMM_Analysis_ID")

  return(PhenotypesNormToDMSODF)
}


# "condition_ID","bio_rep_id","plate_ID", "cell_line_modifications", "cell_number", "treatment_0", "concentration_0",
# "treatment_1", "concentration_1", "split_timepoint", "logistic_growth_phase", "logistic_growth_duration", "second_phase", "second_phase_duration",
# "after_24h_confluency",  "after_48h_confluency","after_72h_confluency","last_confluency","conf_96hr_48hr_ratio",
# "split_first_conf_difference", "phenotype","treatment","CL_treat_conc_id","treat_conc_id","NLMM_Analysis_ID"

#########################################################################################################################################################################################################################################################################
#' Calculating AUC for each Phenotype Normalized to Control Treatment - Quality Control

#'
#' This function calculates the AUC for all modeling-independent phenotypes (i.e. confluency after 24hr of treatment) for KO and its respective WT separately.
#' @param PhenotypesNormToDMSODF A dataframe containing phenotypes normalized to DMSO, condition IDs and analysis group IDs.
#' @param wildtype Variable representing name of wildtype cell line.
#' @return
#'\item {AUCsPerConditionDF} {A dataframe containing the AUC for WT and KO, and the AUC normalized to WT for each biological replicate for each given phenotype.}
#' @author Caroline Barry
#' @import dplyr
#' @import stats
#' @export


getAUCsPerConditionNormalizedToDMSO <- function(PhenotypesNormToDMSODF, wildtype){

  #############  if cell lines are to be normalized to user-sepcified wildtype (i.e. KO/WT)   #############
  if(isTRUE(wildtype!="self")==TRUE){
    treat_conc_combos <- unique(PhenotypesNormToDMSODF$treatment)

    AUCsPerConditionDF <- data.frame(stringsAsFactors = FALSE)

    ListOfPhenotypeRatios <- c("Confluency_after_24h_ratio","Confluency_after_48h_ratio","Confluency_after_72h_ratio","Confluency_after_96h_ratio","Second_Phase_Relative_Change_in_Confluency_ratio",
                               "First_Phase_Duration_ratio","First_Phase_Change_in_Confluency_ratio")

    i = 0

    for(k in 1:length(treat_conc_combos)){

      PhenotypesNormToDMSODFSubset = subset(PhenotypesNormToDMSODF, PhenotypesNormToDMSODF$treatment==treat_conc_combos[k])#k

      # get list of concentrations that wildtype has been treated with
      PhenotypesNormToDMSODFSubsetWT <- PhenotypesNormToDMSODFSubset %>%
        group_by(concentration_0) %>%
        dplyr::filter(cell_line_modification ==wildtype)

      # get list of concentrations that non-wildtype cell lines have been treated with
      PhenotypesNormToDMSODFSubset <- PhenotypesNormToDMSODFSubset %>%
        group_by(concentration_0) %>%
        dplyr::filter(cell_line_modification !=wildtype)

      #filter for concentrations only tested on wildtype
      PhenotypesNormToDMSODFSubset = subset(PhenotypesNormToDMSODFSubset, treatment_0 %in% PhenotypesNormToDMSODFSubsetWT$treatment_0 & concentration_0 %in% PhenotypesNormToDMSODFSubsetWT$concentration_0 & treatment_1 %in% PhenotypesNormToDMSODFSubsetWT$treatment_1 & concentration_1 %in% PhenotypesNormToDMSODFSubsetWT$concentration_1)

      CL_list <- unique(PhenotypesNormToDMSODFSubset$cell_line_modification)

      #loop through each KO cell line
      for(j in 1:length(CL_list)){

        i= i + 1

        PhenotypesNormToDMSODFSubsetCL = subset(PhenotypesNormToDMSODFSubset, PhenotypesNormToDMSODFSubset$cell_line_modification==CL_list[j])#j

        PhenotypesNormToDMSODFSubsetWT2 = subset(PhenotypesNormToDMSODFSubsetWT,
                                                 NLMM_Analysis_ID %in% PhenotypesNormToDMSODFSubsetCL$NLMM_Analysis_ID)

        PlateList <- unique(PhenotypesNormToDMSODFSubsetCL$plate_ID)

        m=0

        #loop through each plate ID
        if(length(PlateList) > 1){

          TestDF <- data.frame()

          for(l in 1:length(PlateList)){

            PhenotypesNormToDMSODFSubsetCL2 = subset(PhenotypesNormToDMSODFSubsetCL, PhenotypesNormToDMSODFSubsetCL$plate_ID==PlateList[l])#l

            PhenotypesNormToDMSODFSubsetWT3 = subset(PhenotypesNormToDMSODFSubsetWT2,
                                                     plate_ID %in% PhenotypesNormToDMSODFSubsetCL2$plate_ID)

            # get concentrations in numerical order for KO df
            KOCompoundConcentrationString <- PhenotypesNormToDMSODFSubsetCL2$concentration_0

            KOCompoundConcentrationNumber <- lapply(KOCompoundConcentrationString, function(x) gsub("[^0-9.-]", "", x))

            PhenotypesNormToDMSODFSubsetCL2 <- PhenotypesNormToDMSODFSubsetCL2%>%
              ungroup()%>%
              mutate(conc_values = unlist(KOCompoundConcentrationNumber))

            PhenotypesNormToDMSODFSubsetCL2$conc_values <- log2(as.numeric(unlist(PhenotypesNormToDMSODFSubsetCL2$conc_values)))

            #loop through each phenotype
            for(n in 1:length(ListOfPhenotypeRatios)){

              PhenotypeRatio <- ListOfPhenotypeRatios[n]

              PhenotypeName1 <- gsub("_ratio", "", PhenotypeRatio)

              PhenotypeName2 <- gsub("_", " ", PhenotypeName1)

              #in the case where there's only one biorep...
              if(isTRUE(length(unique(PhenotypesNormToDMSODFSubsetCL2$NLMM_Analysis_ID))==1)==TRUE){

                ko_auc_val$value <- NaN
                ko_auc_val$abs.error <- "Not enough data"
                ko_min_conc <- NaN
                ko_max_conc <- NaN
              } else{

                #in the case where there's more than one biorep, compute AUC
                ko_min_conc <-  min(as.numeric(unlist(PhenotypesNormToDMSODFSubsetCL2$conc_values)))

                ko_max_conc <- max(as.numeric(unlist(PhenotypesNormToDMSODFSubsetCL2$conc_values)))


                ko_auc_val <- integrate(approxfun(PhenotypesNormToDMSODFSubsetCL2$conc_values,unlist(PhenotypesNormToDMSODFSubsetCL2[,PhenotypeRatio]), ties = mean), lower = ko_min_conc, upper = ko_max_conc)
              }

              # get concentrations in numerical order for WT df
              WTCompoundConcentrationString <- PhenotypesNormToDMSODFSubsetWT3$concentration_0

              WTCompoundConcentrationNumber <- lapply(WTCompoundConcentrationString, function(x) gsub("[^0-9.-]", "", x))

              PhenotypesNormToDMSODFSubsetWT3 <- PhenotypesNormToDMSODFSubsetWT3%>%
                ungroup()%>%
                mutate(conc_values = unlist(WTCompoundConcentrationNumber))

              PhenotypesNormToDMSODFSubsetWT3$conc_values <- log2(as.numeric(unlist(PhenotypesNormToDMSODFSubsetWT3$conc_values)))

              #in the case where there's only one biorep...
              if(isTRUE(length(unique(PhenotypesNormToDMSODFSubsetWT3$NLMM_Analysis_ID))==1)==TRUE){
                wt_auc_val$value <- NaN
                wt_auc_val$abs.error <- "Not enough data"
                wt_min_conc <- NaN
                wt_max_conc <- NaN
              } else{

                #in the case where there's more than one biorep, compute AUC
                wt_min_conc <-  min(as.numeric(unlist(PhenotypesNormToDMSODFSubsetWT3$conc_values)))

                wt_max_conc <- max(as.numeric(unlist(PhenotypesNormToDMSODFSubsetWT3$conc_values)))

                wt_auc_val <- integrate(approxfun(PhenotypesNormToDMSODFSubsetWT3$conc_values,unlist(PhenotypesNormToDMSODFSubsetWT3[,PhenotypeRatio]), ties = mean), lower = wt_min_conc, upper = wt_max_conc)
              }
              #populate TestDF with AUC values

              m =m + 1

              TestDF[m,1] <- unique(PhenotypesNormToDMSODFSubsetCL2$treatment)
              TestDF[m,2] <- unique(PhenotypesNormToDMSODFSubsetCL2$cell_line_modification)
              TestDF[m,3] <- unique(PhenotypesNormToDMSODFSubsetCL2$plate_ID)
              TestDF[m,4] <-log2(ko_auc_val$value)
              TestDF[m,5] <-log2(wt_auc_val$value)
              TestDF[m,6] <-log2(ko_auc_val$value/wt_auc_val$value)
              TestDF[m,7] <- PhenotypeName1
            }
          }

          #add TestDF to empty dataframe
          AUCsPerConditionDF <- rbind(AUCsPerConditionDF, TestDF)
        }

      }
    }

    colnames(AUCsPerConditionDF) <- c("treatment","cell_line_modification","plate_ID", "KO_AUCs", "WT_AUCs", "log_KO/WT_AUC_Ratio","Phenotype")

    #create CL_group
    AUCsPerConditionDF <- AUCsPerConditionDF %>% mutate(CL_group= paste(cell_line_modification,wildtype,sep="|"))

    #split dataframe into WT dataframe and KO dataframe
    AUCsPerConditionWTDF <- AUCsPerConditionDF[,-4]

    AUCsPerConditionWTDF <- AUCsPerConditionWTDF%>% dplyr::rename(log_AUC_values = WT_AUCs) %>%
      mutate(cell_line_modification = wildtype)

    AUCsPerConditionDF <- AUCsPerConditionDF[,-5]

    AUCsPerConditionDF <- AUCsPerConditionDF%>% dplyr::rename(log_AUC_values = KO_AUCs)

    #merge the WT and KO dataframes
    AUCsPerConditionDF <- merge(AUCsPerConditionDF,AUCsPerConditionWTDF, by = c("treatment", "cell_line_modification", "plate_ID","log_AUC_values","log_KO/WT_AUC_Ratio","Phenotype","CL_group"), all=TRUE)
  }



  ############# in case where user desires cell lines are normalized to themselves (i.e. KO/KO)   #############
  else if(isTRUE(wildtype=="self")==TRUE){
    treat_conc_combos <- unique(PhenotypesNormToDMSODF$treatment_0)

    AUCsPerConditionDF <- data.frame(stringsAsFactors = FALSE)

    ListOfPhenotypeRatios <- c("Confluency_after_24h_ratio","Confluency_after_48h_ratio","Confluency_after_72h_ratio","Confluency_after_96h_ratio","Second_Phase_Relative_Change_in_Confluency_ratio",
                               "First_Phase_Duration_ratio","First_Phase_Change_in_Confluency_ratio")

    i = 0

    for(k in 1:length(treat_conc_combos)){

      #subset data for one treatment_0 (allows for treatment_0 AND treatment_0 + treatment_1 to be left in dataframe!)
      PhenotypesNormToDMSODFSubset = subset(PhenotypesNormToDMSODF, PhenotypesNormToDMSODF$treatment_0==treat_conc_combos[k])#k

      CL_list <- unique(PhenotypesNormToDMSODFSubset$cell_line_modification)

      #loop through cell lines
      for(j in 1:length(CL_list)){

        #subset data for one cell line
        PhenotypesNormToDMSODFSubsetCL = subset(PhenotypesNormToDMSODFSubset, PhenotypesNormToDMSODFSubset$cell_line_modification==CL_list[j])#j

        # get list of concentrations that condition with only treatment_0 have been treated with
        PhenotypesNormToDMSODFSubsetWT <- PhenotypesNormToDMSODFSubsetCL %>%
          group_by(concentration_0) %>%
          dplyr::filter(is.na(treatment_1))

        # get list of concentrations that condition with treatment_0 & treatment_1 have been treated with
        PhenotypesNormToDMSODFSubsetCL <- PhenotypesNormToDMSODFSubsetCL %>%
          group_by(concentration_0) %>%
          dplyr::filter(!(is.na(treatment_1)))

        #filter to make sure the same concentrations between the two conditions is kept
        PhenotypesNormToDMSODFSubsetCL = subset(PhenotypesNormToDMSODFSubsetCL, treatment_0 %in% PhenotypesNormToDMSODFSubsetWT$treatment_0 & concentration_0 %in% PhenotypesNormToDMSODFSubsetWT$concentration_0)

        PhenotypesNormToDMSODFSubsetWT = subset(PhenotypesNormToDMSODFSubsetWT, treatment_0 %in% PhenotypesNormToDMSODFSubsetCL$treatment_0 & concentration_0 %in% PhenotypesNormToDMSODFSubsetCL$concentration_0)

        i= i + 1

        PhenotypesNormToDMSODFSubsetWT2 = subset(PhenotypesNormToDMSODFSubsetWT,plate_ID %in% PhenotypesNormToDMSODFSubsetCL$plate_ID)

        PlateList <- unique(PhenotypesNormToDMSODFSubsetCL$plate_ID)

        m=0

        if(length(PlateList) > 1){

          TestDF <- data.frame()

          #loop through plate IDs
          for(l in 1:length(PlateList)){

            PhenotypesNormToDMSODFSubsetCL2 = subset(PhenotypesNormToDMSODFSubsetCL, PhenotypesNormToDMSODFSubsetCL$plate_ID==PlateList[l])#l

            PhenotypesNormToDMSODFSubsetWT3 = subset(PhenotypesNormToDMSODFSubsetWT2,plate_ID %in% PhenotypesNormToDMSODFSubsetCL2$plate_ID)

            # get concentrations in numerical order for condition with treatment_0 & treatment_1 df
            KOCompoundConcentrationString <- PhenotypesNormToDMSODFSubsetCL2$concentration_0

            KOCompoundConcentrationNumber <- lapply(KOCompoundConcentrationString, function(x) gsub("[^0-9.-]", "", x))

            PhenotypesNormToDMSODFSubsetCL2 <- PhenotypesNormToDMSODFSubsetCL2%>%
              ungroup()%>%
              mutate(conc_values = unlist(KOCompoundConcentrationNumber))

            PhenotypesNormToDMSODFSubsetCL2$conc_values <- log2(as.numeric(unlist(PhenotypesNormToDMSODFSubsetCL2$conc_values)))

            #loop through phenotypes
            for(n in 1:length(ListOfPhenotypeRatios)){

              PhenotypeRatio <- ListOfPhenotypeRatios[n]

              PhenotypeName1 <- gsub("_ratio", "", PhenotypeRatio)

              PhenotypeName2 <- gsub("_", " ", PhenotypeName1)

              #in the case there's only one biorep...
              if(isTRUE(length(unique(PhenotypesNormToDMSODFSubsetCL2$NLMM_Analysis_ID))==1)==TRUE){

                ko_auc_val$value <- NaN
                ko_auc_val$abs.error <- "Not enough data"
                ko_min_conc <- NaN
                ko_max_conc <- NaN
              } else{

                #in the case there's more than one biorep, compute AUC for condition with treatment_0 & treatment_1
                ko_min_conc <-  min(as.numeric(unlist(PhenotypesNormToDMSODFSubsetCL2$conc_values)))

                ko_max_conc <- max(as.numeric(unlist(PhenotypesNormToDMSODFSubsetCL2$conc_values)))


                ko_auc_val <- integrate(approxfun(PhenotypesNormToDMSODFSubsetCL2$conc_values,unlist(PhenotypesNormToDMSODFSubsetCL2[,PhenotypeRatio]), ties = mean), lower = ko_min_conc, upper = ko_max_conc)
              }

              # get concentrations in numerical order for ondition with only treatment_0 df
              WTCompoundConcentrationString <- PhenotypesNormToDMSODFSubsetWT3$concentration_0

              WTCompoundConcentrationNumber <- lapply(WTCompoundConcentrationString, function(x) gsub("[^0-9.-]", "", x))

              PhenotypesNormToDMSODFSubsetWT3 <- PhenotypesNormToDMSODFSubsetWT3%>%
                ungroup()%>%
                mutate(conc_values = unlist(WTCompoundConcentrationNumber))

              PhenotypesNormToDMSODFSubsetWT3$conc_values <- log2(as.numeric(unlist(PhenotypesNormToDMSODFSubsetWT3$conc_values)))

              #in the case there's only one biorep...
              if(isTRUE(length(unique(PhenotypesNormToDMSODFSubsetWT3$NLMM_Analysis_ID))==1)==TRUE){
                wt_auc_val$value <- NaN
                wt_auc_val$abs.error <- "Not enough data"
                wt_min_conc <- NaN
                wt_max_conc <- NaN
              } else{

                #in the case there's more than one biorep, compute AUC for condition with treatment_0 & treatment_1
                wt_min_conc <-  min(as.numeric(unlist(PhenotypesNormToDMSODFSubsetWT3$conc_values)))

                wt_max_conc <- max(as.numeric(unlist(PhenotypesNormToDMSODFSubsetWT3$conc_values)))

                wt_auc_val <- integrate(approxfun(PhenotypesNormToDMSODFSubsetWT3$conc_values,unlist(PhenotypesNormToDMSODFSubsetWT3[,PhenotypeRatio]), ties = mean), lower = wt_min_conc, upper = wt_max_conc)
              }

              #populate TestDF with AUC values
              m =m + 1

              TestDF[m,1] <- unique(PhenotypesNormToDMSODFSubsetCL2$treatment)
              TestDF[m,2] <- unique(PhenotypesNormToDMSODFSubsetCL2$cell_line_modification)
              TestDF[m,3] <- unique(PhenotypesNormToDMSODFSubsetCL2$plate_ID)
              TestDF[m,4] <-log2(ko_auc_val$value)
              TestDF[m,5] <-log2(wt_auc_val$value)
              TestDF[m,6] <-log2(ko_auc_val$value/wt_auc_val$value)
              TestDF[m,7] <- PhenotypeName1
            }
          }
          #add TestDF to empty dataframe

          AUCsPerConditionDF <- rbind(AUCsPerConditionDF, TestDF)
        }

      }
    }

    colnames(AUCsPerConditionDF) <- c("treatment","cell_line_modification","plate_ID", "KO_AUCs", "WT_AUCs", "log_KO/WT_AUC_Ratio","Phenotype")

    #create CL_group
    AUCsPerConditionDF <- AUCsPerConditionDF %>% mutate(CL_group=paste(cell_line_modification,wildtype,sep="|"))

    #split dataframe into condition with only treatment_0 and both treatment_0 & treatment_1
    AUCsPerConditionWTDF <- AUCsPerConditionDF[,-4]

    AUCsPerConditionWTDF <- AUCsPerConditionWTDF%>% dplyr::rename(log_AUC_values = WT_AUCs) %>%
      mutate(cell_line_modification = "self")

    AUCsPerConditionDF <- AUCsPerConditionDF[,-5]

    AUCsPerConditionDF <- AUCsPerConditionDF%>% dplyr::rename(log_AUC_values = KO_AUCs)

    #merge two dataframes back together

    AUCsPerConditionDF <- merge(AUCsPerConditionDF,AUCsPerConditionWTDF, by = c("treatment", "cell_line_modification", "plate_ID","log_AUC_values","log_KO/WT_AUC_Ratio","Phenotype","CL_group"), all=TRUE)
  }

  return(AUCsPerConditionDF)

}

#########################################################################################################################################################################################################################################################################

#' Calculate Significance of AUC values for Modeling-independent Phenotype

#'
#' This function calculates the significance of the AUC values for each modeling-independent phenotype of KO compared to WT treated with the same compound for each biological replicate by using the Bayes-moderated paired t-test and applying FDR for multiple testing correction.
#' @param AUCsPerConditionDF A dataframe containing area under the curve (AUC) of the dose-response curves per treatment per phenotype normalized to control treatment, condition IDs and analysis group IDs.
#' @param wildtype Variable representing name of the wildtype cell line.

#' @return
#'\item{AUCSignificanceDF}{ A dataframe containing treatment, cell line modification, mean log2 fold-change, t-statistic and p-value for each phenotype and condition.}
#' @author Caroline Barry
#' @import dplyr
#' @import rstatix
#' @import stats
#' @export

getAUCSignificance <- function(AUCsPerConditionDF, wildtype){


  AUCsPerConditionDF <- subset(AUCsPerConditionDF, cell_line_modification!=wildtype)

  # run t-test and p-value correction per phenotype

  ### Confluency_after_24h #######################################################

  # subset dataframe for desired phenotype
  df_24h <- AUCsPerConditionDF %>%
    filter(Phenotype == "Confluency_after_24h")

  # determine maximum number of biological replicates used
  max_n_col <- df_24h %>%
    filter(cell_line_modification != wildtype) %>%
    group_by(treatment, cell_line_modification) %>%
    dplyr::summarise(n = n()) %>%
    pull(n) %>%
    max()

  dm_24h <- df_24h %>%
    dplyr::filter(cell_line_modification != wildtype) %>%  # remove WT, because we won't need them for the analysis
    dplyr::group_by(treatment, cell_line_modification) %>%  # cast data to a wide format
    reframe(data.frame(t(c(`log_KO/WT_AUC_Ratio`, rep(NA, max_n_col - length(`log_KO/WT_AUC_Ratio`)))))) %>%  # Since max replicate number is 6, we cast to a matrix with 6 columns
    unite(col = "condition", treatment, cell_line_modification, sep = "_") %>%
    column_to_rownames("condition") %>%
    as.matrix()

  # fit linear model to log2 fold change values for each KO vs WT comparison for each compound
  fm_24h <- lmFit(object = dm_24h)

  # apply empirical Bayes moderation to shrink the standard errors and get moderated statistics from the linear model fits
  fm_24h <- eBayes(fm_24h)

  # retrieve Benjamini-Hochberg corrected p-values
  rslt_24h <- topTable(fm_24h, number = "all", sort.by = "none")


  ### Confluency_after_48h #######################################################

  # subset dataframe for desired phenotype
  df_48h <- AUCsPerConditionDF %>%
    filter(Phenotype == "Confluency_after_48h")

  # determine maximum number of biological replicates used
  max_n_col <- df_48h %>%
    filter(cell_line_modification != wildtype) %>%
    group_by(treatment, cell_line_modification) %>%
    dplyr::summarise(n = n()) %>%
    pull(n) %>%
    max()

  dm_48h <- df_48h %>%
    dplyr::filter(cell_line_modification != wildtype) %>%  # remove WT, because we won't need them for the analysis
    dplyr::group_by(treatment, cell_line_modification) %>%  # cast data to a wide format
    reframe(data.frame(t(c(`log_KO/WT_AUC_Ratio`, rep(NA, max_n_col - length(`log_KO/WT_AUC_Ratio`)))))) %>%  # Since max replicate number is 6, we cast to a matrix with 6 columns
    unite(col = "condition", treatment, cell_line_modification, sep = "_") %>%
    column_to_rownames("condition") %>%
    as.matrix()

  # fit linear model to log2 fold change values for each KO vs WT comparison for each compound
  fm_48h <- lmFit(object = dm_48h)

  # apply empirical Bayes moderation to shrink the standard errors and get moderated statistics from the linear model fits
  fm_48h <- eBayes(fm_48h)

  # retrieve Benjamini-Hochberg corrected p-values
  rslt_48h <- topTable(fm_48h, number = "all", sort.by = "none")


  ### Confluency_after_72h #######################################################

  # subset dataframe for desired phenotype
  df_72h <- AUCsPerConditionDF %>%
    filter(Phenotype == "Confluency_after_72h")

  # determine maximum number of biological replicates used
  max_n_col <- df_72h %>%
    filter(cell_line_modification != wildtype) %>%
    group_by(treatment, cell_line_modification) %>%
    dplyr::summarise(n = n()) %>%
    pull(n) %>%
    max()

  dm_72h <- df_72h %>%
    dplyr::filter(cell_line_modification != wildtype) %>%  # remove WT, because we won't need them for the analysis
    dplyr::group_by(treatment, cell_line_modification) %>%  # cast data to a wide format
    reframe(data.frame(t(c(`log_KO/WT_AUC_Ratio`, rep(NA, max_n_col - length(`log_KO/WT_AUC_Ratio`)))))) %>%  # Since max replicate number is 6, we cast to a matrix with 6 columns
    unite(col = "condition", treatment, cell_line_modification, sep = "_") %>%
    column_to_rownames("condition") %>%
    as.matrix()

  # fit linear model to log2 fold change values for each KO vs WT comparison for each compound
  fm_72h <- lmFit(object = dm_72h)

  # apply empirical Bayes moderation to shrink the standard errors and get moderated statistics from the linear model fits
  fm_72h <- eBayes(fm_72h)

  # retrieve Benjamini-Hochberg corrected p-values
  rslt_72h <- topTable(fm_72h, number = "all", sort.by = "none")


  ### Confluency_after_96h #######################################################

  # subset dataframe for desired phenotype
  df_96h <- AUCsPerConditionDF %>%
    filter(Phenotype == "Confluency_after_96h")

  # determine maximum number of biological replicates used
  max_n_col <- df_96h %>%
    filter(cell_line_modification != wildtype) %>%
    group_by(treatment, cell_line_modification) %>%
    dplyr::summarise(n = n()) %>%
    pull(n) %>%
    max()

  dm_96h <- df_96h %>%
    dplyr::filter(cell_line_modification != wildtype) %>%  # remove WT, because we won't need them for the analysis
    dplyr::group_by(treatment, cell_line_modification) %>%  # cast data to a wide format
    reframe(data.frame(t(c(`log_KO/WT_AUC_Ratio`, rep(NA, max_n_col - length(`log_KO/WT_AUC_Ratio`)))))) %>%  # Since max replicate number is 6, we cast to a matrix with 6 columns
    unite(col = "condition", treatment, cell_line_modification, sep = "_") %>%
    column_to_rownames("condition") %>%
    as.matrix()

  # fit linear model to log2 fold change values for each KO vs WT comparison for each compound
  fm_96h <- lmFit(object = dm_96h)

  # apply empirical Bayes moderation to shrink the standard errors and get moderated statistics from the linear model fits
  fm_96h <- eBayes(fm_96h)

  # retrieve Benjamini-Hochberg corrected p-values
  rslt_96h <- topTable(fm_96h, number = "all", sort.by = "none")



  ### Second_Phase_Relative_Change_in_Confluency #################################

  # subset dataframe for desired phenotype
  df_2nd_phase <- AUCsPerConditionDF %>%
    filter(Phenotype == "Second_Phase_Relative_Change_in_Confluency")

  # determine maximum number of biological replicates used
  max_n_col <- df_2nd_phase %>%
    filter(cell_line_modification !=wildtype) %>%
    group_by(treatment, cell_line_modification) %>%
    dplyr::summarise(n = n()) %>%
    pull(n) %>%
    max()

  dm_2nd_phase <- df_2nd_phase %>%
    dplyr::filter(cell_line_modification != wildtype) %>%  # remove WT, because we won't need them for the analysis
    dplyr::group_by(treatment, cell_line_modification) %>%  # cast data to a wide format
    reframe(data.frame(t(c(`log_KO/WT_AUC_Ratio`, rep(NA, max_n_col - length(`log_KO/WT_AUC_Ratio`)))))) %>%  # Since max replicate number is 6, we cast to a matrix with 6 columns
    unite(col = "condition", treatment, cell_line_modification, sep = "_") %>%
    column_to_rownames("condition") %>%
    as.matrix()

  # fit linear model to log2 fold change values for each KO vs WT comparison for each compound
  fm_2nd_phase <- lmFit(object = dm_2nd_phase)

  # apply empirical Bayes moderation to shrink the standard errors and get moderated statistics from the linear model fits
  fm_2nd_phase <- eBayes(fm_2nd_phase)

  # retrieve Benjamini-Hochberg corrected p-values
  rslt_2nd_phase <- topTable(fm_2nd_phase, number = "all", sort.by = "none")


  ### First_Phase_Duration #################################

  # subset dataframe for desired phenotype
  df_1st_phase_dur <- AUCsPerConditionDF %>%
    filter(Phenotype == "First_Phase_Duration")

  # determine maximum number of biological replicates used
  max_n_col <- df_1st_phase_dur %>%
    filter(cell_line_modification != wildtype) %>%
    group_by(treatment, cell_line_modification) %>%
    dplyr::summarise(n = n()) %>%
    pull(n) %>%
    max()

  dm_1st_phase_dur <- df_1st_phase_dur %>%
    dplyr::filter(cell_line_modification != wildtype) %>%  # remove WT, because we won't need them for the analysis
    dplyr::group_by(treatment, cell_line_modification) %>%  # cast data to a wide format
    reframe(data.frame(t(c(`log_KO/WT_AUC_Ratio`, rep(NA, max_n_col - length(`log_KO/WT_AUC_Ratio`)))))) %>%  # Since max replicate number is 6, we cast to a matrix with 6 columns
    unite(col = "condition", treatment, cell_line_modification, sep = "_") %>%
    column_to_rownames("condition") %>%
    as.matrix()

  # fit linear model to log2 fold change values for each KO vs WT comparison for each compound
  fm_1st_phase_dur <- lmFit(object = dm_1st_phase_dur)

  # apply empirical Bayes moderation to shrink the standard errors and get moderated statistics from the linear model fits
  fm_1st_phase_dur <- eBayes(fm_1st_phase_dur)

  # retrieve Benjamini-Hochberg corrected p-values
  rslt_1st_phase_dur <- topTable(fm_1st_phase_dur, number = "all", sort.by = "none")


  ### First_Phase_Change_in_Confluency #################################

  df_1st_phase <- AUCsPerConditionDF %>%
    filter(Phenotype == "First_Phase_Change_in_Confluency")

  # determine maximum number of biological replicates used
  max_n_col <- df_1st_phase %>%
    filter(cell_line_modification != wildtype) %>%
    group_by(treatment, cell_line_modification) %>%
    dplyr::summarise(n = n()) %>%
    pull(n) %>%
    max()

  dm_1st_phase <- df_1st_phase %>%
    dplyr::filter(cell_line_modification != wildtype) %>%  # remove WT, because we won't need them for the analysis
    dplyr::group_by(treatment, cell_line_modification) %>%  # cast data to a wide format
    reframe(data.frame(t(c(`log_KO/WT_AUC_Ratio`, rep(NA, max_n_col - length(`log_KO/WT_AUC_Ratio`)))))) %>%  # Since max replicate number is 6, we cast to a matrix with 6 columns
    unite(col = "condition", treatment, cell_line_modification, sep = "_") %>%
    column_to_rownames("condition") %>%
    as.matrix()

  # fit linear model to log2 fold change values for each KO vs WT comparison for each compound
  fm_1st_phase <- lmFit(object = dm_1st_phase)

  # apply empirical Bayes moderation to shrink the standard errors and get moderated statistics from the linear model fits
  fm_1st_phase <- eBayes(fm_1st_phase)

  # retrieve Benjamini-Hochberg corrected p-values
  rslt_1st_phase <- topTable(fm_1st_phase, number = "all", sort.by = "none")


  ### Merge the various phenotype-specific dataframes and compute mean log2 fold-change for each condition (KO vs WT for a given compound) #################################
  Significance24h <- tibble::rownames_to_column(rslt_24h, "Condition")
  Significance24h <- Significance24h %>% mutate(Condition = paste(Condition, "Confluency_after_24h", sep = "/"))

  Significance48h <- tibble::rownames_to_column(rslt_48h, "Condition")
  Significance48h <- Significance48h %>% mutate(Condition = paste(Condition, "Confluency_after_48h", sep = "/"))

  Significance72h <- tibble::rownames_to_column(rslt_72h, "Condition")
  Significance72h <- Significance72h %>% mutate(Condition = paste(Condition, "Confluency_after_72h", sep = "/"))

  Significance96h <- tibble::rownames_to_column(rslt_96h, "Condition")
  Significance96h <- Significance96h %>% mutate(Condition = paste(Condition, "Confluency_after_96h", sep = "/"))

  Significance2ndPhase_relChg <- tibble::rownames_to_column(rslt_2nd_phase, "Condition")
  Significance2ndPhase_relChg <- Significance2ndPhase_relChg %>% mutate(Condition = paste(Condition, "Second_Phase_Relative_Change_in_Confluency", sep = "/"))

  Significance1stPhase_Dur <- tibble::rownames_to_column(rslt_1st_phase_dur, "Condition")
  Significance1stPhase_Dur <- Significance1stPhase_Dur %>% mutate(Condition = paste(Condition, "First_Phase_Duration", sep = "/"))

  Significance1stPhase_Chg <- tibble::rownames_to_column(rslt_1st_phase, "Condition")
  Significance1stPhase_Chg <- Significance1stPhase_Chg %>% mutate(Condition = paste(Condition, "First_Phase_Change_in_Confluency", sep = "/"))


  data_frame_merge <- rbind(Significance24h, Significance48h)
  data_frame_merge <- rbind(data_frame_merge, Significance72h)
  data_frame_merge <- rbind(data_frame_merge, Significance96h)
  data_frame_merge <- rbind(data_frame_merge, Significance2ndPhase_relChg)
  data_frame_merge <- rbind(data_frame_merge, Significance1stPhase_Dur)
  data_frame_merge <- rbind(data_frame_merge, Significance1stPhase_Chg)

  df2 <- AUCsPerConditionDF[!(AUCsPerConditionDF$cell_line_modification == wildtype),]

  # create condition column label
  df2 <- df2 %>%
    mutate(Condition = paste(treatment, cell_line_modification, sep = "_"))

  # compute mean log2 fold-change per condition
  df2 <- df2 %>%
    mutate(Condition = paste(Condition,Phenotype, sep = "/")) %>% dplyr::group_by(Condition,treatment, CL_group, Phenotype,cell_line_modification) %>%
    dplyr::summarise_at("log_KO/WT_AUC_Ratio",list(mean = mean, sd = sd,se = ~ sd(.) / sqrt(length(.))),na.rm=TRUE)

  # merge dataframe containing p-values to dataframe containing mean log2 fold-change per condition
  CompleteAUCDF <- full_join(data_frame_merge, df2,by= 'Condition')

  # rename certain columns
  CompleteAUCDF <-CompleteAUCDF%>% dplyr::rename(mean_log_AUC_value_ratio = mean,
                                                 p_value= adj.P.Val)


  return(CompleteAUCDF)

}



#########################################################################################################################################################################################################################################################################
#' Create Summary Matrix of AUC Values

#'
#' This function creates a matrix of AUC values where the 5 phenotypes ("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency") are the column names and the conditions (cell line & compound) are the rownames.
#' @param CompleteAUCDF A dataframe containing treatment, cell line modification, mean log2 fold-change KO vs WT AUC ratio for this given phenotype. The KO is compared to the WT from the same biological replicate. These ratios are then averaged across the biological replicates.
#' @return
#'\item{AUCDF5PhenotypesDF}{ A matrix containing the mean AUC values normalized to WT, cell line modification, and treatment for the 5 phenotypes.}

#' @author Caroline Barry
#' @import dplyr
#' @export

getSummaryAUCMatrix <- function(CompleteAUCDF){

  # subset dataframe for desired phenotypes
  AUCDFSubset <- subset(CompleteAUCDF, Phenotype %in% c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"))

  # then subset for each phenotype individually
  Confluency_after_24h_DF <- subset(AUCDFSubset, Phenotype=="Confluency_after_24h")
  Confluency_after_24h_DF <- Confluency_after_24h_DF %>% dplyr::rename(Confluency_after_24h=mean_log_AUC_value_ratio)

  Confluency_after_48h_DF <- subset(AUCDFSubset, Phenotype=="Confluency_after_48h")
  Confluency_after_48h_DF <- Confluency_after_48h_DF %>% dplyr::rename(Confluency_after_48h=mean_log_AUC_value_ratio)

  Confluency_after_72h_DF <- subset(AUCDFSubset, Phenotype=="Confluency_after_72h")
  Confluency_after_72h_DF <- Confluency_after_72h_DF %>% dplyr::rename(Confluency_after_72h=mean_log_AUC_value_ratio)

  Confluency_after_96h_DF <- subset(AUCDFSubset, Phenotype=="Confluency_after_96h")
  Confluency_after_96h_DF <- Confluency_after_96h_DF %>% dplyr::rename(Confluency_after_96h=mean_log_AUC_value_ratio)

  Second_Phase_Relative_Change_in_Confluency_DF <- subset(AUCDFSubset, Phenotype=="Second_Phase_Relative_Change_in_Confluency")
  Second_Phase_Relative_Change_in_Confluency_DF <- Second_Phase_Relative_Change_in_Confluency_DF %>% dplyr::rename(Second_Phase_Relative_Change_in_Confluency=mean_log_AUC_value_ratio)

  # merge the AUC value for each cell line and compound across the phenotypes
  AUCDF5PhenotypesDF<- merge(Confluency_after_24h_DF[,c(8,11,12)], Confluency_after_48h_DF[,c(8,11,12)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  AUCDF5PhenotypesDF <- merge(AUCDF5PhenotypesDF, Confluency_after_72h_DF[,c(8,11,12)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  AUCDF5PhenotypesDF <- merge(AUCDF5PhenotypesDF, Confluency_after_96h_DF[,c(8,11,12)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  AUCDF5PhenotypesDF <- merge(AUCDF5PhenotypesDF, Second_Phase_Relative_Change_in_Confluency_DF[,c(8,11,12)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  AUCDF5PhenotypesDF <- AUCDF5PhenotypesDF %>% mutate(CL_Treat = paste(cell_line_modification,"|",treatment)) %>% mutate(CL_Treat=gsub(" ", "", CL_Treat))

  # assign concatenation of cell line modification and compound as the rownames
  rownames(AUCDF5PhenotypesDF) <- AUCDF5PhenotypesDF$CL_Treat

  # filter for desired columns
  AUCDF5PhenotypesDF <- AUCDF5PhenotypesDF[c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency","cell_line_modification","treatment")]

  return(AUCDF5PhenotypesDF)
}

#########################################################################################################################################################################################################################################################################

#' Create Summary Matrix of P-Values

#'
#' This function gets p-values for 5 phenotypes: Confluency_after_24h, Confluency_after_48h, Confluency_after_72h, Confluency_after_96h, Second_Phase_Relative_Change_in_Confluency. These p-values correspond
#' with the t-test conducted on the AUC values of WT and KO across their biological replicates.
#' @param AUCSignificanceDF A dataframe containing treatment, cell line modification, Welsh t-statistic and Holmes corrected p-value for each given phenotype.
#' @return
#'\item{PValue5PhenotypesDF}{ A matrix containing the p-values for the 5 phenotypes for each condition (cell line and compound).}

#' @author Caroline Barry
#' @import dplyr
#' @export

getSummaryPValueMatrix <- function(AUCSignificanceDF){

  # subset dataframe for desired phenotypes
  AUCSignificanceDFSubset <- subset(AUCSignificanceDF, Phenotype %in% c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"))

  # then subset for each phenotype individually
  Confluency_after_24h_DF <- subset(AUCSignificanceDFSubset, Phenotype=="Confluency_after_24h")
  Confluency_after_24h_DF <- Confluency_after_24h_DF %>% dplyr::rename(Confluency_after_24h=p_value)

  Confluency_after_48h_DF <- subset(AUCSignificanceDFSubset, Phenotype=="Confluency_after_48h")
  Confluency_after_48h_DF <- Confluency_after_48h_DF %>% dplyr::rename(Confluency_after_48h=p_value)

  Confluency_after_72h_DF <- subset(AUCSignificanceDFSubset, Phenotype=="Confluency_after_72h")
  Confluency_after_72h_DF <- Confluency_after_72h_DF %>% dplyr::rename(Confluency_after_72h=p_value)

  Confluency_after_96h_DF <- subset(AUCSignificanceDFSubset, Phenotype=="Confluency_after_96h")
  Confluency_after_96h_DF <- Confluency_after_96h_DF %>% dplyr::rename(Confluency_after_96h=p_value)

  Second_Phase_Relative_Change_in_Confluency_DF <- subset(AUCSignificanceDFSubset, Phenotype=="Second_Phase_Relative_Change_in_Confluency")
  Second_Phase_Relative_Change_in_Confluency_DF <- Second_Phase_Relative_Change_in_Confluency_DF %>% dplyr::rename(Second_Phase_Relative_Change_in_Confluency=p_value)

  # merge the AUC value for each cell line and compound across the phenotypes
  PValue5PhenotypesDF<- merge(Confluency_after_24h_DF[,c(8,11,6)], Confluency_after_48h_DF[,c(8,11,6)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  PValue5PhenotypesDF <- merge(PValue5PhenotypesDF, Confluency_after_72h_DF[,c(8,11,6)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  PValue5PhenotypesDF <- merge(PValue5PhenotypesDF, Confluency_after_96h_DF[,c(8,11,6)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  PValue5PhenotypesDF <- merge(PValue5PhenotypesDF, Second_Phase_Relative_Change_in_Confluency_DF[,c(8,11,6)], by = c("treatment","cell_line_modification"), all.x = TRUE, all.y = TRUE)

  PValue5PhenotypesDF <- PValue5PhenotypesDF %>% mutate(CL_Treat = paste(cell_line_modification,"|",treatment)) %>% mutate(CL_Treat=gsub(" ", "", CL_Treat))

  # assign concatenation of cell line modification and compound as the rownames
  rownames(PValue5PhenotypesDF) <- PValue5PhenotypesDF$CL_Treat

  # filter for desired columns
  PValue5PhenotypesDF <- PValue5PhenotypesDF[c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency","cell_line_modification","treatment")]

  return(PValue5PhenotypesDF)
}




