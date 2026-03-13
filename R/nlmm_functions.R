## Extracting p-values

#' Fitting Data to Non-linear Mixed Effects Model (NLMM) and Plotting the Fits and Saving the Model Estimates to MySQL I - Untreated Interaction Parameter Test

#'
#' This function calls a function that fits only untreated data in a pairwise fashion (KO untreated & WT untreated) to NLMM with logistic growth function, resulting in 3 model estimates: initial confluency, growth rate, and maximum capacity.
#'  Additionally, this function calls a function that saves these estimates in the MySQL database. This function also calls a function that plots the model fits.
#' @param IncucyteDataAndMetaDataAndConditionsLabeledDF Dataframe that combines Incucyte Confluency Dataset, Meta Data, Condition IDs and analysis IDs.
#' @param treatmentendtime Numerical variable represents the time at which treatment ends.
#' @param ListOfControlCL List mapping each cell line to its respective control cell line.
#' @param ListOfNLMMFitFailed List containing IDs of analysis groups that failed to fit to the model.
#' @param ListOfFailedParameterExtract  List containing IDs of analysis groups where model parameters weren't able of being extracted.
#' @param ListOfFailedPlots  List containing IDs of analysis groups where model fits were unable to be plotted.
#' @param TodaysDate The date.
#' @param AnalysisName A variable representing user-defined name of analysis.


#' @return
#'\item{AllFailedPairedUntreatIDs}{ List of all NLMM Analysis IDs of analysis groups that didn't fit NLMM at the first attempt.}
#' @author Caroline Barry
#' @export

runNLMM<- function(IncucyteDatasetSubset,treatmentendtime,ListOfControlCL, ListOfNLMMFitFailed,ListOfFailedParameterExtract,ListOfFailedPlots,TodaysDate,AnalysisName){

  NumberOfBioReps<-length(as.list(unique(IncucyteDatasetSubset$bio_rep_id)))

  NumbersOfBioReps <- NumberOfBioReps + 1


  IncucyteDatasetSubsetControlCL <- subset(IncucyteDatasetSubset, cell_line_modifications %in% ListOfControlCL)

  ControlCellLine <- unique(IncucyteDatasetSubsetControlCL$cell_line_modifications)

  #find a way to make it also work for clone_nlmm_analysis
  StringPortion <- paste(AnalysisName,"_nlmm_analysis_")

  StringPortion <- gsub(" ", "",StringPortion)

  ConditionValue <- as.numeric(gsub(StringPortion, "", unique(IncucyteDatasetSubset$NLMM_Analysis_ID)))

  KOConditionID <- paste("condition_", ConditionValue,sep="")

  KOCompoundDFSubsetOriginal <- subset(IncucyteDatasetSubset, condition_ID == KOConditionID)

  CompoundName <-unique(KOCompoundDFSubsetOriginal$treatment_0)


  KOCellLine <-unique(KOCompoundDFSubsetOriginal$cell_line_modifications)

  IncucyteDatasetSubset <- IncucyteDatasetSubset %>%
    dplyr::group_by(bio_rep_id) %>%
    dplyr::filter(timepoint<=treatmentendtime)


  IncucyteDatasetSubset <- IncucyteDatasetSubset %>%
    mutate(
      cell_line_modifications = factor(cell_line_modifications, levels = c(ControlCellLine, KOCellLine)),
      treatment_0 = factor(treatment_0),
      plate_ID = factor(plate_ID),
      bio_rep_id = factor(bio_rep_id)
    )


  unique(IncucyteDatasetSubset$condition_ID)

  # Create 'groupedData' object
  #model per cell line
  DF <- groupedData(mean_conf ~ timepoint| bio_rep_id, #
                    data = IncucyteDatasetSubset,
                    labels = list(x = "Time", y = "Confluence"),
                    units = list(x = "(h)", y = "(%)"))


  #get number of conditions that are currently being fitted

  NumberOfConditions <- length(unique(DF$condition_ID))

  NumbersOfConditions <- NumberOfConditions -1

  ############################### Fitting Paired Untreated Data to NLMM ###############################

  possibleError <- tryCatch({withCallingHandlers({

    fmm <- runNLMMFitting(DF)

  }, error=function(e){


    e

  },
  warning = function(w) {
    w
  }

  )},error=function(e){

    print("Attempting second model fitting...")

    possibleError <- tryCatch({withCallingHandlers({

      plot(DF)
      fmm <- runNLMMFittingFortTest(DF)


    }, error=function(e){
      e
    },

    warning = function(w) {

      w

    }

    )},error=function(e){
      print(unique(DF$NLMM_Analysis_ID))
      cat("ERROR :",conditionMessage(e), "\n")
      e

      ListOfNLMMFitFailed <<- c(ListOfNLMMFitFailed, list(unique(DF$NLMM_Analysis_ID)))
    }

    )

  })

  if(unique(DF$NLMM_Analysis_ID)%in% unlist(possibleError)==FALSE){

    ############################### Extracting p-values & Uploading to MySQL ###############################

    possibleError2 <- tryCatch({withCallingHandlers({

      InteractionEstimateTable <- createNLMMTable(fmm, DF,TodaysDate)

    }, error=function(e){

      print(unique(DF$NLMM_Analysis_ID))
      {cat("ERROR :",conditionMessage(e), "\n")
        e

        ListOfFailedParameterExtract <<- c(ListOfFailedParameterExtract, list(unique(DF$NLMM_Analysis_ID)))
      }

    },

    warning = function(w) {

      w

    }
    )},error=function(e){
      print(unique(DF$NLMM_Analysis_ID))
      cat("ERROR :",conditionMessage(e), "\n")
      e

      ListOfFailedParameterExtract <<- c(ListOfFailedParameterExtract, list(unique(DF$NLMM_Analysis_ID)))
    }
    )
  }

  if(isTRUE(!(unique(DF$NLMM_Analysis_ID)%in% unlist(possibleError) & unique(DF$NLMM_Analysis_ID)%in% unlist(possibleError2)))==TRUE){
    ############################### Fitting Paired Untreated Data to NLMM ###############################

    possibleError3 <- tryCatch({withCallingHandlers({
      fitted_plot<- plotNLMMModelFits(fmm, DF, NumberOfBioReps, NumberOfConditions, IncucyteDatasetSubset,ControlCellLine)

    }, error=function(e){

      print(unique(DF$NLMM_Analysis_ID))
      {cat("ERROR :",conditionMessage(e), "\n")
        e

        ListOfFailedPlots <<- c(ListOfFailedPlots, list(unique(DF$NLMM_Analysis_ID)))
      }

    },

    warning = function(w) {
      w
    }
    )},error=function(e){
      print(unique(DF$NLMM_Analysis_ID))
      cat("ERROR :",conditionMessage(e), "\n")
      e

      ListOfFailedPlots <<- c(ListOfFailedPlots, list(unique(DF$NLMM_Analysis_ID)))
    }
    )

  }

  ListOfNLMMtTestResults <- list(NLMMOutput = fmm,ListOfNLMMFitFailed = ListOfNLMMFitFailed,
                                 NLMMEstimatesTable = InteractionEstimateTable, ListOfFailedParameterExtract = ListOfFailedParameterExtract,
                                 ModelFitPlot = fitted_plot, ListOfFailedPlots=ListOfFailedPlots)

  return(ListOfNLMMtTestResults)
}

#########################################################################################################################################################################################################################################################################

#' Fits Data to Non-linear Mixed Effects Model Using Logistic Growth Function I - Interaction Parameter Test

#'
#' This function fits data to the Interaction Parameter Test - NLMM in an analysis group-wise manner (KO + compound, WT + compound, KO + control, WT + control) and returns model estimates:
#' initial confluency (y0), maximum capacity (K), and growth rate (r). The parameters K and r are on a log scale. This function computes random effects per biological replicate. The confidence intervals
#'  will be used to deduce significance of conditions when comparing KO to WT in the heatmap for the NLMM model estimates.
#' @param DF A subset of a dataframe containing confluency data, condition IDs, NLMM Analysis IDs, and Meta Data. This subset only contains one unique nlmm analysis group.
#' @details This function uses dependent random effects: correlation structure across parameters with regards to the random effects; Fridolin saw negative correlation between low y0 and high growth rate

#' @return
#'\item{fmm}{ A matrix containing the model estimates from the NLMM Interaction Parameter Test.}
#' @author Caroline Barry & Fridolin Kielisch

runNLMMFitting <- function(DF){

  fm_lis <- nlsList(mean_conf ~ logSSlogisGrowth(timepoint, l2K, y0, l2r), data = DF)

  ctrl <- lmeControl(maxIter=10000, niterEM=10000)

  fmm <- nlme(fm_lis, control = ctrl)

  strt <- c(fixef(fmm)[1], rep(0, 1),
            fixef(fmm)[2], rep(0, 1),
            fixef(fmm)[3], rep(0, 1))

  fmm <- update(fmm,
                fixed = l2K + y0 + l2r ~ cell_line_modifications,
                random = l2K + y0 + l2r ~ 1,
                start = strt)

  return(fmm)

}

#########################################################################################################################################################################################################################################################################
#' Creates NLMM Estimates Table - Untreated NLMM Interaction Parameter Test

#'
#' This function creates a dataframe of the model estimates per analysis group that are extracted from fmm and uploads this to MySQL.
#' @param fmm A matrix containing the model estimates from the NLMM Interaction Parameter Test.
#' @param DF A subset of a dataframe containing confluency data, condition IDs, NLMM Analysis IDs, and Meta Data. This subset only contains one unique nlmm analysis group.
#' @param TodaysDate The date.

#' @return
#'\item{condition}{ Unique condition ID from DF dataframe.}
#' @author Caroline Barry

createNLMMTable <- function(fmm, DF, TodaysDate){

  print(paste0("Will extract parameters of: ",unique(DF$NLMM_Analysis_ID)))

  boo <- summary(fmm)

  goo<- boo$tTable

  l2k_pval <- t(as.data.frame(goo[2,c(1,2,5)]))

  colnames(l2k_pval) <- c("L2FC_K", "K_SE", "K_pvalue")

  # assigning the rownames to null
  rownames(l2k_pval) <- NULL

  y0_pval <- t(as.data.frame(goo[4,c(1,2,5)]))

  colnames(y0_pval) <- c("y0", "y0_SE", "y0_pvalue")

  # assigning the rownames to null
  rownames(y0_pval) <- NULL

  l2r_pval <- t(as.data.frame(goo[6,c(1,2,5)]))

  colnames(l2r_pval) <- c("L2FC_r", "r_SE", "r_pvalue")

  # assigning the rownames to null
  rownames(l2r_pval) <- NULL

  #delete irrelevant columns (degrees of freedom)
  InteractionEstimateTable <- data.frame()

  InteractionEstimateTable <- rbind(InteractionEstimateTable, y0_pval)

  InteractionEstimateTable <- cbind(InteractionEstimateTable, l2r_pval)

  InteractionEstimateTable <- cbind(InteractionEstimateTable, l2k_pval)

  # colnames(InteractionEstimateTable) <- c("y0_interaction_pvalue", "r_interaction_pvalue", "K_interaction_pvalue")


  #wanna take nlmm_id and paste it to each row
  InteractionEstimateTable <- InteractionEstimateTable %>%
    mutate(NLMM_Analysis_ID = unique(DF$NLMM_Analysis_ID)) %>%
    mutate(date = TodaysDate)

  InteractionEstimateTable <- InteractionEstimateTable[,c(10,11,1:9)]

  return(InteractionEstimateTable)

}

#########################################################################################################################################################################################################################################################################

#' #' Plot the NLMM Model Fits - Untreated Interaction Parameter Test
#'
#' #'
#' #' This function plots the NLMM Model Fits for each analysis group (KO + compound, KO + control, WT + compound, WT + control) fitted to the model along with the confluency data for the biological replicates that were used to generate these model estimates.
#' #' @param fmm A matrix containing the model estimates from the NLMM Interaction Parameter Test.
#' #' @param DF A subset of a dataframe containing confluency data, condition IDs, NLMM Analysis IDs, and Meta Data. This subset only contains one unique nlmm analysis group.
#' #' @param NumberOfBioReps A value representing the number of biological replicates in DF.
#' #' @param NumberOfConditions A value representing the number of unique conditions in DF.
#' #' @param IncucyteDatasetSubsetOriginal Dataframe containing confluency data, condition IDs, NLMM Analysis IDs, and Meta Data where all timepoints are included.
#' #' @param IncucyteDatasetSubset Dataframe containing confluency data, condition IDs, NLMM Analysis IDs, and Meta Data where only timepoints before the split point are included.
#' #' @param ControlCellLine Cell line the knockout is being compared to.
#'
#' #' @author Caroline Barry & Fridolin Kielisch
#' #' @export
#'
plotNLMMModelFits <- function(fmm, DF, NumberOfBioReps, NumberOfConditions, IncucyteDatasetSubset,ControlCellLine){

  mm_K <- as.data.frame(emmeans(fmm, specs = c("cell_line_modifications"), param = "l2K", data = DF))
  mm_r <- as.data.frame(emmeans(fmm, specs = c("cell_line_modifications"), param = "l2r", data = DF))
  mm_y0 <- as.data.frame(emmeans(fmm, specs = c("cell_line_modifications"), param = "y0", data = DF))#

  print(paste0("Will plot fitted curves of: ",unique(DF$NLMM_Analysis_ID)))


  IncucyteDatasetSubsetWT <- subset(IncucyteDatasetSubset, cell_line_modifications == ControlCellLine)

  IncucyteDatasetSubsetKO <- subset(IncucyteDatasetSubset, cell_line_modifications != ControlCellLine)

  CondList <- c(unique(IncucyteDatasetSubsetWT$condition_ID), unique(IncucyteDatasetSubsetKO$condition_ID))

  #delete irrelevant columns (degrees of freedom)
  fixedK <- mm_K[,2]
  names(fixedK) <-paste(CondList)
  fixedr <- mm_r[,2]
  names(fixedr) <- paste(CondList)
  fixedy0 <- mm_y0[,2]
  names(fixedy0) <- paste(CondList)

  n<-120

  re <- ranef(fmm)

  re <- re[order(match(rownames(re), unique(IncucyteDatasetSubset$bio_rep_id ))), , drop = FALSE]

  df_plot <- data.frame(t = rep(seq(0, n, length.out = n), as.numeric(NumberOfBioReps)))#bioreps #lines for each biorep of each condition
  df_plot_cnd <- data.frame(t = rep(seq(0, n, length.out = n), as.numeric(NumberOfConditions)))#12) #model estimate line

  strsplit(rownames(re), split = "|", fixed = TRUE) |>
    unlist() |>
    matrix(ncol = 3, byrow = T) -> cndt
  cndt <- cndt[,1] # condtions in the order they appear in ranef


  df_plot$cond <- rep(cndt, each = n)
  df_plot_cnd$cond <- rep(unique(cndt), each = n)

  df_plot$biorep <- rep(rownames(re), each = n)


  K <- fixedK[cndt]
  r <- fixedr[cndt]
  y0 <- fixedy0[cndt]

  K <- K + re[,1]
  r <- r + re[,3]
  y0 <- y0 + re[,2]


  #was 45 but changed to 51 for DMSO
  t <- seq(0, n, length.out = n)
  conf <- as.numeric(NumberOfBioReps)*n
  for(i in 1:as.numeric(NumberOfBioReps)) {
    conf[((i-1)*n + 1):(i*n)] <- logSSlogisGrowth(t, l2K = K[i], y0 = y0[i], l2r = r[i])
  }

  df_plot$conf <- conf

  fixedK_chgOrder <- fixedK[unique(cndt)]
  fixedr_chgOrder <- fixedr[unique(cndt)]
  fixedy0_chgOrder <- fixedy0[unique(cndt)]


  conf <- as.numeric(NumberOfConditions)*n
  for(i in 1:as.numeric(NumberOfConditions)) {
    conf[((i-1)*n + 1):(i*n)] <- logSSlogisGrowth(t, l2K = fixedK_chgOrder[i],
                                                  y0 = fixedy0_chgOrder[i],
                                                  l2r = fixedr_chgOrder[i])
  }
  df_plot_cnd$conf <- conf


  df_avg_conf_plot <- IncucyteDatasetSubset
  colnames(df_avg_conf_plot)[3] <- "cond"
  if(nrow(IncucyteDatasetSubset)>0){
    colnames(IncucyteDatasetSubset)[3] <- "cond"}

  #wanna add condition name labels!

  cndt_names_df <- DF %>%
    dplyr::group_by(condition_ID) %>%
    slice(n())

  cndt_names_df <- cndt_names_df %>%
    dplyr::rename(cond= condition_ID)


  cndt2 <- cndt

  cndt_vals <-sort(as.numeric(gsub("condition_","", cndt2)))

  cndt_vals <- as.list(unique(cndt_vals))

  cndt3<- list()

  for(v in 1:length(cndt_vals)){
    cndt3[v]<- paste("condition_", cndt_vals[v], sep="")

  }

  cndt_names_df$cond <- factor(cndt_names_df$cond,
                               levels=cndt3)
  cndt_names <- cndt_names_df$CL_treat_conc_id

  names(cndt_names) <- cndt3

  cndt_names <- as.list(cndt_names)

  cndt_labeller_fun <- function(variable,value){
    return(cndt_names[value])
  }


  cndt_names <- unlist(cndt_names)


  #create one big df

  df_plot <- df_plot %>%
    dplyr::rename(df_plot_t = t,
                  df_plot_conf = conf,
                  bio_rep_id = biorep
    )


  df_plot_cnd <- df_plot_cnd %>%
    dplyr::rename(df_plot_cnd_t = t,
                  df_plot_cnd_conf = conf)



  #plot each condition separately

  ConditionList <- unique(df_plot$cond)

  fitted_plot <-
    ggplot() +
    geom_line(data=df_plot,aes(x = df_plot_t, y = df_plot_conf, color = bio_rep_id),size = 3) +
    geom_line(data= df_plot_cnd,aes(x = df_plot_cnd_t, y = df_plot_cnd_conf), color = "black") +
    geom_point(data = df_avg_conf_plot,aes(x = timepoint, y = mean_conf, color = bio_rep_id),size = 5) +
    labs(x="Timepoint (hr)", y="Confluency (%)")+
    facet_wrap(facets = vars(cond),labeller=labeller(
      cond = cndt_labeller_fun))+
    theme(strip.text.x = element_text(size = 22),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          axis.text=element_text(size=20),
          axis.title=element_text(size=24,face="bold"),
          plot.title = element_text(size=24,face="bold"))

  fitted_plot<- ggarrange(fitted_plot, legend="bottom")

  return(fitted_plot)

}
