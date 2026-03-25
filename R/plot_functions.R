#' Plot Growth Curves for each Analysis Group

#' This function plots the growth curves without split points of one given analysis group (KO + treatment, KO + control, WT + treatment, and WT + control) across all biological replicates for one given condition (treatment & concentration) together on the same graph.
#' @param AvgIncucyteDataAndMetaDataAndConditionsLabeledDF Dataframe that combines Incucyte Confluency Dataset, Meta Data, Condition IDs and analysis IDs. The confluency is averaged across technical replicates.
#' @param CLTRTList List containing user-specified cell lines and compounds to plot. These are separated with "|" (i.e. cell line | compound).
#' @param ColorListTreat List of colors for treated conditions. These will be used for the line coloring in the plots.
#' @param ColorListControl List of colors for control conditions. These will be used for the line coloring in the plots.
#' @param wildtype Variable representing name of wildtype cell line.
#' @param treatmentendtime Numerical variable represents the time at which treatment ends.
#' @param AnalysisValue The date.
#' @param LinetypeListTreat List of line types for treated conditions. Default is NULL.
#' @param LinetypeListControl List of line types for control conditions. Default is NULL.
#'
#' @return
#'\item{PlotObject}{ Saved ggplot object of the growth curves.}
#'
#' @author Caroline Barry
#' @import dplyr
#' @import ggplot2
#' @export

plotGrowthCurves<- function(IncucyteDataAndMetaDataAndConditionsSubset,CLTRTList,ColorListTreat,ColorListControl,wildtype,treatmentendtime, AnalysisValue,LinetypeListTreat=NULL,LinetypeListControl=NULL){

  IncucyteDataAndMetaDataAndConditionsSubsetValue <- subset(IncucyteDataAndMetaDataAndConditionsSubset, IncucyteDataAndMetaDataAndConditionsSubset$condition_ID == paste("condition_",AnalysisValue,sep=""))

  CompoundName <- unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$treatment)

  ControlNameSubset <- IncucyteDataAndMetaDataAndConditionsSubset %>% dplyr::filter(treatment!=CompoundName)

  ControlName<- unique(ControlNameSubset$treatment)

  CompoundNameConcentration <-gsub("_"," ",unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$treat_conc_id))

  ControlNameConcentration <-gsub("_"," ",unique(ControlNameSubset$treat_conc_id))

  CellLineModification <- unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$cell_line_modifications)

  ListOfCLCmpnd <- paste(CLTRTList, CompoundName, sep=" | ")

  names(ColorListTreat) <- ListOfCLCmpnd

  ListOfCLCntrl <- paste(CLTRTList, ControlName, sep=" | ")

  names(ColorListControl) <- ListOfCLCntrl

  ColorKey <- c(ColorListTreat,ColorListControl)

  if(isTRUE(is.null(LinetypeListTreat) &is.null(LinetypeListControl))==FALSE){
    names(LinetypeListTreat) <- ListOfCLCmpnd
    names(LinetypeListControl) <- ListOfCLCntrl
    LineKey <- c(LinetypeListTreat,LinetypeListControl)
  }


  ListOfPlateIDs = unique(IncucyteDataAndMetaDataAndConditionsSubset$plate_ID)

  # create a unique identifier for observation
  IncucyteDataAndMetaDataAndConditionsSubset <- IncucyteDataAndMetaDataAndConditionsSubset %>%
    mutate(uniq_id=paste(CL_treat_conc_id,plate_ID,sep="|"))

  # create mapping for unique id to color, depends on RSOperator
  MAPPING <- IncucyteDataAndMetaDataAndConditionsSubset %>% distinct(CL_treat_conc_id,plate_ID,uniq_id)
  RS_COLS = brewer.pal(8, "Set2")
  RS_COLS = RS_COLS[1:n_distinct(MAPPING$CL_treat_conc_id)]
  names(RS_COLS) = unique(MAPPING$CL_treat_conc_id)
  PLOT_COLS = RS_COLS[MAPPING$CL_treat_conc_id]
  names(PLOT_COLS) = MAPPING$uniq_id

  IncucyteDataAndMetaDataAndConditionsSubsetPlate <- IncucyteDataAndMetaDataAndConditionsSubset %>%
    dplyr::group_by(bio_rep_id) %>%
    dplyr::filter(timepoint<=treatmentendtime)

  if(isTRUE(is.null(LinetypeListTreat) &is.null(LinetypeListControl))==FALSE){
    PlotObject <- ggplot(IncucyteDataAndMetaDataAndConditionsSubsetPlate,mapping=aes(x=timepoint, y=mean,group=interaction(plate_ID,CL_treat_conc_id),col=CL_TRT, linetype = CL_TRT)) +
      theme_classic() +
      geom_line(aes(x=timepoint,y=mean,  group=interaction(plate_ID,CL_TRT),col = CL_TRT),size = 7)  +
      scale_linetype_manual("Cell line | Treatment",breaks = c(paste(CellLineModification,gsub("_"," ",CompoundName),sep=" | "),paste(wildtype,gsub("_"," ",CompoundName),sep=" | "),
                                                               paste(CellLineModification,ControlName,sep=" | "),paste(wildtype,ControlName,sep=" | ")),values = LineKey)+
      ylim(0.0,100)+
      scale_color_manual("Cell line | Treatment", breaks = c(paste(CellLineModification,gsub("_"," ",CompoundName),sep=" | "),paste(wildtype,gsub("_"," ",CompoundName),sep=" | "),paste(CellLineModification,gsub("_"," ",ControlName),sep=" | "),paste(wildtype,gsub("_"," ",ControlName),sep=" | ")), values=ColorKey)+
      scale_fill_manual("Cell line | Treatment", breaks =c(paste(CellLineModification,gsub("_"," ",CompoundName),sep=" | "),paste(wildtype,gsub("_"," ",CompoundName),sep=" | "),paste(CellLineModification,gsub("_"," ",ControlName),sep=" | "),paste(wildtype,gsub("_"," ",ControlName),sep=" | ")), values=ColorKey)+
      labs(title = paste(IncucyteDataAndMetaDataAndConditionsSubsetPlate$cell_line, CellLineModification, IncucyteDataAndMetaDataAndConditionsSubsetPlate$cell_number, '-', CompoundNameConcentration, '\n'), col='Treatment_Concentration_Bio_rep') +
      labs(x="\nTimepoint (hr)", y="Confluency (%)\n")+
      theme(axis.text=element_text(size=80),
            plot.margin = margin(1,1,1,1, "cm"),
            legend.text = element_text(size = 58),
            axis.title=element_text(size=90),
            legend.title=element_text(size=60),
            plot.title = element_text(size=85),
            axis.line = element_line(colour = 'black', linewidth = 4),
            axis.ticks = element_line(colour = "black", linewidth = 4),
            legend.box.background = element_rect(colour = "black", linewidth = 4))
  } else{
    PlotObject <- ggplot(IncucyteDataAndMetaDataAndConditionsSubsetPlate,mapping=aes(x=timepoint, y=mean,group=interaction(plate_ID,CL_treat_conc_id),col=CL_TRT, linetype = CL_TRT)) +
      theme_classic() +
      geom_line(aes(x=timepoint,y=mean,  group=interaction(plate_ID,CL_TRT),col = CL_TRT),size = 7)  +
      ylim(0.0,100)+
      scale_color_manual("Cell line | Treatment", breaks = c(paste(CellLineModification,gsub("_"," ",CompoundName),sep=" | "),paste(wildtype,gsub("_"," ",CompoundName),sep=" | "),paste(CellLineModification,gsub("_"," ",ControlName),sep=" | "),paste(wildtype,gsub("_"," ",ControlName),sep=" | ")), values=ColorKey)+
      scale_fill_manual("Cell line | Treatment", breaks =c(paste(CellLineModification,gsub("_"," ",CompoundName),sep=" | "),paste(wildtype,gsub("_"," ",CompoundName),sep=" | "),paste(CellLineModification,gsub("_"," ",ControlName),sep=" | "),paste(wildtype,gsub("_"," ",ControlName),sep=" | ")), values=ColorKey)+
      labs(title = paste(IncucyteDataAndMetaDataAndConditionsSubsetPlate$cell_line, CellLineModification, IncucyteDataAndMetaDataAndConditionsSubsetPlate$cell_number, '-', CompoundNameConcentration, '\n'), col='Treatment_Concentration_Bio_rep') +
      labs(x="\nTimepoint (hr)", y="Confluency (%)\n")+
      theme(axis.text=element_text(size=80),
            plot.margin = margin(1,1,1,1, "cm"),
            legend.text = element_text(size = 58),
            axis.title=element_text(size=90),
            legend.title=element_text(size=60),
            plot.title = element_text(size=85),
            axis.line = element_line(colour = 'black', linewidth = 4),
            axis.ticks = element_line(colour = "black", linewidth = 4),
            legend.box.background = element_rect(colour = "black", linewidth = 4))

  }

  PlotObject = PlotObject + guides(colour = guide_legend(override.aes = list(size = 16),ncol = 1))

  return(PlotObject)
}

#########################################################################################################################################################################################################################################################################

#' Plot Dose-response Curves for each Phenotype

#' This function plots the dose-response curves of all modeling-independent phenotypes normalized to control treatment (i.e. DMSO). Subunits from the same BAF subtype are plotted on the same graph.
#' @param PhenotypesNormToDMSODF Dataframe containing phenotypes normalized to control treatment, condition IDs and analysis group IDs.
#' @param BAFSubtypeDF A dataframe containing cell line subtype information (defined by the user).
#' @param all.group.colors A user-defined list specifying which colors should represent each of the cell lines in the dose-response curves.
#' @param all.group.lines A user-defined list specifying which line types should represent each of the cell lines in the dose-response curves.
#' @param wildtype A variable representing the wildtype cell line.
#' @param compound A variable representing the name of the compound.
#' @param treatmentlabel A variable representing how compound should be labeled in the plot.
#' @param controllabel A variable representing how the control compound should be labeled in the plot.
#' @param phenotype A variable specifying which phenotype is being plotted.
#' @param subtype A variable specifying which BAF subtype is being plotted.
#' @param scaling Specifies whether y-axis range should be the same each heatmap produced for each phenotype ("scaled") or not ("unscaled"; the default, where simply the minimum & maximum values are used).
#'
#' @return
#'\item{PlotObject}{ Saved ggplot object of the dose-response curve.}
#'
#' @author Caroline Barry
#' @import dplyr
#' @import ggplot2
#' @export

plotDoseResponseCurves <- function(PhenotypesNormToDMSODF, BAFSubtypeDF,all.group.colors,all.group.lines,wildtype,compound,treatmentlabel,controllabel, phenotype,subtype,scaling = "unscaled"){

    PhenotypesNormToDMSODFSubset = subset(PhenotypesNormToDMSODF, PhenotypesNormToDMSODF$treatment==compound)

    PhenotypeMean <- paste(phenotype,"_ratio_mean")

    PhenotypeMean <- gsub(" ","",PhenotypeMean)

    PhenotypeName1 <- phenotype

    PhenotypeName2 <- gsub("_", " ", PhenotypeName1)

    PhenotypeStdDev <- paste(phenotype,"_ratio_sd")

    PhenotypeStdDev <- gsub(" ","",PhenotypeStdDev)

    BAFSubtypeDFSubset <- subset(BAFSubtypeDFNoWT2, Subtype %in% c(subtype,"control"))

    SubtypeLabel <- gsub(" ", "_",subtype)

    group.lines <-all.group.lines[names(all.group.lines)%in%c(BAFSubtypeDFSubset$CL_name)|names(all.group.lines)==gsub("_"," ",wildtype)]

    BAFSubtypeDFSubset <- BAFSubtypeDFSubset %>%mutate(KO_partner_label=paste(CL,"|",wildtype)) %>% mutate(KO_partner_label=gsub(" ", "", KO_partner_label))

    group.colors <-all.group.colors[names(all.group.colors)%in%c(gsub("_", " ", BAFSubtypeDFSubset$CL_name),wildtype)]

    PhenotypesNormToDMSODFSubtypeSubset <-subset(PhenotypesNormToDMSODFSubset, KO_partner %in% BAFSubtypeDFSubset$KO_partner_label)

    PhenotypesNormToDMSODFSubtypeSubsetWT <- PhenotypesNormToDMSODFSubtypeSubset %>%
      group_by(concentration_0) %>%
      dplyr::filter(cell_line_modification ==wildtype)

    PhenotypesNormToDMSODFSubtypeSubset = subset(PhenotypesNormToDMSODFSubtypeSubset, treatment %in% PhenotypesNormToDMSODFSubtypeSubsetWT$treatment & concentration_0 %in% PhenotypesNormToDMSODFSubtypeSubsetWT$concentration_0)

    # get concentrations in numerical order
    CompoundConcentrationString <- PhenotypesNormToDMSODFSubtypeSubset$concentration_0

    CompoundConcentrationNumber <- lapply(CompoundConcentrationString, function(x) gsub("[^0-9.-]", "", x))

    CompoundConcentrationUnits <-lapply(CompoundConcentrationString, function(x) gsub('[0-9.]+', '', x))

    PhenotypesNormToDMSODFSubtypeSubset <- PhenotypesNormToDMSODFSubtypeSubset%>%
      mutate(conc_values = CompoundConcentrationNumber)

    PhenotypesNormToDMSODFSubtypeSubset$conc_values <- log2(as.numeric(PhenotypesNormToDMSODFSubtypeSubset$conc_values))

    PhenotypesNormToDMSODFSubtypeSubset <- PhenotypesNormToDMSODFSubtypeSubset%>%
      mutate(cell_line_modification = gsub("_", " ", cell_line_modification))%>%
      mutate(KO_partner = gsub("_", " ", KO_partner))

    print(gsub(" ", "_",subtype))

    print(PhenotypeName1)

    print(unique(PhenotypesNormToDMSODFSubtypeSubset$treatment_0))

    PlotObject<-PhenotypesNormToDMSODFSubtypeSubset %>%
      ggplot(aes(x=conc_values,  y=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean], group=interaction(KO_partner,cell_line_modification), col = cell_line_modification ))+
      geom_line(aes(x=conc_values,y=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean],group = interaction(cell_line_modification,KO_partner),col = cell_line_modification),size = 7)  +

      {if(isTRUE(scaling == "scaled")==TRUE){
        ylim(min(PhenotypesNormToDMSODF[,PhenotypeMean]-PhenotypesNormToDMSODF[,PhenotypeStdDev]),max(PhenotypesNormToDMSODF[,PhenotypeMean]+PhenotypesNormToDMSODF[,PhenotypeStdDev]))
      }
        if(isTRUE(scaling == "unscaled")==TRUE){
          ylim(min(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]-PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev]),max(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]+PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev]))}
      }+
      geom_hline(yintercept=1)+
      geom_errorbar(aes(x=conc_values,ymin=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]-PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev], ymax=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]+PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev], group = interaction(cell_line_modification,KO_partner), color =cell_line_modification), width = 0.35,  lwd=4) +#width = 0.05,lwd=0.7

      {if(isTRUE(scaling == "scaled")==TRUE){
        coord_cartesian(ylim = c(min(PhenotypesNormToDMSODF[,PhenotypeMean]-PhenotypesNormToDMSODF[,PhenotypeStdDev]),max(PhenotypesNormToDMSODF[,PhenotypeMean]+PhenotypesNormToDMSODF[,PhenotypeStdDev])))
      }

        if(isTRUE(scaling == "unscaled")==TRUE){coord_cartesian(ylim = c(min(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]-PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev]),max(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]+PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev])))}
      }+
      scale_color_manual("Cell line", values=group.colors)+#
      scale_fill_manual("Cell line", values=group.colors)+
      geom_point(size = 13.0) +
      theme_classic() +
      labs(col='Cell line modification') +
      labs(x=paste("\nlog2", "[",CompoundConcentrationUnits[1],"]"),
           y=paste(gsub("_"," + ",treatmentlabel),"/",controllabel, "\n"))+
      theme(axis.text=element_text(size=80),
            legend.text = element_text(size = 54),
            legend.title = element_text(size = 54),
            axis.text.x = element_text(size = 80),
            axis.text.y = element_text(size = 80),
            axis.title=element_text(size=90),
            plot.margin = margin(1.25,0.90,0.90,0.90, "cm"),
            legend.position=c(.219,.236),
            legend.key.size=unit(3,"lines"),
            legend.background=element_blank(),
            axis.line.x.bottom=element_line(linewidth=4),
            axis.line.y.left=element_line(linewidth=4)
      )

  return(PlotObject)
}

#########################################################################################################################################################################################################################################################################

#' Plot Heatmap for each Phenotype

#' This function plots an unclustered heatmap of the AUC value from the KO normalized to the WT across all biological replicates for a given modeling-independent phenotype treated with the same compound
#' (i.e. the confluency after 96 hrs of Adavosertib). The asterisks denote significance level of the treatment for KO compared to WT from the same biological replicate via their AUC values.
#' @param CompleteAUCDFSubset A dataframe containing treatment, cell line modification, mean KO/WT AUC ratio and standard deviation for a given phenotype. The KO is normalized to the WT from the same biological replicate. These ratios are then averaged across the biological replicates. Includes significances calculated using Bayes t-test and corrected using FDR.
#' @param CompoundMOADF A dataframe containing compound mechanism of action information (defined by the user).
#' @param MOAColors A list of colors for each compound mechanism of action (defined by the user) that will be used to plot the heatmap.
#' @param BAFSubtypeDF A dataframe containing cell line subtype information (defined by the user).
#' @param SubtypeColors A list of colors for each cell line subtype (defined by the user) that will be used to plot the heatmap.
#' @param ListOfCompoundOrder A list providing the desired order the compounds should be in on the heatmap (defined by the user).
#' @param ListOfCLOrder A list providing the desired order the cell lines should be in on the heatmap (defined by the user).
#' @param kolabel Variable specifying how knockout cell line should be labeled.
#' @param wtlabel Variable specifying how wildtype cell line should be labeled.
#' @param scaling Specifies whether y-axis range should be the same each heatmap produced for each phenotype ("scaled") or not ("unscaled"; the default).
#' @param minColorRange For the case when scaling is set to "scaled", variable represents user-defined minimum value for y-axis. Default is NULL.
#' @param maxColorRange For the case when scaling is set to "scaled", variable represents user-defined minimum value for y-axis. Default is NULL.
#'
#' @return
#'\item{comprehensive_heatmap}{ Saved ComplexHeatmap object.}
#'
#' @author Caroline Barry
#' @import ComplexHeatmap
#' @import dplyr
#' @import ggplot2
#' @export


plotAUCHeatmapUnclustered <- function(CompleteAUCDFSubset,CompoundMOADF,MOAColors,BAFSubtypeDF,SubtypeColors,ListOfCompoundOrder,ListOfCLOrder, kolabel,wtlabel, scaling = "unscaled",minColorRange=NULL,maxColorRange=NULL){

  # prepare matrix of mean log2 fold-change of AUC values
  HeatmapDF <- data.frame(stringsAsFactors = FALSE,matrix(nrow = length(unique(CompleteAUCDFSubset$treatment)), ncol = length(unique(CompleteAUCDFSubset$cell_line_modification))))

  CompleteAUCDFSubset <- CompleteAUCDFSubset%>%
    mutate(cell_line_modification = gsub("_", " ", cell_line_modification))

  colnames(HeatmapDF) <- unique(CompleteAUCDFSubset$cell_line_modification)

  rownames(HeatmapDF) <- unique(CompleteAUCDFSubset$treatment)

  CL_list <- unique(CompleteAUCDFSubset$cell_line_modification)

  treat_list <- unique(CompleteAUCDFSubset$treatment)

  for(k in 1:length(CL_list)){

    CompleteAUCDFCLSubset <- subset(CompleteAUCDFSubset, cell_line_modification==CL_list[k])

    for(j in 1:length(treat_list)){

      CompleteAUCDFCLSubsetTreat <- subset(CompleteAUCDFCLSubset, treatment ==treat_list[j])

      if (nrow(CompleteAUCDFCLSubsetTreat) == 0){
        next} else

          HeatmapDF[j,k] <- CompleteAUCDFCLSubsetTreat$mean_log_AUC_value_ratio
    }
  }

  HeatmapDF <- as.matrix(HeatmapDF)

  # prepare matrix of adjusted p-values
  PValueDF <- data.frame(stringsAsFactors = FALSE,matrix(nrow = length(unique(CompleteAUCDFSubset$treatment)), ncol = length(unique(CompleteAUCDFSubset$cell_line_modification))))

  colnames(PValueDF) <- unique(CompleteAUCDFSubset$cell_line_modification)

  rownames(PValueDF) <- unique(CompleteAUCDFSubset$treatment)

  CL_list <- unique(CompleteAUCDFSubset$cell_line_modification)

  treat_list <- unique(CompleteAUCDFSubset$treatment)

  for(k in 1:length(CL_list)){

    CompleteAUCDFCLSubset <- subset(CompleteAUCDFSubset, cell_line_modification==CL_list[k])

    for(j in 1:length(treat_list)){

      CompleteAUCDFCLSubsetTreat <- subset(CompleteAUCDFCLSubset, treatment ==treat_list[j])

      if (nrow(CompleteAUCDFCLSubsetTreat) == 0){
        next} else

          PValueDF[j,k] <- CompleteAUCDFCLSubsetTreat$p_value
    }
  }

  PValueDF <- as.matrix(PValueDF)

  CompoundOrder <- rownames(HeatmapDF)

  print(paste("rownames of heatmap: ", rownames(HeatmapDF)))

  print(paste("colnames of heatmap: ", colnames(HeatmapDF)))

  CompoundMOADF <- CompoundMOADF %>% arrange(match(Compound,CompoundOrder))

  ha = HeatmapAnnotation(MOA = CompoundMOADF$MOA, col =MOAColors,
                         which="row",
                         gap = unit(1, 'mm'),
                         annotation_legend_param = list(
                           title_gp = gpar(
                             fontface = "bold")
                         ),
                         annotation_name_gp= gpar(
                           fontface = "bold")
  )

  CLOrder <- colnames(HeatmapDF)

  BAFSubtypeDFNoWT <- BAFSubtypeDF %>% arrange(match(CL,CLOrder))

  BAFAnn = HeatmapAnnotation(Subtype = BAFSubtypeDF$Subtype, col =SubtypeColors,
                             which="col",
                             gap = unit(1, 'mm'),
                             annotation_legend_param = list(
                               title_gp = gpar(
                                 fontface = "bold")),
                             annotation_name_gp= gpar(
                               fontface = "bold")
  )


  if(isTRUE(scaling == "unscaled")==TRUE){
    minval <- min(CompleteAUCDFSubset$mean_log_AUC_value_ratio)
    maxval <- max(CompleteAUCDFSubset$mean_log_AUC_value_ratio)
  } else if(isTRUE(scaling == "scaled")==TRUE){
    minval <- minColorRange
    maxval <- maxColorRange
  }

  comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
    if(isTRUE(PValueDF[i,j] < 0.001)==TRUE) {
      grid.text("****",x,y)
    } else if(isTRUE(PValueDF[i,j] < 0.01)==TRUE) {
      grid.text("***",x,y)
    }
    else if(isTRUE(PValueDF[i,j]  < 0.05)==TRUE) {
      grid.text("**",x,y)
    }
    else if(isTRUE(PValueDF[i,j]  < 0.1)==TRUE) {
      grid.text("*",x,y)
    }
    if(isTRUE(is.na(PValueDF[i,j]))==TRUE) {
      grid.text("/",x,y)}

  },right_annotation=ha,

  name = paste("Log2",kolabel,"/",wtlabel,"AUC"),use_raster = TRUE,column_title = paste(PhenotypeName, "\n"),
  column_order = ListOfCLOrder,
  row_order= ListOfCompoundOrder,
  col = colorRamp2(c(minval, 0, maxval),c("blue", "white", "red")),
  top_annotation=BAFAnn)

  return(comprehensive_heatmap)
}

#########################################################################################################################################################################################################################################################################
#' Plot Clustered Summary Heatmap of AUC Normalized to WT for 5 Phenotypes

#'
#' Using the best clustering method and K number of clusters, this function plots a clustered heatmap of the mean log2 fold-change AUC value KO vs WT across the 5 phenotypes (Confluency_after_24h, Confluency_after_48h, Confluency_after_72h, Confluency_after_96h, Second_Phase_Relative_Change_in_Confluency) treated with the same compoundd (i.e. Adavosertib). The asterisks denote
#' significance level of the treatment for KO compared to WT from the same biological replicate via their AUC values.
#' @param AUCDF5PhenotypesDF A dataframe containing treatment, cell line modification, mean KO/WT AUC ratio and standard deviation for each given phenotype. The KO is normalized to the WT from the same biological replicate. These ratios are
#' then averaged across the biological replicates.
#' @param PValue5PhenotypesDF A dataframe containing treatment, cell line modification, Welsh t-statistic and Holmes corrected p-value for each given phenotype.
#' @param bestClustMethod A dataframe containing the best clustering algorithm (i.e. K-means) and K clusters based on the CH score.
#' @param BAFSubtypeDF A dataframe containing cell line subtype information (defined by the user).
#' @param CompoundMOADF A dataframe containing compound mechanism of action information (defined by the user).
#' @param SubtypeColors A list of colors for each cell line subtype (defined by the user) that will be used for the heatmap.
#' @param MOAColors A list of colors for each compound mechanism of action (defined by the user) that will be used for the heatmap.
#' @param kolabel Variable denoting how KO cell line should be written on the heatmap.
#' @param wtlabel Variable denoting how WT cell line should be written on the heatmap.
#' @param columnk Numerical variable representing the user-specified k number of clusters for the columns of the heatmap. If columnk = 0, then phenotypes will be in order of Confluency_after_24h, Confluency_after_48h, Confluency_after_72h, Confluency_after_96h, Second_Phase_Relative_Change_in_Confluency.
#' @param rowtextsize Numerical value representing text size of row labels.
#' @param scaling Variable denoting whether minimum and maximum values for color range should be user-defined values ("scaled") or simply the minimum and maximum values of the input data ("unscaled"). Unscaled is the default.
#' @param minColorRange Minimum value for color range. Default is NULL.
#' @param maxColorRange Maximum value for color range. Default is NULL.
#'
#' @return
#'\item{Returns list containing saved ComplexHeatmap object, PCA plot, & visualization of silhouette width for each cluster.}
#'
#' @author Caroline Barry
#' @import ComplexHeatmap
#' @import dplyr
#' @import ggplot2
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export


getAUCSummaryHeatmapClustered <- function(InputDF,AUCDF5PhenotypesDF,PValue5PhenotypesDF,bestClustMethod,BAFSubtypeDF,CompoundMOADF,SubtypeColors,MOAColors,kolabel,wtlabel,columnk,rowtextsize,scaling="unscaled",minColorRange=NULL,maxColorRange=NULL){

  HeatmapDF <- as.matrix(AUCDF5PhenotypesDF[1:5])

  PValueDF <- as.matrix(PValue5PhenotypesDF[1:5])

  CompoundOrder <- rownames(HeatmapDF)

  BAFSubtypeDF <- BAFSubtypeDF %>% mutate(CL=gsub(" ","_",CL))

  ColumnNamesDF <- data.frame(CompoundOrder,str_extract(CompoundOrder, "[^|]+"),sub('.*[|]', '', CompoundOrder))

  colnames(ColumnNamesDF) <- c("CL_Treat", "CL", "Compound")

  BAFSubtypeDF <- BAFSubtypeDF %>% arrange(match(CL,ColumnNamesDF$CL))

  CompoundMOADF <- CompoundMOADF %>% arrange(match(Compound,ColumnNamesDF$Compound))

  anno_df1 = merge(ColumnNamesDF, BAFSubtypeDF, by="CL", all = T)
  anno_df1 <- na.omit(anno_df1)
  anno_df1 <- anno_df1 %>% arrange(match(CL_Treat,ColumnNamesDF$CL_Treat))

  anno_df2 = merge(ColumnNamesDF, CompoundMOADF, by="Compound", all = T)
  anno_df2 <- na.omit(anno_df2)
  anno_df2 <- anno_df2 %>% arrange(match(CL_Treat,ColumnNamesDF$CL_Treat))

  ha1 = HeatmapAnnotation(Subtype=anno_df1$Subtype,
                          col =SubtypeColors,
                          which="row",
                          gap = unit(1, 'mm'),
                          annotation_name_gp= gpar(fontsize = 9)
  )

  ha2 = HeatmapAnnotation(MOA=anno_df2$MOA,
                          col =MOAColors,
                          which="row",
                          gap = unit(1, 'mm'),
                          annotation_name_gp= gpar(fontsize = 9)
  )

  ###### if best clustering method is K-means ########
  if(isTRUE(bestClustMethod$Method=="K-means")==TRUE){

    set.seed(123)
    k2m_data <- kmeans(InputDF, as.numeric(as.character(bestClustMethod$Cluster_Num)), nstart=25)

    if(isTRUE(columnk==0)==TRUE){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "K-means clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        row_split =cbind(k2m_data$cluster),
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "K-means clustering of conditions across phenotypes",
        row_split = k2m_data$cluster,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )

      }
    }

    if(isTRUE(columnk!=0)==TRUE){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "K-means clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        row_split = k2m_data$cluster,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "K-means clustering of conditions across phenotypes",
        row_split = k2m_data$cluster,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      }
    }


    ClusterPlotObject<-factoextra::fviz_cluster(k2m_data, data = InputDF,
                                                ellipse.type = "convex",
                                                star.plot = TRUE,
                                                # repel = TRUE,
                                                geom="point",
                                                ggtheme = theme_minimal(),
                                                main = "K-means Clustering"
    )  + theme(axis.text.y = element_text(size = 13),axis.text.x = element_text(size = 13),
               axis.title=element_text(size=14,face="bold"),
               legend.text = element_text(size=13),legend.title = element_text(size=14),
               plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
               axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
               axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))



    k2m_data <- factoextra::eclust(InputDF, "kmeans", k = as.numeric(as.character(bestClustMethod$Cluster_Num)), nstart = 25, graph = F)
    AvgSilPlotObject <-factoextra::fviz_silhouette(k2m_data, palette = "jco",
                                                   ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                                     axis.title=element_text(size=14,face="bold"),
                                                                                     legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                                     plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                                     axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                                     axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

  }

  ###### if best clustering method is PAM ########
  if(isTRUE(bestClustMethod$Method=="PAM")==TRUE){
    set.seed(123)

    pam_data <- pam(InputDF,as.numeric(as.character(bestClustMethod$Cluster_Num)))

    if(columnk==0){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "PAM clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        row_split = pam_data$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "PAM clustering of conditions across phenotypes",
        row_split = pam_data$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )

      }
    }

    if(columnk!=0){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "PAM clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        row_split = pam_data$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "PAM clustering of conditions across phenotypes",
        row_split = pam_data$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      }
    }
    ClusterPlotObject<-fviz_cluster(pam_data, data = InputDF,
                                    ellipse.type = "convex",
                                    star.plot = TRUE,
                                    repel = TRUE,
                                    geom="point",
                                    ggtheme = theme_minimal(),
                                    main ="PAM Clustering"
    )  + theme(axis.text.y = element_text(size = 13),axis.text.x = element_text(size = 13),
               axis.title=element_text(size=14,face="bold"),
               legend.text = element_text(size=13),legend.title = element_text(size=14),
               plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
               axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
               axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

    pam_data <- eclust(InputDF, "pam", k = as.numeric(as.character(bestClustMethod$Cluster_Num)), nstart = 25, graph = F)
    AvgSilPlotObject <-fviz_silhouette(pam_data, palette = "jco",
                                       ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                         axis.title=element_text(size=14,face="bold"),
                                                                         legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                         plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                         axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                         axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

  }


  ###### if best clustering method is C-means ########
  if(isTRUE(bestClustMethod$Method=="C-means")==TRUE){
    set.seed(123)
    fannyx <- fanny(InputDF, as.numeric(as.character(bestClustMethod$Cluster_Num)),metric = "euclidean", stand = FALSE,memb.exp = 1.5)

    if(columnk==0){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "C-means clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        row_split = fannyx$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "C-means clustering of conditions across phenotypes",
        row_split = fannyx$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )

      }
    }

    if(columnk!=0){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "C-means clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        row_split = fannyx$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "C-means clustering of conditions across phenotypes",
        row_split = fannyx$clustering,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      }
    }

    ClusterPlotObject<-factoextra::fviz_cluster(fannyx, data = InputDF,
                                                ellipse.type = "convex",
                                                palette = "jco",
                                                geom="point",
                                                ggtheme = theme_minimal(),
                                                repel = TRUE,
                                                main = "C-means Clustering")+ theme(axis.text.y = element_text(size = 13),axis.text.x = element_text(size = 13),
                                                                                    axis.title=element_text(size=14,face="bold"),
                                                                                    legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                                    plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                                    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                                    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


    AvgSilPlotObject<-fviz_silhouette(fannyx, palette = "jco",
                                      ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                        axis.title=element_text(size=14,face="bold"),
                                                                        legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                        plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))



  }

  ###### if best clustering method is hierarchical ########
  if(isTRUE(bestClustMethod$Method=="Hierarchical")==TRUE){

    set.seed(123)

    dist_man <- dist(InputDF2, method="manhattan")

    hc_m2 <- hclust(d=dist_man, method="average")

    groupward6 <- cutree(hc_m2, k = as.numeric(as.character(bestClustMethod$Cluster_Num)))

    if(columnk==0){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "Hierarchical clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        row_split = groupward6,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "Hierarchical clustering of conditions across phenotypes",
        row_split = groupward6,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_order= c("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency"),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )

      }
    }

    if(columnk!=0){

      if(isTRUE(scaling!="unscaled")==TRUE){

        minval <- minColorRange
        maxval <- maxColorRange

        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "Hierarchical clustering of conditions across phenotypes",
        col = colorRamp2(c(minval, 0, maxval), c("blue", "white", "red")),
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        row_split = groupward6,
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      } else{
        comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
          if(PValueDF[i,j] < 0.001) {
            grid.text("****",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          } else if(PValueDF[i,j] < 0.01) {
            grid.text("***",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.05) {
            grid.text("**",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
          else if(PValueDF[i,j]  < 0.1) {
            grid.text("*",x,y,0.008, 0.008,gpar(fontsize = 0.75))
          }
        },right_annotation=c(ha1,ha2),
        name = paste("Log2", kolabel, "/",wtlabel,"AUC"),use_raster = TRUE,column_title = "Hierarchical clustering of conditions across phenotypes",
        row_split = groupward6,
        column_names_gp = grid::gpar(fontsize = 6),
        row_names_gp = grid::gpar(fontsize = rowtextsize),
        column_names_rot = 25,
        heatmap_height = unit(200, "mm"),
        column_dend_height = unit(3, "mm")
        )
      }
    }

    ClusterPlotObject<-fviz_dend(hc_m2, k = as.numeric(as.character(bestClustMethod$Cluster_Num)),
                                 cex = 0.5,
                                 color_labels_by_k = TRUE,
                                 rect = TRUE ,
                                 main="Hierarchical Clustering")


    average <- eclust(InputDF, "hclust", k = as.numeric(as.character(bestClustMethod$Cluster_Num)), hc_metric = "manhattan",hc_method = "average", graph = F)
    AvgSilPlotObject<-fviz_silhouette(average, palette = "jco",
                                      ggtheme = theme_classic())+ theme(axis.text.y = element_text(size = 13),
                                                                        axis.title=element_text(size=14,face="bold"),
                                                                        legend.text = element_text(size=13),legend.title = element_text(size=14),
                                                                        plot.title=element_text(size=25),plot.subtitle=element_text(size=18),
                                                                        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                                        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))
  }


  return(list(Heatmap = comprehensive_heatmap, ClusterViz = ClusterPlotObject,AvgSil = AvgSilPlotObject ))
}

