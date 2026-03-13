#' Plot Growth Curves without Split Points for each Analysis Group

#' This function plots the growth curves without split points of one given analysis group (KO + treatment, KO + control, WT + treatment, and WT + control) across all biological replicates for one given condition (treatment & concentration) together on the same graph.
#' @param AvgIncucyteDataAndMetaDataAndConditionsLabeledDF Dataframe that combines Incucyte Confluency Dataset, Meta Data, Condition IDs and analysis IDs. The confluency is averaged across technical replicates.
#' @param CLTRTList List containing user-specified cell lines and compounds to plot. These are separated with "|" (i.e. cell line | compound).
#' @param ColorListTreat List of colors for treated conditions from CLTRTList. These will be used for the line coloring in the plots.
#' @param ColorListControl List of colors for control conditions from CLTRTList. These will be used for the line coloring in the plots.
#' @param wildtype Variable representing name of wildtype cell line.
#' @param treatmentendtime Numerical variable represents the time at which treatment ends.
#' @param TodaysDate The date.
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

  # CompoundConcentration <-unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$concentration_0)

  CellLineModification <- unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$cell_line_modifications)

  # ListOfCLCmpndConc <- paste(CLTRTList, CompoundNameConcentration, sep=" | ")#CompoundName

  ListOfCLCmpnd <- paste(CLTRTList, CompoundName, sep=" | ")

  names(ColorListTreat) <- ListOfCLCmpnd

  # ListOfCLCntrlConc <- paste(CLTRTList, ControlNameConcentration, sep=" | ")

  ListOfCLCntrl <- paste(CLTRTList, ControlName, sep=" | ")

  names(ColorListControl) <- ListOfCLCntrl

  ColorKey <- c(ColorListTreat,ColorListControl)

  if(isTRUE(is.null(LinetypeListTreat) &is.null(LinetypeListControl))==FALSE){
    names(LinetypeListTreat) <- ListOfCLCmpnd
    names(LinetypeListControl) <- ListOfCLCntrl
    LineKey <- c(LinetypeListTreat,LinetypeListControl)
  }



  # ControlNameConcentration <-gsub("_"," ",unique(IncucyteDataAndMetaDataAndConditionsSubsetValue$treat_conc_id))

  ListOfPlateIDs = unique(IncucyteDataAndMetaDataAndConditionsSubset$plate_ID)

  # create a unique identifier for observation
  IncucyteDataAndMetaDataAndConditionsSubset <- IncucyteDataAndMetaDataAndConditionsSubset %>%
    # mutate(concentration_0=case_when(is.na(treatment_0) ~ "none",
    #                                  TRUE ~ concentration_0)) %>%
    # mutate(treatment_0 = case_when(is.na(treatment_0) ~ "none",
    #                                TRUE ~ treatment_0)) %>%
    mutate(uniq_id=paste(CL_treat_conc_id,plate_ID,sep="|"))

  # create mapping for unique id to color, depends on RSOperator
  MAPPING <- IncucyteDataAndMetaDataAndConditionsSubset %>% distinct(CL_treat_conc_id,plate_ID,uniq_id)
  RS_COLS = brewer.pal(8, "Set2")
  RS_COLS = RS_COLS[1:n_distinct(MAPPING$CL_treat_conc_id)]
  names(RS_COLS) = unique(MAPPING$CL_treat_conc_id)
  PLOT_COLS = RS_COLS[MAPPING$CL_treat_conc_id]
  names(PLOT_COLS) = MAPPING$uniq_id

  # # get points after split point
  # IncucyteDataAndMetaDataAndConditionsSubsetPlateOriginal <- IncucyteDataAndMetaDataAndConditionsSubset %>%
  #   dplyr::group_by(bio_rep_id) %>%
  #   dplyr::filter(timepoint>split_timepoint & timepoint<=treatmentendtime)

  IncucyteDataAndMetaDataAndConditionsSubsetPlate <- IncucyteDataAndMetaDataAndConditionsSubset %>%
    dplyr::group_by(bio_rep_id) %>%
    dplyr::filter(timepoint<=treatmentendtime)
  # dplyr::filter(case_when(second_phase==1 & timepoint timepoint<=split_timepoint, second_phase ==0 ~ <=treatmentendtime))


  PlotObject <- ggplot(IncucyteDataAndMetaDataAndConditionsSubsetPlate,mapping=aes(x=timepoint, y=mean,group=interaction(plate_ID,CL_treat_conc_id),col=CL_TRT, linetype = CL_TRT)) + #uniq_id #CL_treat_conc_id
    # geom_line(linewidth=1.5) +
    theme_classic() +
    # geom_point(data = IncucyteDataAndMetaDataAndConditionsSubsetPlateOriginal,shape="x", size=9, aes(x = timepoint, y = mean, color = uniq_id),size = 6) +
    geom_line(aes(x=timepoint,y=mean,  group=interaction(plate_ID,CL_TRT),col = CL_TRT),size = 7)  + #CL_treat_conc_id


    # geom_point(size = 4.0) +#4.0
    scale_linetype_manual("Cell line | Treatment",breaks = c(paste(CellLineModification,gsub("_"," ",CompoundName),sep=" | "),paste(wildtype,gsub("_"," ",CompoundName),sep=" | "),
                                                             paste(CellLineModification,ControlName,sep=" | "),paste(wildtype,ControlName,sep=" | ")),values = LineKey)+


    ylim(0.0,100)+
    # scale_color_manual(values = PLOT_COLS)+
    # scale_color_manual("Cell line | Treatment", breaks = c(paste(CellLineModification,gsub("_"," ",CompoundNameConcentration),sep=" | "),paste(wildtype,gsub("_"," ",CompoundNameConcentration),sep=" | "),paste(CellLineModification,gsub("_"," ",ControlNameConcentration),sep=" | "),paste(wildtype,gsub("_"," ",ControlNameConcentration),sep=" | ")), values=ColorKey)+
    # scale_fill_manual("Cell line | Treatment", breaks =c(paste(CellLineModification,gsub("_"," ",CompoundNameConcentration),sep=" | "),paste(wildtype,gsub("_"," ",CompoundNameConcentration),sep=" | "),paste(CellLineModification,gsub("_"," ",ControlNameConcentration),sep=" | "),paste(wildtype,gsub("_"," ",ControlNameConcentration),sep=" | ")), values=ColorKey)+
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

  PlotObject = PlotObject + guides(colour = guide_legend(override.aes = list(size = 16),ncol = 1))

  return(PlotObject)
}

#########################################################################################################################################################################################################################################################################

#' Plot Each Modeling-Independent Phenotype - Treatment Normalized to Control Treatment

#'
#' This function plots the dose-response curves of all modeling-independent phenotypes normalized to control treatment (i.e. DMSO). Subunits from the same BAF subtype are plotted on the same graph.
#' @param PhenotypesNormToDMSODF Dataframe containing phenotypes normalized to control treatment, condition IDs and analysis group IDs.
#' @param BAFSubtypeDF A dataframe containing cell line subtype information (defined by the user).
#' @param all.group.colors A user-defined list specifying which colors should represent each of the cell lines in the dose-response curves.
#' @param all.group.lines A user-defined list specifying which lines should represent each of the cell lines in the dose-response curves.
#' @param wildtype A variable representing the wildtype cell line.
#' @param control A variable representing the cell line used for quality control.
#' @param TodaysDate The date.
#' @param AnalysisName User-defined unique name of the analysis being conducted.
#' @param scaling Specifies whether y-axis range should be the same each heatmap produced for each phenotype ("scaled") or not ("unscaled"; the default).
#' @param Phenotype Specifies which phenotype should be plotted (i.e. "Confluency_after_24h").
#' @author Caroline Barry
#' @import dplyr
#' @import ggplot2
#' @export


plotDoseResponseCurves <- function(PhenotypesNormToDMSODF, BAFSubtypeDF,all.group.colors,all.group.lines,wildtype, compound,treatmentlabel,controllabel, phenotype,subtype,scaling = "unscaled"){#control , TodaysDate, AnalysisName,scaling = "unscaled", phenotypes = "all"

  if(isTRUE(wildtype!="self")==TRUE){
    PhenotypesNormToDMSODFSubset = subset(PhenotypesNormToDMSODF, PhenotypesNormToDMSODF$treatment==compound)

    # PhenotypeMean <- ListOfPhenotypeMeans[n]
    #
    # PhenotypeName1 <- gsub("_ratio_mean", "", PhenotypeMean)
    #
    # PhenotypeName2 <- gsub("_", " ", PhenotypeName1)
    #
    # PhenotypeStdDev <- ListOfPhenotypeStdDevs[n]

    PhenotypeMean <- paste(phenotype,"_ratio_mean")

    PhenotypeMean <- gsub(" ","",PhenotypeMean)

    PhenotypeName1 <- phenotype

    PhenotypeName2 <- gsub("_", " ", PhenotypeName1)

    PhenotypeStdDev <- paste(phenotype,"_ratio_sd")

    PhenotypeStdDev <- gsub(" ","",PhenotypeStdDev)


    BAFSubtypeDFSubset <- subset(BAFSubtypeDFNoWT2, Subtype %in% c(subtype,"control"))

    # print("BAFSubtypeDFSubset: ")
    # print(BAFSubtypeDFSubset)

    SubtypeLabel <- gsub(" ", "_",subtype)#l

    group.lines <-all.group.lines[names(all.group.lines)%in%c(BAFSubtypeDFSubset$CL_name)|names(all.group.lines)==gsub("_"," ",wildtype)]

    BAFSubtypeDFSubset <- BAFSubtypeDFSubset %>%mutate(KO_partner_label=paste(CL,"|",wildtype)) %>% mutate(KO_partner_label=gsub(" ", "", KO_partner_label))

    # print("BAFSubtypeDFSubset: ")
    # print(BAFSubtypeDFSubset)

    # group.colors <-all.group.colors[names(all.group.colors)%in%c(gsub("_", " ", BAFSubtypeDFSubset$KO_partner_label))]

    group.colors <-all.group.colors[names(all.group.colors)%in%c(gsub("_", " ", BAFSubtypeDFSubset$CL_name),wildtype)]

    # print("PhenotypesNormToDMSODFSubset: ")
    # print(head(PhenotypesNormToDMSODFSubset))

    PhenotypesNormToDMSODFSubtypeSubset <-subset(PhenotypesNormToDMSODFSubset, KO_partner %in% BAFSubtypeDFSubset$KO_partner_label)

    # print("PhenotypesNormToDMSODFSubtypeSubset: ")
    # print(head(PhenotypesNormToDMSODFSubtypeSubset))

    PhenotypesNormToDMSODFSubtypeSubsetWT <- PhenotypesNormToDMSODFSubtypeSubset %>%
      group_by(concentration_0) %>%
      dplyr::filter(cell_line_modification ==wildtype)

    # print("PhenotypesNormToDMSODFSubtypeSubsetWT: ")
    # print(head(PhenotypesNormToDMSODFSubtypeSubsetWT))


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
      ggplot(aes(x=conc_values,  y=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean], group=interaction(KO_partner,cell_line_modification), col = cell_line_modification ))+# geom_point() + #KO_partner
      # geom_line(aes(x=conc_values,y=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean], group=interaction(KO_partner, cell_line_modification),col = cell_line_modification, linetype ="solid"),size = 7)  +#size = 0.9 #KO_partner
      geom_line(aes(x=conc_values,y=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean],group = interaction(cell_line_modification,KO_partner),col = cell_line_modification),size = 7)  +#linetype ="solid"

      {if(isTRUE(scaling == "scaled")==TRUE){
        ylim(min(PhenotypesNormToDMSODF[,PhenotypeMean]-PhenotypesNormToDMSODF[,PhenotypeStdDev]),max(PhenotypesNormToDMSODF[,PhenotypeMean]+PhenotypesNormToDMSODF[,PhenotypeStdDev]))
      }
        if(isTRUE(scaling == "unscaled")==TRUE){
          ylim(min(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]-PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev]),max(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]+PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev]))}
      }+

      # xlim(min(PhenotypesNormToDMSODFSubtypeSubset$conc_values), max(PhenotypesNormToDMSODFSubtypeSubset$conc_values))+
      geom_hline(yintercept=1)+
      geom_errorbar(aes(x=conc_values,ymin=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]-PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev], ymax=PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]+PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev], group = interaction(cell_line_modification,KO_partner), color =cell_line_modification), width = 0.35,  lwd=4) +#width = 0.05,lwd=0.7
      # scale_x_log10()+

      {if(isTRUE(scaling == "scaled")==TRUE){
        coord_cartesian(ylim = c(min(PhenotypesNormToDMSODF[,PhenotypeMean]-PhenotypesNormToDMSODF[,PhenotypeStdDev]),max(PhenotypesNormToDMSODF[,PhenotypeMean]+PhenotypesNormToDMSODF[,PhenotypeStdDev])))
      }

        if(isTRUE(scaling == "unscaled")==TRUE){coord_cartesian(ylim = c(min(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]-PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev]),max(PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]+PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeStdDev])))}
      }+
      # scale_color_manual("KO partner", breaks = gsub("_"," ",BAFSubtypeDFSubset$KO_partner_label), values=group.colors)+#
      # scale_fill_manual("KO partner", breaks = gsub("_"," ",BAFSubtypeDFSubset$KO_partner_label), values=group.colors)+
      scale_color_manual("Cell line", values=group.colors)+#
      scale_fill_manual("Cell line", values=group.colors)+
      geom_point(size = 13.0) +#4.0
      # scale_linetype_manual("cell line",values = group.lines)+
      theme_classic() +
      labs(#title = paste(PhenotypeName2,":",subtype, "subunits", "\n"),
        col='Cell line modification') +
      # labs(x=paste("\nlog2(",gsub("_"," ",PhenotypesNormToDMSODFSubtypeSubset$treatment_0), "[",CompoundConcentrationUnits[1],"] )"), #y=paste(treatmentlabel,"/",controllabel, "\n"))+
      labs(x=paste("\nlog2", "[",CompoundConcentrationUnits[1],"]"),
           y=paste(gsub("_"," + ",treatmentlabel),"/",controllabel, "\n"))+
      theme(axis.text=element_text(size=80),#30
            legend.text = element_text(size = 54),#25
            legend.title = element_text(size = 54),
            axis.text.x = element_text(size = 80),#23
            axis.text.y = element_text(size = 80),#63,
            axis.title=element_text(size=90),#63
            plot.margin = margin(1.25,0.90,0.90,0.90, "cm"),
            legend.position=c(.219,.236),
            legend.key.size=unit(3,"lines"),
            legend.background=element_blank(),
            axis.line.x.bottom=element_line(linewidth=4),
            axis.line.y.left=element_line(linewidth=4)
            #plot.title = element_text(size=40)
      )
  }


  #special case where wildtype is the CL itself
  if(isTRUE(wildtype=="self")==TRUE){
    PhenotypesNormToDMSODFSubset = subset(PhenotypesNormToDMSODF, PhenotypesNormToDMSODF$treatment_0==compound)

    # PhenotypeMean <- ListOfPhenotypeMeans[n]
    #
    # PhenotypeName1 <- gsub("_ratio_mean", "", PhenotypeMean)
    #
    # PhenotypeName2 <- gsub("_", " ", PhenotypeName1)
    #
    # PhenotypeStdDev <- ListOfPhenotypeStdDevs[n]

    PhenotypeMean <- paste(phenotype,"_ratio_mean")



    PhenotypeStdDev <- paste(phenotype,"_ratio_sd")


    PhenotypeStdDev <- gsub(" ","",PhenotypeStdDev)
    PhenotypeMean <- gsub(" ","",PhenotypeMean)

    BAFSubtypeDFSubset <- subset(BAFSubtypeDFNoWT2, Subtype %in% c(subtype,"control"))

    # print("BAFSubtypeDFSubset: ")
    # print(BAFSubtypeDFSubset)

    DesiredColumnsList <- c("condition_ID"  ,"cell_line_modification" , "treatment_0" ,  "concentration_0" ,"treatment" ,  "CL_treat_conc_id", "treat_conc_id" , "NLMM_Analysis_ID","treat_type",PhenotypeMean,PhenotypeStdDev )

    PhenotypesNormToDMSODFSubset <- subset(PhenotypesNormToDMSODFSubset, select = names(PhenotypesNormToDMSODFSubset) %in% DesiredColumnsList)

    names(PhenotypesNormToDMSODFSubset)[names(PhenotypesNormToDMSODFSubset) == PhenotypeMean] <- "Phenotype_Mean"

    names(PhenotypesNormToDMSODFSubset)[names(PhenotypesNormToDMSODFSubset) == PhenotypeStdDev] <- "Phenotype_Std_Dev"


    PhenotypeName1 <- phenotype

    PhenotypeName2 <- gsub("_", " ", PhenotypeName1)

    SubtypeLabel <- gsub(" ", "_",subtype)#l

    # group.lines <-all.group.lines[names(all.group.lines)%in%c(PhenotypesNormToDMSODFSubset$treat_type_CL)|names(all.group.lines)==gsub("_"," ",wildtype)]

    # BAFSubtypeDFSubset <- BAFSubtypeDFSubset %>%mutate(KO_partner_label=paste(CL,"|",wildtype)) %>% mutate(KO_partner_label=gsub(" ", "", KO_partner_label))

    # print("BAFSubtypeDFSubset: ")
    # print(BAFSubtypeDFSubset)

    group.colors <-all.group.colors[names(all.group.colors)%in%c(gsub("_", " ", BAFSubtypeDFSubset$CL))]

    # print("PhenotypesNormToDMSODFSubset: ")
    # print(head(PhenotypesNormToDMSODFSubset))

    PhenotypesNormToDMSODFSubtypeSubset <-subset(PhenotypesNormToDMSODFSubset, cell_line_modification %in% BAFSubtypeDFSubset$CL)

    # print("PhenotypesNormToDMSODFSubtypeSubset: ")
    # print(head(PhenotypesNormToDMSODFSubtypeSubset))

    PhenotypesNormToDMSODFSubtypeSubsetWT <- PhenotypesNormToDMSODFSubtypeSubset %>%
      group_by(concentration_0) %>%
      dplyr::filter(treat_type =="Compound")

    # print("PhenotypesNormToDMSODFSubtypeSubsetWT: ")
    # print(head(PhenotypesNormToDMSODFSubtypeSubsetWT))


    PhenotypesNormToDMSODFSubtypeSubset = subset(PhenotypesNormToDMSODFSubtypeSubset, treatment_0 %in% PhenotypesNormToDMSODFSubtypeSubsetWT$treatment_0
                                                 & concentration_0 %in% PhenotypesNormToDMSODFSubtypeSubsetWT$concentration_0
    )


    # #only keep WT from same plates
    # controlcl <- subset(BAFSubtypeDFNoWT2, Subtype=="control")
    #
    # PhenotypesNormToDMSODFSubtypeSubsetKO <- PhenotypesNormToDMSODFSubtypeSubset %>%
    #   dplyr::filter(cell_line_modification !=as.character(controlcl$CL))
    #
    # PhenotypesNormToDMSODFSubtypeSubset = subset(PhenotypesNormToDMSODFSubtypeSubset,  NLMM_Analysis_ID %in% PhenotypesNormToDMSODFSubtypeSubsetKO$NLMM_Analysis_ID)
    #

    # get concentrations in numerical order
    CompoundConcentrationString <- PhenotypesNormToDMSODFSubtypeSubset$concentration_0

    CompoundConcentrationNumber <- lapply(CompoundConcentrationString, function(x) gsub("[^0-9.-]", "", x))

    CompoundConcentrationUnits <-lapply(CompoundConcentrationString, function(x) gsub('[0-9.]+', '', x))


    PhenotypesNormToDMSODFSubtypeSubset <- PhenotypesNormToDMSODFSubtypeSubset%>%
      mutate(conc_values = lapply(concentration_0, function(x) gsub("[^0-9.-]", "", x)))

    PhenotypesNormToDMSODFSubtypeSubset$conc_values <- log2(as.numeric(PhenotypesNormToDMSODFSubtypeSubset$conc_values))

    PhenotypesNormToDMSODFSubtypeSubset <- PhenotypesNormToDMSODFSubtypeSubset%>%
      mutate(cell_line_modification = gsub("_", " ", cell_line_modification))%>%
      mutate(treat_type = gsub("_", " ", treat_type))


    print(gsub(" ", "_",subtype))

    print(PhenotypeName1)

    print(unique(PhenotypesNormToDMSODFSubtypeSubset$treatment_0))

    PhenotypesNormToDMSODFSubtypeSubset$cell_line_modification <- as.factor(PhenotypesNormToDMSODFSubtypeSubset$cell_line_modification)

    PhenotypesNormToDMSODFSubtypeSubset$treat_type <- as.factor(PhenotypesNormToDMSODFSubtypeSubset$treat_type)

    #PhenotypesNormToDMSODFSubtypeSubset[,PhenotypeMean]
    PlotObject<-PhenotypesNormToDMSODFSubtypeSubset %>%
      ggplot(aes(x=conc_values,  y=Phenotype_Mean, group=cell_line_modification, col = cell_line_modification))+
      # geom_point() +
      geom_line(aes(x=conc_values,y=Phenotype_Mean, group=interaction(cell_line_modification, treat_type),col = cell_line_modification, linetype =treat_type),size = 0.9)  +#

      {if(isTRUE(scaling == "scaled")==TRUE){
        #PhenotypesNormToDMSODF[,PhenotypeMean]
        minval <- min(PhenotypesNormToDMSODF[,PhenotypeMean]-PhenotypesNormToDMSODF[,PhenotypeStdDev])
        maxval <- max(PhenotypesNormToDMSODF[,PhenotypeMean]+PhenotypesNormToDMSODF[,PhenotypeStdDev])

        ylim(as.numeric(minval),as.numeric(maxval))
      }
        if(isTRUE(scaling == "unscaled")==TRUE){
          minval <- min(PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Mean-PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Std_Dev)
          maxval <- max(PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Mean+PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Std_Dev)

          ylim(minval,maxval)}
      }+

      # xlim(min(PhenotypesNormToDMSODFSubtypeSubset$conc_values), max(PhenotypesNormToDMSODFSubtypeSubset$conc_values))+
      geom_hline(yintercept=1)+
      geom_errorbar(aes(x=conc_values,ymin=Phenotype_Mean-Phenotype_Std_Dev, ymax=Phenotype_Mean+Phenotype_Std_Dev, color =cell_line_modification), width = 0.05,  lwd=0.7) +
      # scale_x_log10()+

      {if(isTRUE(scaling == "scaled")==TRUE){
        coord_cartesian(ylim = c(min(PhenotypesNormToDMSODF[,PhenotypeMean]-PhenotypesNormToDMSODF[,PhenotypeStdDev]),max(PhenotypesNormToDMSODF[,PhenotypeMean]+PhenotypesNormToDMSODF[,PhenotypeStdDev])))
      }

        if(isTRUE(scaling == "unscaled")==TRUE){coord_cartesian(ylim = c(min(PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Mean-PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Std_Dev),max(PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Mean+PhenotypesNormToDMSODFSubtypeSubset$Phenotype_Std_Dev)))}
      }+
      scale_color_manual(breaks = BAFSubtypeDFSubset$CL_name, values=group.colors)+#breaks = gsub("_"," ",BAFSubtypeDFSubset$CL),
      # scale_fill_manual(breaks = BAFSubtypeDFSubset$CL_name,values=group.colors)+
      geom_point(size = 4.0) +
      scale_linetype_manual("Treatment type",values = all.group.lines)+
      theme_classic() +
      labs(title = paste(PhenotypeName2,":",subtype, "subunits", "\n"), col='Cell line modification') +
      labs(x=paste("\nlog2(",gsub("_"," ",PhenotypesNormToDMSODFSubtypeSubset$treatment_0), "[",CompoundConcentrationUnits[1],"] )"), y=paste(treatmentlabel,"/",controllabel, "\n"))+
      theme(axis.text=element_text(size=30),
            legend.text = element_text(size = 25),
            axis.text.x = element_text(size = 23),
            axis.title=element_text(size=25,face="bold"),
            plot.margin = margin(1.25,0.90,0.90,0.90, "cm"),
            plot.title = element_text(size=40))
  }



  return(PlotObject)
}

#########################################################################################################################################################################################################################################################################

## Individual Heatmaps

#' Plot Unclustered Heatmap of AUC

#'
#' This function plots an unclustered heatmap of the AUC value from the KO normalized to the WT across all biological replicates for a given modeling-independent phenotype treated with the same compound
#' (i.e. the confluency after 96 hrs of Adavosertib). The asterisks denote significance level of the treatment for KO compared to WT from the same biological replicate via their AUC values.
#' @param CompleteAUCDF A dataframe containing treatment, cell line modification, mean KO/WT AUC ratio and standard deviation for this given phenotype. The KO is normalized to the WT from the same biological replicate. These ratios are then averaged across the biological replicates. Includes significances calculated using Bayes t-test and corrected using FDR.
#' @param CompoundMOADF A dataframe containing compound mechanism of action information (defined by the user).
#' @param MOAColors A list of colors for each compound mechanism of action (defined by the user) that will be used to plot the heatmap.
#' @param BAFSubtypeDF A dataframe containing cell line subtype information (defined by the user).
#' @param SubtypeColors A list of colors for each cell line subtype (defined by the user) that will be used to plot the heatmap.
#' @param ListOfCompoundOrder A list providing the desired order the compounds should be in on the heatmap (defined by the user).
#' @param ListOfCLOrder A list providing the desired order the cell lines should be in on the heatmap (defined by the user).
#' @param wildtype Variable representing name of wildtype cell line.
#' @param TodaysDate The date.
#' @param AnalysisName User-defined unique name of the analysis being conducted.
#' @param scaling Specifies whether y-axis range should be the same each heatmap produced for each phenotype ("scaled") or not ("unscaled"; the default).
#' @param phenotypes Specifies whether a heatmap for each phenotype should be generated ("all") or if only select phenotypes should be plotted (user-specified; i.e. c("Confluency_after_24h","Confluency_after_48h")).
#' @author Caroline Barry
#' @import ComplexHeatmap
#' @import dplyr
#' @import ggplot2
#' @export
#'


plotAUCHeatmapUnclustered <- function(CompleteAUCDFSubset,CompoundMOADF,MOAColors,BAFSubtypeDF,SubtypeColors,ListOfCompoundOrder,ListOfCLOrder, kolabel,wtlabel, scaling = "unscaled",minColorRange=NULL,maxColorRange=NULL){#phenotypes = "all".

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
                         # annotation_name_gp= gpar(fontsize = 25)
                         annotation_legend_param = list(

                           title_gp = gpar(#fontsize = 9,
                             fontface = "bold")#,
                           #labels_gp = gpar(fontsize = 9)
                         ),
                         annotation_name_gp= gpar(#fontsize = 11,
                           fontface = "bold")
                         # labels_gp= gpar(fontsize = 25)
  )




  CLOrder <- colnames(HeatmapDF)

  # BAFSubtypeDFNoWT <- subset(BAFSubtypeDF,CL!=wildtype)

  BAFSubtypeDFNoWT <- BAFSubtypeDF %>% arrange(match(CL,CLOrder))

  # BAFAnn = HeatmapAnnotation(Subtype = BAFSubtypeDF$Subtype, col =SubtypeColors,
  #                            which="col",
  #                            gap = unit(1, 'mm'),
  #                            annotation_name_gp= gpar(fontsize = 25)
  #                            # labels_gp= gpar(fontsize = 25)
  #                            )

  BAFAnn = HeatmapAnnotation(Subtype = BAFSubtypeDF$Subtype, col =SubtypeColors,
                             which="col",
                             gap = unit(1, 'mm'),
                             # heatmap_legend_param = list(legend_gp = gpar(fontsize = 19)),
                             annotation_legend_param = list(

                               title_gp = gpar(#fontsize = 9,
                                 fontface = "bold")#,
                               #labels_gp = gpar(fontsize = 9)
                             ),
                             annotation_name_gp= gpar(#fontsize =11,
                               fontface = "bold")
                             # annotation_name_gp= gpar(fontsize = 25)
                             # labels_gp= gpar(fontsize = 25)
  )


  if(isTRUE(scaling == "unscaled")==TRUE){
    minval <- min(CompleteAUCDFSubset$mean_log_AUC_value_ratio)
    maxval <- max(CompleteAUCDFSubset$mean_log_AUC_value_ratio)
  } else if(isTRUE(scaling == "scaled")==TRUE){
    minval <- minColorRange
    maxval <- maxColorRange
  }

  # comprehensive_heatmap <- Heatmap(HeatmapDF,cell_fun = function(j, i, x, y, w, h, fill) {
  #
  #   if(isTRUE(PValueDF[i,j] < 0.001)==TRUE) {
  #     grid.text("****",x,y,gp = gpar(fontsize = 30))
  #   } else if(isTRUE(PValueDF[i,j] < 0.01)==TRUE) {
  #     grid.text("***",x,y,gp = gpar(fontsize = 30))
  #   }
  #   else if(isTRUE(PValueDF[i,j]  < 0.05)==TRUE) {
  #     grid.text("**",x,y,gp = gpar(fontsize = 30))
  #   }
  #   else if(isTRUE(PValueDF[i,j]  < 0.1)==TRUE) {
  #     grid.text("*",x,y,gp = gpar(fontsize = 30))
  #   }
  # },
  # # right_annotation=ha,
  # row_order= ListOfCompoundOrder,
  # column_names_gp = grid::gpar(fontsize = 7),
  # row_names_gp = grid::gpar(fontsize = 7),
  # column_order = ListOfCLOrder, name = paste("\nLog2", kolabel, "/",wtlabel,"AUC\n"),
  # use_raster = TRUE,column_title = paste(PhenotypeName,"\n"),column_title_gp = gpar(fontsize = 31,fontface = "bold"),
  # # heatmap_legend_param = list(legend_gp = gpar(fontsize = 19)),
  # col = colorRamp2(c(minval, 0, maxval),c("blue", "white", "red")),
  # row_title_gp = gpar(fontsize = 14),
  # heatmap_legend_param = list(
  #
  #   title_gp = gpar(fontsize = 12,
  #                   fontface = "bold"),
  #   labels_gp = gpar(fontsize = 10),
  #   # legend_height = unit(4, "cm"),
  #   # legend_width = unit(7, "cm"),
  #   direction = "horizontal"),
  # width = unit(8, "cm"), height = unit(2.5, "cm")
  #
  # # top_annotation=BAFAnn
  # )

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

  # comprehensive_heatmap <- draw(comprehensive_heatmap)#, heatmap_legend_side = "bottom"

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
#' @param columnk Numerical variable representing the user-specified k number of clusters for the columns of the heatmap.
#' @param rowtextsize Numerical value representing text size of row labels.
#' @param scaling Variable denoting whether minimum and maximum values for color range should be user-defined values ("scaled") or simply the minimum and maximum values of the input data ("unscaled"). Unscaled is the default.
#' @param minColorRange Minimum value for color range. Default is NULL.
#' @param maxColorRange Maximum value for color range. Default is NULL.
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

  # minval <- min(AUCDF5PhenotypesDF[1:5])
  # maxval <- max(AUCDF5PhenotypesDF[1:5])

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
        # row_split = k2m_data$cluster,
        row_split =cbind(k2m_data$cluster),
        # row_km = as.numeric(as.character(bestClustMethod$Cluster_Num)),
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
        # column_km = columnk,
        # column_km_repeats = 100,
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
        # column_km = columnk,
        # column_km_repeats = 100,
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

