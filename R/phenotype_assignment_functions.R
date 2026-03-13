#' Computing the First Derivative of Growth Curves

#'
#' This function computes the first derivative of the growth curve with the intention of being uploaded to MySQL.
#' @param IncucyteDataAndMetaDataAndConditionsDF Dataframe that combines Incucyte Confluency Dataset, Meta Data and Condition ID table.

#' @return
#'\item{IncucyteDataFirstDerivativeDF}{ Dataframe containing the Incucyte Confluency Dataset, Meta Data, Condition ID table and first derivative of the growth curve(s).}
#' @author Caroline Barry
#' @import dplyr
#' @export


getFirstDerivative <- function(IncucyteDataAndMetaDataAndConditionsDF){
  # (d1+d2+d3/3)-(d2+d3+d4/3)
  # need to filter out measurement intervals not equal to 6!
  # group by bio rep and subtract lag and lead timepoint from one another
  # filter out those with difference != 6
  IncucyteDataAndMetaDataAndConditionsDFRounded <- IncucyteDataAndMetaDataAndConditionsDF %>%
    group_by(bio_rep) %>%
    mutate(simple_time = round(timepoint))

  IncucyteDataAndMetaDataAndConditionsDFRounded <- subset(IncucyteDataAndMetaDataAndConditionsDFRounded,IncucyteDataAndMetaDataAndConditionsDFRounded$simple_time %% 6 == 0)

  IncucyteDataAndMetaDataAndConditionsDFRounded <-IncucyteDataAndMetaDataAndConditionsDFRounded[,-28]


  # need to avg tech_reps across bio_reps
  # calculate average of 3 tech reps per bio rep
  IncucyteAveragedDataAndMetaDataAndConditionsDF<- IncucyteDataAndMetaDataAndConditionsDFRounded %>%
    group_by(bio_rep,project_ID, plate_ID,timepoint,condition_ID) %>%
    dplyr::summarise(mean_conf=mean(confluency))

  # averaging over three conf data points but (conf1, conf2, conf3)/3 will belong to timepoint 1 and so on....
  # will have to leave out last two timepoints
  ConfluencyAveraged3DF <- do.call(rbind,
                                   sapply(seq(1, nrow(IncucyteAveragedDataAndMetaDataAndConditionsDF), 1), function(i){
                                     x <- IncucyteAveragedDataAndMetaDataAndConditionsDF[ i:(i + 2),6 , drop = FALSE]#5
                                     res <- rbind(x, colSums(x)/3) #colSums
                                     rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                                     res
                                   }))

  # get well ID labels back!
  # replaced desired_df with conf_avg_df
  IncucyteAveraged3DataAndMetaDataAndConditionsDF <- cbind(IncucyteAveragedDataAndMetaDataAndConditionsDF, avg_conf_3 = ConfluencyAveraged3DF[,4])

  # #these avg_conf_3 need to be assigned to the second/mid timepoint
  IncucyteAveraged3DataAndMetaDataAndConditionsDF <- IncucyteAveraged3DataAndMetaDataAndConditionsDF %>%
    group_by(bio_rep) %>%
    dplyr::mutate(avg_conf_3 = replace(avg_conf_3, tail(row_number(), 2), NA))

  IncucyteAveraged3DataAndMetaDataAndConditionsDF <- IncucyteAveraged3DataAndMetaDataAndConditionsDF %>%
    group_by(bio_rep) %>%
    mutate(avg_conf_3= dplyr::lag(avg_conf_3))

  IncucyteDataFirstDerivativeDF <- IncucyteAveraged3DataAndMetaDataAndConditionsDF %>%
    group_by(bio_rep) %>%
    mutate(conf_diff = avg_conf_3 -  dplyr::lag(avg_conf_3, default =  dplyr::first(avg_conf_3)))

  # these conf_diff need to be assigned to a newly created timepoint in between the two timepoints being subtracted
  FirstDerivativeTimepointsDF <- do.call(rbind,
                                         sapply(seq(1, nrow(IncucyteAveragedDataAndMetaDataAndConditionsDF), 1), function(i){
                                           x <- IncucyteDataFirstDerivativeDF[ i:(i + 1),4 , drop = FALSE]#5
                                           res <- rbind(x, colSums(x)/2) #colSums
                                           rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                                           res
                                         }))
  #get well ID labels back!
  #replaced desired_df with conf_avg_df
  IncucyteDataFirstDerivativeDF <- cbind(IncucyteDataFirstDerivativeDF, mid_timepoint = FirstDerivativeTimepointsDF[,3])

  IncucyteDataFirstDerivativeDF <- IncucyteDataFirstDerivativeDF %>%
    group_by(bio_rep) %>%
    mutate(conf_diff= dplyr::lead(conf_diff))

  return(IncucyteDataFirstDerivativeDF)
}

#########################################################################################################################################################################################################################################################################

#' Computing the Second Derivative of Growth Curves

#'
#' This function computes the second derivative of the growth curve with the intention of being uploaded to MySQL.
#' @param IncucyteDataFirstDerivativeDF Dataframe containing Meta Data, Condition IDs, Confluency data, and first derivative of the growth curves.
#' @param treatmentstarttime Numerical variable represents the time at which treatment begins.
#' @param treatmentendtime Numerical variable represents the time at which treatment ends.

#' @return
#'\item{IncucyteFirstAndSecondDerivativeWithSplitPointsDF}{ Dataframe containing Meta Data, Condition IDs, Confluency data, and both first and second derivatives of the growth curve(s).}
#' @author Caroline Barry
#' @import dplyr
#' @export

getSecondDerivativeAndSplitPoint <- function(IncucyteDataFirstDerivativeDF,treatmentstarttime, treatmentendtime){

  IncucyteDataAndMetaDataAndConditionsFirstDerivativeDF24hr <- subset(IncucyteDataFirstDerivativeDF, timepoint>=treatmentstarttime & timepoint<=treatmentendtime)

  IncucyteDataAndMetaDataAndConditionsFirstAndSecondDerivativeDF <- IncucyteDataAndMetaDataAndConditionsFirstDerivativeDF24hr %>%
    group_by(bio_rep) %>%
    mutate(DiffConf_diff = conf_diff - dplyr::lag(conf_diff, default = dplyr::first(conf_diff)))

  # need to assign DiffConf_diff to an original timepoint inbetween the two timepoints being subtracted
  # identifying min diff of conf_diff
  MinimumSecondDerivativeDF <- IncucyteDataAndMetaDataAndConditionsFirstAndSecondDerivativeDF %>%
    group_by(bio_rep) %>%
    slice(which.min(DiffConf_diff))

  # merge split_df with conf_diff_df
  MinimumSecondDerivativeDF <-MinimumSecondDerivativeDF %>%
    mutate(min_diff_diff_timepoint = 1)

  IncucyteFirstAndSecondMimimumDerivativeDF <- merge(IncucyteDataAndMetaDataAndConditionsFirstAndSecondDerivativeDF, MinimumSecondDerivativeDF,
                                                     by = c("bio_rep", "timepoint","condition_ID","project_ID","plate_ID","mean_conf","avg_conf_3","mid_timepoint", "conf_diff"), all.x = TRUE, all.y = TRUE)

  IncucyteFirstAndSecondMimimumDerivativeDF<-IncucyteFirstAndSecondMimimumDerivativeDF[,-12]

  IncucyteFirstAndSecondMimimumDerivativeDF <-IncucyteFirstAndSecondMimimumDerivativeDF %>%
    dplyr::rename(
      DiffConf_diff = DiffConf_diff.x,
      min_diff_diff = DiffConf_diff.y)

  # filtering out sudden increases in conf_diff after reaching max conf_diff
  IncucyteFirstAndSecondMimimumDerivativeDF <- IncucyteFirstAndSecondMimimumDerivativeDF %>%
    group_by(bio_rep) %>%
    dplyr::filter(cumany(DiffConf_diff==min_diff_diff)) %>%
    slice(seq_len(match(0, DiffConf_diff, nomatch = n())))%>%
    ungroup

  IncucyteRecoveryOrDeathCasesDF <- IncucyteFirstAndSecondMimimumDerivativeDF %>%
    group_by(bio_rep) %>%
    dplyr::filter(conf_diff<0 |DiffConf_diff >= 0) %>%
    dplyr::filter(row_number()== 1) %>%
    dplyr::mutate(last_timepoint = 1)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF <- merge(IncucyteDataAndMetaDataAndConditionsFirstAndSecondDerivativeDF, IncucyteRecoveryOrDeathCasesDF, by = c("bio_rep","timepoint", "mean_conf",
                                                                                                                                                                    "condition_ID","project_ID", "plate_ID","mid_timepoint", "conf_diff", "avg_conf_3","DiffConf_diff"), all.x = TRUE, all.y = TRUE)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF[,-11]

  return(IncucyteFirstAndSecondDerivativeWithSplitPointsDF)

}

#########################################################################################################################################################################################################################################################################

#' Creating a Split Point Table

#'
#' This function computes modeling-independent phenotypes: confluency after 24hr, 48hr, 72hr, and 96 hr of treatment, along with the change of confluency in and duration of the first phase, and ratio of change in confluency after 96hr/48hr
#'  of treatment. This function also creates and returns a Split Point Table.
#' @param IncucyteFirstAndSecondDerivativeWithSplitPointsDF Dataframe containing Meta Data, Condition IDs, Confluency data, and both first and second derivatives of the growth curve(s).
#' @param treatmentstarttime Numerical variable represents the time at which treatment begins.
#' @param treatmentendtime Numerical variable represents the time at which treatment ends.
#' @return
#'\item{IncucyteSplitPointTable}{ Dataframe containing computed desired phenotypes:
#'project_ID
#'plate_ID
#'condition_ID
#'logistic_growth?
#'logistic_growth_duration
#'second_phase?
#'second_phase_duration
#'split_timepoint
#'split_first_conf_difference
#'after_24h_confluency
#'after_48h_confluency
#'after_72h_confluency
#'last_confluency
#'conf_96hr_48hr_ratio
#'phenotype}
#' @author Caroline Barry
#' @import dplyr
#' @export


createSplitPointTable <- function(IncucyteFirstAndSecondDerivativeWithSplitPointsDF,treatmentstarttime, treatmentendtime){


  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr <- subset(IncucyteFirstAndSecondDerivativeWithSplitPointsDF, timepoint>=treatmentstarttime & timepoint<=treatmentendtime)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstTP <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    group_by(bio_rep) %>%
    mutate(logistic_growth_phase = 1) %>%
    slice(1)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstTP <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstTP %>%
    dplyr::rename(first_confluency = mean_conf,
                  first_tp = timepoint)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrLastTP <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    group_by(bio_rep) %>%
    mutate(logistic_growth_phase = 1) %>%
    slice(n())

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrLastTP <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrLastTP %>%
    dplyr::rename(last_confluency = mean_conf,
                  last_tp = timepoint)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstAndLastTP <- merge(IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstTP,IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrLastTP,
                                                                               by = c("bio_rep", "condition_ID", "project_ID", "plate_ID", "logistic_growth_phase"))


  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstAndLastTP <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstAndLastTP[,-c(8:12,15:19)]


  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>% left_join(IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hrFirstAndLastTP,
                                                                                                                               by = c("bio_rep", "condition_ID", "project_ID", "plate_ID"))


  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter24hrConf <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    group_by(bio_rep) %>%
    mutate(after_24h_confluency = case_when(timepoint >= first_tp + 24 ~ mean_conf))

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter24hrConf <-IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter24hrConf%>%
    group_by(bio_rep) %>%
    dplyr::filter(!is.na(after_24h_confluency)) %>%
    slice(1)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter48hrConf <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    group_by(bio_rep) %>%
    mutate(after_48h_confluency = case_when(timepoint >= first_tp + 48 ~ mean_conf))

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter48hrConf <-IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter48hrConf%>%
    group_by(bio_rep) %>%
    dplyr::filter(!is.na(after_48h_confluency)) %>%
    slice(1)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter72hrConf <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    group_by(bio_rep) %>%
    mutate(after_72h_confluency = case_when(timepoint >= first_tp + 72 ~ mean_conf))

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter72hrConf <-IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter72hrConf%>%
    group_by(bio_rep) %>%
    dplyr::filter(!is.na(after_72h_confluency)) %>%
    slice(1)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter24hrConf <- IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter24hrConf[, c(1,4:6, 17)]

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter48hrConf <- IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter48hrConf[, c(1,4:6, 17)]

  IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter72hrConf <- IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter72hrConf[, c(1,4:6, 17)]

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>% left_join(IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter24hrConf,
                                                                                                                               by = c("bio_rep", "condition_ID", "project_ID","plate_ID"))

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>% left_join(IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter48hrConf,
                                                                                                                               by = c("bio_rep", "condition_ID", "project_ID","plate_ID"))

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>% left_join(IncucyteFirstAndSecondDerivativeWithSplitPointsDFAfter72hrConf,
                                                                                                                               by = c("bio_rep", "condition_ID", "project_ID","plate_ID"))

  # logistic growth phase? yes/no (if duration = 0h) and second phase? yes/no
  SplitPointTableInfo <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    group_by(bio_rep) %>%
    dplyr::filter(last_timepoint==1)

  SplitPointTableInfo <- SplitPointTableInfo %>%
    dplyr::rename(split_conf = mean_conf,
                  split_timepoint = timepoint) %>%
    mutate(second_phase = 1)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    group_by(bio_rep) %>%
    dplyr::distinct()%>%
    dplyr::filter(row_number()==1)

  IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr <- IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr %>%
    dplyr::filter(!(bio_rep %in% SplitPointTableInfo$bio_rep))

  IncucyteSplitPointTable<- merge(IncucyteFirstAndSecondDerivativeWithSplitPointsDF24hr, SplitPointTableInfo,
                                  by = c("bio_rep", "condition_ID","project_ID",
                                         "plate_ID","mid_timepoint", "conf_diff","avg_conf_3", "DiffConf_diff",
                                         "last_timepoint", "logistic_growth_phase", "last_confluency","first_confluency",
                                         "first_tp", "last_tp", "after_24h_confluency", "after_48h_confluency", "after_72h_confluency"), all.x = TRUE, all.y = TRUE)

  IncucyteSplitPointTable <-IncucyteSplitPointTable %>%
    group_by(bio_rep) %>%
    mutate(logistic_growth_duration = case_when(second_phase==1 ~ split_timepoint-first_tp,
                                                is.na(second_phase) ~ last_tp-first_tp)) %>%
    mutate(second_phase_duration = case_when(second_phase==1 ~ last_tp-split_timepoint,
                                             is.na(second_phase) ~ 0))%>%
    mutate(split_first_conf_difference = case_when(second_phase==1 ~ split_conf-first_confluency,
                                                   is.na(second_phase) ~ last_confluency-first_confluency)) %>%
    mutate(conf_96hr_48hr_ratio = last_confluency/after_48h_confluency) %>%
    mutate(phenotype = if_else(conf_96hr_48hr_ratio >=1.05, "recovering",
                               if_else(conf_96hr_48hr_ratio <=0.95, "dying", "arrested")))%>%
    mutate(logistic_growth_phase = case_when(is.na(logistic_growth_phase) ~ 1,
                                             TRUE ~ logistic_growth_phase)) %>%
    mutate(second_phase = case_when(is.na(second_phase) ~ 0,
                                    TRUE ~ second_phase))

  IncucyteSplitPointTable <-IncucyteSplitPointTable[,c(3,4,2,10,23,22,24,20,25,15,16,17,11,26,27)]

  # reformat table so it's in order of:
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

  return(IncucyteSplitPointTable)

}
