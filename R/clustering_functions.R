#' Run C-Means Clustering

#'
#' This function runs c-means clustering on either the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes ("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency") or the mean
#'  AUC values normalized to WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied; and returns the performance of the clustering algorithm applied with user-specified number of clusters.
#' @param InputDF2 A dataframe that either contains the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  or the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied.
#' @param AUCDF5PhenotypesDF A matrix containing the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes.
#' @param LabelTypesList A list of label types (defined by the user) to cluster with and validate.
#' @param K Number of clusters to test.
#' @param UnsupervisedResultsDF A dataframe containing performance of other clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and external cluster validation
#' (Meila’s Variation of Information (MVI),corrected Rand Index (CRI)). Default is null unless object is otherwise passed in the function's argument.
#' @return
#'\item{UnsupervisedResultsDF}{ An updated UnsupervisedResultsDF containing the performance of hierarchical clustering described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and how well they group
#' by user-defined label types (i.e. MOA of compounds) via external cluster validation (Meila’s Variation of Information (MVI), corrected Rand Index (CRI)).}

#' @author Caroline Barry
#' @import ComplexHeatmap
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export


runCMeans <- function(InputDF2,K,AUCDF5PhenotypesDF=NULL,LabelTypesList=NULL,UnsupervisedResultsDF=NULL,phenotype=NULL,type=NULL){

  fannyx <- fanny(InputDF2, K,metric = "euclidean", stand = FALSE,memb.exp = 1.5)

  ######################## cluster validation - Dunn Index ########################
  # Dunn Index = 1 --> indicates clusters are perfectly separated & compact

  fcm_stats <- cluster.stats(dist(InputDF2), fannyx$clustering)


  ######################## cluster validation  - Connectivity ########################
  # low Connectivity --> clusters are homogenous

  ConnectivityScore <-connectivity(distance = NULL, fannyx$clustering, Data = InputDF2, neighbSize = 20,method = "euclidean")

  ######################## External Cluster Validation ########################

  if(isTRUE(is.null(phenotype))==FALSE){

    UnsupervisedResultsDF <- data.frame(matrix(ncol=7,nrow=0))

    RowInput <- c("C-means", K, fannyx$silinfo$avg.width, fcm_stats$dunn,fcm_stats$within.cluster.ss,fcm_stats$ch, ConnectivityScore,phenotype,type)

    UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)

    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Phenotype","Type")

  } else{

    UnsupervisedResultsDF <- data.frame(matrix(ncol=10,nrow=0))

    for(LabelType in LabelTypesList){

      # corrected Rand Index
      CRIValue <-cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], fannyx$clustering)$corrected.rand

      # Meila’s Variation of Information
      MVIValue <- cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], fannyx$clustering)$vi


      RowInput <- c("C-means", K, fannyx$silinfo$avg.width, fcm_stats$dunn,fcm_stats$within.cluster.ss,fcm_stats$ch, ConnectivityScore, LabelType, CRIValue, MVIValue)

      UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)

    }

    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Label_Type","CRI","MVI")

  }

  return(UnsupervisedResultsDF)

}

#########################################################################################################################################################################################################################################################################

#' Run Hierarchical Clustering

#'
#' This function runs hierarchical clustering on either the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes ("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency") or the
#'  mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied; and returns the performance of the clustering algorithm applied with user-specified number of clusters.
#' @param InputDF2 A dataframe that either contains the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  or the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied.
#' @param AUCDF5PhenotypesDF A matrix containing the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes.
#' @param LabelTypesList A list of label types (defined by the user) to cluster with and validate.
#' @param K Number of clusters to test.
#' @param UnsupervisedResultsDF A dataframe containing performance of other clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and external cluster validation
#' (Meila’s Variation of Information (MVI),corrected Rand Index (CRI)). Default is null unless object is otherwise passed in the function's argument.

#' @return
#'\item{UnsupervisedResultsDF}{ An updated UnsupervisedResultsDF containing the performance of hierarchical clustering described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and how well they group
#' by user-defined label types (i.e. MOA of compounds) via external cluster validation (Meila’s Variation of Information (MVI),corrected Rand Index (CRI)).}
#' @author Caroline Barry
#' @import ComplexHeatmap
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export


runHierarchicalClustering <- function(InputDF2,K,AUCDF5PhenotypesDF=NULL,LabelTypesList=NULL,UnsupervisedResultsDF=NULL,phenotype=NULL,type=NULL){

  ########### Hierachical Clustering - Average Linkage Method ############

  dist_man <- dist(InputDF2, method="manhattan")

  hc_m2 <- hclust(d=dist_man, method="average")

  coph_m2 <- cophenetic(hc_m2)

  groupward6 <- cutree(hc_m2, k = K)

  ######################## cluster validation for Hierachical Clustering - Average Linkage Method - Silhouette Analysis ########################

  average <- eclust(InputDF2, "hclust", k = K, hc_metric = "manhattan",hc_method = "average", graph = F)

  ######################## cluster validation for Hierachical Clustering - Average Linkage Method - Dunn Index ########################

  averagea <- cluster.stats(dist(InputDF2), average$cluster)

  ######################## cluster validation for Hierachical Clustering - Average Linkage Method - Connectivity ########################

  ConnectivityScore <-connectivity(distance = NULL, average$cluster, Data = InputDF2, neighbSize = 20,method = "euclidean")

  ######################## External Cluster Validation ########################

  if(isTRUE(is.null(phenotype))==FALSE){

    UnsupervisedResultsDF <- data.frame(matrix(ncol=7,nrow=0))


    RowInput <- c("Hierarchical", K, averagea$avg.silwidth, averagea$dunn,averagea$within.cluster.ss,averagea$ch, ConnectivityScore,phenotype,type)

    UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)

    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Phenotype","Type")

  } else{

    UnsupervisedResultsDF <- data.frame(matrix(ncol=10,nrow=0))


    for(LabelType in LabelTypesList){

      # corrected Rand Index
      CRIValue <-cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], average$cluster)$corrected.rand


      # Meila’s Variation of Information
      MVIValue <- cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], average$cluster)$vi


      RowInput <- c("Hierarchical", K, averagea$avg.silwidth, averagea$dunn,averagea$within.cluster.ss,averagea$ch, ConnectivityScore, LabelType, CRIValue, MVIValue)

      UnsupervisedResultsDF <- rbind(UnsupervisedResultsDF, RowInput)

    }

    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Label_Type","CRI","MVI")

  }



  return(UnsupervisedResultsDF)

}

#########################################################################################################################################################################################################################################################################

#' Run K-Means Clustering

#'
#' This function runs k-means clustering on either the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes ("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency") or the
#' mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied; and returns the performance of the clustering algorithm applied with user-specified number of clusters.
#' @param InputDF2 A dataframe that either contains the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  or the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied.
#' @param AUCDF5PhenotypesDF A matrix containing the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes.
#' @param LabelTypesList A list of label types (defined by the user) to cluster with and validate.
#' @param K Number of clusters to test.
#' @param UnsupervisedResultsDF A dataframe containing performance of other clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and external cluster validation
#' (Meila’s Variation of Information (MVI),corrected Rand Index (CRI)). Default is null unless object is otherwise passed in the function's argument.

#' @return
#'\item{UnsupervisedResultsDF}{ A dataframe containing the performance of k-means clustering described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and external cluster validation (Meila’s Variation of Information (MVI),corrected Rand Index (CRI)).}
#' @author Caroline Barry
#' @import ComplexHeatmap
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export

runKMeans <- function(InputDF2,K,AUCDF5PhenotypesDF=NULL,LabelTypesList=NULL,UnsupervisedResultsDF=NULL,phenotype=NULL,type=NULL){

  k2m_data <- kmeans(InputDF2, K, nstart=25)

  ######################## cluster validation for K-means - Silhouette Analysis ########################

  k2m_data <- factoextra::eclust(InputDF2, "kmeans", k = K, nstart = 25, graph = F)

  ######################## cluster validation for K-means - Dunn Index ########################
  km_stats <- cluster.stats(dist(InputDF2), k2m_data$cluster)
  km_stats$dunn

  ######################## cluster validation for K-means - CH Index ########################
  km_stats$ch

  ######################## cluster validation for K-means - Connectivity ########################
  ConnectivityScore <- connectivity(distance = NULL, k2m_data$cluster, Data = InputDF2, neighbSize = 20,
                                    method = "euclidean")

  ######################## cluster validation for K-means - WSS ########################
  #Elbow method
  km_stats$within.cluster.ss

  ######################## external cluster validation ########################

  if(isTRUE(is.null(phenotype))==FALSE){

    UnsupervisedResultsDF <- data.frame(matrix(ncol=7,nrow=0))


    RowInput <- c("K-means", K, km_stats$avg.silwidth, km_stats$dunn,km_stats$within.cluster.ss,km_stats$ch, ConnectivityScore,phenotype,type)

    UnsupervisedResultsDF <- rbind(UnsupervisedResultsDF, RowInput)

    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Phenotype","Type")

  } else{

    UnsupervisedResultsDF <- data.frame(matrix(ncol=10,nrow=0))


    for(LabelType in LabelTypesList){

      # Corrected Rand Index
      CRIValue <- cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], k2m_data$cluster)$corrected.rand

      # Meila’s Variation of Information
      MVIValue <-cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], k2m_data$cluster)$vi

      RowInput <- c("K-means", K, km_stats$avg.silwidth, km_stats$dunn,km_stats$within.cluster.ss,km_stats$ch, ConnectivityScore, LabelType, CRIValue, MVIValue)

      UnsupervisedResultsDF <- rbind(UnsupervisedResultsDF, RowInput)

    }
    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Label_Type","CRI","MVI")
  }


  return(UnsupervisedResultsDF)

}

#########################################################################################################################################################################################################################################################################

#' Run PAM Clustering

#'
#' This function runs PAM (partitioning around medoids (AKA: k-medoids)) clustering on either the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes ("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h",
#' "Second_Phase_Relative_Change_in_Confluency")  or the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied;
#' and returns the performance of the clustering algorithm applied with user-specified number of clusters.
#' @param InputDF2 A dataframe that either contains the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  or the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied.
#' @param AUCDF5PhenotypesDF A matrix containing the mean log2 fold-change AUC values KO vs WT, cell line modification, and treatment for the 5 phenotypes.
#' @param LabelTypesList A list of label types (defined by the user) to cluster with and validate.
#' @param K Number of clusters to test.
#' @param UnsupervisedResultsDF A dataframe containing performance of other clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and external cluster validation
#' (Meila’s Variation of Information (MVI),corrected Rand Index (CRI)). Default is null unless object is otherwise passed in the function's argument.

#' @return
#'\item{UnsupervisedResultsDF}{ An updated UnsupervisedResultsDF containing the performance of PAM clustering described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and external cluster validation (MVI,CRI).}

#' @author Caroline Barry
#' @import ComplexHeatmap
#' @import cluster
#' @import factoextra
#' @import dplyr
#' @import stats
#' @import clValid
#' @export

runPAM <- function(InputDF2,K, AUCDF5PhenotypesDF=NULL,LabelTypesList=NULL,UnsupervisedResultsDF=NULL,phenotype=NULL,type=NULL){


  ########### PAM - partitioning around medoids (AKA: k-medoids) ############

  pam_data <- pam(InputDF2,K)

  ######################## cluster validation for PAM - Silhouette Analysis ########################

  pam_data <- eclust(InputDF2, "pam", k = K, nstart = 25, graph = F)

  ######################## cluster validation for PAM - Dunn Index ########################

  pm_stats <- cluster.stats(dist(InputDF2), pam_data$cluster)

  ######################## cluster validation for K-means - Connectivity ########################

  ConnectivityScore <- connectivity(distance = NULL, pam_data$cluster, Data = InputDF2, neighbSize = 20,method = "euclidean")

  ######################## External Cluster Validation ########################

  if(isTRUE(is.null(phenotype))==FALSE){

    UnsupervisedResultsDF <- data.frame(matrix(ncol=8,nrow=0))


    RowInput <- c("PAM", K, pm_stats$avg.silwidth, pm_stats$dunn, pm_stats$within.cluster.ss, pm_stats$ch, ConnectivityScore, phenotype,type)

    UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)

    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Phenotype","Type")

  } else{

    UnsupervisedResultsDF <- data.frame(matrix(ncol=10,nrow=0))


    for(LabelType in LabelTypesList){

      # corrected Rand Index
      CRIValue <- cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], pam_data$cluster)$corrected.rand

      # Meila’s Variation of Information
      MVIValue <-cluster.stats(d = dist(InputDF2),AUCDF5PhenotypesDF[,c(LabelType)], pam_data$cluster)$vi


      RowInput <- c("PAM", K, pm_stats$avg.silwidth, pm_stats$dunn, pm_stats$within.cluster.ss, pm_stats$ch, ConnectivityScore, LabelType, CRIValue, MVIValue)

      UnsupervisedResultsDF<- rbind(UnsupervisedResultsDF, RowInput)

    }

    colnames(UnsupervisedResultsDF) <- c("Method","Cluster_Num","Avg_Sil_Width", "Dunn_Score","WSS","CH_Index", "Connectivity_Score","Label_Type","CRI","MVI")
  }

  return(UnsupervisedResultsDF)

}


#########################################################################################################################################################################################################################################################################

#' Comparing the Performance of Clustering Algorithms

#'
#' This function plots the performance of various clustering algorithms on either the mean AUC values normalized to WT, cell line modification, and treatment for the 5 phenotypes ("Confluency_after_24h","Confluency_after_48h","Confluency_after_72h","Confluency_after_96h","Second_Phase_Relative_Change_in_Confluency") or the mean AUC
#' values normalized to WT, cell line modification, and treatment for the 5 phenotypes after the Principal Components Analysis has been applied.
#' @param InputDF2 A dataframe that either contains the mean AUC values normalized to WT, cell line modification, and treatment for the 5 phenotypes  or the mean AUC values normalized to WT, cell line modification, and treatment for the 5 phenotypes  after Principal Components Analysis has been applied.
#' @param UnsupervisedResultsDF A dataframe containing the performance of k-means, PAM (partitioning around the medoids), hierarchical, and c-means clustering algorithms described by how well the clusters are separated (average Silhuouette Width, Dunn Score, Connectivity Score) and external cluster validation (Meila’s Variation of Information (MVI),corrected Rand Index (CRI)).
#' @param Label A variable representing a label type (defined by the user) to cluster with and use for external validation of how clusters group.

#' @author Caroline Barry
#' @import dplyr
#' @import ggplot2
#' @export

runUnsupervisedAlgorithmComparison <- function(InputDF2, UnsupervisedResultsDF, Label){

  ################ Determining the Best Clustering Algorithm Based on Clustering Analysis Results #####################

  UnsupervisedResultsDF%<>% mutate(Cluster_Size= str_pad(Cluster_Num, width = 2,side ="left",'0'))

  UnsupervisedResultsDF%<>% mutate(Method_Clust = paste(Cluster_Size,Method, sep="_"))

  UnsupervisedResultsDFUnique <- subset(UnsupervisedResultsDF, Label_Type==unique(UnsupervisedResultsDF$Label_Type)[1])

  # plot the average silhouette scores for all clustering algorithms
  AllAvgSilPlotObject<-ggplot(UnsupervisedResultsDFUnique, aes(y = as.numeric(as.character(Avg_Sil_Width)), x = Method_Clust , group = Cluster_Size,fill=Cluster_Size)) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title =  "Average Silhouette Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Average Silhouette Score") +
    theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),legend.title = element_text(size=12),
          plot.title=element_text(size=25),plot.subtitle=element_text(size=18))



  # plot the Dunn Indexes for all clustering algorithms
  AllDunnPlotObject<-ggplot(UnsupervisedResultsDFUnique, aes(y = as.numeric(as.character(Dunn_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Dunn Indexes", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Dunn Index")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))


  # plot the connectivity scores for all clustering algorithms
  AllConnPlotObject<-ggplot(UnsupervisedResultsDFUnique, aes(y = as.numeric(as.character(Connectivity_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Connectivity Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Connectivity")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))

  # plot the WSS for all clustering algorithms
  AllWSSPlotObject<-ggplot(UnsupervisedResultsDFUnique, aes(y = as.numeric(as.character(WSS)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Total Within-Cluster Sum of Squares (WSS)", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("WSS")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))

  # plot the CH index for all clustering algorithms
  AllCHPlotObject<-ggplot(UnsupervisedResultsDFUnique, aes(y = as.numeric(as.character(CH_Index)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Calinski — Harabasz (CH) Index", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("CH Index")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))

  # select best clustering algorithm and K clusters based on the method with the best CH score (how well the clusters are separated)
  bestClustMethod <- UnsupervisedResultsDFUnique %>%
    arrange(desc(as.numeric(as.character(CH_Index))))


  bestClustMethod <- bestClustMethod[1,]


  ########### now plotting the clustering accuracy for each label type (i.e. Drug MOA) ##############

  UnsupervisedResultsDFSubset <- subset(UnsupervisedResultsDF, Label_Type==Label)

  UnsupervisedResultsDFSubset %<>% mutate(Label2 = gsub("_", " ", Label))

  # plotting Corrected Rand Indexes for all clustering algorithms
  CRIPlotObject<-ggplot(UnsupervisedResultsDFSubset, aes(y = as.numeric(as.character(CRI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Corrected Rand Indexes:",unique(UnsupervisedResultsDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Corrected Rand Index")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))


  # plotting Meila's Variation of Information for all clustering algorithms
  MVIPlotObject<-ggplot(UnsupervisedResultsDFSubset, aes(y =as.numeric(as.character(MVI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Meila's Variation of Information:",unique(UnsupervisedResultsDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("MVI")+
    theme(axis.text.x = element_text(size = 8.5,angle = 90, vjust = 0.5, hjust=1))



  ################# comparing only the hierarchical clustering ###################

  HCDF <- subset(UnsupervisedResultsDF, Method=="Hierarchical")

  HCDFUnique <- subset(HCDF, Label_Type==unique(HCDF$Label_Type)[1])

  # plotting average silhouette scores
  HCAvgSilPlotObject<-ggplot(HCDFUnique, aes(y = as.numeric(as.character(Avg_Sil_Width)), x = Method_Clust , group = Cluster_Size,fill=Cluster_Size)) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title =  "Average Silhouette Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Average Silhouette Score")



  # plotting Dunn Indexes
  HCDunnPlotObject<-ggplot(HCDFUnique, aes(y = as.numeric(as.character(Dunn_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Dunn Indexes", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Dunn Index")



  # plotting connectivity scores
  HCConnPlotObject<-ggplot(HCDFUnique, aes(y = as.numeric(as.character(Connectivity_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Connectivity Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Connectivity")

  # plotting WSS
  HCWSSPlotObject<-ggplot(HCDFUnique, aes(y = as.numeric(as.character(WSS)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Total Within-Cluster Sum of Squares (WSS)", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("WSS")

  # plotting CH Index
  HCCHPlotObject<-ggplot(HCDFUnique, aes(y = as.numeric(as.character(CH_Index)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Calinski — Harabasz (CH) Index", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("CH Index")



  # now plotting the clustering accuracy for each label type (i.e. Drug MOA)

  HCDFSubset <- subset(HCDF, Label_Type==Label)

  HCDFSubset %<>% mutate(Label2 = gsub("_", " ", Label))

  # plotting Corrected Rand Indexes
  HCCRIPlotObject<-ggplot(HCDFSubset, aes(y = as.numeric(as.character(CRI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Corrected Rand Indexes:",unique(HCDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Corrected Rand Index")



  # plotting Meila's Variation of Information
  HCMVIPlotObject<-ggplot(HCDFSubset, aes(y = as.numeric(as.character(MVI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Meila's Variation of Information:",unique(HCDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("MVI")



  ################# comparing only the k-means ###################

  KMDF <- subset(UnsupervisedResultsDF, Method=="K-means")

  KMDFUnique <- subset(KMDF, Label_Type==unique(KMDF$Label_Type)[1])

  # plotting average silhouette scores
  KMAvgSilPlotObject<-ggplot(KMDFUnique, aes(y = as.numeric(as.character(Avg_Sil_Width)), x = Method_Clust , group = Cluster_Size,fill=Cluster_Size)) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title =  "Average Silhouette Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Average Silhouette Score")



  # plotting Dunn Indexes
  KMDunnPlotObject<-ggplot(KMDFUnique, aes(y = as.numeric(as.character(Dunn_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Dunn Indexes", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Dunn Index")



  # plotting connectivity scores
  KMConnPlotObject<-ggplot(KMDFUnique, aes(y = as.numeric(as.character(Connectivity_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Connectivity Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Connectivity")


  # plotting CH Index
  KMCHPlotObject<-ggplot(KMDFUnique, aes(y = as.numeric(as.character(CH_Index)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Calinski — Harabasz (CH) Index", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("CH Index")


  # plotting WSS
  KMWSSPlotObject<-ggplot(KMDFUnique, aes(y = as.numeric(as.character(WSS)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Total Within-Cluster Sum of Squares (WSS)", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("WSS")


  # now plotting the clustering accuracy for each label type (i.e. Drug MOA)

  KMDFSubset <- subset(KMDF, Label_Type==Label)

  KMDFSubset %<>% mutate(Label2 = gsub("_", " ", Label))

  # plotting Corrected Rand Indexes
  KMCRIPlotObject<-ggplot(KMDFSubset, aes(y = as.numeric(as.character(CRI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Corrected Rand Indexes:",unique(KMDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Corrected Rand Index")



  # plotting Meila's Variation of Information
  KMMVIPlotObject<-ggplot(KMDFSubset, aes(y = as.numeric(as.character(MVI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Meila's Variation of Information:",unique(KMDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("MVI")



  ################# comparing only the PAM ###################

  PAMDF <- subset(UnsupervisedResultsDF, Method=="PAM")

  PAMDFUnique <- subset(PAMDF, Label_Type==unique(PAMDF$Label_Type)[1])

  # plotting average silhouette scores
  PAMAvgSilPlotObject<-ggplot(PAMDFUnique, aes(y = as.numeric(as.character(Avg_Sil_Width)), x = Method_Clust , group = Cluster_Size,fill=Cluster_Size)) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title =  "Average Silhouette Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Average Silhouette Score")



  # plotting Dunn Indexes
  PAMDunnPlotObject<-ggplot(PAMDFUnique, aes(y = as.numeric(as.character(Dunn_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Dunn Indexes", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Dunn Index")



  # plotting connectivity scores
  PAMConnPlotObject<-ggplot(PAMDFUnique, aes(y = as.numeric(as.character(Connectivity_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Connectivity Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Connectivity")


  # plotting CH Index
  PAMCHPlotObject<-ggplot(PAMDFUnique, aes(y = as.numeric(as.character(CH_Index)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Calinski — Harabasz (CH) Index", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("CH Index")

  # plotting WSS
  PAMWSSPlotObject<-ggplot(PAMDFUnique, aes(y = as.numeric(as.character(WSS)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Total Within-Cluster Sum of Squares (WSS)", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("WSS")

  # now plotting the clustering accuracy for each label type (i.e. Drug MOA)

  PAMDFSubset <- subset(PAMDF, Label_Type==Label)

  PAMDFSubset %<>% mutate(Label2 = gsub("_", " ", Label))

  # plotting Corrected Rand Indexes
  PAMCRIPlotObject<-ggplot(PAMDFSubset, aes(y = as.numeric(as.character(CRI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Corrected Rand Indexes:",unique(PAMDF$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Corrected Rand Index")


  # plotting Meila's Variation of Information
  PAMMVIPlotObject<-ggplot(PAMDFSubset, aes(y = MVI, x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Meila's Variation of Information:",unique(PAMDF$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("MVI")



  ################# comparing only the c-means ###################

  CMDF <- subset(UnsupervisedResultsDF, Method=="C-means")

  CMDFUnique <- subset(CMDF, Label_Type==unique(CMDF$Label_Type)[1])

  # plotting average silhouette scores
  CMAvgSilPlotObject<-ggplot(CMDFUnique, aes(y = as.numeric(as.character(Avg_Sil_Width)), x = Method_Clust , group = Cluster_Size,fill=Cluster_Size)) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title =  "Average Silhouette Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Average Silhouette Score")



  # plotting Dunn Indexes
  CMDunnPlotObject<-ggplot(CMDFUnique, aes(y = as.numeric(as.character(Dunn_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Dunn Indexes", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Dunn Index")



  # plotting connectivity scores
  CMConnPlotObject<-ggplot(CMDFUnique, aes(y = as.numeric(as.character(Connectivity_Score)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Connectivity Scores", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Connectivity")

  # plotting CH Index
  CMCHPlotObject<-ggplot(CMDFUnique, aes(y = as.numeric(as.character(CH_Index)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Calinski — Harabasz (CH) Index", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("CH Index")

  # plotting WSS
  CMWSSPlotObject<-ggplot(CMDFUnique, aes(y = as.numeric(as.character(WSS)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = "Total Within-Cluster Sum of Squares (WSS)", subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("WSS")


  # now plotting the clustering accuracy for each label type (i.e. Drug MOA)

  CMDFSubset <- subset(CMDF, Label_Type==Label)

  CMDFSubset %<>% mutate(Label2 = gsub("_", " ", Label))

  # plotting Corrected Rand Indexes
  CMCRIPlotObject<-ggplot(CMDFSubset, aes(y = as.numeric(as.character(CRI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Corrected Rand Indexes:",unique(CMDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("Corrected Rand Index")


  # plotting Meila's Variation of Information
  CMMVIPlotObject<-ggplot(CMDFSubset, aes(y = as.numeric(as.character(MVI)), x = Method_Clust, group = Cluster_Size,fill=Cluster_Size )) +
    geom_bar(stat= "identity") +
    theme_minimal()+
    labs(title = paste("Meila's Variation of Information:",unique(CMDFSubset$Label2)), subtitle = "for each Clustering Algorithm") +
    xlab("Clustering Algorithm") +
    ylab("MVI")



  return(list(AllMethods = list(AvgSil = AllAvgSilPlotObject, Dunn = AllDunnPlotObject, Connectivity = AllConnPlotObject,
                                WSS = AllWSSPlotObject, CH = AllCHPlotObject,
                                LabelAccuracy = list(CRI = CRIPlotObject, MVI = MVIPlotObject)),
              HC = list(AvgSil = HCAvgSilPlotObject, Dunn = HCDunnPlotObject, Connectivity = HCConnPlotObject,
                        WSS = HCWSSPlotObject, CH = HCCHPlotObject,
                        LabelAccuracy = list(CRI = HCCRIPlotObject, MVI = HCMVIPlotObject)),
              KM = list(AvgSil = KMAvgSilPlotObject, Dunn = KMDunnPlotObject, Connectivity = KMConnPlotObject,
                        WSS = KMWSSPlotObject, CH = KMCHPlotObject,
                        LabelAccuracy = list(CRI = KMCRIPlotObject, MVI = KMMVIPlotObject)),
              PAM = list(AvgSil = PAMAvgSilPlotObject, Dunn = PAMDunnPlotObject, Connectivity = PAMConnPlotObject,
                         WSS = PAMWSSPlotObject, CH = PAMCHPlotObject,
                         LabelAccuracy = list(CRI = PAMCRIPlotObject, MVI = PAMMVIPlotObject)),
              CM = list(AvgSil = CMAvgSilPlotObject, Dunn = CMDunnPlotObject, Connectivity = CMConnPlotObject,
                        WSS = CMWSSPlotObject, CH = CMCHPlotObject,
                        LabelAccuracy = list(CRI = CMCRIPlotObject, MVI = CMMVIPlotObject)),
              bestClustMethod = bestClustMethod))

}
