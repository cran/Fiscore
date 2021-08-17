#' @title cluster_ID
#'
#' @description Function to select an optimal number of clusters and a model to be fitted during the EM phase of clustering for Gaussian Mixture Models. The function provides summaries and helps to visualise clusters based on Fi-score using scatter plotting and dimension reduction plots.
#'
#' @param pdb_df data frame containing processed PDB file with Fi-score values
#' @param max_range number of clusters to consider during model selection; default 20 clusters
#' @param secondary_structures include information on secondary structure elements from PDB when plotting, default value is TRUE
#' @param clusters number of clusters to test not based on the best BIC output, user also needs to supply modelNames
#' @param modelNames can only be supplied when clusters are also specified, this option will model based on the user parameters
#'
#' @return A data frame object that contains a summary of clusters
#' @ImportFrom ggplot2 aes
#' @ImportFrom ggplot2 ggplot
#' @ImportFrom ggplot2 geom_point
#' @ImportFrom ggplot2 ggtitle
#' @ImportFrom plotly  ggplotly
#' @import mclust
#' @ImportFrom mclust mclustBIC
#' @ImportFrom mclust Mclust
#' @ImportFrom mclust MclustDR
#' @ImportFrom methods show
#' @export
#' @examples
#' path_to_processed_PDB<- system.file("extdata", "pdb_df.tabular", package="Fiscore")
#' # basic usage of cluster_ID
#' pdb_df<-read.table(path_to_processed_PDB)
#' head(cluster_ID(pdb_df))

cluster_ID<-function(pdb_df, max_range=20, secondary_structures=TRUE, clusters=NULL, modelNames=NULL){



  #prepare data frame
  df<-as.data.frame(pdb_df$"df_resno")
  df<-cbind(df,pdb_df$"Fi_score")
  colnames(df)<-c("Residue_number","Fi_score")
  rownames(df)<-df$"Residue_number"

  #calculate Bayesian information criterion and plot different GMM models
  BIC <- mclust::mclustBIC(df, G=1:max_range)
  plot(BIC,with="BIC")
  #report the best cluster value
  print(summary(BIC))

  #To select other BIC values based on the report
  if(is.null(clusters)){
  model <- mclust::Mclust(df, x = BIC)}
  if(!is.null(clusters)&&!is.null(modelNames)){
    model <- mclust::Mclust(df,G=clusters, modelNames=modelNames)}


  #prepare a model for reporting and plotting by extracting relevant information
  model_report<-as.data.frame(model$"data")
  model_report$"Cluster"<-as.factor(model$"classification")
  model_report$"Type"<-pdb_df$"Type"

if(secondary_structures==TRUE){

  #to avoid namescape conflicts
  Residue_number_val<-model_report$"Residue_number"
  Cluster_val<-model_report$"Cluster"
  Type_val<-model_report$"Type"
  Fi_score_val<-model_report$"Fi_score"

  plot<-ggplot2::ggplot(model_report, ggplot2::aes(x=Residue_number_val, y=Fi_score_val, color=Cluster_val, shape=Type_val))+ggplot2::geom_point(alpha=0.5,size=2)+ggplot2::ggtitle(label="Cluster distribution across the protein")
  plot<-plotly::ggplotly(plot)
  methods::show(plot)  }else if(secondary_structures==FALSE){
    #to avoid namescape conflicts
    Residue_number_val<-model_report$"Residue_number"
    Cluster_val<-model_report$"Cluster"
    Fi_score_val<-model_report$"Fi_score"

    plot<-ggplot2::ggplot(model_report, ggplot2::aes(x=Residue_number_val, y=Fi_score_val, color=Cluster_val))+ggplot2::geom_point(alpha=0.5,size=2)+ggplot2::ggtitle(label="Cluster distribution across the protein")

  plot<-plotly::ggplotly(plot)
  methods::show(plot)  }



  #Dimension reduction based clustering visualisation

  model_dir <- mclust::MclustDR(model)
  print(summary(model_dir))
  plot(model_dir, what = "scatterplot", main="Distribution of structure feature clusters ")


  return(model_report)
}

