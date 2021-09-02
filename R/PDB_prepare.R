#' @title PDB_prepare
#'
#' @description Function to prepare a PDB file after it was pre-processed to generate Fi-score and normalised B factor values as well as secondary structure designations
#'
#' @param file_name PDB file name to load that was split into chains, e.g. '6KZ5_A.pdb'

#' @return  returns a processed data frame with Fi-score 'Fi_score', normalised B factor values 'B_normalised' and secondary structure designations
#' @ImportFrom bio3d read.pdb
#' @ImportFrom bio3d clean.pdb
#' @ImportFrom bio3d torsion.pdb
#' @ImportFrom stringr str_extract
#' @ImportFrom stats sd
#' @ImportFrom stats complete.cases
#' @ImportFrom dplyr case_when
#' @ImportFrom dplyr mutate
#' @export
#' @examples
#' path_to_processed_PDB<- system.file("extdata", "3nf5_A.pdb", package="Fiscore")
#' # you can call PDB_prepare with the set path
#' head(PDB_prepare(path_to_processed_PDB))
PDB_prepare<-function(file_name){



  #Helper functions for the analysis
 MINMAX_normalisation_func<-function(array){

    #MIN-MAX normalisation based on the input array
    #input numeric array
    #returns normalised array values

    #check for cases where all B-factor values are 0
    if((max(array)-min(array)==0)&&(min(array)==0)&&(max(array)==0)){return (0)}

    return ((array-min(array))/(max(array)-min(array)))
  }

  # Fi-score scoring and PDB file processing


  #prepare PDB file ------
  pdb_file_temp<-bio3d::read.pdb(file_name)
  #clean file to remove terminal residues and ligand data
  pdb_file_temp<-bio3d::clean.pdb(pdb_file_temp, consecutive = TRUE, force.renumber = FALSE, fix.chain = FALSE, fix.aa = TRUE, rm.wat = TRUE, rm.lig = TRUE, rm.h = FALSE, verbose = FALSE)

  #warnings
  if((length(bio3d::pdbseq(pdb_file_temp)))==0){stop("The file has no amino acid residues")}
  if((length(bio3d::pdbseq(pdb_file_temp)))<=5){stop("This is a peptide and not protein")}

  #extract specific features------------

  #prepare helix dataframe
  #
  #  Reference: https://www.wwpdb.org/documentation/file-format-content/format23/sect5.html

  #                TYPE OF HELIX          CLASS NUMBER
  #                                       (COLUMNS 39 - 40)
  #      ---------------------------------------------------
  #      Right-handed alpha (default)       1
  #      Right-handed omega                 2
  #      Right-handed pi                    3
  #      Right-handed gamma                 4
  #      Right-handed 310                   5
  #      Left-handed alpha                  6
  #      Left-handed omega                  7
  #      Left-handed gamma                  8
  #      27 ribbon/helix                    9
  #       Polyproline                       10
  #
  #
  #TYPE OF SHEET
  #
  #The sense indicates whether strand n is parallel (sense = 1) or anti-parallel (sense = -1) to strand n-1. Sense is equal to zero (0) for the #first strand of a sheet.
  #
  #TURNS
  #
  ##Turns include those sets of residues which form beta turns, i.e., have a hydrogen bond linking (C- O)i to (N-H)i +3. Turns which link residue i to i+2 (gamma-bends) may also be included. Others may be also be classified as turns.

  feature_list<-list()

  if("helix" %in% attributes(pdb_file_temp)$names){
    #test if attribute is not NULL
    if(!is.null(pdb_file_temp$helix$start)){
    helix_df<-as.data.frame(pdb_file_temp$helix)
    type<-as.vector(helix_df$type)
    helix_df$Type<-dplyr::case_when(type==1~'Right-handed alpha helix',type==2 ~ 'Right-handed omega helix',type==3 ~ 'Right-handed pi helix',type==4 ~ 'Right-handed gamma helix', type==5 ~ 'Right-handed 310 helix',type==6~'Left-handed alpha helix', type==7 ~ 'Left-handed omega helix',type==8 ~ 'Left-handed gamma helix',type==9 ~ '27 ribbon/helix helix',type==10 ~ 'Polyproline helix', TRUE ~ as.character(type))

    feature_list[["helix"]]<-helix_df}
  }
  if("sheet" %in% attributes(pdb_file_temp)$names){
    #prepare sheet data frame
    #test if attribute is not NULL
    if(!is.null(pdb_file_temp$sheet$start)){
    sheet_df<-as.data.frame(pdb_file_temp$sheet)
    type<-as.vector(sheet_df$sense)
    sheet_df$Type<-dplyr::case_when(type==0~'Parallel sheet',type==1~'Parallel sheet',type==-1 ~ 'Antiparalel sheet', TRUE ~ as.character(type))
    feature_list[["sheet"]]<-sheet_df}
  }
  if("turn" %in% attributes(pdb_file_temp)$names){
    #prepare turn data frame
    #test if attribute is not NULL
    if(!is.null(pdb_file_temp$turn$start)){
    turn_df<-as.data.frame(pdb_file_temp$turn)
    type<- as.vector(turn_df$turnId)
    turn_df$Type<-type
    feature_list[["turn"]]<-turn_df}
  }

   #calculate torsion angles-------------------------
  torsion_angles<-bio3d::torsion.pdb(pdb_file_temp)
  #extract torsion angle table
  pdb_df<-torsion_angles$tbl
  #leave rows that contain full dihedral angle visualization
  #NOTE: terminal residues do not contain all of the angles
  pdb_df<-pdb_df[stats::complete.cases(pdb_df[,c("phi","psi")]),]
  #extract residue numbers
  df_resno<-as.numeric(sapply(rownames(pdb_df), function(x){stringr::str_extract(x,"[0-9]{1,4}")}))
  #extract residue names
  df_res<-as.vector(sapply(rownames(pdb_df), function(x){stringr::str_extract(x,"[A-Z]{3}")}))

  #construct the dataframe to contain residue names and numbers
  pdb_df<-cbind.data.frame(pdb_df,df_resno)
  pdb_df<-cbind.data.frame(pdb_df, df_res)

  #B-factor extraction ----
  #Preparing full data frame that includes dihedral angles coordinates and residue info
  #Extracting B factor information for C alpha atom to match dihedral angles
  pdb_b<-pdb_file_temp$atom[which((pdb_file_temp$atom$elety=="CA")&(pdb_file_temp$atom$resno%in%pdb_df$"df_resno")),c("resno","resid","b")]
  #Adding B factor information
  pdb_df$B_factor<-pdb_b$b
  if(all(pdb_df$B_factor==0)){warning("All B-factors are 0 and the analysis will be limited")}
  #***B-factor normalization and adding norm column ----

  pdb_df$B_normalised<-MINMAX_normalisation_func(pdb_df$B_factor)

  #Fi-score ------
  #Generate Fi-score per residue and store in the dataframe
  #calculate fi_score for the whole protein and individual aa


  fi_score<-c()


  #phi/psi normalization

  #normalization is only scaled by SD
  psi_SD<-stats::sd(pdb_df$psi)
  phi_SD<-stats::sd(pdb_df$phi)
  for(index in seq(1,nrow(pdb_df))){
    fi_score<-c(fi_score, ((pdb_df[index,'phi'])*(pdb_df[index,'psi'])*pdb_df[index,'B_normalised']/(psi_SD*phi_SD))) }

  pdb_df$Fi_score<-fi_score

  #Incorporate information for Ca to indicate what secondary structure element it belongs to

  Type_vals<-rep("NA",nrow(pdb_df))
  pdb_df$Type<-Type_vals

  if(length(feature_list)!=0){

    for(feature in feature_list){
      #feature is a data frame that contains available information on helix, sheet or turn

      for(i in 1:nrow(feature)){

        #set feature range
        range<-feature[i,"start"]:feature[i,"end"]
        #mutate column to add what type of the secondary structure the residue belongs to
        pdb_df<-dplyr::mutate(pdb_df,Type=ifelse((pdb_df$df_resno%in%range),feature[i,"Type"],pdb_df$Type))    }
    }

  }




 return(pdb_df)
}
