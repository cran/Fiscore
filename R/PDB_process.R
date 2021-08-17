#' @title PDB_process
#'
#' @description Function to preprocess and inspect a PDB file
#'
#' @param file_name PDB file name to load, e.g. '6KZ5.pdb'
#' @param path location where to transfer split PDB files, default will create a new directory in your working directory
#' @return   generates split chain PDB files in the default or selected directory and then returns the names of the files
#' @import bio3d
#' @export
#' @examples
#' path_to_PDB_file<- system.file("extdata", "3nf5.pdb", package="Fiscore")
#' # basic usage of PDB_process calls the selected path to load a large file
PDB_process<-function(file_name, path = "split_PDB"){



  pdb_file <- bio3d::read.pdb(file_name) #read a PDB file

  #split PDB file into separate chains
  split_pdb<-bio3d::pdbsplit(file_name, path=path)

  #extract chain names using the convention of _ to indicate chain
  files<-list.files(path, pattern = "_",full.names=FALSE)

  return (files) #return list of split files
}



