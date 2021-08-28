


get_geneID2GO_custom <- function() {
  
  annotations_file_path <- file.path(annotations_base_path, "topgo_viseago_custom_annotations.txt")
  GTOGO_custom <- read.table(annotations_file_path, 
                             sep="\t", comment="!", quote = "\"", 
                             header = TRUE)
  
  # convert from table format to list format
  geneID2GO_custom <- by(GTOGO_custom$GOID, GTOGO_custom$gene_id, function(x) as.character(x))
  #examine result
  #head (geneID2GO)
  
  return(geneID2GO_custom)
  
  
}