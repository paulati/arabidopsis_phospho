# We build a custom annotations file, function build_custom_annotations, that can be used both by 
# * topGo enrichment analysis
# * Viseago semantic similarity analysis
#---------------------------------------------------

base_path <- '/u01/home/paula/arabidopsis/github/arabidopsis_phospho'
source(file.path(base_path, 'common', 'config.R'))

library(ViSEAGO)
library(biomaRt)

#--------------------------------------------------
# ViSEAGO annotations

convert_to_data_frame <- function(gene_id, data_list) {
  
  #print(gene_id)
  
  terms <- data_list[[gene_id]]
  
  gene_term_df <- data.frame("gene_id" = gene_id,
                             "evidence_code" = names(terms),
                             "go_term_id" = terms,
                             stringsAsFactors = FALSE)  
  
  return(gene_term_df)
}


build_ViSEAGO_annotations <- function(){

  ensembl <- ViSEAGO::Ensembl2GO(biomart = "plants_mart", host = 'plants.ensembl.org', version = NULL)
  ensembl_id <- "athaliana_eg_gene"
  
  
  genes_to_go <- ViSEAGO::annotate(
    ensembl_id,
    ensembl
  )
  
  gene_term_mf_df_list <- lapply(names(genes_to_go@MF), function(x) convert_to_data_frame(x, genes_to_go@MF))
  #length(gene_term_mf_df_list)
  
  gene_terms_mf_df <- rbindlist(gene_term_mf_df_list, use.names = TRUE)
  gene_terms_mf_df$ontology <- "MF"
  
  gene_term_bp_df_list <- lapply(names(genes_to_go@BP), function(x) convert_to_data_frame(x, genes_to_go@BP))
  #length(gene_term_bp_df_list)
  
  gene_terms_bp_df <- rbindlist(gene_term_bp_df_list, use.names = TRUE)
  gene_terms_bp_df$ontology <- "BP"
  
  gene_term_cc_df_list <- lapply(names(genes_to_go@CC), function(x) convert_to_data_frame(x, genes_to_go@CC))
  #length(gene_term_cc_df_list)
  
  gene_terms_cc_df <- rbindlist(gene_term_cc_df_list, use.names = TRUE)
  gene_terms_cc_df$ontology <- "CC"
  
  gene_terms_df <- rbind(gene_terms_mf_df, gene_terms_bp_df, gene_terms_cc_df)
  #length(unique(gene_terms_df$gene_id))
  
  out_file_path <- file.path(annotations_base_path, "ViSEAGO_Ensembl2GO_annotations.csv")
  write.table(gene_terms_df, out_file_path, sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)                           

}
#--------------------------------------------------

compare_biomart_viseago_annotations <- function() {

  # compare annotations created by biomart and filtered by from ViSEAGO:

  annotations_viseago_file_path <- file.path(annotations_base_path, "ViSEAGO_Ensembl2GO_annotations.csv")
  annotations_viseago <- read.table(annotations_viseago_file_path, sep = "\t", header = TRUE)
  
  mart <- biomaRt::useMart(biomart = "plants_mart",
                           dataset = "athaliana_eg_gene",
                           host = 'plants.ensembl.org')
  
  # Get ensembl gene ids and GO terms
  annot_attributes <- c( "ensembl_gene_id", "tair_symbol", "go_id", "go_linkage_type", "namespace_1003")
  GTOGO <- biomaRt::getBM(attributes = annot_attributes, 
                          filters = "with_go", 
                          values = TRUE,
                          mart = mart)
  
    
  #mf:
  annotations_viseago_mf <- unique(annotations_viseago[annotations_viseago$ontology == 'MF', c('gene_id', 'go_term_id')])
  nrow(annotations_viseago_mf)
  GTOGO_mf <- unique(GTOGO[GTOGO$namespace_1003 == 'molecular_function', c('ensembl_gene_id', 'go_id')])
  nrow(GTOGO_mf)
  annotations_mf <- merge(annotations_viseago_mf, GTOGO_mf,
                         by.x = c('gene_id', 'go_term_id'),
                         by.y= c('ensembl_gene_id', 'go_id'),
                         all.x = FALSE,
                         all.y = FALSE)
  nrow(annotations_mf)

  #cc:
  annotations_viseago_cc <- unique(annotations_viseago[annotations_viseago$ontology == 'CC', c('gene_id', 'go_term_id')])
  nrow(annotations_viseago_cc)
  GTOGO_cc <- unique(GTOGO[GTOGO$namespace_1003 == 'cellular_component', c('ensembl_gene_id', 'go_id')])
  nrow(GTOGO_cc)
  annotations_cc <- merge(annotations_viseago_cc, GTOGO_cc,
                         by.x = c('gene_id', 'go_term_id'),
                         by.y= c('ensembl_gene_id', 'go_id'),
                         all.x = FALSE,
                         all.y = FALSE)
  nrow(annotations_cc)

  #bp:
  annotations_viseago_bp <- unique(annotations_viseago[annotations_viseago$ontology == 'BP', c('gene_id', 'go_term_id')])
  nrow(annotations_viseago_bp)
  GTOGO_bp <- unique(GTOGO[GTOGO$namespace_1003 == 'biological_process', c('ensembl_gene_id', 'go_id')])
  nrow(GTOGO_bp)
  annotations_bp <- merge(annotations_viseago_bp, GTOGO_bp,
                         by.x = c('gene_id', 'go_term_id'),
                         by.y= c('ensembl_gene_id', 'go_id'),
                         all.x = FALSE,
                         all.y = FALSE)
  nrow(annotations_bp)

  # both versions are equivalent, we save GTOGO as custom annotations because these have more info (fields)

}
  
#--------------------------------------------------

# biomart annotations

build_custom_annotations <- function() {

  mart <- biomaRt::useMart(biomart = "plants_mart",
                           dataset = "athaliana_eg_gene",
                           host = 'plants.ensembl.org')
  
  # Get ensembl gene ids and GO terms
  annot_attributes <- c( "ensembl_gene_id", "tair_symbol", "go_id", "go_linkage_type", "namespace_1003")
  GTOGO <- biomaRt::getBM(attributes = annot_attributes, 
                          filters = "with_go", 
                          values = TRUE,
                          mart = mart)
  
  # fields required by ViSEAGO:
  #"ensembl_gene_id","go_id","go_linkage_type","namespace_1003"),
  
  GTOGO <- GTOGO[GTOGO$go_id != '',]
  
  # ViSEAGO requires  taxid = athaliana_eg_gene
  GTOGO$taxid <- 'athaliana_eg_gene'
 
  
  # remove terms not associated to any ontology:
  ontology_mask <- GTOGO$namespace_1003 != ''
  GTOGO_with_ontology <- GTOGO[ontology_mask, ]
  
  # remove terms not included in GO.db (obsolet):
  all_terms <- keys(GO.db)
  obsolete_terms_mask <- !(GTOGO_with_ontology$go_id %in% all_terms)
  sum(obsolete_terms_mask)
  
  GTOGO_with_ontology_not_obsolete <- GTOGO_with_ontology[! obsolete_terms_mask, ]
  
  # rename categories:
  mf_mask <- GTOGO_with_ontology_not_obsolete$namespace_1003 == 'molecular_function'
  cc_mask <- GTOGO_with_ontology_not_obsolete$namespace_1003 == 'cellular_component'
  bp_mask <- GTOGO_with_ontology_not_obsolete$namespace_1003 == 'biological_process'
  GTOGO_with_ontology_not_obsolete[mf_mask, 'namespace_1003'] <- 'MF'
  GTOGO_with_ontology_not_obsolete[cc_mask, 'namespace_1003'] <- 'CC'
  GTOGO_with_ontology_not_obsolete[bp_mask, 'namespace_1003'] <- 'BP'
  
  
  # viseago: custom annotation file columns names required: "taxid","gene_id","gene_symbol","GOID","evidence"
  colnames(GTOGO_with_ontology_not_obsolete)<-c("gene_id", "gene_symbol", "GOID", "evidence", "category", "taxid")
  
  #remove category column because viseago adds this column and it raises an error if already exists:
  category_column_index <- 5
  GTOGO_with_ontology_not_obsolete <- GTOGO_with_ontology_not_obsolete[, -category_column_index]
  
  out_file_path <- file.path(annotations_base_path, 'topgo_viseago_custom_annotations.txt')
  write.table(GTOGO_with_ontology_not_obsolete, out_file_path, sep="\t", 
              col.names = TRUE, row.names = FALSE,
              quote = FALSE)

}


#--------------------------------------------------

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