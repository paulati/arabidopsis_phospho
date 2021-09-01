current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
base_path <- dirname(current_working_dir)
source(file.path(base_path, 'common', 'config.R'))

library(GO.db)

load_long_format_genes <- function(data) {
  
  
  partial_results <- apply(data, MARGIN = 1,  function(row) {
    genes <- row["genes"]
    gene <- unlist(str_split(genes, pattern = ", "))
    
    gene_df <- data.frame(gene)
    gene_df["group"] <- row["group"]
    gene_df["ontology"] <- row["ontology"]
    gene_df["go_id"] <- row["go_id"]
    gene_df["genes_count"] <- row["genes_count"]
    
    result <- gene_df[, c("group", "ontology", "go_id", "gene", "genes_count")]
    return(result)
  })
  
  data_long_format <- do.call("rbind", partial_results)

  return(data_long_format)
  
}

build_perc_file <- function(et_data, group_id, group_file_name) {
  
  group_file_path <- file.path(go_terms_lists_base_path, group_file_name)
  group_data <- read.table(group_file_path, sep="\t", header= TRUE)
  group_data_uniques <- unique(group_data[, c( "group", "ontology", "go_id", "genes", "genes_count")])
  group_data_long <- load_long_format_genes(group_data_uniques)
  
  data_merge <- merge(group_data_long, et_data, by.x = "gene", by.y = "V1", all.x = TRUE, all.y = FALSE)
  
  go_db <- as.list(GOTERM)
  
  ontologies <- c("BP", "CC", "MF")
  
  for (ontology in ontologies) {
    
    result <- data.frame(group_id = character(),
                         ontology = character(),
                         go_id = character(),
                         go_term = character(),
                         genome_genes_count = numeric(),
                         et_genes_count = numeric(),
                         stringsAsFactors = FALSE)
    
    row_index <- 1
    
    mask <- data_merge$ontology == ontology
    
    data_merge_ontology <- data_merge[mask, ]
    
    unique_go_ids <- unique(data_merge_ontology$go_id)
    
    for(unique_go_id in unique_go_ids) {
      
      go_term <- go_db[[unique_go_id]]@Term
      
      #print(unique_go_id)
      mask <- data_merge_ontology$go_id == unique_go_id
      
      genes <- data_merge_ontology$gene[mask]
      et_genes <- data_merge_ontology$et_gene[mask]
      
      genes_unique <- unique(genes)
      et_genes_unique <- unique(na.omit(et_genes))
      
      genes_unique_count <- length(genes_unique)
      et_genes_unique_count <- length(et_genes_unique)
      
      row <- c(group_id, ontology, unique_go_id, go_term, genes_unique_count, et_genes_unique_count)
      result[row_index, ] <- row
      row_index <- row_index + 1
    }
    
    result$proportion <- as.numeric(result$et_genes_count) / as.numeric(result$genome_genes_count)
    result$percentage <- round(result$proportion *100, 1)
    
    file_name <- paste0("genes_prop_", group_id, "_", ontology, ".csv")
    if(! dir.exists(genes_proportion_tables_base_path)) {
      dir.create(genes_proportion_tables_base_path, showWarnings = FALSE)
    }
    file_path <- file.path(genes_proportion_tables_base_path, file_name)
    write.table(result, file_path, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
  }
  
}

build_perc_files <- function() {

  # etiolated data
  et_file_name <- sample_file_names[comparison_etiolated_vs_all]
  et_file_path <- file.path(sample_base_path, et_file_name)
  et_data <- read.table(et_file_path, header=FALSE)
  et_data$et_gene <- et_data$V1

  # et CG2 y CG3
  group_id <- "CG2_CG3"
  group_2_3_et_file_name <- "topGo_genes_by_term_groups_2_3_et.csv"
  build_perc_file(et_data, group_id, group_2_3_et_file_name)

  # et CG1
  group_id <- "CG1"
  group_1_et_file_name <- "topGo_genes_by_term_groups_1_et.csv"
  build_perc_file(et_data, group_id, group_1_et_file_name)
  
}


# build_perc_files()
