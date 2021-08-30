base_path <- '/u01/home/paula/arabidopsis/github/arabidopsis_phospho'
source(file.path(base_path, 'common', 'config.R'))
source(file.path(base_path, 'topgo', 'topgo_alternative.R'))
source(file.path(base_path, 'topgo', 'annotations.R'))

library(topGO)
library(utils)

# unzip compressed files for required data
prepare_data <- function(comparisons) {
  
  for (comparison in comparisons) {
    
    sample_file_name <- sample_file_names[comparison]
    sample_file_path <- file.path(data_base_path, sample_file_name)
    
    if(! file.exists(sample_file_path)) {
      zip_file_name <- zip_sample_file_names[comparison]
      zip_file_path <- file.path(data_base_path, zip_file_name)
      sample_file_path <- unzip(zip_file_path, exdir = data_base_path)
    }
    
    population_file_name <- population_file_names[comparison]
    population_file_path <- file.path(data_base_path, population_file_name)
    
    if(! file.exists(population_file_path)) {
      zip_file_name <- zip_population_file_names[comparison]
      zip_file_path <- file.path(data_base_path, zip_file_name)
      population_file_path <- unzip(zip_file_path, exdir = data_base_path)
    }
    
  }
  
  
}

run_test_topGO <- function(topGo_data, result_out_file_path, alternative='greater') {
  
  # perform TopGO test using parentchild algorithm
  # ontology_parentchild <- topGO::runTest(
  #   topGo_data,
  #   algorithm ="parentchild",
  #   statistic = "fisher"
  # )
  
  
  if(alternative == 'greater') {
    test_stat_greater <- new("parentChild", testStatistic = GOFisher_greater, name ="parentchild fisher greater")
    ontology_parentchild <- getSigGroups(topGo_data, test_stat_greater)
  } else if(alternative == 'less') {
    test_stat_less <- new("parentChild", testStatistic = GOFisher_less, name ="parentchild fisher less")
    ontology_parentchild <- getSigGroups(topGo_data, test_stat_less)
  } else {
    ontology_parentchild <- NULL
  }
  
  terms_count <- length(score(ontology_parentchild))
  
  results_tab <- GenTable(object = topGo_data, 
                          #classicFisher = results_classic, 
                          parentchildFisher = ontology_parentchild,
                          #weight01Fisher = results_weight01,  
                          topNodes = terms_count)
  
  
  to_fix_parentchild_mask <- results_tab$parentchildFisher == '< 1e-30'
  results_tab$parentchildFisher[to_fix_parentchild_mask] <- '1e-30'
  
  adjusted_parentchildFisher <- p.adjust(results_tab$parentchildFisher, method = 'fdr', n = nrow(results_tab))
  results_tab$parentchildFisher_fdr <- adjusted_parentchildFisher
  
  
  if (! is.null(result_out_file_path)) {
    # save partial results:
    write.table(results_tab, result_out_file_path, sep='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  
  # fdr correction:
  scores <- score(ontology_parentchild)
  adj_scores <- p.adjust(scores, method = 'fdr', n = terms_count)
  score(ontology_parentchild) <- adj_scores
  
  result <- list("ontology_parentchild" = ontology_parentchild, 
                 "ontology_parentchild_df" = results_tab)
  return(result)
  
}


prepare_topGO_data <- function(ontology, comparison, description, data_file_path, background_file_path) {
  
  data <- read.table(data_file_path, stringsAsFactors = FALSE)
  data$V1 <- toupper(data$V1)
  
  gene_names_set <- as.character(data$V1)
  
  population_data <- read.table(background_file_path)
  gene_names_all <- population_data$V1
  
  #gene_names_all <- get_gene_names_all(comparison, background_base_path)
  
  # si no es factor no anda la creacion de topGOdata
  gene_list <- factor(as.integer(gene_names_all %in% gene_names_set))
  names(gene_list) <- gene_names_all
  
  #significant_terms <- get_go_terms(gene_list, geneID2GO, out_file_path, ontology)
  #significant_terms
  
  #geneID2GO <- get_geneID2GO()
  geneID2GO <- get_geneID2GO_custom()
  
  set_GOdata <- new("topGOdata", 
                    ontology = ontology,
                    allGenes = gene_list, 
                    annot = annFUN.gene2GO, 
                    gene2GO = geneID2GO, 
                    description= description)  
  
  
  return(set_GOdata)
  
}

build_go_terms_matrix <- function(significant_terms_list) {
  
  significant_terms_names <- names(significant_terms_list)
  
  all_terms <- c()
  matrix_columns <- c()
  
  for(name in significant_terms_names) {
    
    comparison_terms <- significant_terms_list[[name]]
    if(length(comparison_terms) > 0) {
      all_terms <- c(all_terms, comparison_terms)
      matrix_columns <- c(matrix_columns, name)
    }
    
  }
  all_terms_unique <- unique(all_terms)
  
  result <- matrix(0, nrow = length(all_terms_unique), ncol = length(matrix_columns))
  colnames(result) <- matrix_columns
  rownames(result) <- all_terms_unique
  
  for(name in significant_terms_names) {
    comparison_terms <- significant_terms_list[[name]]
    for(term in comparison_terms) {
      result[term, name] <- 1
    }
  }
  
  return(result)
  
  
}


go_terms_enrichment <- function(ontologies, comparisons, out_file_prefix = '') {
  
  prepare_data(comparisons)
  
  description <- "Ensembl athaliana_eg_gene plants.ensembl.org Ensembl Plants Genes 50"
  
  alternative <- "greater"
  
  significant_terms <- list()
  
  for (ontology in ontologies) {
  
    for (comparison in comparisons) {
      
      sample_file_path <- file.path(data_base_path, sample_file_names[comparison])
      population_file_path <- file.path(data_base_path, population_file_names[comparison])
      ontology_comparison <- prepare_topGO_data(ontology, comparison, description, sample_file_path, population_file_path)
      
      result_out_file_name <- paste0(out_file_prefix, '_', tolower(ontology), '_', comparison, '.csv')
      result_out_file_path <- file.path(go_terms_out_base_path, result_out_file_name)
      ontology_parentchild_comparison_obj <- run_test_topGO(ontology_comparison, result_out_file_path, alternative)
      ontology_parentchild_comparison <- ontology_parentchild_comparison$ontology_parentchild_df
    
      significant_terms_mask <- ontology_parentchild_comparison$parentchildFisher_fdr < 0.01
      significant_terms_data <- ontology_parentchild_comparison[significant_terms_mask, c('GO.ID', 'Term')]
      
      significant_terms_label <- paste(significant_terms_data$GO.ID, significant_terms_data$Term, sep='_')
      
      item_name <- paste0(tolower(ontology), '_', comparison)    
      significant_terms[[item_name]] <- significant_terms_label
      
    }
    
  }

  go_terms_matrix <- build_go_terms_matrix(significant_terms)
 
  go_terms_matrix_file_name <- paste0(out_file_prefix, '_', go_terms_matrix_base_file_name, '.csv')
  go_terms_matrix_file_path <- file.path(go_terms_out_base_path, go_terms_matrix_file_name)
  write.table(x = go_terms_matrix, file = go_terms_matrix_file_path, sep = "\t",
            col.names = TRUE, row.names = TRUE, quote = FALSE)
   
}


fig2_enrichment <- function() {

# fig 2 enrichment

  ontologies <- c('MF', 'CC', 'BP')
  comparisons <- c(comparison_experimental_vs_all, comparison_predOr_vs_all, comparison_union_experimental_predOr_vs_all)
  out_file_prefix <- "fig2"
  go_terms_enrichment(ontologies, comparisons, out_file_prefix)

}

fig3_enrichment <- function() {

  # fig 3 enrichment

  ontologies <- c('MF', 'CC', 'BP')
  comparisons <- c(comparison_etiolated_vs_all, comparison_etiolated_vs_union11, comparison_union11_vs_all)
  out_file_prefix <- "fig3"
  go_terms_enrichment(ontologies, comparisons, out_file_prefix)

}


