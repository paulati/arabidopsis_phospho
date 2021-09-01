current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
base_path <- dirname(current_working_dir)
source(file.path(base_path, 'common', 'config.R'))
source(file.path(base_path, 'topgo', 'go_terms_analysis.R'))


library(stringr)

# comparison_pattern <- g2_pattern
process_group <- function(matrix_file_path, comparisons, comparison_pattern) {
  
  matrix_data <- read.table(matrix_file_path, sep="\t")
  
  colnames(matrix_data)  
  
  ontologies <- c('BP', 'MF', 'CC')
  
  result <- list()
  
  for(ontology in ontologies) {
    
    mask <- rep(TRUE, nrow(matrix_data))
    
    for (i in 1:length(comparisons)) {
      
      comparison <- comparisons[i]
      comparison_value <- comparison_pattern[i]
      
      col_name <- paste0(tolower(ontology), "_", comparison)
      
      mask <- mask & matrix_data[col_name] == comparison_value
      
    }
    
    group_pattern_ontology <- rownames(matrix_data[mask, ])
    result[[ontology]] <- group_pattern_ontology
    
  }
  
  return(result)
   
}

process_fig_common <- function(groups_comb, groups, ontologies, comparisons) {

  for (group_comb in groups_comb) {
    
    result <- data.frame("group" = character(),
                         "ontology" = character(),
                         "comparison" = character(),
                         "go_id" = character(),
                         "genes" = character(),
                         "genes_count" = numeric(),
                         stringsAsFactors = FALSE)
    
    row_index <- 1
    
    group_ids <- group_comb[["group_ids"]]
    
    for(group_id in group_ids) {
      
      for(ontology in ontologies) {
        
        group <- groups[[group_id]][[ontology]]
        
        go_ids_named <- sapply(group, function(x) unlist(strsplit(x,split = "_"))[1])
        go_ids_values <- unname(go_ids_named)
        go_ids <- data.frame(go_ids_values)
        
        
        for(comparison in comparisons) {
          
          description <- ''
          data_file_path <- file.path(sample_base_path, sample_file_names[comparison])
          background_file_path <- file.path(population_base_path, population_file_names[comparison])
          
          set_GOdata <- prepare_topGO_data(ontology, comparison, description, data_file_path, background_file_path)
          
          all_terms <- set_GOdata@graph@nodes
          
          go_terms_in_all <- sapply(go_ids_values, function(x) x %in% all_terms)
          if(length(go_terms_in_all) > 0) {
            print(sum(go_terms_in_all) == length(go_ids_values))  
          } else {
            print(0 == length(go_ids_values))  
          }
          
          for(go_id_value in go_ids_values) {
            
            genes_go_id <- unlist(genesInTerm(set_GOdata, go_id_value))
            genes_count <- length(genes_go_id)
            genes_go_id_str <- paste(genes_go_id, collapse = ", ")
            #genes_go_id_str_clean <- gsub(pattern = "\"", replacement = "", x = genes_go_id_str)
            #genes_go_id_str_clean <- gsub(pattern = "\n", replacement = "", x = genes_go_id_str_clean)
            
            result[row_index, ] <- c(group_id, ontology, comparison, go_id_value, genes_go_id_str, genes_count)
            row_index <- row_index + 1
            
          }
          
          
        }
        
      }
    }    
    
    out_file_name <- group_comb[["out_file_name"]]
    out_file_path <- file.path(base_path, 'data', 'lists_by_groups', out_file_name)
    write.table(result, out_file_path, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)    
    
  }  
    
}

process_fig_2_groups <- function(group_ids_to_include) {
  
  #matrix_file_path <- file.path(base_path, "results/topGo_may2021/alternative/plots_fig2", "go_terms_matrix.csv")
  matrix_file_path <- file.path(base_path, "data", "results_topgo", "fig2_go_terms_matrix.csv")  
  
  comparisons <- c(comparison_predOr_vs_all, comparison_experimental_vs_all, comparison_union_experimental_predOr_vs_all)
  
  g1_pattern <- c(1, 1, 1)
  g2_pattern <- c(0, 1, 0)
  g3_pattern <- c(1, 0, 1)
  g4_pattern <- c(1, 0, 0)
  g5_pattern <- c(0, 1, 1)
  g6_pattern <- c(1, 1, 0)
  g7_pattern <- c(0, 0, 1)
  
  g_patterns <- list(g1_pattern, g2_pattern, g3_pattern,
                  g4_pattern, g5_pattern, g6_pattern,
                  g7_pattern)
  
  
  # groups <- list("2" = process_group(matrix_file_path, comparisons, g2_pattern),
  #                "3"  = process_group(matrix_file_path, comparisons, g3_pattern),
  #                "4"  = process_group(matrix_file_path, comparisons, g4_pattern),
  #                "5"  = process_group(matrix_file_path, comparisons, g5_pattern)
  # )
  

  groups <- list()
  
  for(group_id in group_ids_to_include) {
    
    key <- as.character(group_id)
    
    value <- process_group(matrix_file_path, comparisons, g_patterns[[group_id]])
    
    groups[[key]] <- value
    
  }

  return(groups)
  
  
}

process_fig_2 <- function() {
  
  group_ids_to_include <- c(2, 3, 4, 5)
  
  groups <- process_fig_2_groups(group_ids_to_include)

  groups_comb <- list(
    "comb_1" = list(
      "group_ids" = c("2", "3"),
      "out_file_name" = "topGo_genes_by_term_groups_2_3.csv"
    ), 
    "comb_2" = list(
      "group_ids" = c("4", "5"),
      "out_file_name" = paste0("topGo_genes_by_term_groups_4_5.csv")
    )
  )
  
  ontologies <- c("BP", "MF", "CC")
  
  comparisons <- c(comparison_predOr_vs_all, comparison_experimental_vs_all, comparison_union_experimental_predOr_vs_all)

  process_fig_common(groups_comb, groups, ontologies, comparisons)

  
}

process_fig_3_groups <- function(group_ids_to_include) {
  
  #matrix_file_path <- file.path(base_path, "results/topGo_may2021/alternative/plots_fig4", "go_terms_matrix_et.csv")
  matrix_file_path <- file.path(base_path, "data", "results_topgo", "fig3_go_terms_matrix.csv")
  
  comparisons <- c(comparison_etiolated_vs_all, comparison_etiolated_vs_union11, comparison_union11_vs_all)
  
  g1_pattern <- c(0, 0, 1)
  g2_pattern <- c(1, 1, 0)
  g3_pattern <- c(1, 1, 1)
  g4_pattern <- c(1, 0, 1)
  g5_pattern <- c(1, 0, 0)
  g6_pattern <- c(0, 1, 0)

  g_patterns <- list(g1_pattern, g2_pattern, g3_pattern,
                     g4_pattern, g5_pattern, g6_pattern)
  
  
  # groups <- list("2" = process_group(matrix_file_path, comparisons, g2_pattern),
  #                "3"  = process_group(matrix_file_path, comparisons, g3_pattern),
  #                "4"  = process_group(matrix_file_path, comparisons, g4_pattern),
  #                "5"  = process_group(matrix_file_path, comparisons, g5_pattern)
  # )
  
  
  groups <- list()
  
  for(group_id in group_ids_to_include) {
    
    key <- as.character(group_id)
    
    value <- process_group(matrix_file_path, comparisons, g_patterns[[group_id]])
    
    groups[[key]] <- value
    
  }
  
  return(groups)
  
}  

process_fig_3 <- function() {
  
  # groups_ids_to_include <- c(2, 3, 4, 5)
  # 
  # groups <- process_fig_2_groups(groups_ids_to_include)
  # 
  # groups_comb <- list(
  #   "comb_1" = list(
  #     "group_ids" = c("2", "3"),
  #     "out_file_name" = "topGo_genes_by_term_groups_2_3.csv"
  #   ), 
  #   "comb_2" = list(
  #     "group_ids" = c("4", "5"),
  #     "out_file_name" = paste0("topGo_genes_by_term_groups_4_5.csv")
  #   )
  # )
  # 
  # ontologies <- c("BP", "MF", "CC")
  # 
  # comparisons <- c(comparison_predOr_vs_all, comparison_experimental_vs_all, comparison_union_experimental_predOr_vs_all)
  # 
  # process_fig_common(groups_comb, groups, ontologies, comparisons)
  
  group_ids_to_include <- c(1:3)
  
  groups <- process_fig_3_groups(group_ids_to_include)
  
  groups_comb <- list(
    "comb_1" = list(
      "group_ids" = c("1"),
      "out_file_name" = "topGo_genes_by_term_groups_1_et.csv"
    ), 
    "comb_2" = list(
      "group_ids" = c("2", "3"),
      "out_file_name" = paste0("topGo_genes_by_term_groups_2_3_et.csv")
    )
  )
  
  ontologies <- c("BP", "MF", "CC")
  
  comparisons <- c(comparison_etiolated_vs_all, comparison_etiolated_vs_union11,  comparison_union11_vs_all)
  
  process_fig_common(groups_comb, groups, ontologies, comparisons)
  
}

# grupos que no son etiolados
process_lists_groups <- function() {
  
  group_ids_to_include <- c(1:7)
  
  groups <- process_fig_2_groups(group_ids_to_include)
  
  groups_comb <- list(
    "comb_1" = list(
      "group_ids" = c("1"),
      "out_file_name" = "topGo_genes_by_term_groups_1.csv"
    ), 
    "comb_2" = list(
      "group_ids" = c("2"),
      "out_file_name" = paste0("topGo_genes_by_term_groups_2.csv")
    ),
    "comb_3" = list(
      "group_ids" = c("3"),
      "out_file_name" = "topGo_genes_by_term_groups_3.csv"
    ), 
    "comb_4" = list(
      "group_ids" = c("4"),
      "out_file_name" = "topGo_genes_by_term_groups_4.csv"
    ), 
    "comb_5" = list(
      "group_ids" = c("5"),
      "out_file_name" = "topGo_genes_by_term_groups_5.csv"
    ), 
    "comb_6" = list(
      "group_ids" = c("6"),
      "out_file_name" = "topGo_genes_by_term_groups_6.csv"
    ), 
    "comb_7" = list(
      "group_ids" = c("7"),
      "out_file_name" = "topGo_genes_by_term_groups_7.csv"
    ) 
  )
  
  ontologies <- c("BP", "MF", "CC")
  
  comparisons <- c(comparison_predOr_vs_all, comparison_experimental_vs_all, comparison_union_experimental_predOr_vs_all)
  
  process_fig_common(groups_comb, groups, ontologies, comparisons)
    
  
}

# grupos que son etiolados
process_et_groups <- function() {
  
  group_ids_to_include <- c(1:6)
  
  groups <- process_fig_3_groups(group_ids_to_include)
  
  groups_comb <- list(
    "comb_1" = list(
      "group_ids" = c("1"),
      "out_file_name" = "topGo_genes_by_term_groups_1_et.csv"
    ), 
    "comb_2" = list(
      "group_ids" = c("2"),
      "out_file_name" = paste0("topGo_genes_by_term_groups_2_et.csv")
    ),
    "comb_3" = list(
      "group_ids" = c("3"),
      "out_file_name" = "topGo_genes_by_term_groups_3_et.csv"
    ), 
    "comb_4" = list(
      "group_ids" = c("4"),
      "out_file_name" = "topGo_genes_by_term_groups_4_et.csv"
    ), 
    "comb_5" = list(
      "group_ids" = c("5"),
      "out_file_name" = "topGo_genes_by_term_groups_5_et.csv"
    ), 
    "comb_6" = list(
      "group_ids" = c("6"),
      "out_file_name" = "topGo_genes_by_term_groups_6_et.csv"
    ) 
  )
  
  ontologies <- c("BP", "MF", "CC")
  
  comparisons <- c(comparison_etiolated_vs_all, comparison_etiolated_vs_union11,  comparison_union11_vs_all)
  
  process_fig_common(groups_comb, groups, ontologies, comparisons)
  
  
  

}



check_quantities <- function() {

  groups <- data.frame("BP"= integer(7),
                       "MF"= integer(7),
                       "CC"= integer(7))


  groups_et <- data.frame("BP"= integer(6),
                       "MF"= integer(6),
                       "CC"= integer(6))

  ontologies <- c("BP", "MF", "CC")


  for(group_id in c(1:nrow(groups))) {


    file_name <- paste0("topGo_genes_by_term_groups_", as.character(group_id), ".csv")
    file_path <- file.path(base_path, 'data', 'lists_by_groups', file_name)
    check_quantities
    data <- read.table(file_path, sep="\t", header = TRUE, stringsAsFactors = FALSE)

    for (ontology in ontologies) {

      mask <- data$ontology == ontology
      go_ids <- unique(data[mask, "go_id"])
      go_ids_count <- length(go_ids)

      groups[group_id, ontology] <- go_ids_count

    }


  }

  for(group_id in c(1:nrow(groups_et))) {


    file_name <- paste0("topGo_genes_by_term_groups_", as.character(group_id), "_et.csv")
    file_path <- file.path(base_path, 'data', 'lists_by_groups', file_name)

    data <- read.table(file_path, sep="\t", header = TRUE, stringsAsFactors = FALSE)

    for (ontology in ontologies) {

      mask <- data$ontology == ontology
      go_ids <- unique(data[mask, "go_id"])
      go_ids_count <- length(go_ids)

      groups_et[group_id, ontology] <- go_ids_count

    }


  }

  return(list("groups"=groups, "groups_et"=groups_et))
  
}












