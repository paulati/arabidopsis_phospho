base_path <- '/u01/home/paula/arabidopsis/github/arabidopsis_phospho'
source(file.path(base_path, 'common', 'config.R'))
source(file.path(base_path, 'topgo', 'go_terms_analysis.R'))
source(file.path(base_path, 'semantic_similarity', 'semantic_similarity_MDSplot.R'))

library(topGO)
library(stringr)
library(GOSemSim)


prepare_viseago_data <- function(ontology, comparisons, alternative, in_base_path, out_base_path) {
  
  description <- "Ensembl athaliana_eg_gene plants.ensembl.org Ensembl Plants Genes 50"

  custom_annotations <- ViSEAGO::Custom2GO(annotations_custom_file_path)
  genes_to_go_custom <- ViSEAGO::annotate(
    id="athaliana_eg_gene",
    object = custom_annotations
  )
  
  result <- list(
    "genes_to_go" = genes_to_go_custom
    
  )

  for (comparison in comparisons) {
    
    sample_file_path <- file.path(data_base_path, sample_file_names[comparison])
    population_file_path <- file.path(data_base_path, population_file_names[comparison])
    ontology_comparison <- prepare_topGO_data(ontology, comparison, description, sample_file_path, population_file_path)
  
    result_out_file_path <- file.path(out_base_path, paste0(ontology, '_', comparison, '_mix.csv'))
    ontology_parentchild_comparison_obj <- run_test_topGO(ontology_comparison, result_out_file_path, alternative)
    ontology_parentchild_comparison <- ontology_parentchild_comparison_obj$ontology_parentchild
    
    key_1 <- paste0('ontology_', comparison)
    result[[key_1]] <- ontology_comparison
    
    key_2 <- paste0('ontology_parentchild_', comparison)
    result[[key_2]] <- ontology_parentchild_comparison    
  }
  
  return(result)
  
}

compute_ss_distances_terms <- function(object, distance, valid_go_terms) {
  
  # check object
  if(!is(object,"GO_SS") & !is(object,"GO_clusters")){
    
    stop("object must be a GO_SS or GO_clusters class from ViSEAGO::build_GO_SS() or ViSEAGO::GOterms_heatmap(), respectively")
  }
  
  if(is(object,"GO_SS")){
    
    # extract GO terms
    raw_terms<-slot(slot(object,"enrich_GOs"),"data")$GO.ID
    
    # revisar aca que devuelve esto:
    #term=slot(slot(object,"enrich_GOs"),"data")$term
    
    valid_terms_mask <- unlist(sapply(raw_terms, function(t) t %in% valid_go_terms))
    
    Terms <- raw_terms[valid_terms_mask]
    
    # for each distance
    for (x in distance){
      
      print(x)
      
      # only wang is valid because IC is computed for the complete set and this value is different from that for go terms subset
      # check measure
      x=match.arg(x,c("Wang"))
      
      # convert GO similarity matrix to distance
      values=list(as.dist(1-mgoSim(Terms,Terms,semData=object,measure=x,combine=NULL)))
      
      # add distance name to values
      names(values)=x
      
      # return values
      slot(object,"terms_dist")<-c(slot(object,"terms_dist"),values)
    }
  }
  
  return(object)
  
}


prepare_ss_viseago_data <- function(viseago_data, comparisons, go_terms_subset=NULL) {

  # reflection:
  # https://stackoverflow.com/questions/6702202/does-the-r-programming-language-have-reflection

  
  genes_to_go <- viseago_data$genes_to_go

  param_list <- list()
  for (comparison in comparisons) {
    
    ontology_comparison_variable_name <- paste0('ontology_', comparison)
    ontology_comparison_variable_value <- viseago_data[[ontology_comparison_variable_name]]
    assign(ontology_comparison_variable_name, ontology_comparison_variable_value)
           
    ontology_parentchild_comparison_variable_name <- paste0('ontology_parentchild_', comparison)
    ontology_parentchild_comparison_variable_value <- viseago_data[[ontology_parentchild_comparison_variable_name]]
    assign(ontology_parentchild_comparison_variable_name, ontology_parentchild_comparison_variable_value)
    
    comparison_element <- c(ontology_comparison_variable_name, ontology_parentchild_comparison_variable_name)
    
    param_list[[comparison]] <- comparison_element
    
  }
  
  
  ontology_sResults <- ViSEAGO::merge_enrich_terms(
    Input = param_list,
    cutoff = 0.01,
    envir = environment()
  )
  
  myGOs <- ViSEAGO::build_GO_SS(
    gene2GO = genes_to_go,
    enrich_GO_terms = ontology_sResults
  )
  
  if(is.null(go_terms_subset)) {
    
    # compute all available Semantic Similarity (SS) measures
    myGOs_distances <- ViSEAGO::compute_SS_distances(
      myGOs,
      #distance="Resnik"
      distance = c("Wang", "Resnik", "Lin","Rel","Jiang")
    )
    
    result <- myGOs_distances
    
  } else {
    
      # alternative: only Wang distance makes sense with go terms subset
      myGOs_subset_distances <- compute_ss_distances_terms(
        myGOs,
        distance = c("Wang"), 
        valid_go_terms = go_terms_subset
      )
    
    result <- myGOs_subset_distances
    
  }
  
  return(result)
  
}

main <- function(ontology, alternative, comparisons, 
                 go_terms_subset=NULL, groups_label=NULL, go_term_size=NULL,
                 include_title = TRUE, show_legend = TRUE,
                 show_measures = TRUE,
                 x_tick_values=NULL, y_tick_values= NULL,
                 x_range=NULL, y_range=NULL) {

  out_base_path <- file.path(ss_base_path, tolower(ontology))
  dir.create(out_base_path, showWarnings = FALSE)
  
  viseago_data <- prepare_viseago_data(ontology, comparisons, alternative, data_base_path, out_base_path)
  
  myGOs_distances <- prepare_ss_viseago_data(viseago_data, comparisons, go_terms_subset)
  
  if(include_title) {
    
    plot_title_tail <- paste(comparisons, collapse = ", ")
    plot_title <- paste0("MDS. ", str_to_title(alternative), ". ")
    
    if(! is.null(groups_label)) {
      plot_title <- paste0(plot_title, toupper(groups_label), '. ')
    }
    
    plot_title <- paste0(plot_title, toupper(ontology), " ",  plot_title_tail)
  
  } else {
    
    plot_title <- ''
  
    }
  
  p <- MDSplot_custom(myGOs_distances, "GOterms",
                     title = plot_title,
                     file=NULL,
                     go_term_size, 
                     show_legend = show_legend,
                     show_measures = show_measures,
                     x_tick_values, y_tick_values,
                     x_range, y_range)
  
  
  result <- list("plot" = p, "title" = plot_title)
  
  
  return (result)



}

cg1_go_terms_subset <- function(ontology) {

  group_id <- 1
  
  file_name <- paste0("topGo_genes_by_term_groups_", as.character(group_id), "_et.csv")  
  file_path <- file.path(go_terms_lists_base_path, file_name)
  
  data <- read.table(file_path, sep="\t", header = TRUE, stringsAsFactors = FALSE)
  
  mask <- data$ontology == ontology
  
  result <- unique(data[mask, 'go_id'])
  
  return(result)
}

cg2_cg3_go_terms_subset <- function(ontology) {
  
  group_ids <- c(2, 3)
  
  result <- c()
  
  for(group_id in group_ids) {
  
    file_name <- paste0("topGo_genes_by_term_groups_", as.character(group_id), "_et.csv")  
    file_path <- file.path(go_terms_lists_base_path, file_name)
  
    data <- read.table(file_path, sep="\t", header = TRUE, stringsAsFactors = FALSE)
  
    mask <- data$ontology == ontology
    
    group_result <- unique(data[mask, 'go_id'])
    
    result <- c(result, group_result)
  
  }
  
  return(result)
}

get_go_term_sizes <- function(ontology, group_id) {
  
  file_name <- paste0("genes_prop_", group_id, "_", ontology, ".csv")
  file_path <- file.path(genes_proportion_tables_base_path, file_name)
  data <- read.table(file_path, sep="\t", header = TRUE)
  result <- data[, c("go_id", "percentage")]
  colnames(result) <- c("GO.ID", "go_term_size")
  return(result)
}
