current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
base_path <- dirname(current_working_dir)
source(file.path(base_path, 'common', 'config.R'))
source(file.path(base_path, 'semantic_similarity', 'viseago_topGo_combinatory.R'))

plot_ontology <- function(ontology, show_legend) {

  alternative <- 'greater'
  #alternative <- 'less'
  
  comparisons <- c(comparison_etiolated_vs_all, comparison_etiolated_vs_union11, comparison_union11_vs_all)
  
  #cg1
  #go_terms_subset_1 <- cg1_go_terms_subset(ontology)
  
  #cg2, cg33
  go_terms_subset_2_3 <- cg2_cg3_go_terms_subset(ontology)
  #length(go_terms_subset_2_3)
  
  #groups_label <- "CG1"
  groups_label <- "CG2 CG3"
  
  group_id <- "CG2_CG3"
  
  go_term_size <- get_go_term_sizes(ontology, group_id)
  
  plot_and_title <- main(ontology, alternative, comparisons, go_terms_subset_2_3, groups_label, go_term_size,
                         include_title = FALSE, show_legend = show_legend, show_measures = FALSE, 
                         y_tick_values = seq(from = -0.5, to = 0.5, by = 0.1),
                         y_range = c(-0.5, 0.6))
  
  p <- plot_and_title$plot
  
  return(p)
  
}


fig4_body <- function() {

  plot_BP <- plot_ontology("BP", show_legend = FALSE)
  plot_BP
  
  # Sys.setenv("plotly_username"="...")
  # Sys.setenv("plotly_api_key"="...")
  #file_name <- paste0("202107-sizes-fig4-body")
  #api_create(plot_BP, filename = file_name)

}

fig4_supp <- function() {

  plot_MF <- plot_ontology("MF", show_legend = FALSE)
  # plot_MF
  
  plot_CC <- plot_ontology("CC", show_legend = FALSE)
  # plot_CC
  
  fig <- subplot(plot_MF, plot_CC,
                 nrows = 2,
                 #margin = 0.05,
                 shareX = TRUE,
                 shareY = FALSE,
                 titleX = TRUE,
                 titleY = TRUE,
                 heights = c(0.45, 0.45))
  
  
  fig
  
  # Sys.setenv("plotly_username"="...")
  # Sys.setenv("plotly_api_key"="...")
  #file_name <- paste0("202107-sizes-fig4-supp")
  #api_create(fig, filename = file_name)
  
}


# fig4_body()
# fig4_supp()

