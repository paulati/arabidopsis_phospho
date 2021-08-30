# pasarlo a yaml

# CONFIGURATION

base_path <- '/u01/home/paula/arabidopsis/github/arabidopsis_phospho'

annotations_base_path <- file.path(base_path, 'data') 

data_base_path <- file.path(base_path, 'data/results')

sample_base_path <- file.path(base_path, 'data/results')

population_base_path <- file.path(base_path, 'data/results')

significance_base_path <- file.path(base_path, 'results/topGo_may2021/alternative')

go_terms_out_base_path <- file.path(base_path, 'data/results_topgo')

go_terms_lists_base_path <- file.path(base_path, 'data', 'lists_by_groups')

genes_proportion_tables_base_path <- file.path(base_path, 'data', 'lists_by_groups')

plots_base_path <- file.path(base_path, 'data', 'results_semantic_similarity')

ss_base_path <- file.path(base_path, 'data', 'results_semantic_similarity')

annotations_custom_file_path <- file.path(base_path, "data/topgo_viseago_custom_annotations.txt")

# CONSTANTS

# figure 2
comparison_experimental_vs_all <- 'experimental_vs_all'
comparison_predOr_vs_all <- 'predOr_vs_all'
comparison_union_experimental_predOr_vs_all <- 'union_exp_predOr_vs_all'

# figure 4
comparison_etiolated_vs_all <- 'etiolated_vs_all'
comparison_etiolated_vs_union11 <- 'etiolated_vs_union11'
comparison_union11_vs_all <- 'union11_vs_all'

# FILE NAMES

go_terms_matrix_base_file_name <- "go_terms_matrix"

sample_file_names <- list()
sample_file_names[comparison_experimental_vs_all] <- 'experimental_ids.txt'
sample_file_names[comparison_predOr_vs_all] <- 'predOr_phosphat_musiteDeep.csv'
sample_file_names[comparison_union_experimental_predOr_vs_all] <- 'union_experimental_predOr_phosphat_musiteDeep_ids.csv'
sample_file_names[comparison_etiolated_vs_all] <- 'total_fosfoproteinas_unicas_cuantificadas.csv'
sample_file_names[comparison_etiolated_vs_union11] <- 'total_fosfoproteinas_unicas_cuantificadas.csv'
sample_file_names[comparison_union11_vs_all] <- 'union_experimental_predOr_phosphat_musiteDeep_ids_plus_11.csv'
 
population_file_names <- list()
population_file_names[comparison_experimental_vs_all] = 'all_ids.txt'
population_file_names[comparison_predOr_vs_all] = 'all_ids.txt'
population_file_names[comparison_union_experimental_predOr_vs_all] = 'all_ids.txt'
population_file_names[comparison_etiolated_vs_all] = 'all_ids.txt'
population_file_names[comparison_etiolated_vs_union11] = 'union_experimental_predOr_phosphat_musiteDeep_ids_plus_11.csv'
population_file_names[comparison_union11_vs_all] = 'all_ids.txt'

zip_sample_file_names <- list()
zip_sample_file_names[comparison_experimental_vs_all] <- 'ExpRS.zip'
zip_sample_file_names[comparison_predOr_vs_all] <- 'PredRS.zip'
zip_sample_file_names[comparison_union_experimental_predOr_vs_all] <- 'UnRS.zip'
zip_sample_file_names[comparison_etiolated_vs_all] <- 'total_fosfoproteinas_unicas_cuantificadas.zip'
zip_sample_file_names[comparison_etiolated_vs_union11] <- 'total_fosfoproteinas_unicas_cuantificadas.zip'
zip_sample_file_names[comparison_union11_vs_all] <- 'union_experimental_predOr_phosphat_musiteDeep_ids_plus_11.zip'

zip_population_file_names <- list()
zip_population_file_names[comparison_experimental_vs_all] = 'all_ids.zip'
zip_population_file_names[comparison_predOr_vs_all] = 'all_ids.zip'
zip_population_file_names[comparison_union_experimental_predOr_vs_all] = 'all_ids.zip'
zip_population_file_names[comparison_etiolated_vs_all] = 'all_ids.zip'
zip_population_file_names[comparison_etiolated_vs_union11] = 'union_experimental_predOr_phosphat_musiteDeep_ids_plus_11.csv'
zip_population_file_names[comparison_union11_vs_all] = 'all_ids.zip'



#under-representation




# data_file_name <- 'predOr_phosphat_musiteDeep.csv'
# data_file_name <- 'experimental_ids.txt'
# data_file_name <- 'denise_experimental_ids.txt'
#data_file_name <- 'total_fosfoproteinas_unicas_cuantificadas.csv' # etiolated
#data_file_name <- 'pred_vot2_phosphat_musiteDeep_netphos.csv'
#data_file_name <- 'nature_dataset.txt'
#data_file_name <- 'union_experimenta_predOr_phosphat_musiteDeep_ids.csv'
#data_file_name <- 'union_experimental_predOr_phosphat_musiteDeep_ids_plus_11.csv'
#data_file_name <- 'LIGHT_list.csv'







#comparison <- 'etiolated_vs_all'
# comparison <- 'denise_vs_predOr'
# comparison <- 'denise_vs_predVot2'
# comparison <- 'denise_vs_experimental'
# comparison <- 'predOr_vs_all'
# comparison <- 'experimental_vs_all'
# comparison <- 'predVot2_vs_all'
# comparison <- 'atlas_vs_all'
#comparison <- 'union_exp_predOr_vs_all'
# comparison <- 'etiolated_vs_union11'
# comparison <- 'union11_vs_all'
#comparison <- 'light_vs_all'
#comparison <- 'light_vs_union11'
#comparison <- 'dark_vs_all'
#comparison <- 'dark_vs_union11'
