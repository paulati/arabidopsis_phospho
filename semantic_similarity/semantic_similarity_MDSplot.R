library(RColorBrewer)
library(colorspace)


build_combinatory <- function(comparisons) {
  
  max_m <- length(comparisons)
  
  combinatory_nums <- sapply(c(1:max_m), function(x) choose(length(comparisons), x))
  combinations_count <- sum(combinatory_nums)
  
  result <- matrix(0,
                   ncol= length(comparisons) + 1, 
                   nrow = combinations_count,
                   dimnames = list(NULL, c(comparisons, 'sum')))
  
  result_row_index <- 1
  
  
  for (m in c(1:max_m)) {
    
    combinations <- combn(comparisons, m, simplify = TRUE)
    
    
    # las columnas definen los elementos
    # las filas son cada uno de las comparaciones que forman parte de ese elemento
    
    for (i in c(1:ncol(combinations))) {
      
      
      for (required_1 in combinations[, i]) {
        
        result[result_row_index, required_1] <- 1
        
      }
     
      result[result_row_index, 'sum'] <- length(combinations[, i])
      
      result_row_index <-  result_row_index + 1
    }
    
  }
  
  return(result)
   
}


# combinatory_data_matrix <- combinatory_comparisons
load_color_scale <- function(combinatory_data_matrix) {
  
  result <- as.data.frame(combinatory_data_matrix)
  
  columns_count <- dim(combinatory_data_matrix)[2]
  comparisons_columns_indexes <- c(1:(columns_count-1))
  comparisons <- colnames(combinatory_data_matrix)[comparisons_columns_indexes]
  
  basic_colors_count <- sum(combinatory_data_matrix[, 'sum'] == 1)
  
  #display.brewer.all()
  #palette()
  
  #name <- "Set1"
  # 1. Return the hexadecimal color specification 
  #basic_colors <- brewer.pal(basic_colors_count, name)
  if(basic_colors_count == 2) {
  
    basic_colors <- c("#0000ff", "#ff0000")
    
  } else if(basic_colors_count == 3) {
    # basic_colors <- c("#4d4dff", "#f64c4d", "#79e9de")  
    basic_colors <- c("#0000ff", "#ff0000", "#4aedea")
    
    
  } else {
    name <- "Set1"
    # 1. Return the hexadecimal color specification 
    basic_colors <- brewer.pal(basic_colors_count, name)
  }
  
  names(basic_colors) <- comparisons
  
  #('predOr_vs_all', 'exp_vs_all', 'union_vs_all')
  
  rows_colors_to_combine <- list()
  
  for (row_index in c(1:nrow(combinatory_data_matrix))) {
    
    row_colors <- c()
    
    for(column_index in comparisons_columns_indexes) {
      
      comparison <- comparisons_columns_indexes[column_index]
      
      comparison_color <- basic_colors[column_index]
    
      if(combinatory_data_matrix[row_index, column_index] == 1) {
        
        row_colors <- c(row_colors, comparison_color)
        
      }
     
        
      
        
    }
    
    rows_colors_to_combine[[row_index]] <- row_colors
    
  }
  
  
  rows_mixed_colors <- c()
  
  for (row_colors_to_combine in rows_colors_to_combine) {
    
    mixed_color <- NULL
    
    for(row_color in row_colors_to_combine) {

      #hex(hex2RGB(row_color))
      
      #row_color_RGB_vector <- c(col2rgb(row_color))
      row_color_RGB_vector <- hex2RGB(row_color)
      row_color_RGB <- RGB(row_color_RGB_vector@coords[1], 
                           row_color_RGB_vector@coords[2], 
                           row_color_RGB_vector@coords[3])
      
      if(is.null(mixed_color)) {
        
        mixed_color <- row_color_RGB
        
      } else {

        alpha <- 0.5
        mixed_color <- mixcolor(alpha, mixed_color, row_color_RGB)
        
      }
      
      
    }
    
    mixed_color_hex <- hex(mixed_color, fixup = TRUE)
    
    # mixed_color@coords
    # mixed_color_int <- sapply(mixed_color@coords, function(x) as.integer(x))
    # mixed_color_int_RGB <- RGB(mixed_color_int[1], mixed_color_int[2], mixed_color_int[3])
    # 
    # hex(sapply(mixed_color@coords, function(x) as.integer(x)), gamma = NULL, fixup = TRUE)
    # 
    # #mixed_color_hex <- hex(mixed_color_int_RGB, gamma = NULL, fixup = FALSE)
    # #rows_mixed_colors <- c(rows_mixed_colors, mixed_color_hex)
    rows_mixed_colors <- c(rows_mixed_colors, mixed_color_hex)
    
  }
  
  result[, 'color'] <- rows_mixed_colors
  
  return(result)
}


load_legends <- function(combinatory_data_matrix) {
  
  result <- as.data.frame(combinatory_data_matrix)
  
  columns_count <- dim(combinatory_data_matrix)[2]
  
  # PAULA TODO: hacerlo dinamico excluyendo por nombre color y sum:
  comparisons_columns_indexes <- c(1:(columns_count-2))
  comparisons <- colnames(combinatory_data_matrix)[comparisons_columns_indexes]

  
  rows_legends_to_combine <- list()
  
  for (row_index in c(1:nrow(combinatory_data_matrix))) {
    
    row_legends <- c()
    
    for(column_index in comparisons_columns_indexes) {
      
      comparison <- comparisons_columns_indexes[column_index]
      
      comparison_legend <- comparisons[column_index]
      
      if(combinatory_data_matrix[row_index, column_index] == 1) {
        
        row_legends <- c(row_legends, comparison_legend)
        
      }
    
    }
    
    rows_legends_to_combine[[row_index]] <- row_legends
    
  }
  
  
  rows_mixed_legends <- c()
  
  for (row_legends_to_combine in rows_legends_to_combine) {
    
    mixed_legend <- ''
    
    for(row_legend in row_legends_to_combine) {
      
      if(nchar(mixed_legend) == 0) {
        
        mixed_legend <- row_legend
        
      } else {
      
        mixed_legend <- paste0(mixed_legend, ' AND ', row_legend)
        
      }
      
    }
    
    rows_mixed_legends <- c(rows_mixed_legends, mixed_legend)
    
  }
  
  result[, 'legend'] <- rows_mixed_legends
  
  return(result)  

    
  
}

get_terms_categories <- function(object) {

  # paula TODO revisar
  threshold <- object@enrich_GOs@cutoff[[1]][1]
    
  comparisons <- names(object@enrich_GOs@summary)

  categories <- data.frame("GO.ID" = object@enrich_GOs@data$GO.ID,
                           stringsAsFactors = FALSE)    
  
  for (comparison in comparisons) {
  
    category_key <- paste0(comparison, ".pvalue")
    category_terms_mask <- object@enrich_GOs@data[[category_key]] < threshold
    categories[comparison] <- as.integer(category_terms_mask)
  
  }
  

  # # load category:
  # cat_1_terms_mask <- object@enrich_GOs@data[["predOr_vs_all.pvalue"]] < threshold
  # #cat1_terms <- object@enrich_GOs@data$GO.ID[cat_1_terms_mask]
  # 
  # cat_2_terms_mask <- object@enrich_GOs@data[["exp_vs_all.pvalue"]] < threshold
  # #cat2_terms <- object@enrich_GOs@data$GO.ID[cat_2_terms_mask]
  # 
  # #cat_3_terms_mask <- object@enrich_GOs@data[["denise1000_vs_all.pvalue"]] < threshold
  # 
  # cat_4_terms_mask <- object@enrich_GOs@data[["union_vs_all.pvalue"]] < threshold
  # 
  # categories <- data.frame("GO.ID" = object@enrich_GOs@data$GO.ID,
  #                          "predOr_vs_all" = as.integer(cat_1_terms_mask), 
  #                          "exp_vs_all" = as.integer(cat_2_terms_mask),
  #                          #"denise1000_vs_all" = as.integer(cat_3_terms_mask),
  #                          "union_vs_all" = as.integer(cat_4_terms_mask),
  #                          stringsAsFactors = FALSE)
  
  columns_count <- length(colnames(categories))
  #categories$sum <- categories$predOr_vs_all + categories$exp_vs_all + categories$denise1000_vs_all + categories$union_vs_all
  
  tmp_sum <- rep(0, nrow(categories))
  for (i in c(2:columns_count)) {
    tmp_sum <- tmp_sum + categories[, i]
  }
  categories$sum <- tmp_sum
  
  
  combinatory_comparisons <- build_combinatory(comparisons)
  
  combinatory_comparisons_with_colors <- load_color_scale(combinatory_comparisons)
  
  combinatory_comparisons_with_colors_and_legends <- load_legends(combinatory_comparisons_with_colors)
  
  categories['category'] <- rep(NA, nrow(categories))
  categories['legend'] <- rep(NA, nrow(categories))
  
  
  #agrego a categories los colores y leyendas
  
  for (i in c(1:nrow(categories))) {
   
    row_condition_columns <- c(comparisons, 'sum')
    
    row_condition <- categories[i, row_condition_columns]
    
    mask <- NULL
    
    partial_masks <- list()
    
    for(colname in row_condition_columns) {
      
      value <- row_condition[, colname]
      
      mask <- combinatory_comparisons_with_colors_and_legends[, colname] == value
      
      partial_masks[[colname]] <- mask
      
    }
    
    
    logical_result <- NULL
    
    for(m in partial_masks) {
      
      if(is.null(logical_result)) {
        logical_result <- m
      } else {
        logical_result <- logical_result & m
      }
      
    }
    
    categories$category[i] <- combinatory_comparisons_with_colors_and_legends[logical_result, 'color']
    
    categories$legend[i] <- combinatory_comparisons_with_colors_and_legends[logical_result, 'legend']
    
                                                                 
  }
  
  categories_order <- combinatory_comparisons_with_colors_and_legends$color
  legend_order <- combinatory_comparisons_with_colors_and_legends$legend
  categories$category <- factor(categories$category, categories_order)
  categories$legend <- factor(categories$legend, legend_order)
  
  
  
  # only 1
  # denise1000_cat_mask <- (categories$sum == 1) & (categories$denise1000_vs_all == 1) 
  # predOr_cat_mask <- (categories$sum == 1) & (categories$predOr_vs_all == 1)
  # exp_cat_mask <- (categories$sum == 1) & (categories$exp_vs_all == 1)
  # union_cat_mask <- (categories$sum == 1) & (categories$union_vs_all == 1)
  # 
  # # only 2
  # intersection_2_exp_predOr_cat_mask <- categories$sum == 2  & (categories$exp_vs_all == 1) & (categories$predOr_vs_all == 1)
  # intersection_2_exp_denise1000_cat_mask <- categories$sum == 2  & (categories$exp_vs_all == 1) & (categories$denise1000_vs_all == 1)
  # intersection_2_exp_union_cat_mask <- categories$sum == 2  & (categories$exp_vs_all == 1) & (categories$union_vs_all == 1)
  # intersection_2_predOr_denise1000_cat_mask <- categories$sum == 2  & (categories$predOr_vs_all == 1) & (categories$denise1000_vs_all == 1)
  # intersection_2_predOr_union_cat_mask <- categories$sum == 2  & (categories$predOr_vs_all == 1) & (categories$union_vs_all == 1)
  # intersection_2_denise1000_union_cat_mask <- categories$sum == 2  & (categories$denise1000_vs_all == 1) & (categories$union_vs_all == 1)
  # 
  # # only 3
  # # no exp:
  # intersection_3_denise1000_union_predOr_cat_mask <- categories$sum == 3  & (categories$denise1000_vs_all == 1) & (categories$union_vs_all == 1) & (categories$predOr_vs_all == 1)
  # # no predOr
  # intersection_3_denise1000_union_exp_cat_mask <- categories$sum == 3  & (categories$denise1000_vs_all == 1) & (categories$union_vs_all == 1) & (categories$exp_vs_all == 1)
  # #no denise1000
  # intersection_3_predOr_exp_union_cat_mask <- categories$sum == 3  & (categories$predOr_vs_all == 1) & (categories$exp_vs_all == 1) & (categories$union_vs_all == 1)
  # # no union
  # intersection_3_predOr_exp_denise1000_cat_mask <- categories$sum == 3  & (categories$predOr_vs_all == 1) & (categories$exp_vs_all == 1) & (categories$denise1000_vs_all == 1)
  # 
  # # all - only 4
  # intersection_cat_mask <- categories$sum == 4
  
  # categories$category[denise1000_cat_mask] <- "yellow"
  # categories$legend[denise1000_cat_mask] <- "et vs genome"
  # 
  # categories$category[predOr_cat_mask] <- "blue" #"predOr_vs_all"
  # categories$legend[predOr_cat_mask] <- "predOr vs genome" #"predOr_vs_all"
  # 
  # categories$category[exp_cat_mask] <- "red" #"exp_vs_all"
  # categories$legend[exp_cat_mask] <- "exp vs genome" #"predOr_vs_all"
  # 
  # categories$category[union_cat_mask] <- "turquoise" #"exp_vs_all"
  # categories$legend[union_cat_mask] <- "union vs genome" #"predOr_vs_all"
  # 
  # #---
  # 
  # categories$category[intersection_2_exp_predOr_cat_mask] <- "darkviolet"
  # categories$legend[intersection_2_exp_predOr_cat_mask] <- "exp vs genome AND predOr vs genome" 
  # 
  # categories$category[intersection_2_exp_denise1000_cat_mask] <- "darkorange" 
  # categories$legend[intersection_2_exp_denise1000_cat_mask] <- "exp vs genome AND et vs genome" 
  # 
  # categories$category[intersection_2_exp_union_cat_mask] <- "palevioletred" 
  # categories$legend[intersection_2_exp_union_cat_mask] <- "exp vs genome AND union vs genome" 
  # 
  # categories$category[intersection_2_predOr_denise1000_cat_mask] <- "seagreen" 
  # categories$legend[intersection_2_predOr_denise1000_cat_mask] <- "predOr vs genome AND et vs genome"
  # 
  # categories$category[intersection_2_predOr_union_cat_mask] <- "dodgerblue" 
  # categories$legend[intersection_2_predOr_union_cat_mask] <- "predOr vs genome AND union vs genome"
  # 
  # categories$category[intersection_2_denise1000_union_cat_mask] <- "palegreen" 
  # categories$legend[intersection_2_denise1000_union_cat_mask] <- "et vs genome AND union vs genome"
  # 
  # #---
  # 
  # categories$category[intersection_3_denise1000_union_predOr_cat_mask] <- "lavender"
  # categories$legend[intersection_3_denise1000_union_predOr_cat_mask] <- "et vs genome AND union vs genome AND predOr vs genome"
  # 
  # categories$category[intersection_3_denise1000_union_exp_cat_mask] <- "lavenderblush"
  # categories$legend[intersection_3_denise1000_union_exp_cat_mask] <- "et vs genome AND union vs genome AND exp vs genome"
  # 
  # categories$category[intersection_3_predOr_exp_union_cat_mask] <- "gainsboro"
  # categories$legend[intersection_3_predOr_exp_union_cat_mask] <- "preOr vs genome AND exp vs genome AND union vs genome"
  # 
  # categories$category[intersection_3_predOr_exp_denise1000_cat_mask] <- "thistle" 
  # categories$legend[intersection_3_predOr_exp_denise1000_cat_mask] <- "preOr vs genome AND exp vs genome AND et vs genome"
  # 
  # 
  # #---  
  # 
  # categories$category[intersection_cat_mask] <- "darkgray" #intersection all"
  # categories$legend[intersection_cat_mask] <- "exp vs genome AND predOr vs genome AND et vs genome AND union vs genome"
  # 
  
    
  # # TODO ver en que orden es mas conveniente mostrarlos
  # #categories$category <- as.factor(categories$category)
  # categories$category <- factor(categories$category, levels = c("yellow", "blue", "red", "turquoise",
  #                                                               "darkviolet", "darkorange", "palevioletred", "seagreen", "dodgerblue", "palegreen",
  #                                                               "lavender", "lavenderblush", "gainsboro", "thistle" ,
  #                                                               "darkgray"))
  # 
  # #categories$legend <- as.factor(categories$legend)
  # categories$legend <- factor(categories$legend, levels = c("et vs genome", "predOr vs genome", "exp vs genome", "union vs genome",
  #                                                           # 2:
  #                                                           "exp vs genome AND predOr vs genome" , "exp vs genome AND et vs genome" , "exp vs genome AND union vs genome" ,
  #                                                           "predOr vs genome AND et vs genome", "predOr vs genome AND union vs genome", "et vs genome AND union vs genome",
  #                                                           # 3:
  #                                                           "et vs genome AND union vs genome AND predOr vs genome", "et vs genome AND union vs genome AND exp vs genome",
  #                                                           "predOr vs genome AND exp vs genome AND union vs genome", "preOr vs genome AND exp vs genome AND et vs genome",
  #                                                           # 4:
  #                                                           "exp vs genome AND predOr vs genome AND et vs genome AND union vs genome"))

  
  return(categories)
  
}

#no hace falta pasar go_terms_subset porque estan d=slot(object,"terms_dist")
get_mds_coords <- function(object, type, go_term_size=NULL) {
  
  # check class
  if(!is(object,"GO_SS") & !is(object,"GO_clusters")){
    stop("object must come from ViSEAGO::compute_SS_distances() or ViSEAGO::GOterms_heatmap()")
  }
  
  # check class
  if(is(object,"GO_SS") & type=="GOclusters"){
    stop("show_clusters is only available for GO_clusters class after clusters SS distance calculation with ViSEAGO::compute_SS_distances()")
  }
  
  # check type argument
  type=match.arg(type,c("GOterms","GOclusters"))
  
  # for SS_dist from object
  if(type=="GOterms"){
    
    # import SS_dist from object
    d=slot(object,"terms_dist")
    
    # if empty
    if(length(d)==0){
      stop("Please use GO_terms class object with computed SS distances using ViSEAGO::compute_SS_distance()")
    }
    
    # measures
    measures=names(d)
    
    ## MDS
    
    # run MDS
    res.mds=lapply(measures,function(x){
      
      # run MDS
      res.mds <-cmdscale(d[[x]], eig =T, k = 2)
      
      # extract point values
      res.mds<-res.mds$points
      
      # convert to data.table
      data.table(
        GO.ID=attr(res.mds,"dimnames")[[1]],
        Dim.1=res.mds[,1],
        Dim.2=res.mds[,2],
        measure=x
      )
    })
    
    # bind columns
    res.mds=rbindlist(res.mds)
    
    # for GO_SS
    if(is(object,"GO_SS")){

      #find terms description from go ids in res.mds:
      res.mds_go_ids <- res.mds$GO.ID
      enrich_GOs_data <- slot(slot(object,"enrich_GOs"),"data")
      term_description_mask <- sapply(enrich_GOs_data$GO.ID, function(i) i %in% res.mds_go_ids)
      term_description <- enrich_GOs_data$term[term_description_mask]
      
        # GO names
        res.mds=data.table(
          res.mds,
          #term=slot(slot(object,"enrich_GOs"),"data")$term
          term=term_description
        )
      
      
      # GO categories (paula)
      categories_raw <- get_terms_categories(object)
      categories <- categories_raw[, c("GO.ID", "category", "legend")]
      res.mds = merge(res.mds, categories, by = "GO.ID", all.x = TRUE)
      
      # GO sizes (paula)
      if(! is.null(go_term_size)) {
        res.mds = merge(res.mds, go_term_size, by = "GO.ID", all.x = TRUE)
      } else {
        # default value:
        res.mds$go_term_size <- 20
      }
      
      
      # add levels to measures
      res.mds$measure<-factor(
        res.mds$measure,
        levels=unique(res.mds$measure)
      )
      
    }else{
      
      #TODO
    }
    
  }else{
    
    # # import SS_dist from object
    # d=slot(object,"clusters_dist")
    # 
    # # if empty
    # if(length(d)==0){
    #   stop("Please use GO_clusters class object with computed clusters distances using ViSEAGO::compute_SS_distance()")
    # }
    # 
    # # measures
    # measures=names(d)
    # 
    # ## MDS
    # 
    # # run MDS
    # res.mds=lapply(measures,function(x){
    #   
    #   # run MDS
    #   res.mds <-cmdscale(d[[x]], eig = T, k = 2)
    #   
    #   # extract point values
    #   res.mds<-res.mds$points
    #   
    #   # convert to data.table
    #   data.table(
    #     GO.cluster=attr(res.mds,"dimnames")[[1]],
    #     Dim.1=res.mds[,1],
    #     Dim.2=res.mds[,2],
    #     measure=x
    #   )
    # })
    # 
    # # bind columns
    # res.mds=rbindlist(res.mds)
    # 
    # # custom text
    # res.mds[,
    #         `:=`(
    #           "text"=res.mds$GO.cluster,
    #           GO.cluster=gsub("_.+$","",res.mds$GO.cluster)
    #         )
    # ]
    # 
    # # add levels to measures
    # res.mds$measure<-factor(
    #   res.mds$measure,
    #   levels=unique(res.mds$measure)
    # )
    # 
    # # add levels to GO.cluster 
    # res.mds$GO.cluster<-factor(
    #   res.mds$GO.cluster,
    #   levels=unique(res.mds$GO.cluster)
    # )
    # 
    # # add GO.ID for GO
    # res.mds[,"text":=paste("GO.cluster:",text)]
    # 
    # # add GO.ID for GO
    # res.mds[,"text":=gsub("_GO","<br>GO.ID: GO",text)]
    # 
    # # add GO.name for term definition
    # res.mds[,"text":=gsub("_"," <br>GO.name: ",text)]
  }
  
  ## plot
  
  # mds coordinates
  #values=unlist(res.mds$Dim.1,res.mds$Dim.2)
  
  
  #result <- list("res.mds" = res.mds, "values" = values)
  #return(result)
  
  
  result <- list("measures" = measures, "res.mds" = res.mds)
  return(result)
  
}


# object <- myGOs_distances
# type <- "GOterms"
# file <- NULL
# TODO: hacer que reciba una lista de objects en lugar de dos objects:
MDSplot_custom <- function(object, type, title = "MultiDimensional Scaling plot", 
                          file, go_term_size=NULL, show_legend = TRUE,
                          show_measures = TRUE, 
                          x_tick_values=NULL, y_tick_values= NULL,
                          x_range=NULL, y_range=NULL){
    
  
    mds_data <- get_mds_coords(object, type, go_term_size)
    res.mds <- mds_data$res.mds
    measures <- mds_data$measures
    
  
    # init graph
    p <- plot_ly()
    
    # for GO_SS
    if(is(object,"GO_SS")){
      
      # add trace to plot by measure
      for(x in seq_along(measures)){
        
        # default visualization
        
        if(x==1){visible=TRUE}else{visible=FALSE}
        
        data_mds <- res.mds[measures[x], on="measure"]
        legend_categories <- levels(data_mds$legend)
        
        for (legend_category in legend_categories) {
          
          data_mds_legend_category_mask <- data_mds$legend == legend_category
          data_mds_legend <- data_mds[data_mds_legend_category_mask, ]
          
          # create trace
          p<-add_markers(
            p,
            #data=res.mds[measures[x],on="measure"],
            data = data_mds_legend,
            x = ~Dim.1,
            y = ~Dim.2,
            #name=measures[x],
            name = legend_category,
            text=~paste('GO.ID:',GO.ID,'<br>GO.name:',term),
            showlegend=TRUE,
            marker =list(
              # fix paula:
              size = ~go_term_size,
              opacity = 0.7,
              #color = data_mds_legend$category # "royalblue"
              color = ~category,  # "royalblue"
              line = list(width=2, color=~category)
            ),
            visible=visible
          )
          
        }
        
      }
      
      xaxis_values = list(title="Dimension 1")
      if(! is.null(x_tick_values)) {
        xaxis_values[["tickvals"]] = x_tick_values
      }
      if(! is.null(x_range)) {
        xaxis_values[["range"]] = x_range
      }
        
      yaxis_values = list(title="Dimension 2") 
      if(! is.null(y_tick_values)) {
        yaxis_values[["tickvals"]] = y_tick_values
      }
      if(! is.null(y_range)) {
        yaxis_values[["range"]] = y_range
      }
      
      
      # add custom layout with dropdown menu
      p<-layout(
        
        p,
        
        # add title
        title = title,
        
        # increase font size
        font = list(size=14),
        
        #paula aca 20210714
        
        # add axis legends
        xaxis = xaxis_values,
        yaxis = yaxis_values,
        # param yticks
        
        # plot_ly(x=~x, y=~y) %>%
        #   layout(
        #     xaxis = list(
        #       range=c(20,40)
        #     )
        #   )
        
        # a <- list(
        #   autotick = FALSE,
        #   ticks = "outside",
        #   tick0 = 0,
        #   dtick = 0.25,
        #   ticklen = 5,
        #   tickwidth = 2,
        #   tickcolor = toRGB("blue")
        # )
        
        
        # add margin
        margin = list(t = 100),

        legend = list(orientation = "v"),   # show entries horizontally
                     #xanchor = "right"),  # use center of legend as anchor
#                     x = 0.5),
        
                
        # legend = list(orientation = "h",   # show entries horizontally
        #               xanchor = "center",  # use center of legend as anchor
        #               x = 0.5),
        
        # legend = list(orientation = "v",   # show entries horizontally
        #               xanchor = "left",  # use center of legend as anchor
        #               legend.position = c(1.1, 0.2)
        # ),
        
        
        showlegend = show_legend,
        
        if(show_measures) {
        
          # add scrolling menu for availables measures
          updatemenus = list(
            
            # measure dropdown menu
            list(
              
              # button position
              x = 0.1,
              y = 1.1,
              
              buttons = lapply(seq_along(measures), function(x){
                
                # init visibility to all FASLE
                values = rep(FALSE, length(measures) * length(legend_categories))
                
                # turn to true for each measures
                min_index <- length(legend_categories) * (x-1) + 1
                max_index <- min_index + length(legend_categories) - 1 #7 categorias para cada medida, length(legend_categories)
                values[min_index:max_index] <- TRUE
                #values[x] <- TRUE
                
                # output button
                list(
                  method = "restyle",
                  args = list("visible", as.list(values)),
                  label = measures[x]
                )
              })
            )
            
          )
        }          
      )
    }
    
    # for GO_SS
    # if(is(object,"GO_clusters")){
    #   
    #   # color labels of cutting tree
    #   colors= unique(
    #     get_leaves_attr(
    #       slot(object,"dendrograms")$GO,
    #       "edgePar"
    #     )
    #   )
    #   
    #   # extract the GO. cluster column
    #   count<-slot(slot( object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE]
    #   
    #   # count the number of terms by clusters
    #   count<-count[,list("nb"=.N),by="GO.cluster"]
    #   
    #   # add count to res.mds
    #   res.mds=merge(
    #     res.mds,
    #     count,
    #     by="GO.cluster",
    #     sort=FALSE
    #   )
    #   
    #   # ordering clusters
    #   res.mds[,
    #           GO.cluster:=factor(GO.cluster,levels=sort(unique(as.numeric(GO.cluster))))
    #   ]
    #   
    #   # for terms
    #   if(type=="GOterms"){
    #     
    #     # create trace
    #     p<-add_markers(
    #       p,
    #       data=res.mds,
    #       x=~Dim.1,
    #       y=~Dim.2,
    #       color=~GO.cluster,
    #       text = ~paste("cluster:",GO.cluster,"<br>GO.ID:",GO.ID,"<br>GO.name:",term),
    #       showlegend=TRUE,
    #       colors=colors,
    #       marker =list(
    #         size =20,
    #         opacity = 0.7
    #       )
    #     )
    #     
    #     # add custom layout with dropdown menu
    #     p<-layout(
    #       p,
    #       
    #       # add title
    #       title=paste(measures,"distance MultiDimensional Scaling plot"),
    #       
    #       # increase font size
    #       font=list(size=14),
    #       
    #       # add axis legends
    #       xaxis=list(title="Dimension 1"),
    #       yaxis=list(title="Dimension 2"),
    #       
    #       # add margin
    #       margin=list(t=100)
    #     )
    #   }else{
    #     
    #     # add count to text
    #     res.mds[,"text":=paste(sub("<br>.+$","",text),"<br>GO.count:",res.mds$nb,sub("^.+<br>GO.ID","<br>GO.ID",text))]
    #     
    #     # add trace to plot by measure
    #     for(x in seq_along(measures)){
    #       
    #       # default visualization
    #       if(x==1){visible=TRUE}else{visible=FALSE}
    #       
    #       # create trace
    #       p<-add_markers(
    #         p,
    #         data=res.mds[measures[x],on="measure"],
    #         x=~Dim.1,
    #         y=~Dim.2,
    #         name=measures[x],
    #         text=~text,
    #         showlegend=FALSE,
    #         sizes=c(20,50),
    #         size=~nb,
    #         marker =list(
    #           sizemode = 'diameter',
    #           opacity = 0.4,
    #           line=list(color=colors),
    #           color=colors
    #         ),
    #         visible=visible
    #       )
    #     }
    #     
    #     # add custom layout with dropdown menu
    #     p<-plotly::layout(p,
    #                       
    #                       # add title
    #                       title="MultiDimensional Scaling plot",
    #                       
    #                       # increase font size
    #                       font=list(size=14),
    #                       
    #                       # add axis legends
    #                       xaxis=list(title="Dimension 1"),
    #                       yaxis=list(title="Dimension 2"),
    #                       
    #                       # add margin
    #                       margin=list(t=100),
    #                       
    #                       # add scrolling menu for availables measures
    #                       updatemenus = list(
    #                         
    #                         # measure dropdown menu
    #                         list(
    #                           
    #                           # button position
    #                           x = 0.1,
    #                           y = 1.1,
    #                           
    #                           buttons = lapply(seq_along(measures),function(x){
    #                             
    #                             # init visibility to all FASLE
    #                             values=rep(FALSE,length(measures))
    #                             
    #                             # turn to true for each measures
    #                             values[x]<-TRUE
    #                             
    #                             # output button
    #                             list(method ="restyle",
    #                                  args = list("visible",as.list(values)),
    #                                  label=measures[x])
    #                           })
    #                         )
    #                       )
    #     )
    #   }
    # }
    
    # return or print
    if(is.null(file)){
      
      # return the plot
      p
      
    }else{
      
      # print heatmap
      orca(p,file)
    }
    
    return(p)
  }

