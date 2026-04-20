#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Correlation analysis, clustering, network prep
# 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

source("R/utils.R")
source("R/config_plotting.R")

calculate_correlations <- function(melt_species, cell_type){

  # Remove whole cells and rearrange to wide 
  tableau_data <- melt_species %>%
    # Remove Cholesterol and recalculate molpct
    filter(cell %in% cell_type,
           organelle != "Whole cells") %>%
    mutate(organelle = droplevels(organelle)) %>%
    filter(!class == "Chol") %>%
    group_by(organelle, exp, id) %>%
    mutate(total_molpct = sum(molpct, na.rm = T)) %>%
    ungroup() %>%
    mutate(molpct_no_chol = molpct/total_molpct*100)%>%
    # Rearrange to wide
    mutate(organelle_exp = interaction(cell, organelle, exp,id, sep = "_")) %>%
    dplyr::arrange(organelle, exp, id, species)%>%
    dplyr::select(species, organelle_exp, molpct_no_chol) %>%
    pivot_wider(values_from = molpct_no_chol,
                names_from = organelle_exp)

  # Sort alphabetically (for consistency)
  tableau_data <- tableau_data[order(tableau_data$species),]
    
  # Convert tablau data to matrix
  matrix_data <- tableau_data %>%
    column_to_rownames("species")  %>%
    as.matrix()
  
  # Replace NAs with 0
  matrix_data[is.na(matrix_data)] <- 0
  
  # Correlate species
  cor_matrix <- cor(t(matrix_data))
  
  return(cor_matrix)

}

cluster_cor_matrix <- function(cor_matrix,
                               lipidQuan_output_folder, 
                               exp_id,
                               cell_type, 
                               n_cluster, 
                               cluster_names) {
  

  
  # Calculate distances based on the correlation space
  dist_mat <- dist(cor_matrix, method = "euclidean")
    
  # Cluster 
  hclust_rows <- hclust(d = dist_mat, method = "complete")

  clusters <- cutree(hclust_rows, k = n_cluster)
  
  clusters <- tibble(species = names(clusters), cluster = clusters) %>%
    mutate(cluster_name = cluster_names[match(cluster, 1:n_cluster)])
  
  # Save heatmap with clusters
  annotation_df <- clusters %>%
    column_to_rownames(var = "species")
  
  heatmap_filepath <- file.path(lipidQuan_output_folder,"Results", paste0("PAPSL_", exp_id,"_",paste0(cell_type, collapse = "_"), ".pdf"))
  
  pheatmap::pheatmap(
    cor_matrix,
    color = color2,
    breaks = seq(-1, 1, length.out = 102),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize = 8,
    cutree_rows = n_cluster,
    cutree_cols = n_cluster,
    cellwidth = 1,
    cellheight = 1,
    border_color = NA,
    show_rownames = FALSE,
    annotation_row = annotation_df["cluster_name"],
    show_colnames = FALSE,
    filename = heatmap_filepath
  )
  
  cat("* Cluster heatmap saved as", heatmap_filepath, "\n")
  
  return(clusters)
  
}

# Generate correlation table for network visualization
summarise_cor <- function(cor_matrix, clusters) {  
  cor_table <- cor_matrix %>%
    as_tibble(rownames = "species") %>%
    pivot_longer(cols = -species, values_to = "cor", names_to = "species_target") %>%
    filter(species != species_target) %>%
    left_join(., clusters,by = join_by(species))%>%
    dplyr::bind_cols(., 
                     lipid_info(.$species), 
                     lipid_info(.$species_target, suffix = "_target")) %>%
    dplyr::mutate(same_class = class == class_target, 
                  same_db = db == db_target,
                  same_length = length == length_target,
                  same_cat = cat == cat_target)
  
  cor_table_summarise <- cor_table %>%
    dplyr::select(same_class, same_db, same_length, same_cat) %>%
    pivot_longer(cols = everything()) %>%
    group_by(name,value) %>%
    summarise(n = n(), .groups = "drop")
  
  cor_table <- cor_table %>%
    mutate(same_classlabel = ifelse(same_class,
                                    paste0("Same-class comparison (n = ",
                                           filter(cor_table_summarise, name == "same_class", value == TRUE)$n,
                                           ")"),
                                    paste0("Different-class comparison (n = ",
                                           filter(cor_table_summarise, name == "same_class", value == FALSE)$n
                                           ,")")
    ),
    same_dblabel = ifelse(same_class,
                          paste0("Same-db comparison (n = ",
                                 filter(cor_table_summarise, name == "same_db", value == TRUE)$n,
                                 ")"),
                          paste0("Different-db comparison (n = ",
                                 filter(cor_table_summarise, name == "same_db", value == FALSE)$n
                                 ,")")
    ),
    same_lengthlabel = ifelse(same_class,
                              paste0("Same-length comparison (n = ",
                                     filter(cor_table_summarise, name == "same_length", value == TRUE)$n,
                                     ")"),
                              paste0("Different-length comparison (n = ",
                                     filter(cor_table_summarise, name == "same_length", value == FALSE)$n
                                     ,")")
    ),
    same_catlabel = ifelse(same_class,
                           paste0("Same-cat comparison (n = ",
                                  filter(cor_table_summarise, name == "same_cat", value == TRUE)$n,
                                  ")"),
                           paste0("Different-cat comparison (n = ",
                                  filter(cor_table_summarise, name == "same_cat", value == FALSE)$n
                                  ,")")
    )
    )
}

# Main ----

correlation_and_clustering <- function(lipidQuan_output_folder, 
                                       cell_type, 
                                       output_file,
                                       exp_id,
                                       n_cluster,
                                       cluster_name) {
  
  ## load quantified and filtered data
  load(file.path(lipidQuan_output_folder,
                 "Results", 
                 paste0("PAPSL_",exp_id, ".rdata")))
  
  # Calculate correlations and clusters
  cor_matrix <- calculate_correlations(melt_species, cell_type) 
  
  clusters <- cluster_cor_matrix(cor_matrix, lipidQuan_output_folder, 
                                 exp_id = exp_id, 
                                 cell_type = cell_type, 
                                 n_cluster = n_cluster,
                                 cluster_names = cluster_name)

  cor_table <- summarise_cor(cor_matrix, clusters)
  
  # Save data for Cytoscape and figures
  save_to <- file.path(lipidQuan_output_folder, "Results",output_file)
  
  save(list = c("cor_table"),
       file = save_to
  )
  
  cat("* Data saved as:", save_to, "\n")
  
  write_csv(cor_table[abs(cor_table$cor)>0.7,], file = file.path(lipidQuan_output_folder, 
                                                           "Results",
                                                           paste0("Network_cluster_", paste0(cell_type, collapse = "_"),"_cor0.7_noChol.csv")))

  
}




