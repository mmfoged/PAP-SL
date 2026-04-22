#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Merge and filter data after lipidQuan 
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Packages and helpers

source(here("R/utils.R"))
source(here("R/config_plotting.R"))

## 1: Merging tableauOutput files ----------------------------------------------

#' Merge species data
#'
#' This function merges several LipidQuan output datasets into one. Use for 
#' data recorded and processed on different days. Adds batch number to column
#' names based on the input order of file names.
#' Writes to file named tableauOut_merged'.
#'
#' @param files character vector containing 'tableauOut' file names of datasets to be 
#' merged.
#' @return A tibble with combined datasets. First column contains species names. 
#' @export
mergeSpecies <- function(files,dir){ 
  
  list_datasets <- list()
  
  for (file_name in files){
    
    batch_no <- which(files == file_name)
    
    # LipidQuant adds classes in tableau output, remove these, added later
    # Classes have no ":" in the name
    file_read <- read_csv(file.path(dir, file_name), show_col_types = FALSE) %>%
      filter(grepl(":", pmol))
    
    # Data control, check for infinite values
    inf_columns <- file_read %>%
      summarise(across(where(is.numeric), ~ any(is.infinite(.x)))) %>%
      pivot_longer(everything(), names_to = "column", values_to = "has_inf") %>%
      filter(has_inf) %>%
      dplyr::pull(column)
    
    if(length(inf_columns) > 0) {
      cat("Files with inf-values (inf set to 0):\n")
      cat("File:", file_name, "\nColumn(s):",inf_columns, "\n")
    }
    
    # change column names to include batch number
    sample_cols <- seq_len(ncol(file_read))[-1]
    colnames(file_read) <- c("species",
                             paste0("sample_",
                                    seq_along(sample_cols),
                                    "_",
                                    batch_no
                             ))
    
    list_datasets[[file_name]] <- file_read
  }
  
  
  # Join and convert NAs and Inf values to 0
  merged <- list_datasets %>%
    purrr::reduce(full_join, by = "species") %>%
    dplyr::mutate(across(where(is.numeric), ~ifelse(!is.finite(.x), 0, .x)))
  
  # write to file
  write_csv(merged, file = file.path(dir, "tableauOut_merged.csv")) 
  
  return(merged)   
}

#' Merge tableau output files
#' 
#' Merges all files ending with "tableauOutput.csv" in the input folder.
#' 
#'@param dir folder containing lipidQuan output files
merge_files <- function(dir, samples_path, order_organelles, removed_raws = NULL) {
  
  tableau_filenames <- list.files(path = dir, pattern = "tableauOutput.csv$") 
  
  # Merge all out-files
  merged_species <- mergeSpecies(tableau_filenames, dir) 
  
  # read sample info, add necessary columns if missing
  samples <- read_csv(samples_path, col_types = cols(.default = "c")) %>%
    mutate(
      exp = if ("exp" %in% names(.)) exp else 1,
      no = if ("no" %in% names(.)) no else 1,
      organelle = if ("organelle" %in% names(.)) organelle else "whole_cells",
      file = if ("file" %in% names(.)) file else paste0("raw_",1:nrow(.)),
      id = if ("id" %in% names(.)) id else sprintf("%03d", 1:nrow(.)),
      treatment = if ("treatment" %in% names(.)) treatment else "none",
      cell = ifelse(organelle == "blank","blank",cell)) %>%
    CorrectOrganelle(., order_organelles = order_organelles)

  # Data check
  if(!nrow(samples)==(ncol(merged_species)-1)) {
    stop("Samples file does not fit the merged samples:")
  }
  
  # Add corresponding column names to sample info
  samples$colname <- colnames(merged_species)[-1]
  
  # Long format
  merged_species_long <- merged_species %>%
    pivot_longer(cols = -species,
                 names_to = "colname",
                 values_to = "pmol") %>%
    left_join(samples, ., by = "colname")%>%
    bind_cols(., lipid_info(.$species)) 
  
  # Remove poor quality raw files if defined in process_all.Rmd
  removed_raws <- intersect(removed_raws, merged_species_long$file)
  
  if (!is.null(removed_raws) && length(removed_raws)>0) {
      
    cat("removed poor quality raw files:\n")
    print(removed_raws)
    
    merged_species_long <- merged_species_long %>%
      filter(!file %in% removed_raws)
  
  }
  
  return(merged_species_long)
  
}



## 2: Aggregate on means of replicates -----------------------------------------

# make data frame with aggregated data

aggregate_species <- function(df) {
  
  df %>%
    group_by(across(-c(file, colname, pmol, no, id))) %>%
    summarise(pmol = mean(pmol, na.rm =T), .groups = "drop",
              file = paste0(file, collapse ="|"), 
              colname = paste0(colname, collapse = "|"),
              id = paste0(id, collapse = "|"))
  
}

## 3: Blank filter -------------------------------------------------------------

blank_filter <- function(df, blank_factor, min_pmol, missing_exp_allowed,min_detected) {
  
  # Calculate blank filter
  blank_filter_df <- df %>%
    filter(cell == "blank") %>%
    group_by(exp, species) %>%
    summarise(pmol_blank = mean(pmol, na.rm = T), .groups = "drop") %>%
    mutate(blank_filter = ifelse(pmol_blank * blank_factor > min_pmol, 
                                 pmol_blank * blank_factor, 
                                 min_pmol)) 
  
  # Subtract blanks if pmol is higher than blank filter for any organelle in 
  # each replicate (else set to 0)
  blank_filter_applied <- df %>%
    left_join(., blank_filter_df,by = join_by(exp, species))%>%
    mutate(pmol_blank_filter_pass = ifelse(pmol > blank_filter, 1, 0)) %>%
    group_by(exp, species) %>%
    summarise(keep = sum(pmol_blank_filter_pass)>=min_detected, .groups = "drop")
  
  species_blank_filter <- df %>%
    left_join(., blank_filter_df,by = join_by(exp, species))  %>%
    left_join(., blank_filter_applied,by = join_by(exp, species)) %>%
    filter(keep) %>%
    mutate(pmol = ifelse(pmol>pmol_blank,pmol-pmol_blank, 0))
  
  
  # Replicate filter
  # Remove species that are missing in more than 'missing_exp_allowed' experiments
  
  no_exp <- length(unique(df$exp)) - missing_exp_allowed
  
  rep_filter <- species_blank_filter %>%
    group_by(species, exp) %>%
    summarise(n_above_threshold = sum(pmol>0), .groups = "drop") %>%
    group_by(species) %>%
    summarise(n = sum(n_above_threshold>0), keep = n>= no_exp, .groups = "drop")
  
  keep_species_replicate_filter <- rep_filter %>%
    filter(keep) %>%
    pull(species) %>%
    unique()
  
  species_replicate_filter <- species_blank_filter %>%
    filter(species %in% keep_species_replicate_filter & pmol > 0,
           cell != "blank") %>%
    dplyr::select(-c(pmol_blank, blank_filter, keep))
  
  return(species_replicate_filter)
}

## 4: Database filter ####
# Output: merged.species.pmol, samples.noblank, with NAs


database_filter <- function(df, database_path) {
  
  # Import database
  database <- read_csv(database_path, show_col_types = FALSE) %>%
    dplyr::pull(species)
  
  # Rename D4-labelled species to match format
  df$species[grepl("^D4", df$species)] <- 
  gsub("O-","O",gsub("^D4","",gsub("[[:space:]]", "D4 ", df$species[grepl("^D4", df$species)])))
  
  # add BMP to the database
  database <- c(database, str_replace(database[grepl("^PG ", database)], 
                                      pattern = "PG", 
                                      replacement = "BMP"))
  
  # Filter
  df %>%
    filter(species %in% database & pmol > 0)
  
}





## 5: Calculate molpct, classpct etc ----

calculate_molpct <- function(df, filtered_classes) {
  
  df %>%
    filter(!lipid_class_name(species) %in% filtered_classes) %>%
    mutate(organelle = factor(organelle, 
                              levels = order_organelles[order_organelles %in% organelle]))%>%
    group_by(organelle, exp, id) %>%
    mutate(molpct = pmol / sum(pmol, na.rm = T) * 100) %>%
    group_by(organelle, exp, class, id) %>%
    mutate(molpct_class = sum(molpct, na.rm = T), 
           classpct = pmol / sum(pmol, na.rm = T) * 100) %>%
    group_by(organelle, exp, cat, id) %>%
    mutate(molpct_cat = sum(molpct, na.rm = T), 
           catpct = pmol / sum(pmol, na.rm = T) * 100) %>%
    ungroup() %>%
    mutate(log2_molpct = log2(molpct))
  
}


## 6: Make class-level data and 'tableau' format ----

long_to_tableau_format <- function(df) {

  col_to_tableau <- function(df, col_name, level = "species") {
  
    df %>%
      mutate(organelle_exp = interaction(cell, organelle, exp,id, sep = "_")) %>%
      dplyr::arrange(organelle, exp, id,!!sym(level))%>%
      dplyr::select(!!sym(level), organelle_exp, !!sym(col_name)) %>%
      pivot_wider(values_from = !!sym(col_name),
                  names_from = organelle_exp)
  
  }
  
  # Make column with molpct relative to whole cells
  df <- df %>%
    filter(organelle == "Whole cells") %>%
    group_by(exp, cell, species) %>%
    summarise(wc_molpct = mean(molpct, na.rm = T), .groups = "drop") %>%
    left_join(df,., by = c("exp", "cell", "species")) %>%
    mutate(molpct_rel_to_wc = molpct / wc_molpct)

  # add n experiments per organelle  
  df <- df %>%
    group_by(organelle, species) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(organelle) %>%
    summarise(n = max(n), .groups = "drop") %>%
    left_join(df, ., by = "organelle")
  
  melt_classes_pex <- df %>%
    group_by(id, organelle, exp, cell, treatment, class, cat, file, colname, n) %>%
    summarise(molpct_class = unique(molpct_class), .groups = "drop") %>%
    mutate(log2_molpct_class = log2(molpct_class))
  
  melt_classes <- melt_classes_pex %>%
    filter(!organelle %in% c("LMF", "Peroxisomes", "Ctrl"))
  
  melt_species <- df %>%
    filter(!organelle %in% c("LMF", "Peroxisomes", "Ctrl"))

  tableau_pmol <-  col_to_tableau(melt_species, "pmol", level = "species")
  tableau_molpct <-  col_to_tableau(melt_species, "molpct", level = "species")
  tableau_log2 <-  col_to_tableau(melt_species, "log2_molpct", level = "species")
  tableau_classpct <- col_to_tableau(melt_species, "classpct", level = "species")
  tableau_class <- col_to_tableau(melt_classes, "molpct_class", level = "class")
  
  ## Summary of mean lengths and double bonds per organelle - separate classes
  mean_length_db_classes <- melt_species%>%
    group_by(cell, organelle, exp, db, length, class, cat, id) %>%
    summarise(.groups = "drop", double_bond_contribution = sum(classpct, na.rm = T) / 100 * db, 
              length_contribution      = sum(classpct, na.rm = T)  / 100 * length) %>%
    group_by(cell, organelle, exp, class, cat) %>%
    summarise(.groups = "drop", mean_db = sum(double_bond_contribution), 
              mean_length = sum(length_contribution))
  
  ## Summary of mean lengths and double bonds per organelle - merge classes
  # Summarise mol% without cholesterol
  mean_length_db_all <- melt_species %>%
    filter(!class == "Chol") %>%
    group_by(organelle, exp, id) %>%
    mutate(total_molpct = sum(molpct)) %>%
    ungroup() %>%
    mutate(molpct_no_chol = molpct/total_molpct*100,
           no_chains = ifelse(cat %in% c("LGPL", "LGPLO-"), yes = 1, 
                              ifelse(class %in% "CL", 4, 
                                     ifelse(class == "Chol", NA, 
                                            2)))) %>%
    group_by(id, cell, organelle, exp, db, length, class, no_chains) %>%
    summarise(double_bond_contribution = sum(molpct_no_chol, na.rm = T) / 100 * db, 
              length_contribution = sum(molpct_no_chol, na.rm = T)  / 100 * length,
              .groups = "drop", ) %>%
    group_by(id, cell, organelle, exp) %>%
    summarise(mean_db = sum(double_bond_contribution/no_chains), 
              mean_length = sum(length_contribution/no_chains),
              .groups = "drop", )
  
  
  # Sample info
  tableau_samples <- melt_species %>%
    mutate(organelle_exp = interaction(cell, organelle, exp,id, sep = "_")) %>%
    dplyr::arrange(organelle, exp, cell)%>%
    dplyr::select(organelle_exp, exp, organelle, cell, id, file, treatment) %>%
    unique()

  
  return(list(tableau_pmol = tableau_pmol, 
              tableau_molpct = tableau_molpct, 
              tableau_log2 = tableau_log2,
              tableau_classpct = tableau_classpct,
              tableau_samples = tableau_samples, 
              melt_classes = melt_classes,
              melt_species = melt_species,
              melt_classes_pex = melt_classes_pex,
              melt_species_pex = df,
              tableau_class = tableau_class, 
              mean_length_db_classes = mean_length_db_classes, 
              mean_length_db_all = mean_length_db_all))
}










## 7. filter_plot ----

arrange_filter_plot <- function(merged_species_long, 
                                species_replicate_filter, 
                                species_database_filter,
                                filtered_classes) {
  

filter_plot <- merged_species_long %>%
  filter(!lipid_class_name(species) %in% filtered_classes) %>%  
  filter(pmol > 0) %>%
  group_by(species, length, db, class) %>%
  summarise("merged_species" = n()>0, .groups = "drop")

filter_plot <- species_replicate_filter %>%
  filter(pmol > 0) %>%
  group_by(species, length, db, class) %>%
  summarise("rep_filter" = n()>0, .groups = "drop") %>%
  left_join(filter_plot, ., by = join_by(species, length, db, class)) %>%
  mutate(rep_filter = ifelse(is.na(rep_filter), FALSE, rep_filter))

filter_plot <- species_database_filter %>%
  filter(pmol > 0) %>%
  group_by(species, length, db, class) %>%
  summarise("database_filter" = n()>0, .groups = "drop") %>%
  left_join(filter_plot, ., by = join_by(species, length, db, class)) %>%
  mutate(database_filter = ifelse(is.na(database_filter), FALSE, database_filter)) 

filter_plot <- filter_plot %>%
  mutate(step = ifelse(database_filter, "Included",
                       ifelse(rep_filter, "Remove by database",
                              "Remove by blank"))) %>%
  mutate(lipid = paste(class, length), 
         class = factor(class, levels = order_class[order_class %in% class]))  


lipid_order <- as.character(filter_plot$lipid)[order(filter_plot$class, 
                                                     filter_plot$length, 
                                                     filter_plot$db )]

first_row_lipids <- lipid_order[seq_len(length(lipid_order) %/% 2)]

filter_plot <- filter_plot %>%
  mutate(row = ifelse(lipid %in% first_row_lipids & !grepl("^PG ", lipid), 1, 2), 
         lipid = factor(lipid, levels = unique(lipid_order)))

}

# 8. Summarize datasets and processing ----

#' Save processing parameters
save_processing_parameters <- function(blank_factor,
                                       min_pmol,
                                       missing_exp_allowed,
                                       filtered_classes,
                                       blank_filter_type,
                                       lipidQuan_output_folder){
  
  blank_filter_type <- paste(
    "Each species may only be missing in ",
    missing_exp_allowed,
    "replicate(s). A species is considered missing from a replicate if its pmol value is lower than mean(blank) *",
    blank_factor,
    "OR",
    min_pmol,
    "pmol in all organelles (including whole cells)."
  )
  
  processing_parameters <- list(merge_filter_date = date(),
                                blank_factor = blank_factor, 
                                min_pmol = min_pmol, 
                                missing_exp_allowed = missing_exp_allowed, 
                                manually_filtered_classes = filtered_classes, 
                                blank_filter_type = blank_filter_type)
  
  write.csv(do.call(c, processing_parameters), 
            file = file.path(lipidQuan_output_folder, 
                             "Results/processing_parameters.csv"))
  save(processing_parameters, 
       file=file.path(lipidQuan_output_folder, 
                      "Results/processing_parameters.RData"))
  
}


summarise_dataset_wide <- function(df, title)
{
  
  # Leave out "BMP" to avoid counting PG and BMP as two different classes
  df <- df %>%
    filter(!grepl("BMP", species)) 
  
  n_species <- nrow(df)
  n_classes <- length(unique(LipidClassName(df$species)))
  n_samples <- ncol(df)-1
  
  cat(rep("-", 80), "\n", sep ="")
  cat(title,": \n\n", sep = "")
  cat("number of species", n_species, "\n")
  cat("number of classes", n_classes, "\n")
  cat("number of samples", n_samples, "\n")
}


summarise_dataset_long <- function(df, title){
  
  df <- df %>%
    filter(!class == "BMP") 
  
  n_species <- length(unique(df$species))
  n_classes <- length(unique(df$class))
  n_samples <- length(unique(df$colname))
  
  cat(rep("-", 80), "\n", sep ="")
  cat(title,": \n\n", sep = "")
  cat("number of species", n_species, "\n")
  cat("number of classes", n_classes, "\n")
  cat("number of samples", n_samples, "\n")
  
}

# Main ----

# Merge datasets

merge_and_filter <- function(lipidQuan_output_folder, 
                             samples_path,
                             order_organelles,
                             filtered_classes,
                             removed_rows = NULL, 
                             blank_factor, 
                             min_pmol,
                             missing_exp_allowed,
                             min_detected,
                             database_path,
                             exp_ID){
  
  save_processing_parameters(blank_factor = blank_factor, 
                             min_pmol = min_pmol, 
                             missing_exp_allowed = missing_exp_allowed, 
                             filtered_classes = filtered_classes,
                             lipidQuan_output_folder = lipidQuan_output_folder)
  
  merged_species_long <- merge_files(dir = lipidQuan_output_folder, 
                                     samples_path = samples_path, 
                                     order_organelles = order_organelles, 
                                     removed_raws = removed_rows)
  
  summarise_dataset_long(merged_species_long, title = "Input data")
  
  # Aggregate technical replicates
  species_aggregate <- aggregate_species(merged_species_long)

  summarise_dataset_long(species_aggregate, title = "Aggregating samples")

  # blank and replicate filter
  species_replicate_filter <- blank_filter(df = species_aggregate, 
                                           blank_factor = blank_factor, 
                                           min_pmol = min_pmol, 
                                           missing_exp_allowed = missing_exp_allowed,
                                           min_detected = min_detected)

  
  summarise_dataset_long(species_replicate_filter, title = "Blank filter")
  
  # Database filter
  species_database_filter <- database_filter(species_replicate_filter, 
                                             database_path = database_path)
  
  summarise_dataset_long(species_database_filter,title = "Database filter")
  
  # Calculate molpct, classpct, log2 etc
  processed_data <- calculate_molpct(species_database_filter, filtered_classes)
  
  # Tableau data
  tableau_data <- long_to_tableau_format(processed_data)
  
  list2env(tableau_data, envir = .GlobalEnv)
  
  # Arrange plots
  filter_plot <- arrange_filter_plot(merged_species_long, 
                                     species_replicate_filter, 
                                     species_database_filter,
                                     filtered_classes)
  
  # csv files
  write_csv(
    tableau_pmol,
    file.path(lipidQuan_output_folder, "Results/processed_pmol.csv")
  )
  write_csv(
    tableau_molpct,
    file.path(lipidQuan_output_folder, "Results/processed_molpct.csv")
  )
  write_csv(
    tableau_samples,
    file.path(lipidQuan_output_folder, "Results/processed_samples.csv")
  )
  

  # supplementary data
  save(melt_species_pex, 
       melt_classes_pex, 
       species_database_filter,
       filter_plot, 
       file = file.path(lipidQuan_output_folder, 
                        "Results", 
                        paste0("PAPSL_",exp_ID, "_supplementary.rdata")))
  
  # main data
  save(tableau_samples, 
       tableau_log2, 
       tableau_classpct, 
       tableau_class, 
       tableau_pmol, 
       tableau_molpct, 
       melt_classes, 
       melt_species, 
       mean_length_db_classes, 
       mean_length_db_all,
       file = file.path(lipidQuan_output_folder,
                        "Results", 
                        paste0("PAPSL_",exp_ID, ".rdata")))

  
}


