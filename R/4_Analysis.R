
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Analysis
#
# * t-tests
# * ANOVA on class properties
# * PCA
# 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Load data and settings

source("R/utils.R")

# T-tests ----
run_all_t_tests <- function(df,
                            organelles_for_analysis,
                            test_val, 
                            value_col, 
                            compare, 
                            lipid_level,
                            p_adjust_method) {

  if(!compare %in% c("organelle", "cell")){
    stop("'compare' must be either 'organelle' or 'cell'")
  }
    
  df <- df %>%
    filter(organelle %in% organelles_for_analysis) %>%
    mutate(organelle = droplevels(organelle))
  
  # Make data frame with one row per comparison element
  t_test_df <- df %>%
    group_by(!!sym(lipid_level), !!sym(compare)) %>%
    summarise(vals = list(.data[[value_col]][!is.na(.data[[value_col]])]), .groups = "drop")
  
  # data frame with all comparisons to be made

  comparisons <- crossing(
    !!sym(lipid_level) := pull(df, !!sym(lipid_level)),
    name_a = as.character(pull(df, !!sym(compare))),
    name_b = as.character(pull(df, !!sym(compare)))
  ) %>%
    filter(name_a < name_b) %>%
    distinct()

  replace_null <- function(x) {
    if (is.null(x)) {
      NA_real_
      }else{
        x}
  }
  
  # Report number of comparisons
  n_compare <- length(unique(pull(df, !!sym(compare))))
  n_lipids <- length(unique(pull(df, !!sym(lipid_level))))
  lipid_level_name <- gsub("class", "classes", {{lipid_level}})
  cat("* t-tests: ", nrow(comparisons), " comparisons (", n_compare," ",{{compare}}, 
      "s, ",n_lipids, " ", lipid_level_name,").\n", sep = "")
  
  # Run lipid_t_test on all comparisons
  tests <- comparisons %>%
    left_join(t_test_df, by = c({{lipid_level}}, "name_a" = {{compare}})) %>%
    rename(lipid_a = vals) %>%
    left_join(t_test_df, by = c({{lipid_level}}, "name_b" = {{compare}})) %>%
    rename(lipid_b = vals) %>%
    dplyr::mutate(lipid_a = map(lipid_a, ~ replace_null(.x)), 
                  lipid_b = map(lipid_b, ~ replace_null(.x))) %>%
    rowwise() %>%
    # t-test
    summarise(
      lipid_t_test(
        name_a   = name_a,
        name_b   = name_b,
        lipid_a  = lipid_a,
        lipid_b  = lipid_b,
        lipid    = !!sym(lipid_level),
        test_val = test_val
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(padj = p.adjust(p_val, method = p_adjust_method)) %>%
    # Renaming for clarity
    dplyr::rename(!!sym(lipid_level) := species,
                  !!sym(paste0(compare, "_a")) := compare_a,
                  !!sym(paste0(compare, "_b")) := compare_b)
  
  return(tests)
}

#' lipid_t_test
#'


lipid_t_test <- function(name_a, name_b, lipid_a, lipid_b,lipid, test_val) {
  
  lipid_a <- lipid_a[!is.na(lipid_a)]
  lipid_b <- lipid_b[!is.na(lipid_b)]
  
  
  # results - default if no test is done
  results <- tibble(
    species = lipid,
    comparison = paste(name_a, "-", name_b),
    compare_a = name_a,
    compare_b = name_b,
    estimate_a = NA,
    estimate_b = NA,
    difference = NA,
    ci1 = NA,
    ci2 = NA,
    p_val = NA,
    test = NA
  )
  
  # t-tests, five cases:
  
  variance_ab <- var(c(lipid_a, lipid_b), na.rm = T)
  
  # 0): self-comparison: no test. Still included in results.
  
  if(name_a == name_b) {
    
    results$estimate_a <- mean(lipid_a, na.rm = T)
    results$estimate_b <- mean(lipid_b, na.rm = T)
    results$difference <- 0
    
  }
  
  # 1) too few values for testing
  
  else if (is.na(variance_ab) || variance_ab == 0) return(results)
  
  # 2) lipid detected in both cell lines in the organelle - 
  #    perform two-sample t test.
  
  else if(length(lipid_a)>1 & length(lipid_b)>1){
    t.test1 <- t.test(lipid_b,lipid_a) 
    
    results$p_val <- t.test1$p.value
    results$estimate_a = t.test1$estimate[2]
    results$estimate_b = t.test1$estimate[1]
    results$difference = diff(t.test1$estimate)
    results$ci1 = t.test1$conf.int[1]
    results$ci2 = t.test1$conf.int[2]
    results$test = "two.samples"
    
  }
  
  # 3) lipid detected in cell a only 
  #    perform one-sample t-test against test_val
  else if(var(lipid_a) != 0 & length(lipid_a)>1){
    
    # one-sample t-test, one-tailed, greater than "test_val"
    # difference estimate is positive to keep "name_a - name_b format"
    t.test1 <- t.test(lipid_a, alternative = "greater" , 
                      mu = test_val) 
    
    results$p_val = t.test1$p.value
    results$estimate_a = t.test1$estimate
    results$difference = max(t.test1$estimate-test_val,0)
    results$ci1 = -t.test1$conf.int[1]
    results$ci2 = -t.test1$conf.int[2]
    results$test = "one.sample"
    
  }
  
  # 4) lipid detected in cell b only 
  #    perform one-sample t-test against test_val
  #    difference estimate is negative to keep "lipid_a- lipid_b format"
  else if(var(lipid_b) != 0 & length(lipid_b)>1){
    
    # one-sample t-test, one-tailed, greater than "test_val"
    t.test1<-t.test(lipid_b, 
                    alternative = "greater", 
                    mu = test_val) 
    
    results$p_val = t.test1$p.value
    results$estimate_b = t.test1$estimate
    results$difference =  -max(t.test1$estimate-test_val,0)
    results$ci1 = t.test1$conf.int[1]
    results$ci2 = t.test1$conf.int[2]
    results$test = "one.sample"
    
  }
  
  return(results)
  
}

# ANOVA ------------------------------------------------------

# Run ANOVA comparing features across organelles

anova_test <- function(melt_classes, 
                       mean_length_db_classes,
                       organelles_for_analysis,
                       lipidQuan_output_folder,
                       compare,
                       exp_id) {
  
  cat("* ANOVA tests.\n")
  
  # for adding stars to figures later
  add_p_star <- function(x) {
    
    ifelse(x<0.001, "*** ", ifelse(x<0.01, "** ", ifelse(x<0.05, "*", " ")))
    
  }
  
  anova_summarise <- function(anova_df, 
                              compare, 
                              lipidQuan_output_folder,
                              compare_property,
                              exp_id) {
  
  # Stop if one or more classes result in errors
  if(any(unlist(lapply(anova_df$model, \(x) class(x)[1] == "try-error")))) {
    stop("error in ANOVA test")
  }
  
  # summary of all ANOVAs
  anova_df$summary <- lapply(anova_df$model, FUN = function(x) try(summary(x)))
  
  # Get p-values of all ANOVAs
  anova_df$p_vals  <- lapply(anova_df$summary, FUN = function(x) return(x[[1]][["Pr(>F)"]][1]))%>%
    unlist()
  
  # Correct ANOVA tests using the "holm" correction. 
  anova_df$p_val_correct <- anova_df$p_vals%>%
    p.adjust(method = "holm")
  
  # add column with stars for plotting
  anova_df$sign_holm <- add_p_star(anova_df$p_val_correct)
  
  # add column with stars without p value adjustment. 
  anova_df$sign_uncorrected <- add_p_star(anova_df$p_vals)
  
  # write to csv
  write.csv(
    anova_df[c(
      "class",
      "p_vals",
      "p_val_correct",
      "sign_holm",
      "sign_uncorrected"
    )],
    file = file.path(lipidQuan_output_folder, 
                     paste0("Results/ANOVA_",compare,"_",exp_id,".csv")
                     ),
    row.names = F
  )
  
}

### class 

anova_class_df <- melt_classes %>%
  subset(organelle %in% c(organelles_for_analysis, "Whole cells")) %>%
  mutate(organelle = droplevels(organelle)) %>%
  group_by(class) %>%
  do(model = try(aov(log2_molpct_class ~ organelle, data = .), silent = F)) 

anova_summarise(anova_df = anova_class_df, compare_property = "class",
                lipidQuan_output_folder = lipidQuan_output_folder,
                compare = compare,
                exp_id = exp_id)

### Mean double bonds 

anova_db <- mean_length_db_classes %>%
  subset(!organelle %in% c("Whole cells", "LMF"))%>%
  subset(class %in% c("PCO-", "PC","PS", "PEO-", "PI", "PE", "SM")) %>%
  group_by(class) %>%
  do(model = try(aov(mean_db ~ organelle, data = .)))

anova_summarise(anova_df = anova_db, compare_property = "db",
                lipidQuan_output_folder = lipidQuan_output_folder,
                compare = compare,
                exp_id = exp_id)

### Mean length 

anova_l <- mean_length_db_classes %>%
  subset(!organelle %in% c("LMF", "Whole cells")) %>%
  subset(class %in% c("PCO-", "PC","PS", "PEO-", "PI", "PE", "SM")) %>%
  group_by(class) %>%
  do(model = try(aov(mean_length ~ organelle, data = .)))

anova_summarise(anova_df = anova_l, compare_property = "length",
                lipidQuan_output_folder = lipidQuan_output_folder,
                compare = compare,
                exp_id = exp_id)

}

# PCA ----

organelle_pca <- function(tableau_data, 
                          tableau_samples, 
                          na_replace, 
                          organelles_for_analysis){
  
  # Convert tablau data to matrix
  pca_data <- tableau_data %>%
    column_to_rownames("species")  %>%
    as.matrix()
  
  cat("* PCA:",nrow(pca_data), "species,", ncol(pca_data), "samples", "\n")
  
  # Filter columns
  if(!all(tableau_samples$organelle_exp %in% colnames(pca_data))) stop("tableau_samples and tableau_log2 do not match")
  keep_cols <- tableau_samples[tableau_samples$organelle %in% organelles_for_analysis,]$organelle_exp
  pca_data <- pca_data[,colnames(pca_data) %in% keep_cols]
  #  PCA_samples <- tableau_samples[tableau_samples$organelle_exp %in% keep_cols,]
  
  # Replace NAs or low values with na_replace
  pca_data[is.na(pca_data) | pca_data < na_replace] <- na_replace
  
  # Remove 0-var columns
  pca_data <- pca_data[!apply(pca_data,1, \(x) var(x,na.rm = T)==0),]
  
  # Calculate PCA and organize as tibble
  pca <- prcomp(t(pca_data), scale. = T, center = T)
  
  # organize as tibble
  pca_df <- pca$x %>%
    as_tibble(rownames = "organelle_exp") %>%
    right_join(tableau_samples, ., by = join_by(organelle_exp)) %>%
    CorrectOrganelle(order_organelles = order_organelles)
  
  return(pca_df)
  
}

# Main ----

analyse_merged_data <- function(lipidQuan_output_folder, 
                                organelles_for_analysis,
                                output_file,
                                test_val,
                                compare,
                                lipid_level,
                                value_col,
                                p_adjust_method,
                                pca_na_replace,
                                exp_id,
                                anova_test){

  ## load quantified and filtered data
  load(file.path(lipidQuan_output_folder,
                 "Results", 
                 paste0("PAPSL_",exp_id, ".rdata")))
  
  tTest <- run_all_t_tests(df = melt_species,
                           organelles_for_analysis = organelles_for_analysis,
                           compare = compare,
                           lipid_level = lipid_level,
                           value_col = value_col, 
                           test_val = log2(0.1), 
                           p_adjust_method = "BH")
  
  pca_df <- organelle_pca(tableau_data = tableau_log2,
                          tableau_samples = tableau_samples,
                          na_replace = pca_na_replace,
                          organelles_for_analysis = organelles_for_analysis)
  
  if(anova_test){
  anova_test(melt_classes = melt_classes, 
             mean_length_db_classes = mean_length_db_classes,
             organelles_for_analysis = organelles_for_analysis,
             lipidQuan_output_folder = lipidQuan_output_folder,
             compare = compare,
             exp_id = exp_id)
  }
  
  save_to <- file.path(lipidQuan_output_folder, "Results",output_file)
  
  save(list = c("pca_df",
                "tTest"),
       file = save_to
  )
  
  cat("Saved as:", save_to, "\n")
  
}
