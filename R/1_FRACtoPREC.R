#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Moves Fragment columns to precursor columns, enabling quantification in 
# lipidQuan based on fragment ions rather than precursor ions.
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#' Move Fragment Intensities to Precursor Columns
#'
#' For specific lipid classes, copies values from FRAG1 columns to PREC columns.
#' This enables the downstream lipidQuan quantification to be based on precursor
#' values, which is more reliable for some classes,
#' 
moveFragToPrec <- function(lipidx_df, class_list, frag = "FRAG1") { 
  
  # Add internal standards to class list
  is_list <- paste("is", class_list, sep = "")
  class_list_is <- c(class_list, is_list)
  
  for (rows in 1:nrow(lipidx_df)) {
    
    # row class matches one of the specified classes, or their corresponding internal standards
    if (lipidx_df[rows,]$CLASS %in% class_list_is) { 
      
      # Define FRAG columns 
      FRAG_cols <- grepl(paste(frag,"\\:", sep = ""), colnames(lipidx_df))
      
      # move the frag cols of this this row to PREC          
      lipidx_df[rows, grepl("PREC\\:", colnames(lipidx_df))] <- 
        
        as.numeric(lipidx_df[rows, FRAG_cols]) 
      
    }
  }
  return(lipidx_df)
}

#' Calculate sums of unique precursor intensities of target classes
#' 
#' Processes lipidXplorer output df: Converts columns to numeric and aggregates 
#' rows by their species name for target classes. The latter is done to have
#' single unique values in the PREC column after moving FRAG to PREG columns,
#' which is required to avoid problems with lipidQuan quantification.

sumPrecUnique<-function(lipidx_df, class_list) {  
  
  # Conversion function 
  # Converts non-numerical values (indicated by "none" or "keine") to 0
  convert_none_keine <-
    function(x) {
      x <- ifelse(grepl("none|kein", x, ignore.case = T), 0, x)
      return(as.numeric(x))
    }
  
  # Apply conversion
  # Numerical columns are indicated by ":" in colnames
  intensity_columns <- grep(":", colnames(lipidx_df), value = T)
  
  lipidx_df[intensity_columns] <- lapply(lipidx_df[intensity_columns],
                                         convert_none_keine
  )
  
  # if no target classes in dataset return original df
  classes_in_dataset <- intersect(class_list,lipidx_df$CLASS)
  
  empty_rows <- lipidx_df$CLASS == "" | is.na(lipidx_df$CLASS)
  
  non_processed_classes <- 
    setdiff(unique(lipidx_df$CLASS[!empty_rows]), classes_in_dataset)
  
  if(length(classes_in_dataset)==0) {
    return(lipidx_df)
  }
  
  # separate non target classes from target classes
  non_target_df <- lipidx_df %>%
    filter(CLASS %in% non_processed_classes)
  
  # Process target classes
  calculate_unique_sum <- function(x) {
    if (is.numeric(x)) {return(sum(unique(x), na.rm = TRUE))} 
    else {return(x[1])}
  }
  
  target_df <- lapply(classes_in_dataset, FUN = function(class) {
    
    # Subset data for the class
    class_subset <- lipidx_df[lipidx_df$CLASS == class, ]
    
    # Aggregate several species rows into one (if there are more than one)
    class_subset_sum <- class_subset %>%
      aggregate.data.frame(
        .,
        by = list(NAME2 = .$NAME),
        FUN = calculate_unique_sum
      ) %>%
      select(-NAME2)
    
    return(class_subset_sum)
  }) %>%
    do.call(rbind, .)
  
  # Combine processed and unprocessed classes
  return(rbind(non_target_df, target_df))
  
}

#' FRAC to PREC main function
#'

frac_to_prec <- function(dir,class_list, target_folder = "Move FRAG1 to PREC") {
  
  counter_processed <- 0
  counter_unprocessed <- 0
  
  # Files in folder:
  files_in_folder <- list.files(path = file.path(dir, "unprocessed_output"), 
                                pattern = "-out")
  
  # Only use non-processed files
  non_processed <- files_in_folder[!grepl("^all|^process|FRAG1", files_in_folder)]  
  
  # Loop over all files
  for(file in non_processed){ # For each "-out" file 
    
    lipidx_df <- read.table(file=file.path(dir, "unprocessed_output",file), 
                            sep=",",
                            header = T, 
                            as.is = T, 
                            fill = T, 
                            row.names = NULL, 
                            check.names = F) 
    
    # Remove whitespace before class names
    lipidx_df$CLASS <- str_remove(string = lipidx_df$CLASS, pattern = "^ ")
    
    write_frac_to_prec_to <- file.path(dir, target_folder, paste0(file, "_FRAG1.csv"))
    
    # if the file contains any of the classes specified in class_list
    if(sum(class_list %in% lipidx_df$CLASS) > 0) { 
      
      counter_processed <- counter_processed + 1
      
      # Sum up FRAG intensities for species with several rows
      summed_precursor_df <- sumPrecUnique(lipidx_df, class_list = class_list) 
      
      #  Copy FRAG1 columns to PREC coumns 
      frac_to_prec_processed <- moveFragToPrec(summed_precursor_df, 
                                               class_list = class_list, 
                                               frag = "FRAG1")
      
      # Write processed files
      
      write.csv(frac_to_prec_processed, 
                file = write_frac_to_prec_to, 
                row.names = F, 
                na = "")
      
    }
    
    # Unprocessed files directly saved to the target folde
    else{
      
      counter_unprocessed <-   counter_unprocessed +1
      
      
      
      write.csv(lipidx_df, 
                file = write_frac_to_prec_to, 
                row.names = F, 
                na = "")
      
    }
    
    
    
  } 
  
  cat("Processed files:", counter_processed,"\n",
      "Not processed:",counter_unprocessed,"\n",
      "saved to",file.path(dir, target_folder),"\n") 
}
