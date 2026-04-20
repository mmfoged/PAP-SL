# Geometric functions ----

#' Calculate geometric mean
#'
#' Log-transforms data, calculates arithmetic mean, then back-transforms to non-log.
#'
#' @param x Numeric vector of values.
#' @param na.rm Logical. Remove NAs?
#' @return Geometric mean.
geom_mean <- function(x, na.rm = FALSE) {
  exp(mean(log(x), na.rm = na.rm))
}

#' Calculate mean and standard deviation based on log-transformed data
#'
#' First log-transforms data, then calculates mean and standard deviation, 
#' and then back-transforms mean and standard deviation to non-log. 
#' Based on Hmisc::smean.sdl.
#'
#'@param x a numeric vector, see Hmisc::smean.sdl.
#'@param mult see Hmisc::smean.sdl 
#'@param na.rm should NA values be removed? Note that default is na.rm = T
#'@param zero.to.na Should 0-values be converted to NA? 
geom_mean_sdl<-function(x, mult = 1 , na.rm = T, zero.to.na = T){
  
  if(zero.to.na) {
    x[x==0] <- NA
  }
  
  log_x <- log(x)
  
  mean_x <- mean(log_x, na.rm = na.rm)
  
  sd_x = sd(log_x, na.rm = na.rm) 
  
  out1 <- data.frame(
    y = mean_x,
    ymin = mean_x - sd_x * mult,
    ymax = mean_x + sd_x * mult
  )
  
  out1 <- exp(out1) 
  
  return(out1)
  
}

# Lipid info ----

#' Get lipid class name from species name
#'
#' This function removes information on double bonds and OH from species names,
#' thus providing the class names.
#'
#' @param x vector of lipid species names
#' @return vector containing class names
#' @export

lipid_class_name <- function(x){
  
  x %>%
    str_remove(pattern = "^ ") %>%
    str_replace(pattern = " O-", "O- ") %>%
    str_remove(pattern =" \\d.*") %>%
    str_remove(pattern = "[[:space:]]\\:?$")
  
}

#' Get lipid cat names from class names
#'
#' This functions returns cat names from class names
#'
#' @param x vector of lipid class names
#' @return vector containing category names
#' @export
lipid_cat_name <- function(x){
  cat_convert <- list(
    data.frame(class = c("PA", "PE", "PI", "PS", "PG", "CL", "PC", "PGP", "BMP", "PCD13", "PED4", "PCD4", "D4PC", "D4PE", "D13PC", "13C215NPE", "13C315NPS"), 
               cat = "GPL"),
    data.frame(class = c("LPA", "LPE", "LPI", "LPS", "LPG", "LPC", "LPCD13", "LPCD4", "LPED4", "D13LPC", "D4LPC", "D4LPE", "13C215NLPE", "D4LPE", "D4LPC"),
               cat = "LGPL"),
    data.frame(class = c("PAO-", "PEO-", "PIO-", "PSO-", "PGO-","PEP-", "PCP-", "PAP-", "PIP-", "PSP-", "PGP-", "PCO-","PA_O", "PC_O", "PE_O", "PG_O", "PI_O", "PCOD4", "PCOD13", "D13PCO-", "D4PEO-","13C215NPEO-",  "PEOD4", "D4PCO-", "PAF", "PAFD13"),
               cat = "GPLO-"),
    data.frame(class =c("LPAO-", "LPAP-","LPEP-", "LPCP-", "LPIP-", "LPSP-", "LPEO-", "LPIO-", "LPSO-", "LPGO-", "LPCO-", "LPCOD13", "LPCOD4", "LPEOD4", "D4LPEO-", "D13LPCO-", "D4LPCO-"), 
               cat = "LGPLO-"),
    data.frame(class = c("Cer", "D4SM","HexCer", "LHexCer","diHexCer", "triHexCer", "CerP", "SM", "SMD13", "SHexCer", "GM1", "GM2", "GM3", "LCB","LSM", "LCBP", "DihLCB","13C215NtriHexCer",  "D13SM", "13C215NLCB","13C215NCer", "13C215NSM", "13C215NdiHexCer", "13C215NHexCer"),
               cat = "SL"),
    data.frame(class = c("Chol", "CE"),
               cat = "ST"),
    data.frame(class = c("DAG", "TAG","MAG", "TAGO-", "TAGP-", "DAGO-", "DAGP-","MAGO-", "MAGP-"),
               cat = "GL"),
    data.frame(class = c("FA", "FAlc", "FAld", "FAOH"), 
               cat = "FA")
  ) %>%
    do.call(rbind, .)
  
  
  out <- cat_convert$cat[match(x, cat_convert$class)]
  
  if(any(is.na(out))) {
    stop("Classes not recognized:", paste0(x[is.na(out)], sep = " "))
  }
  
  return(out)
}


# Helper function for lipid_info
# returns the info on chain length, double bonds and OH for each species
lipid_chain<-function(species){
  
  # Modify species vector before extracting information:
  species <- as.character(species)
  
  # remove space at beginning
  species <- stringr::str_remove(species, "^ ") 
  
  # remove space at end
  species <- stringr::str_remove(species, " $") 
  
  # Modify ether lipid names to have no space before 'O-'
  species <- stringr::str_replace(string = species,  pattern = " O-", replacement = "O- ") 
  
  # Make data frame
  df <- data.frame(species = species) 
  
  # Check that the species is defined by its sum composition and does not have to chains defined
  df$sumComp<- !grepl("\\d-\\d", species) 
  
  # extract information on chain lengths if the species has two chains defined
  df2 <- tidyr::extract(df, col = species, into = c("length1", "db1", "length2", "db2"), regex = " (.*):(.*)-(.*):(.*)", remove = F, convert = T) 
  
  # add length of chains
  df2$length <- as.numeric(df2$length1)+ as.numeric(df2$length2) 
  
  # add chain double-bonds
  df2$db = as.numeric(df2$db1)+ as.numeric(df2$db2) 
  
  # paste info on length and double bonds together in one column
  df2$sumSpecies <- with(df2, paste(length, ":", db, sep = "")) 
  
  # add info on length and double bonds from the species with sum compositions defined
  df2$sumSpecies[df2$sumComp] <-sub(pattern = ".*[[:space:]]",replacement =  "", x = df2$species[df2$sumComp]) 
  
  # If OH is not defined, add ';0' at end of sumSpecies
  df2$sumSpecies <- ifelse(grepl(";", df2$sumSpecies), df2$sumSpecies, paste(df2$sumSpecies, ";0", sep = "")) 
  
  return(df2$sumSpecies) 
}

#' Get information on lipid species
#'
#' This functions returns a data frame with information on acyl chain length,
#' double bonds, number of OH, class and category. The function accepts lipid
#' species with sum compositions or 2 acyl chains defined.
#'
#' @param species vector of lipid species names. Species must be in the format
#' X Y:Z;O, were X is the class name, Y is the acyl chain length, Z is the
#' number of acyl chain double bonds and O is the number of OH. Alternatively,
#' two acyl chains may be defined as X Y1:X1-Y2-X2. If to chains are defined,
#' the OH number cannot be included.
#'
#' @return data frame with columns 'class', 'length', 'db', 'OH' and 'cat'
#' @export
lipid_info<-function(species, suffix = ""){
  
  class <- lipid_class_name(species) # class name
  cat <- lipid_cat_name(class) # category name
  chain <- lipid_chain(species) # chain info
  
  #extract length, double bonds and OH from chain
  frame1<- tidyr::extract(data.frame(y = chain),y, c("length", "db", "OH"), regex = '(\\d*):(.*);(.*)', remove = T) # separate length, db and OH into columns
  frame2<-cbind(class, frame1, cat) #combine
  frame2$length[frame2$length==""]<-0 # change value to 0 when length is not defined
  frame2$db[frame2$db==""]<-0 # change value to 0 when db is not defined (for example for 'Chol :')
  
  #change columns to factors and numeric
  frame2$cat<-factor(frame2$cat, levels = unique(frame2$cat))
  frame2$length<-as.numeric(as.character(frame2$length))
  frame2$db <-as.numeric(as.character(frame2$db))
  frame2$OH <-as.numeric(as.character(frame2$OH))
  
  colnames(frame2) <- paste0(colnames(frame2), suffix)
  
  return(frame2)
}

# Order organelles / classes ----

#' order_organelles
#' @export
order_organelles <- c(
  "Whole cells",
  "LMF",
  "Mitochondria",
  "ER",
  "Golgi",
  "Cis-Golgi",
  "Trans-Golgi",
  "Plasma membrane",
  "Lysosomes",
  "Peroxisomes",
  "Ctrl",
  "blank"
)


#' order_class
#' @export
order_class<-c("Chol","CE",
               "LCB", "DihLCB","SM", "Cer", "CerP", "HexCer", "diHexCer", "triHexCer", "GM3", "GM2", "GM1", "SHexCer",
               "DAG", "TAG", "PA", "PC", "PI", "PE", "PS", "PG","CL",
               "LPA", "LPC", "LPI", "LPE", "LPS", "LPG",
               "PAO-", "PCO-", "PIO-", "PEO-", "PSO-", "PGO-",
               "LPAO-", "LPCO-", "LPIO-", "LPEO-", "LPSO-", "LPGO-", "FA", "FAlc", "MAGO-", "DAGO-", "TAGO-","ADAG", "PAF")

# Correct organelles /classes / cells ----

#' Correct organelle
#'
#' This function is a shortcut for organising organelles as factors and 
#' changing spelling errors in their names
#'
#' @param samples a data frame containing sample info
#' @export
CorrectOrganelle <- function(samples, order_organelles) {
  samples %>%
    mutate(
      organelle = case_when(
        grepl("cis", organelle, ignore.case = TRUE) ~ "Cis-Golgi",
        grepl("trans", organelle, ignore.case = TRUE) ~ "Trans-Golgi",
        grepl("mito", organelle, ignore.case = TRUE) ~ "Mitochondria",
        grepl("plasma", organelle, ignore.case = TRUE) ~ "Plasma membrane",
        grepl("blank|b0", organelle, ignore.case = TRUE) ~ "blank",
        organelle == "b" ~ "blank",
        grepl("whole", organelle, ignore.case = TRUE) ~ "Whole cells",
        TRUE ~ organelle
      ),
      organelle = fct_relevel(
        organelle,
        order_organelles[order_organelles %in% unique(organelle)]
      )
    )
}

#' Correct cell_organelle
#'
#' Shortcut to organise cell_organelle as factor in HKh-2 vs HCT116 dataset
#'
#' @param samples a data frame containing sample info
#' @export
CorrectCell_organelle<-function(samples){
  samples$cell_organelle<-as.character(with(samples, paste(cell, organelle, sep = "_")))
  levels1<-c(paste("HKH2", order_organelles, sep = "_"), paste("HCT116", order_organelles, sep = "_"), "none_blank")
  levels1<-levels1[levels1 %in% samples$cell_organelle]
  samples$cell_organelle<-factor(samples$cell_organelle, levels = levels1)
  return(samples)
}  

#' Correct class
#'
#' Shortcut to organise class as factor
#'
#' @param samples a data frame containing sample info
#' @export
CorrectClass<-function(samples){
  samples$class<-factor(samples$class, levels = order_class[order_class %in% unique(samples$class)])
  return(samples)
}

