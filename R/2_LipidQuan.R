#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# LipidQuan
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

readFile<-function(dataPath){
  #
  # Auxiliary function that reads a data set by using the destination stored in the dataPath parameter.
  #
  data <- read.table(file=paste0(dataPath), sep=",",header = T, as.is = T, fill = T, row.names = NULL)
  
  return(data)
}

rmSpaceInBeginning<-function(data){
  
  #
  # Auxiliary function that removes the first letter in NAME and SPECIE, if the word begins with at [SPACE]
  #
  
  # for each row, check if either the NAME or SPECIE column begins with a [SPACE] and remove this [SPACE] if true.
  data[,"NAME"] <- ifelse(substring(data$NAME,1,1) == " ", substring(data$NAME, 2), data$NAME)
  data[, "SPECIE"] <- ifelse(substring(data$SPECIE,1,1) == " ", substring(data$SPECIE, 2), data$SPECIE)
  if(sum(is.na(data$NAME))>0){stop("NA-values not allowed in NAME column")}
  
  return(data)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Merge Files
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
mergeDataSets<-function(dataList, database, multiply = NULL, list = NULL){

  # Find all unique columns for all data sets
  uniqueCols <- vector("list", length(dataList))
  for(i in 1:length(dataList)){
    data<-readFile(dataList[i])
    uniqueCols[[i]] <- colnames(data)
  }
  uniqueCols <- unique(unlist(uniqueCols))
  
  
  # take all *FA* columns except SUMFA columns
  FA_cols <- uniqueCols[grep("FA",uniqueCols)]
  # remove SUMFA columns
  FA_cols <- FA_cols[-grep("^SUM",FA_cols)]
  
  
  # take all *FRAG* columns
  FRAG_cols <- uniqueCols[grep("FRAG",uniqueCols)]
  
  
  # take all *NLS* columns
  NLS_cols <- uniqueCols[grep("NLS",uniqueCols)]
  
  
  # merge *FA* *FRAG* and *NLS* together
  FA_FRAG_NLS_cols <- c(FA_cols, FRAG_cols, NLS_cols)
  
  # change long FAxxINTENS* name to FAxxINTENS_xx, where xx is a number for all the different column types (FA, FRAG and NLS)
  FA_FRAG_NLS_cols <- gsub("^(\\w+).*_(\\w+).raw", "\\1_\\2",FA_FRAG_NLS_cols)
  FA_FRAG_NLS_cols <- unique(FA_FRAG_NLS_cols)
  
  #### merge data sets together 
  firstRun <- TRUE
  for(dataPath in dataList){
    # load data from dataList
    data <- readFile(dataPath)
    
    # if a row in the NAME or SPECIE column starts with a [SPACE], remove this [SPACE]
    data <- rmSpaceInBeginning(data)
    
    
    # select specific columns for a given data set (columns that are always present in each dataset).
    selectedCols <- subset(data, select = c("ERROR", "CLASS", "LENGTH", "DB", "NAME", "SPECIE", "MASS"))
    # Insert values from OH col to selectedCols, if present in the data set or set value to NA.
    if("OH" %in% colnames(data)){
      selectedCols$OH <- data$OH
    }else{
      selectedCols$OH <- NA
    }
    
    
    # change long PREC* name to PREC_xx, where xx is a number
    PREC_tmp <- data[,c(colnames(data)[grep("^PREC",colnames(data))])]
    colnames(PREC_tmp) <- gsub("^(\\w+).*_(\\w+).raw", "\\1_\\2",colnames(PREC_tmp))
    
    # insert new PREC_xx col names to selectedCols
    selectedCols <- cbind(selectedCols,PREC_tmp)
    
    # insert mode column describing whether the measurements comes from POS or NEG measurements
    if(grepl("POS", basename(dataPath))){
      mode <- "POS"
    }else{
      if(grepl("NEG", basename(dataPath))){
        mode <- "NEG"
      }else{
        mode <- NA
      }
    }
    
    # insert mode column into selectedCols
    selectedCols$MODE <- mode
    
    
    # convert FA, FRAG and NLS col names in the input data (e.g. FAxxINTENS* name to FAxxINTENS_xx, where xx is a number), so they can be compared with the obtained colnames in FA_FRAG_NLS_cols.
    colnames(data) <- gsub("^(\\w+).*_(\\w+).raw", "\\1_\\2",colnames(data))
    
    
    # insert FA, FRAG and NLS columns into selectedCols. If they do not exist in the given data set, the value = NA
    dataNames <- colnames(data)
    for(col in FA_FRAG_NLS_cols){
      if(col %in% dataNames){
        selectedCols[,paste(col)] <- data[,paste(col)]  
      }else{
        selectedCols[,paste(col)] <- NA
      }
    }
    
    
    # remove all rows where NAME == ""
    
    
    selectedCols$NAME[is.na(selectedCols$NAME)]<-"0"
    
    selectedCols <- selectedCols[selectedCols$NAME != "",]
    
    
    # merge all data sets into one data set (mergedDataSet).
    if(firstRun == TRUE){ # only used on first run, since mergedDataSet is not defined yet
      mergedDataSet <- selectedCols
      
      firstRun <- FALSE
    }else{ # used for all runs except first run.
      mergedDataSet <- rbind(mergedDataSet, selectedCols)
    }
  }
  
  
  
  # convert different zero-symbols (e.g. " none" and NA) to the same symbol ("0").
  mergedDataSet[mergedDataSet == "0.0" | mergedDataSet == "none" | mergedDataSet == " none" | mergedDataSet == "None" | mergedDataSet == " None" | mergedDataSet == "Keine" | is.na(mergedDataSet)] <- "0"
  
  
  # multiply PREC* values if multiply is set to a value
  if(!is.null(multiply) && !is.null(list)){
    for(PREC in colnames(PREC_tmp)){
      mergedDataSet[,PREC] <- ifelse(mergedDataSet[,"NAME"] %in% list, mergedDataSet[,PREC]*multiply, mergedDataSet[,PREC])
    }
  }
  
  #### Filtering based on 1/0 columns in database
  
  # find all class names in database
  classNames <- unique(database[,"NAME"])
  
  
  return(mergedDataSet)
  
}



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# sort_is
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-



sort_is<-function(data){
  
  #
  # This function moves all internal standards (is) in a given data set so that they appear after all
  # normal species.
  #
  
  
  # save all is-rows
  isTmp <- data[grep("^is",data$NAME),]
  
  # remove all is-rows from data
  data <- data[-grep("^is",data$NAME),]
  
  # rbind all saved is-rows at the buttom of data
  data <- rbind(data,isTmp)
  
  return(data)
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# ID & Filtrations
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
filterDataSet<-function(data, database, filter = T){
  
  #
  # This function selects relevants columns from a data set remove species if the name does not exit
  # in the global ref. list.
  #
  
  # if a row in the NAME or SPECIE column in the database starts with a [SPACE], remove this row
  database <- rmSpaceInBeginning(database)
  
  #### select relevant columns
  data_tmp <- subset(data, select = c("NAME","ERROR", "MASS","SPECIE","MODE"))
  # View(data_tmp)
  PREC_tmp <- data[,c(colnames(data)[grep("^PREC",colnames(data))])]
  NLS_tmp <- data[,c(colnames(data)[grep("^NLS",colnames(data))])]
  
  FRAG_tmp <- data[,c(colnames(data)[grep("^FRAG",colnames(data))])]
  
  FA_tmp <- data[,c(colnames(data)[grep("^FA",colnames(data))])]
  FA_tmp <- FA_tmp[,-grep("^FA[0-9]$",colnames(FA_tmp))] # remove FA1, FA2, etc... 
  
  # remove FAO in the FA_tmp cols if it's present
  if(length(FA_tmp[,grep("^FAO$",colnames(FA_tmp))]) > 0 ){
    FA_tmp <- FA_tmp[,-grep("^FAO$",colnames(FA_tmp))]
  }
  
  
  data <- cbind(data_tmp,PREC_tmp, NLS_tmp)
  
  # create QUAN column to data consisting of the QUAN column in database.
  if (filter) {
    data$QUAN <- database$QUAN[match(data$NAME, database$NAME)]
  } else {
    data$QUAN <- "PREC"
  }
  
  # create SPECIE.GLOBAL column to data consisting of the SPECIE column in database.
  if (filter) {
    data$SPECIE.GLOBAL <- database$SPECIE[match(data$NAME, database$NAME)]
  } else {
    data$SPECIE.GLOBAL <- ""
  }
  # remove specie names that are not included in database
  if (filter) {
    GLOBAL.NAME.CHECK <- database$NAME[match(data$NAME, database$NAME)] # transfer NAME col from database to data
    data <- data[!is.na(GLOBAL.NAME.CHECK), ] # remove all rows whose names were not found in database
  }
  
  #### create SPECIE.ALL col: consists of all species within specie name seperated by "|", e.g. DAG 16:1-16:1|DAG 18:1-14:1
  data$SPECIE.ALL <- NA
  nameList <- unique(data$NAME)
  for(name in nameList){ # for each specie name, find all species and insert them into SPECIE.ALL seperated by "|"
    specie_tmp <- subset(data, NAME == name)$SPECIE
    specie_tmp <- paste(specie_tmp, collapse = '|')
    data[which(data$NAME == name),"SPECIE.ALL"]<-specie_tmp
  }
  
  
  #### remove duplicates
  
  # find potential value > 0 for each NAME in each PREC.* column and set this value as the default for this class name (if it exists), so that they appear after removal of duplicates. (PREC.* always have the same value > 0)
  classNames <- unique(data[,"NAME"])
  for(PREC in colnames(PREC_tmp)){
    for(className in classNames){
      data[data$NAME == className, PREC] <- max(data[data$NAME == className, PREC])
    }
  }
  
  # remove all duplicates of NAME.
  data <- data[!duplicated(data$NAME),]
  
  
  
  
  # (MAYBE DEPRECATED NOW THAT IT IS DONE IN THE mergeDataSet FUNCTION. CHECK SOON WITH RUnit) replace all "Keine" and "None" with the corresponding row value from the NAME column
  data$SPECIE <- ifelse(data$SPECIE == "Keine" | data$SPECIE == "None", data$NAME, data$SPECIE)
  
  
  # replace NA and "" with "none" in SPECIE.GLOBAL
  data$SPECIE.GLOBAL <- ifelse(is.na(data$SPECIE.GLOBAL) | data$SPECIE.GLOBAL == "", "none", data$SPECIE.GLOBAL)
  
  return(data)
}



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# pmol & Clean Up
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#PMOL CALC
pmolCalc <- function(data, database, spikeVariable, zeroThresh, blank.subtract=T, level.subt = 1){
  
  #
  # This function calculates pico mol (pmol) of species based on intensity from measurements
  # (target specie + internal standard) and known quantity of internal standard
  #
  
  #### if a row in the NAME or SPECIE column in database starts with a [SPACE], remove this row
  
  database <- rmSpaceInBeginning(database)
  
  # define PREC columns and BLNK column (last PREC column)
  PREC_names <- colnames(data)[grep("PREC",colnames(data))] # names of all PREC.* columns
  BLNK <- PREC_names[length(PREC_names)] # name of BLNK column (last PREC.* column)
  PREC_names <- PREC_names[-length(PREC_names)] # remove last column from PREC_names since this is BLNK
  
  #### zero adapter: check value >= 0. if not, value <- 0.
  for(PREC in PREC_names){ # for PREC columns
    data[,PREC] <- ifelse(data[,PREC] < 0, 0, data[,PREC])
  }
  # for BLNK column (equal to last PREC column)
  data[,BLNK] <- ifelse(data[,BLNK] < 0, 0, data[,BLNK])
  
  # split data set into exData (all experimental classes) and isData (all internal classes)
  exData <- data[-grep("^is",data$NAME),]
  isData <- data[grep("^is",data$NAME),]
  
  # find all class names (without number) in exData
  classNames <- gsub("^(\\w+.)[[:space:]].*", "\\1",exData[,"NAME"])
  isNames1 <- gsub("^(\\w+.)[[:space:]].*", "\\1",isData[,"NAME"])
  isDatabase <-database$NAME[grepl("^is", database$NAME)]
  
  cat("Number of classes:", length(unique(classNames)), "\n")
  
  missing_is_in_data <- setdiff(classNames, gsub("^is","", isNames1))
  
  missing_is_in_database <- setdiff(isData[,"NAME"], isDatabase)
  
  if(length(missing_is_in_data)> 0) {
    stop("Missing internal standards for classes: ", paste0(missing_is_in_data, sep = " "))
  }
  
  if(length(missing_is_in_database)> 0) {
    cat("IsDATA \n")
    print(isData$NAME)
    cat("isDatabase\n")
    print(isDatabase)
    stop("Missing database internal standard entries for classes: ", paste0(missing_is_in_database, sep = " "))
  }
  
  #### pmol calculation ( PREC*(NAME)/PREC*(isNAME)   x   pmol(isSpecie) )
  
  cat("Number of columns: ", length(PREC_names), "\n")
  for(PREC in PREC_names){
    
    for(i in 1:nrow(exData)){
      # find corresponding internal standard for the current class name.
      
      is <- isData[grep(paste0("is",classNames[i]," "),isData$NAME),]
      
      
      # pmol_isSpecie = spikeVariabel(uL) x [isLP]
      
      pmol_isSpecie <- spikeVariable * subset(database, NAME == is[,"NAME"])$isLP 
      pmol_calc <- exData[i,PREC] / is[,PREC] * pmol_isSpecie
      
      data[i,paste0("PMOL_",PREC)] <- pmol_calc 
      
    }
  }
  
  
  #### pmol BLNK calculation ( BLNK(NAME)/BLNK(isNAME)   x   pmol(isSpecie) )
  for(i in 1:nrow(exData)){
    
    # find corresponding internal standard.
    is <- isData[grep(paste0("is",classNames[i]," "),isData$NAME),]
    
    # pmol_isSpecie <- spikeVariabel(uL) x [isLP]
    
    pmol_isSpecie <- spikeVariable * subset(database, NAME == is[,"NAME"])$isLP
    
    # pmol calculation ( PREC:*(NAME)/PREC:*(isNAME) x pmol(isSpecie) )
    pmol_calc <- exData[i,BLNK] / is[,BLNK] * pmol_isSpecie
    data[i,paste0("PMOL_BLNK_",BLNK)] <- pmol_calc
  }
  
  # Possibility of defining a threshold

  #### subtract pmol BLNK from pmol PREC*
  PMOL_PREC_names <- colnames(data)[grep("^PMOL_PREC",colnames(data))]
  PMOL_BLNK <- colnames(data)[grep("^PMOL_BLNK",colnames(data))]
  
  # blank subtract start
  if (blank.subtract) {
    for (PMOL_PREC in PMOL_PREC_names) {
      x1 <- cbind(data[, PMOL_PREC], data[, PMOL_BLNK])
      x1[is.na(x1)] <- 0
      x2 <- blank.threshold(x1, level = level.subt)
      data[, paste0("SUBT_", PMOL_PREC)] <- x2
    }
  } else {
    for (PMOL_PREC in PMOL_PREC_names) {
      data[, paste0("SUBT_", PMOL_PREC)] <- data[, PMOL_PREC]
    }
    data[, paste0("SUBT_PMOL_PREC_", length(PMOL_PREC_names) + 1)] <- data[,
                                                                           PMOL_BLNK
    ]
  }
  
  #### zero adapter: check value >= 0. if not, value <- 0.
  SUBT_PMOL_PREC_names <- colnames(data)[grep("^SUBT_PMOL_PREC",colnames(data))]
  for(SUBT_PMOL_PREC in SUBT_PMOL_PREC_names){ # for PREC's
    data[,SUBT_PMOL_PREC] <- ifelse(data[,SUBT_PMOL_PREC] < 0, 0, data[,SUBT_PMOL_PREC])
  }
  
  # remove a given row if all PREC* (except last BLNK PREC) values contains zeros.
  data <- data[apply(data[PREC_names],1,function(value) any(value != 0)),]
  # update exData and isData with the newly created columns, and removed rows
  exData <- data[-grep("is",data$NAME),]
  isData <- data[grep("is",data$NAME),]
  
  
  
  # molpct calculation
  
  #### mol% specie calucated from PREC* values after BLNK subtraction
  for(SUBT_PMOL_PREC in SUBT_PMOL_PREC_names){
    
    # calculate mol% species for each specie
    sumSpecies <- sum(data[1:nrow(exData),SUBT_PMOL_PREC],na.rm = TRUE)
    data[1:nrow(exData),paste0("MOL_PCT_SPECIES_",SUBT_PMOL_PREC)] <- 100/sumSpecies*exData[1:nrow(exData),SUBT_PMOL_PREC]
    
  }
  
  
  #### zero adapter: check value >= zeroThresh. if not, value <- 0.
  MOL_PCT_SPECIES_SUBT_PMOL_PREC_names <- colnames(data)[grep("^MOL_PCT_SPECIES_SUBT_PMOL_PREC_",colnames(data))]
  
  # set all values that are under an user defined threshold (zeroThresh) to 0 for mol% species in exData.
  for(MOL_PCT_SPECIES_SUBT_PMOL_PREC in MOL_PCT_SPECIES_SUBT_PMOL_PREC_names){ # for PREC columns
    
    data[1:nrow(exData),MOL_PCT_SPECIES_SUBT_PMOL_PREC] <- ifelse(data[1:nrow(exData),MOL_PCT_SPECIES_SUBT_PMOL_PREC] < zeroThresh, 0, data[1:nrow(exData),MOL_PCT_SPECIES_SUBT_PMOL_PREC])
  }
  
  
  
  # re-calculation of the mol% specie after removal of rows below an user defined threshold (zeroThresh)
  for(MOL_PCT_SPECIES_SUBT_PMOL_PREC_re in MOL_PCT_SPECIES_SUBT_PMOL_PREC_names){
    sumSpecies<-sum(data[1:nrow(exData),MOL_PCT_SPECIES_SUBT_PMOL_PREC_re], na.rm = TRUE)
    data[1:nrow(exData),paste0("FILTERED_",MOL_PCT_SPECIES_SUBT_PMOL_PREC_re)] <- 100/sumSpecies*data[1:nrow(exData),MOL_PCT_SPECIES_SUBT_PMOL_PREC_re]
  }
  
  
  #### sum pmol values for each classes in each PREC* after BLNK subtraction
  # find all unique class names (without numbers)
  classNames <- gsub("^(\\w+.)[[:space:]].*", "\\1",data[1:nrow(exData),"NAME"])
  classNames <-unique(classNames)
  
  
  sumClassValueList <- matrix(numeric(), nrow = length(classNames), ncol = length(SUBT_PMOL_PREC_names)) # store all sumClassValue to be used later in mol% class calculation
  for(i in 1:length(classNames)){
    
    for(j in 1:length(SUBT_PMOL_PREC_names)){
      # sum values for each class
      sumClassValue <- sum(data[grep(paste0("^",classNames[i]," "),data$NAME),SUBT_PMOL_PREC_names[j]],na.rm = TRUE)
      sumClassValueList[i,j] <- sumClassValue
      
      # add sum value to all species with the same class
      data[grep(paste0("^",classNames[i]," "),data$NAME),paste0("CLASS_PMOL_",SUBT_PMOL_PREC_names[j])] <- sumClassValue
    }
  }
  
  
  #### mol% class (100/sum(sumClassValueList) * sumClassValue[j,i]) calculation
  # find sum of all class values for each SUBT_PMOL_PREC*
  for(j in 1:length(SUBT_PMOL_PREC_names)){
    totalSumClassValueList<-sum(sumClassValueList[,j],na.rm = TRUE)
    
    # calculate mol% class for each sumClassValueList col, which corresponds to each SUBT_PMOL_PREC
    for(i in 1:nrow(sumClassValueList)){
      mol_pct_class<-100/totalSumClassValueList * sumClassValueList[i,j]
      
      # insert mol% class calculation into the respective class in data 
      data[grep(paste0("^",classNames[i]," "),data$NAME),paste0("MOL_PCT_CLASS_",SUBT_PMOL_PREC_names[j])]<-mol_pct_class
      
    }
  }
  
  
  return(data)
  
}


# COMPACT OUT
compactOutput_pmolCalc <- function(data){
  
  #
  # This function saves a data.frame with only NAME, CLASS_PMOL_SUBT_PMOL_PREC*, MOL_PCT_CLASS_SUBT_PMOL_PREC* and CLASS FILTERED columns
  #
  
  # create new data frame with class names
  classPmol_molPctClass <- data.frame(NAME = sub(" .*", "", data$NAME))
  
  # insert CLASS_PMOL_SUBT_PMOL_PREC*, and MOL_PCT_CLASS_SUBT_PMOL_PREC* 
  # and CLASS_FILTERED* columns
  
  class_pmol <- grep("^CLASS_PMOL_SUBT_PMOL_PREC", colnames(data), value = TRUE)
  
  class_molpct <- grep("^MOL_PCT_CLASS_SUBT_PMOL_PREC", colnames(data), value = TRUE)
  
  classPmol_molPctClass <- cbind( classPmol_molPctClass, 
                                  data[class_pmol], 
                                  data[class_molpct])
  
  # Class columns have identical rows for lipids of the same class. 
  # Since species are left out, remove duplicates of NAME
  classPmol_molPctClass <- unique(classPmol_molPctClass)
  
  #### sum filtered values for each classes for each sample after BLNK subtraction
  
  classNames <- sub("^([^[:space:]]+).*", "\\1", data$NAME)
  
  FILTERED_columns <- data[grep("^FILTERED",colnames(data))] 
  
  colnames(FILTERED_columns) <- paste0("CLASS_", colnames(FILTERED_columns))
  
  class_filtered <- aggregate.data.frame(FILTERED_columns, 
                                         by = list(NAME = classNames), 
                                         FUN = function(x) sum(x, na.rm = TRUE))
  
  # Set the IS values to NA
  class_filtered[grep("^is", class_filtered$NAME),-(colnames(class_filtered) == "NAME")] <- NA
  
  # Add class_filtered to the output data frame
  classPmol_molPctClass <- merge(classPmol_molPctClass, 
                                 class_filtered, 
                                 by = "NAME", 
                                 sort = F)
  
  return(classPmol_molPctClass)
  
}


# MAKE TABLEAU OUTPUT

makeTableauOutput <- function(classPmol_molPctClass, pmolCalculatedDataSet, molpct.calc = T){
  
  #
  # This function creates an output file of the results in a format that can be used by Tableau.
  #
  
  # take all lipid species (NAME col) from pmolCalculatedDataSet without using the is-rows and their respective 
  if (molpct.calc) {
    lipidSpecies <- pmolCalculatedDataSet[,c(1,grep("^FILTERED", colnames(pmolCalculatedDataSet)))]
  } else {
    lipidSpecies <- pmolCalculatedDataSet[,c(1,grep("^SUBT", colnames(pmolCalculatedDataSet)))]
  }
  
  lipidSpecies <- lipidSpecies[-grep("^is",lipidSpecies$NAME),]
  
  # change colnames: NAME -> mol%, FILTERED* -> Sample_01
  sampleNames <- paste0("Sample_",1:(length(colnames(lipidSpecies))-1))
  
  if (molpct.calc) {
    colnames(lipidSpecies) <- c("mol%", sampleNames)
  } else {
    colnames(lipidSpecies) <- c("pmol", sampleNames)
  }
  
  # sum classes:  
  # take NAME from classPmol_molPctClass without using the is-rows
  classes <- classPmol_molPctClass[,c(1,grep("^CLASS_FILTERED*", colnames(classPmol_molPctClass)))]
  classes <- classes[-grep("^is",classes$NAME),]
  
  # change colnames: NAME -> mol%, FILTERED* -> Sample_01
  sampleNames <- paste0("Sample_",1:(length(colnames(classes))-1))
  colnames(classes) <- c("mol%", sampleNames)
  #sum classes end
  
  # merge data sets together
  
  if (molpct.calc) {
    tableauOutput <- rbind(lipidSpecies, classes)
  } else {
    tableauOutput <- lipidSpecies
  }
  
  # change classes/species with "O-": insert [SPACE] before "O-" and remove [SPACE] after "O-"
  namesWithOIndexes <- grep("O-",tableauOutput$"mol%")
  for(i in namesWithOIndexes){
    nameWithO <- tableauOutput$"mol%"[i]
    nameWithO <- strsplit(nameWithO," ")
    
    if (length(nameWithO[[1]]) == 2) {
      # used when name represents a specie (contains number specs.)
      nameWithO <- paste0(gsub("O-", " O-", nameWithO[[1]][1]), nameWithO[[1]][2])
    } else {
      # used when name represents a class (without number specs.)
      nameWithO <- paste0(gsub("O-", " O-", nameWithO[[1]][1]))
    }
    tableauOutput$"mol%"[i] <- nameWithO
  }
  
  return(tableauOutput)
}

# Single-replicate analysis

LipidQuanAnalysis <- function(
    dataList,
    database,
    spike_var,
    exp_prefix,
    out_exclude,
    lipidQuan_output_folder
) {
  
  # 1: Merge data
  mergedDataSets <- mergeDataSets(dataList = dataList, database = database)
  
  # 2: Move internal standards to the end
  sortedDataSets <- sort_is(mergedDataSets)
  
  # 3: Remove duplicates (do not filter with the database):
  filteredDataSet <- filterDataSet(sortedDataSets, database, filter = F)
  
  # 4: Calculate pmol:
  pmolCalculatedDataSet <- pmolCalc(
    data = filteredDataSet,
    database = database,
    spikeVariable = spike_var,
    zeroThresh = 0,
    blank.subtract = F,
    level.subt = 1
  )
  
  # 5: Compact output
  classPmol_molPctClass <- compactOutput_pmolCalc(pmolCalculatedDataSet)
  
  # 6: Make tableau output
  # Do not calculate molpct - do this after species filtering. 
  tableauOutput <- makeTableauOutput(
    classPmol_molPctClass,
    pmolCalculatedDataSet,
    molpct.calc = F
  )
  
  # 7: Save csv
  output_file <- file.path(
    lipidQuan_output_folder,
    paste0(exp_prefix, "_tableauOutput.csv")
  )
  
  write.csv(tableauOutput, file = output_file, quote = FALSE, row.names = FALSE)
  
  cat("saved to", output_file, "\n")
}

#  Main function ----

lipidQuanBatch <- function(database_file, 
                           exp_no, 
                           replicate_no, 
                           spike_var,
                           out_exclude,
                           lipidQuan_output_folder,
                           lipidX_output_folder) {
  
  # Load database
  database <- 
    read.table(database_file, stringsAsFactors = FALSE, header = TRUE, sep = ",")
  
  # vector of replicates
  
  if (replicate_no > 1) {
    replicate_vector <- paste(exp_no, "_", 1:replicate_no, sep = "")
  } else {
    replicate_vector <- exp_no
  }
  
  # loop over replicates
  for(i in replicate_vector){
    
    # print replicate name
    cat("\n=== Quantifying", i,"===", "\n")
    
    # Experiment name 
    exp_prefix <- i
    
    # List of data to be processed
    lipidx_files <- list.files(lipidX_output_folder)
    
    dataList <- lipidx_files[grepl(exp_prefix, lipidx_files) & grepl(".csv", lipidx_files)]
    
    # Print names of out-files that were not included in the analysis
    cat("exluded out-files:", dataList[grepl(out_exclude, dataList)], "\n")
    
    # Remove excluded files from list of data to be processed
    dataList <- dataList[!grepl(out_exclude, dataList)]
    # Run lipidQuan
    LipidQuanAnalysis(dataList = file.path(lipidX_output_folder,dataList), 
                      database = database, 
                      spike_var = spike_var,
                      exp_prefix = exp_prefix,
                      out_exclude = out_exclude, 
                      lipidQuan_output_folder = lipidQuan_output_folder)
    
  } 
  
}

