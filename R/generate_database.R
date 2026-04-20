# Construct a list of possible lipid species:

# Sphingolipids ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Ceramides can have between 32 and 44 Cs, and between 0 and 3 DBs. 
# This is Based on SM from LIPID MAPS


Cers <- paste(sort(c(paste(seq(32,44,2),0, sep = ":"),
                     paste(seq(32,44,2),1, sep = ":"),
                     paste(seq(32,44,2),2, sep = ":"),
                     paste(seq(32,44,2),3, sep = ":"))
), ";2", sep = "")
Cers

Cs <- list()
for(class in c("Cer", "SM", "CerP", "HexCer", "diHexCer", "triHexCer", "GM3", 
               "GM2", "GM1", "SHexCer")){
  Cs[[class]]<-paste(class, Cers)
}

names(Cs) <- NULL

Cs <- do.call("c", Cs)

# FA chain-containing lipids:  ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# this list is based on LIPID MAPS LPC chain lenghts
# 14:0, 16:0-1, 18:0-3, 20:0-5, 22:0-6
LPGLs <- sort(c(
  paste(seq(14, 22, 2), 0, sep = ":"),
  paste(seq(16, 22, 2), 1, sep = ":"),
  paste(seq(18, 22, 2), 2, sep = ":"),
  paste(seq(18, 22, 2), 3, sep = ":"),
  paste(seq(20, 22, 2), 4, sep = ":"),
  paste(seq(20, 22, 2), 5, sep = ":"),
  paste(22, 6, sep = ":")
))

# Lyso lipids ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


lps <- list()
for (class in c("LPC", "LPE", "LPS", "LPG", "LPA", "LPI", "CE")) {
  lps[[class]] <- (paste(class, LPGLs, sep = " "))
  if (!class %in% "CE") {
    classo <- paste(class, "O-", sep = "") # also make ether version
    lps[[classo]] <- (paste(classo, LPGLs, sep = " "))
  }
}

names(lps) <- NULL

lps <- do.call("c", lps)

# Fatty acids ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
FAs <- list()
for (class in c("FA", "FAlc", "MAGO-", "MAG")) {
  FAs[[class]] <- (paste(class, LPGLs, sep = " "))
}
names(FAs) <- NULL
FAs <- do.call("c", FAs)

# Diacyl-lipids ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Calculate all possible sum compositions of to side chains
GPLs_all <- expand.grid(data.frame(a = LPGLs, b = LPGLs)) 

# Calculate sum compositions
GPLs <- GPLs_all%>%
  apply(1, FUN = function(x) {
    x1 <- (as.numeric(strsplit(x[1], split = ":")[[1]][1]) +
      as.numeric(strsplit(x[2], split = ":")[[1]][1]))
    x2 <- (as.numeric(strsplit(x[1], split = ":")[[1]][2]) +
      as.numeric(strsplit(x[2], split = ":")[[1]][2]))
    return(paste(x1, x2, sep = ":"))
  })

GPLs3 <- GPLs

GPLs <- unique(GPLs)

# Side chains defined

GPLs2 <- with(GPLs_all, paste(b, a, sep = "/"))

GPLs2 <- unique(GPLs2)

# Classes 

## Sum compositions
Gpl <- list()
for(class in c("PC", "PE", "PS", "PG", "PA", "PI", "PCO-", "PEO-", "PSO-", 
               "PGO-", "PAO-", "PIO-", "DAG", "PGP", "DAGO-")){
  Gpl[[class]] <- paste(class, GPLs)
}
names(Gpl) <- NULL
Gpl <- do.call("c", Gpl)

## side chains
Gpl2 <- list()
for(class in c("PC", "PE", "PS", "PG", "PA", "PI", "PCO-", "PEO-", "PSO-", 
               "PGO-", "PAO-", "PIO-", "DAG", "PGP")){
  Gpl2[[class]] <- paste(class, GPLs2)
}

names(Gpl2) <- NULL

Gpl2 <- do.call("c", Gpl2)

## all GPLs (non unique) for adding as a column to Gpl2
Gpl3 <- list()
for(class in c("PC", "PE", "PS", "PG", "PA", "PI", "PCO-", "PEO-", "PSO-", 
               "PGO-", "PAO-", "PIO-", "DAG", "PGP")){ 
  Gpl3[[class]] <- paste(class, GPLs3)
}
names(Gpl3) <- NULL
Gpl3 <- do.call("c", Gpl3)


# Triacyl ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
TAGs_all <- expand.grid(data.frame(a = LPGLs, b = LPGLs, c = LPGLs))
  
TAGs <-  apply(TAGs_all,1, 
                   FUN = function(x){
                   x1 <- as.numeric(strsplit(x[1], split = ":")[[1]][1])+as.numeric(strsplit(x[2], split = ":")[[1]][1])+as.numeric(strsplit(x[3], split = ":")[[1]][1])
                   x2 <- as.numeric(strsplit(x[1], split = ":")[[1]][2])+as.numeric(strsplit(x[2], split = ":")[[1]][2])+as.numeric(strsplit(x[3], split = ":")[[1]][2])
              return(paste(x1, x2, sep = ":"))
            })

TAGs3 <- TAGs # for adding as a column to TAGs2

TAGs <- unique(TAGs)

# TAGs Sum composition
TAGs <- c(paste("TAG", sort(TAGs)), paste("TAGO-", (TAGs)))

# TAGs side chains
TAGs2 <- with(TAGs_all, paste(b, a, c, sep = "/"))

TAGs2 <- c(paste("TAG", (TAGs2)), paste("TAGO-", (TAGs2)))

TAGs3 <- c(paste("TAG", (TAGs3)), paste("TAGO-", (TAGs3)))

# Tetracyl (CL) ----
# The order is very important here
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
CLs_all <- expand.grid(data.frame(a = LPGLs, b = LPGLs, c = LPGLs, d = LPGLs))

CLs <- apply(CLs_all,1, 
            FUN = function(x){
              x1 <- as.numeric(strsplit(x[1], split = ":")[[1]][1]) + as.numeric(strsplit(x[2], split = ":")[[1]][1]) + as.numeric(strsplit(x[3], split = ":")[[1]][1]) + as.numeric(strsplit(x[4], split = ":")[[1]][1])
              x2 <- as.numeric(strsplit(x[1], split = ":")[[1]][2]) + as.numeric(strsplit(x[2], split = ":")[[1]][2]) + as.numeric(strsplit(x[3], split = ":")[[1]][2]) + as.numeric(strsplit(x[4], split = ":")[[1]][2])
              return(paste(x1, x2, sep = ":"))
            })

CLs3 <- CLs # for adding as a column to CLs2

CLs <- unique(CLs)
CLs <- CLs[!grepl("^56|58|60|62|82|84|86|88", CLs)]
CLs <- sort(CLs)

CLs2 <- with(CLs_all, paste(b, a, c, d, sep = "/"))

CLs2 <- CLs2[!grepl("^56|58|60|62|82|84|86|88", CLs)]
CLs3 <- CLs3[!grepl("^56|58|60|62|82|84|86|88", CLs)]

CLs <- paste("CL", CLs)
CLs2 <- paste("CL", CLs2)
CLs3 <- paste("CL", CLs3)

# Combine all ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

database0 <- c("Chol :", lps, Cs, Gpl, Gpl2, TAGs, CLs, "PAF 18:0", FAs)

# Add isotopic species
database_isotopes <- unique(c(database0,
                   
                   str_replace(database0[grepl("^PC ", database0)], pattern = "PC", replacement = "PCD13"), #
                   str_replace(database0[grepl("^PC ", database0)], pattern = "PC", replacement = "PCD4"),
                   str_replace(database0[grepl("^PCO- ", database0)], pattern = "PCO-", replacement = "PCOD4"),
                   str_replace(database0[grepl("^PCO- ", database0)], pattern = "PCO-", replacement = "PCOD13"),
                   str_replace(database0[grepl("^PEO- ", database0)], pattern = "PEO-", replacement = "PEOD4"),
                   str_replace(database0[grepl("^PE ", database0)], pattern = "PE", replacement = "PED4"),
                   str_replace(database0[grepl("^LPC ", database0)], pattern = "LPC", replacement = "LPCD13"),
                   str_replace(database0[grepl("^LPC ", database0)], pattern = "LPC", replacement = "LPCD4"),
                   str_replace(database0[grepl("^LPCO- ", database0)], pattern = "LPCO-", replacement = "LPCOD4"),
                   str_replace(database0[grepl("^LPCO- ", database0)], pattern = "LPCO-", replacement = "LPCOD13"),
                   str_replace(database0[grepl("^LPEO- ", database0)], pattern = "LPEO-", replacement = "LPEOD4"),
                   str_replace(database0[grepl("^LPE ", database0)], pattern = "LPE", replacement = "LPED4"),
                   str_replace(database0[grepl("^SM ", database0)], pattern = "SM", replacement = "SMD13"),
                   str_replace(database0[grepl("^PAF ", database0)], pattern = "PAF", replacement = "PAFD13")
))

# write to file:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
write.csv(unique(database_isotopes), "database.LIPIDMAPS_xxxx_xx_xx.csv", row.names = F)
