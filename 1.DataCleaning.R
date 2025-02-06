#### 1. Clean Raw Data from Database and Shared Drive
# S. Graham
# 2/5/2025

#### Required files: ####
# 1. YieldData/BVT02T.dbf
# 2. YieldData/BVT03T.xlsx
# 3. YieldData/BVT04Td.xlsx
# 4. YieldData/BVT05T.dbf
# 5. YieldData/BVT06T.xlsx
# 6. YieldData/BVT08LincolnSidneyColbySummaryFinal.xls
# 7. YieldData/BVT09CO.xls
# 8. YieldData/BVT_Genovix.csv
# 9. YieldData/BDUP_Genovix.csv
# 10. YieldData/BDUP07-23_missing10.csv
# 11. Barley_NewNames.csv
# 12. YieldData/BDUP02T.dbf
# 13. YieldData/bdup03t.dbf
# 14. YieldData/BDUP06.dbf
# 15. YieldData/bdup04t.dbf
# 16. ListofGenos.csv

#### Output files: ####
# 1. BVTData_Cleaned.csv
# 2. BDUPData_Clean.csv
# 3. GenoSummary.csv
# 4. Cleaned_YieldTrials.csv

#Load Packages
library(foreign)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(remotes)
library(visdat)
library(naniar)
library(reshape2)

setwd("/Users/sydneygraham/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Statistics MS/Barley")

##### Define Functions #####
metadata <- function(df, c) { #data frame and column number with string
  df2 <- separate(df,c, into = c("Experiment", "Year", "Location"), sep = c(3,5))
  df2[,'Year'] <- as.numeric(df2[,'Year'])
  i<- 1
  for (i in 1:nrow(df)) {
    df2[i,'Year'] <- (df2[i,'Year'] + 2000)
    if ((df2[i, 'Location'] == "l") || (df2[i, 'Location'] == "L")) {
      df2[i, 'Location'] <- "Lincoln"
    }
    else if ((df2[i, c+2] == "s") || (df2[i, c+2] == "S")) {
      df2[i, 'Location'] <- "Sidney"
    }
    else if ((df2[i, c+2] == "c") || (df2[i, c+2] == "CO")) {
      df2[i, 'Location'] <- "Colby"
    }
    else if (df2[i, c+2] == "mc") {
      df2[i, 'Location'] <- "McCook"
    }
    else if (df2[i, c+2] == "M") {
      df2[i, 'Location'] <- "Mead"
    }
  }
  
  return(df2)
}

metadataDup <- function(df, c) { #data frame and column number with string
  df2 <- separate(df,c, into = c("Experiment", "Year", "Location"), sep = c(4,6))
  df2[,'Year'] <- as.numeric(df2[,'Year'])
  i<- 1
  for (i in 1:nrow(df)) {
    df2[i,'Year'] <- (df2[i,'Year'] + 2000)
    if ((df2[i, 'Location'] == "l") || (df2[i, 'Location'] == "L")|| (df2[i, 'Location'] == "LIN")) {
      df2[i, 'Location'] <- "Lincoln"
    }
    else if ((df2[i, c+2] == "s") || (df2[i, c+2] == "S") || (df2[i, c+2] == "SID")){
      df2[i, 'Location'] <- "Sidney"
    }
    else if ((df2[i, c+2] == "c") || (df2[i, c+2] == "CO")) {
      df2[i, 'Location'] <- "Colby"
    }
    else if (df2[i, c+2] == "mc") {
      df2[i, 'Location'] <- "McCook"
    }
    else if (df2[i, c+2] == "M") {
      df2[i, 'Location'] <- "Mead"
    }
  }
  
  return(df2)
}

standardize <- function(df) {
  columns_to_add <- c('Experiment.Name', 'Year', 'BLOC', 'Location', 'TPLOT', 'Name1', 'HDATEJULIA', 'HEIGHT', 'TESTWT', 'YLDBUA', 'Pedigree', 'IBLOC', 'Lodging', 'WinSur')
  
  # Convert column names to lowercase for case-insensitive comparison
  existing_columns <- tolower(names(df))
  
  for (col in columns_to_add) {
    if (!(tolower(col) %in% existing_columns)) {
      df[[col]] <- NA
    }
  }
  
  return(df)
}

gms2bushels <- function(df, column_name) {
  # Ensure the column exists in the dataframe
  if (!(column_name %in% colnames(df))) {
    stop("The specified column does not exist in the dataframe.")
  }
  
  # Ensure the column is numeric
  if (!is.numeric(df[[column_name]])) {
    stop("The specified column is not numeric.")
  }
  
  # Define constants
  constant1 <- 43560
  constant2 <- 50 #50 sqft plots
  constant3 <- 453.6
  constant4 <- 48 #48lbs/bushel for barley
  
  # Calculate the new values
  new_values <- round((df[[column_name]] * constant1) / (constant2 * constant3 * constant4), 1)
  
  # Update YLDBUA only where it is NA
  df$YLDBUA[is.na(df$YLDBUA)] <- new_values[is.na(df$YLDBUA)]
  
  return(df)
}

#### Read in & Clean BVT Data ######
# Need standard set of data
# Expt, year, location, rep, bloc, name, YLDBUA, Height, Heading Date (Julian), test weight
standard <- c('Experiment.Name', 'Year', 'BLOC', 'Location', 'TPLOT', 'Name1','HDATEJULIA', 'HEIGHT', 'TESTWT', 'YLDBUA', 'Pedigree','IBLOC','Lodging', 'WinSur','WINSUR')

BVT02 <- read.dbf("YieldData/BVT02T.dbf")
BVT02 <- metadata(BVT02, 1)
BVT02$HD <- (as.numeric(BVT02$HD) +121) #convert to julian, not leap year
colnames(BVT02) <- c("Experiment.Name", "Year", "Location", "Plot", "BLOC", "EntryNo", "ID", "Name1", "Pedigree", "PLOTOLD", "HDATEJULIA", "HEIGHT", "TESTWT", "Plot_WT", "YLDBUA")
BVT02$YLDBUA[BVT02$Location %in% c("Lincoln", "McCook")] <- NA #Yld data seems very wrong, will calculate from grams
BVT02 <- gms2bushels(BVT02,"Plot_WT")
BVT02 <- BVT02[,which(names(BVT02) %in% standard)]
BVT02 <- standardize(BVT02)

BVT03 <- read_excel("YieldData/BVT03T.xlsx")
BVT03 <- data.frame(BVT03)
BVT03 <- metadata(BVT03, 1)
BVT03$HD03 <- (as.numeric(BVT03$HD03) +121) #convert to julian, not leap year
colnames(BVT03) <- c("Experiment.Name","Year", "Location", "BLOC", "IBLOC", "Plot", "Marker", "EntryNo", "Previous.Nursery", "ENTRYOLD", "Name1", "Pedigree", "HDATEJULIA", "HEIGHT", "Plot_WT", "YLD_LBS", "TESTWT", "PLOTOLD")
BVT03$YLDBUA <- NA
BVT03 <- gms2bushels(BVT03,"Plot_WT")
BVT03 <- BVT03[,which(names(BVT03) %in% standard)]
BVT03 <- standardize(BVT03)

BVT04 <- read_excel("YieldData/BVT04Td.xlsx")
BVT04 <- data.frame(BVT04)
BVT04 <- metadata(BVT04, 1) #throws errors because of remmnant packets
BVT04$HD <- (as.numeric(BVT04$HD) +122) #convert to julian, LEAP YEAR
colnames(BVT04) <- c("Experiment.Name","Year", "Location", "BLOC", "IBLOC", "Plot", "EntryNo", "Previous.Nursery", "ENTRYOLD", "Name1", "Pedigree", "HDATEJULIA", "HEIGHT", "Plot_WT" , "YLD_LBS","TESTWT", "Comment")
BVT04$YLDBUA <- NA
BVT04 <- gms2bushels(BVT04,"Plot_WT")
BVT04 <- BVT04[,which(names(BVT04) %in% standard)]
BVT04 <- standardize(BVT04)

BVT05 <- read.dbf("YieldData/BVT05T.dbf")
BVT05 <- metadata(BVT05, 1)
BVT05$HD <- (as.numeric(BVT05$HD) +121) #convert to julian, not leap year
colnames(BVT05) <- c("Experiment.Name", "Year", "Location", "BLOC", "IBLOC", "Plot", "EntryNo", "ENTRYOLD", "Previous.Nuresery", "Name1", "Pedigree", "HDATEJULIA", "HEIGHT", "Plot_WT", "YLD_LBSCOR", "YLD_LBSNCOR", "MOIST", "TESTWT", "Lodging", "STAVG04", "STRANK04")
BVT05$YLDBUA <- NA
BVT05 <- gms2bushels(BVT05,"Plot_WT")
BVT05 <- BVT05[,which(names(BVT05) %in% standard)]
BVT05 <- standardize(BVT05)

BVT06 <- read_excel("YieldData/BVT06T.xlsx")
BVT06 <- data.frame(BVT06)
BVT06 <- metadata(BVT06, 1) #throws errors because of remmnant packets
BVT06$HD06 <- (as.numeric(BVT06$HD06) +121) #convert to julian, not leap year
colnames(BVT06) <- c("Experiment.Name", "Year", "Location", "BLOC", "IBLOC", "Plot", "EntryNo", "ENTRYOLD", "Name1", "HDATEJULIA", "HEIGHT", "Plot_WT", "YLD_LBS", "MOIST", "TESTWT", "Pedigree", "Previous.Nursery")
BVT06$YLDBUA <- NA
BVT06 <- gms2bushels(BVT06,"Plot_WT")
BVT06 <- BVT06[,which(names(BVT06) %in% standard)]
BVT06 <- standardize(BVT06)

BVT08 <- read_excel("YieldData/BVT08LincolnSidneyColbySummaryFinal.xls", sheet = "DataAll")
BVT08$EXPT <- "BVT"
colnames(BVT08) <- c("Experiment.Name", "Year", "Location", "BLOC", "IBLOC", "Plot", "EntryNo", "Previous.Nursery", "PlotOld", "Name1","Pedigree","WinSur", "HDATEJULIA", "HEIGHT", "Plot_WT", "YLD_LBS", "MOIST", "TESTWT")
BVT08$HDATEJULIA <- (as.numeric(BVT08$HDATEJULIA) +122) #convert to julian, LEAP YEAR
BVT08$YLDBUA <- NA
BVT08 <- gms2bushels(BVT08,"Plot_WT")
BVT08 <- BVT08[,which(names(BVT08) %in% standard)]
BVT08 <- data.frame(BVT08)
BVT08 <- standardize(BVT08)

BVT09 <- read_excel("YieldData/BVT09CO.xls")
BVT09 <- data.frame(BVT09)
BVT09 <- metadata(BVT09, 1)
colnames(BVT09) <- c("Experiment.Name", "Year", "Location", "BLOC", "IBLOC", "Order",  "EntryNo", "Plot", "Previous.Nursery", "PlotOld", "Name1","Pedigree", "HDATEJULIA", "HEIGHT", "Plot_WT", "MOIST", "TESTWT")
BVT09$YLDBUA <- NA
BVT09 <- gms2bushels(BVT09,"Plot_WT")
BVT09 <- BVT09[,which(names(BVT09) %in% standard)]
BVT09$HDATEJULIA <- (as.numeric(BVT09$HDATEJULIA) +121) #convert to julian, not leap year
BVT09 <- standardize(BVT09)

genovix <- read.csv("YieldData/BVT_Genovix.csv") #Data exported from genovix, unused columns and remmnant removed, NewNames added
BVT22 <- genovix[genovix$Year == '2022',] #2022 YLDBUA looks good
BVT22$Name1 <- BVT22$NEWNAME
BVT22 <- BVT22[,which(names(BVT22) %in% standard)]
BVT22 <- standardize(BVT22)

BVT21 <- genovix[genovix$Year == '2021',]
BVT21 <- gms2bushels(BVT21,"YLD_GM")
BVT21$Name1 <- BVT21$NEWNAME
BVT21 <- BVT21[,which(names(BVT21) %in% standard)]
BVT21 <- standardize(BVT21)

BVT20 <- genovix[genovix$Year == '2020',]
BVT20 <- gms2bushels(BVT20,"YLD_GM")
BVT20$Name1 <- BVT20$NEWNAME
BVT20 <- BVT20[,which(names(BVT20) %in% standard)]
BVT20 <- standardize(BVT20)

BVT19 <- genovix[genovix$Year == '2019',]
BVT19$YLDBUA <- BVT19$YLDLBSA/48
BVT19 <- gms2bushels(BVT19,"YLD_GM")
BVT19$Name1 <- BVT19$NEWNAME
BVT19 <- BVT19[,which(names(BVT19) %in% standard)]
BVT19 <- standardize(BVT19)

BVT18 <- genovix[genovix$Year == '2018',]
BVT18$Name1 <- BVT18$NAMENEW
BVT18$YLDBUA <- BVT18$YLDLBSA/48
BVT18 <- gms2bushels(BVT18,"YLD_GM")
BVT18 <- BVT18[,which(names(BVT18) %in% standard)]
BVT18 <- standardize(BVT18)

BVT17 <- genovix[genovix$Year == '2017',]
BVT17$Name1 <- BVT17$NEWNAME
BVT17 <- gms2bushels(BVT17,"YLD_GM") #the excluded erroneous Colby data
BVT17 <- BVT17[,which(names(BVT17) %in% standard)]
BVT17 <- standardize(BVT17)

BVT16 <- genovix[genovix$Year == '2016',]
BVT16$YLDBUA <- BVT16$YLDLBSA/48
BVT16 <- BVT16[,which(names(BVT16) %in% standard)]
BVT16 <- standardize(BVT16)

BVT15 <- genovix[genovix$Year == '2015',]
BVT15 <- gms2bushels(BVT15,"YLD_GM")
BVT15 <- BVT15[,which(names(BVT15) %in% standard)]
BVT15 <- standardize(BVT15)

BVT14 <- genovix[genovix$Year == '2014',]
BVT14$YLDBUA <- BVT14$YLDLBSA/48
BVT14 <- BVT14[,which(names(BVT14) %in% standard)]
BVT14 <- standardize(BVT14)

BVT13 <- genovix[genovix$Year == '2013',]
BVT13 <- gms2bushels(BVT13,"YLD_GM")
BVT13 <- BVT13[,which(names(BVT13) %in% standard)]
BVT13 <- standardize(BVT13)

BVT12 <- genovix[genovix$Year == '2012',]
BVT12 <- gms2bushels(BVT12,"YLD_GM")
BVT12 <- BVT12[,which(names(BVT12) %in% standard)]
BVT12 <- standardize(BVT12)

BVT11 <- genovix[genovix$Year == '2011',]
BVT11 <- gms2bushels(BVT11,"YIELD")
BVT11 <- BVT11[,which(names(BVT11) %in% standard)]
BVT11 <- standardize(BVT11)

BVT10 <- genovix[genovix$Year == '2010',]
BVT10 <- gms2bushels(BVT10,"YIELD")
BVT10$HDATEJULIA <- (as.numeric(BVT10$HDATE) +121) #convert to julian, not leap year
BVT10 <- BVT10[,which(names(BVT10) %in% standard)]
BVT10 <- standardize(BVT10)

#Do not know why so much data is missing for this year
BVT07 <- genovix[genovix$Year == '2007',]
BVT07$HDATEJULIA <- (as.numeric(BVT07$HDATE) +121) #convert to julian, not leap year
BVT07 <- BVT07[,which(names(BVT07) %in% standard)]
BVT07 <- standardize(BVT07)

#### Merge Data #######
dflist <- list(BVT02, BVT03, BVT04, BVT05, BVT06, BVT07, BVT08, BVT09, BVT10, BVT11, BVT12, BVT13, BVT14, BVT15, BVT16, BVT17, BVT18, BVT19, BVT20, BVT21, BVT22)
df <- dflist[[1]]

for (i in 2:length(dflist)) {
  to_merge <- dflist[[i]]
  common_columns <- intersect(names(df), names(to_merge))
  df <- merge(df, to_merge, by = common_columns, all = T)
}

#### Clean Data ####
summary(df)
df <- df %>% mutate_at(vars(HDATEJULIA:WinSur), ~ifelse(. %in% c(0, -9), NA, .)) 
df <- df %>% filter(YLDBUA < 179 & YLDBUA > 0) #must have yield data or is removed
df <- df %>% mutate(WINSUR = coalesce(WinSur, WINSUR))
df$WinSur <- NULL
df$HEIGHT[df$HEIGHT < 5] <- NA
df$HEIGHT[df$HEIGHT > 75] <- NA
df$TESTWT[df$TESTWT == 1.11] <- NA
df <- df %>% filter(!is.na(YLDBUA) | !is.na(HEIGHT) | !is.na(WINSUR) | !is.na(TESTWT) | !is.na(Lodging) | !is.na(WINSUR))
nrow(df)

df$Name1 <- gsub("NE", "NB", df$Name1)

df <- df %>%
  mutate(Name1 = case_when(
    Name1 == "NE95713" ~ "P-713",
    Name1 == "NE86954" ~ "P-954",
    Name1 == "NB99845" ~ "P-845",
    Name1 == "TAMBAR 50" ~ "TAMBAR 501",
    TRUE ~ Name1  # Keep unchanged if it doesn't match any condition
  ))

length(unique(df$Name1)) #324 genos
summary(df)

write.csv(df, "BVTData_Cleaned.csv")

#### Read in & Clean BDUP Data #####
bdup <- read.csv("YieldData/BDUP_Genovix.csv")
earlybdup <- read.csv("YieldData/BDUP07-23_missing10.csv")
earlybdup <- earlybdup %>% filter(Year >2007 & Year < 2019)
commoncols <- intersect(colnames(bdup), colnames(earlybdup))
BDUP <- merge(bdup, earlybdup, by = commoncols, all = T)

BDUP$Name1 <- ifelse(trimws(BDUP$NEWNAME) != "", BDUP$NEWNAME, BDUP$Name1)
BDUP <- gms2bushels(BDUP,"YLD_GM")
BDUP$YLDBUA[BDUP$YLDBUA < 0] <- NA

summary_by_year <- BDUP %>%
  group_by(Year) %>%
  summarise(
    CountYLD_GM = sum(!is.na(YLD_GM)),
    CountMOIST = sum(!is.na(MOIST)),
    CountHDATEJULIA = sum(!is.na(HDATEJULIA)),
    CountHEIGHT = sum(!is.na(HEIGHT)),
    CountYLDBUA = sum(!is.na(YLDBUA))
  )

print(summary_by_year) #Rm 2008, 2009, 2011, 2025 due to no data
BDUP <- BDUP %>% filter(!(Year %in% c(2008, 2009, 2011, 2025)))
BDUP <- BDUP %>% filter(Year != 2007) #due to data irregularities

#need new names for 2019-2022, 2016-2017
#have complete for 2022, 2019, 2016, 2017, 2020,2021
names <- read.csv("Barley_NewNames.csv")

bdup_merged <- merge(BDUP, names, by.x = "Name1", by.y = "CrossName", all.x = TRUE)
bdup_merged$Name1 <- ifelse(!is.na(bdup_merged$NewName), bdup_merged$NewName, bdup_merged$Name1) #some cross names all in BDUP2
bdup_merged$NewName <- NULL

BDUP <- bdup_merged

BDUP <- subset(BDUP, !grepl("BDUP2H21", Experiment.Name)) #Remove the BDUP2 from 2021 (and others?)
BDUP <- subset(BDUP, !grepl("BDUP2_22T", Experiment.Name)) 
BDUP <- subset(BDUP, !grepl("BDUP2_23", Experiment.Name))
BDUP <- subset(BDUP, !grepl("BDUP2_24", Experiment.Name))
BDUP <- BDUP[!grepl("NIN|RPN", BDUP$Experiment.Name), ]

#Assign BLOC by plot numbers
hist(BDUP$Plot) #Generally 401-445 == BLOC 1; 446-490 == BLOC 2
#in 2022: 201-245 == BLOC 1; 246-290 == BLOC 2
#In 2019: LN is normal, S (675-719) is all BLOC 1, MD 538-582 = BLOC 1, 583-627 = BLOC 2
BDUP$BLOC <- ifelse(
  is.na(BDUP$BLOC) | BDUP$BLOC == "",
  ifelse(
      (BDUP$Plot >= 401 & BDUP$Plot <= 445) |
      (BDUP$Plot >= 201 & BDUP$Plot <= 245) |
      (BDUP$Plot >= 675 & BDUP$Plot <= 719) |
      (BDUP$Plot >= 538 & BDUP$Plot <= 582), 1,
    ifelse(
      (BDUP$Plot >= 446 & BDUP$Plot <= 490) |
      (BDUP$Plot >= 246 & BDUP$Plot <= 290) |
      (BDUP$Plot >= 583 & BDUP$Plot <= 627), 2,
      NA  # Default to NA if no match
    )
  ),
  BDUP$BLOC #Keep BLOC if not empty
)

BDUP <- BDUP[,which(names(BDUP) %in% standard)]
BDUP <- standardize(BDUP)
summary(BDUP)

BDUP02 <- read.dbf("YieldData/BDUP02T.dbf")
BDUP02 <- metadataDup(BDUP02,1)
BDUP02$HD02 <- (as.numeric(BDUP02$HD02) +121) #convert to julian, not leap year
colnames(BDUP02) <- c("Experiment.Name", "Year", "Location", "Plot", "BLOC", "EntryNo", "ID", "Name1","Source", "Pedigree", "HDATEJULIA","HEIGHT", "Plot_WT", "YLDBUA", "HD01","HT01","YLDBUA01","LR01", "SR01")
BDUP02 <- gms2bushels(BDUP02, "Plot_WT")
BDUP02 <- BDUP02[,which(names(BDUP02) %in% standard)]
BDUP02 <- standardize(BDUP02)

BDUP03 <- read.dbf("YieldData/bdup03t.dbf")
BDUP03 <- metadataDup(BDUP03,1)
BDUP03$HD03 <- (as.numeric(BDUP03$HD03) +121) #convert to julian, not leap year
colnames(BDUP03) <- c("Experiment.Name", "Year", "Location", "BLOC","IBLOC","Plot", "EntryNo","Source","PLOTOLD", "Name1", "Pedigree", "HD02", "HT02", "YLDBUA02", "RANK02", "COMMENT02","HDATEJULIA","HEIGHT", "Plot_WT")
BDUP03$YLDBUA <- NA
BDUP03 <- gms2bushels(BDUP03, "Plot_WT")
BDUP03 <- BDUP03[,which(names(BDUP03) %in% standard)]
BDUP03 <- standardize(BDUP03)

BDUP06 <- read.dbf("YieldData/BDUP06.dbf")
BDUP06 <- metadataDup(BDUP06,1)
BDUP06$HD06 <- (as.numeric(BDUP06$HD06) +121) #convert to julian, not leap year
colnames(BDUP06) <- c("Experiment.Name", "Year", "Location", "BLOC","IBLOC","Plot", "EntryNo","Name1", "HDATEJULIA", "HEIGHT", "Plot_WT" ,"YLDOld", "SOURCE", "PLOTOLD", "YldO", "yld", "HD05", "HT05", "source04", "source03", "plot03", "plot04", "source02", "plot02", "cross", "Pedigree", "hdavg", "htavg", "yld02", "buc", "lba")
BDUP06$YLDBUA <- NA
BDUP06 <- gms2bushels(BDUP06, "Plot_WT")
BDUP06 <- BDUP06[,which(names(BDUP06) %in% standard)]
BDUP06 <- standardize(BDUP06)

#BDUP05 data dropped due to errors

BDUP04 <- read.dbf("YieldData/bdup04t.dbf")
BDUP04 <- metadataDup(BDUP04,1)
BDUP04$HD <- (as.numeric(BDUP04$HD) +122) #convert to julian, LEAP YEAR
colnames(BDUP04) <- c("Experiment.Name", "Year", "Location", "BLOC","IBLOC","Plot", "EntryNo","Source", "Name1","PLOTOLD", "Pedigree", "HDATEJULIA", "HEIGHT", "YLDBUA02","COMMENT02","Plot_WT","YLDLBSA", "HT03","YLD03")
BDUP04$YLDBUA <- NA
BDUP04 <- gms2bushels(BDUP04, "Plot_WT")
BDUP04 <- BDUP04[,which(names(BDUP04) %in% standard)]
BDUP04 <- standardize(BDUP04)

dflist <- list(BDUP02, BDUP03, BDUP04, BDUP06, BDUP)
df2 <- dflist[[1]]

for (i in 2:length(dflist)) {
  to_merge <- dflist[[i]]
  common_columns <- intersect(names(df2), names(to_merge))
  df2 <- merge(df2, to_merge, by = common_columns, all = T)
}

BDUP <- df2 %>%
  mutate(Name1 = case_when(
    Name1 == "NE95713" ~ "P-713",
    Name1 == "NE86954" ~ "P-954",
    Name1 == "NB99845" ~ "P-845",
    Name1 == "TAMBAR 50" ~ "TAMBAR 501",
    TRUE ~ Name1  # Keep unchanged if it doesn't match any condition
  ))

BDUP <- BDUP %>% filter(!is.na(YLDBUA) | !is.na(HEIGHT)| !is.na(TESTWT) | !is.na(Lodging))

length(unique(BDUP$Name1)) #705 in BDUP
BDUP$HDATEJULIA[BDUP$HDATEJULIA == -9] <- NA
summary(BDUP)
write.csv(BDUP, "BDUPData_Clean.csv")

#### Merge BVT and BDUP ####
common_columns <- intersect(names(BDUP), names(df))
all <- merge(BDUP, df, by = common_columns, all = T)

all <- all %>% filter(!is.na(YLDBUA) | !is.na(HEIGHT) | !is.na(WINSUR) | !is.na(TESTWT) | !is.na(Lodging) | !is.na(WINSUR))
all$Name1 <- gsub("NE", "NB", all$Name1)
all[all=="HitchcocK"] <- "Hitchcock"

genos <- read.csv("ListofGenos.csv")
all <- merge(all, genos, by = "Name1", all.x = TRUE)
all$env <- paste0(all$Year, all$Location)

summary_by_geno_all2 <- all %>%
  group_by(Name1) %>%
  summarise(
    RepsYield = sum(!is.na(YLDBUA)),
    RepsHeight = sum(!is.na(HEIGHT)),
    RepsHDate = sum(!is.na(HDATEJULIA)),
    RepsWinSur = sum(!is.na(WINSUR)),
    YearIn = min(Year),
    YearOut = max(Year),
    Duration = ((max(Year) - min(Year)) + 1),
    Genotyped = (sum(!is.na(Genotyped)) > 0),
    Envs = length(unique(env))
  ) %>%
  left_join(
    all %>%
      group_by(Name1, env) %>%
      summarise(
        RepsYieldEnv = sum(!is.na(YLDBUA))
      ) %>%
      pivot_wider(
        names_from = env,
        values_from = RepsYieldEnv,
        values_fill = 0
      ),
    by = "Name1"
  ) %>%
  select_if(~ !(is.numeric(.) && all(. == 0)))

write.csv(summary_by_geno_all2, "GenoSummary.csv")

geno_summary <- summary_by_geno_all %>% filter(Genotyped == TRUE) 
summary(geno_summary)

write.csv(file = "Cleaned_YieldTrials.csv", all)
