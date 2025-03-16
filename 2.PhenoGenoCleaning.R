#### 2. Prep phenotypic & genotypic data for Genomic Selection
# S. Graham
# 2/5/2025

#### Required files: ####
# 1. Cleaned_YieldTrials.csv (1.DataCleaning.R)
# 2. BS4R8_24L_Corrected.csv
# 3. UNL2023FrelsAdjusted-Samples.csv
# 4. UNL2024Frels-B3K_01-02-Samples.csv
# 5. UNL2023FrelsAdjusted.MorexV3.IDnum.vcf
# 6. UNL2024Frels-B3K_01-02.MorexV3.IDnum.vcf

#### Output files: ####
# 1. all_phenoV7.csv
# 2. geno_allV2.csv

#Load Packages
library(BGLR)
library(tidyverse)
library(asreml)
library(vcfR)
library(SoyNAM)
library(cowplot)

setwd("/Users/sydneygraham/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Statistics MS/Barley")

##### Define Functions #####
replaceNAwithMean <- function(mat){
  replaceNAwithMeanVec <- function(vec){
    mu <- mean(vec, na.rm=TRUE)
    vec[is.na(vec)] <- mu
    return(vec)
  }
  return(apply(mat, 2, replaceNAwithMeanVec))
}

#### Load Pheno Data ####
pheno <- read.csv("Cleaned_YieldTrials.csv")
pheno_obs <- read.csv("BS4R8_24L_Corrected.csv")
is.na(pheno_obs$YLDBUA) <- pheno_obs$YLDBUA > 140

ggplot(pheno_obs, aes(Plot.Column, Plot.Range, fill= Check)) + 
  geom_tile() 

#remove KANSAS, McCook, Oklahoma, Colby data
pheno <- pheno %>% filter(!Location %in% c("KANSAS", "Oklahoma", "McCook", "Colby"))

pheno <- pheno %>%
  mutate(Experiment.Name = ifelse(grepl("bvt", Experiment.Name, ignore.case = TRUE), 
                                  "BVT", 
                                  ifelse(grepl("bdup", Experiment.Name, ignore.case = TRUE), 
                                         "BDUP", 
                                         Experiment.Name)))

pheno$trial <- paste(pheno$env, pheno$Experiment.Name)
AllPhenos <- pheno

#### YIELD ####
#No yield data for BDUP: 2016LN, 2017S, 2018S, 2020M, 2023L #Limited Data for BVT: 2007M, 2007S, BDUP: 2015L, 2015M, 2016S, 2019L
trials_to_remove <- c("2016Lincoln BDUP", "2017Sidney BDUP", "2018Sidney BDUP", "2020Mead BDUP", "2023Lincoln BDUP", "2007Mead BVT", "2007Sidney BVT", "2015Lincoln BDUP", "2015Mead BDUP", "2016Sidney BDUP", "2019Lincoln BDUP", "2019Sidney BDUP")
pheno <- AllPhenos[!AllPhenos$trial %in% trials_to_remove, ]

#### Calculate phenotypic BLUPs ####
pheno <- pheno %>% filter(Experiment.Name == "BVT")
trials <- unique(pheno$trial)
pheno <- pheno[order(pheno$Experiment.Name, pheno$BLOC),]
factors <- c("Name1","trial", "Year", "Experiment.Name", "BLOC")
pheno[,factors] <- lapply(pheno[,factors] , factor)

#### Visualization of Raw Data ####
pheno$fake1 <- "Within Environment"
pheno$fake2 <- "Across Environment"
pheno$YLDkgha <-  pheno$YLDBUA*53.8009
median(pheno$YLDkgha, na.rm = T)

pheno_obs$YLDkgha <-  pheno_obs$YLDBUA*53.8009
max(pheno_obs$YLDkgha, na.rm = T)


box <- ggplot(data = pheno, aes(x = env, y = YLDkgha, fill = Location)) + 
  geom_boxplot() +
  theme_grey() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.6)) +
  xlab("Environment") +
  ylab("Yield (kg/ha)") +
  scale_fill_manual(values=c("#1f6f6f", "#54a1a1", "#9fc8c8")) +
  facet_wrap(~fake1)

hist <- ggplot(data = pheno, aes(x = YLDkgha)) + 
  geom_histogram(position = "identity", binwidth = 200) + 
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  facet_wrap(~fake2) 

final <- cowplot::plot_grid(box, hist, 
                   ncol = 2, rel_widths = c(3, 1),
                   align = 'h', axis = 'l')  

save_plot("Yield_Dist.png", final, base_width = 7.25, base_height = 3.5)

WS <- pheno %>% filter(!(is.na(WINSUR)))
boxWS <- ggplot(data = WS, aes(x = env, y = WINSUR, fill = Location)) + 
  geom_boxplot() +
  theme_grey() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.6)) +
  xlab("Environment") +
  ylab("Winter Survival (%)") +
  scale_fill_manual(values=c("#1f6f6f", "#54a1a1", "#9fc8c8")) +
  facet_wrap(~fake1)

histWS <- ggplot(data = WS, aes(x = WINSUR)) + 
  geom_histogram(position = "identity", binwidth = 5) + 
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  facet_wrap(~fake2) 

finalWS <- cowplot::plot_grid(boxWS, histWS, 
                            ncol = 2, rel_widths = c(3, 1),
                            align = 'h', axis = 'l')  

save_plot("WinSur_Dist.png", finalWS, base_width = 7.25, base_height = 3.5)

#Create a list of DFs that will contain BLUPs
BLUP_list <- c()

for (i in 1:length(trials)) {
  temp <- pheno[pheno$trial == trials[i],]
  df_name <- paste(trials[i], "BLUP")
  if (length(unique(temp$BLOC)) > 1) {
      norm_bloc <- asreml(fixed = YLDBUA ~ BLOC + Name1, 
                          data = temp,
                          maxiter = 25,
                          na.action = na.method(y = "include", x = "include"))
      norm_bloc <- update.asreml(norm_bloc)
      
      BLUP <- predict(norm_bloc, classify = "Name1")$pvals 
      prds <- predict(norm_bloc, classify = "Name1", vcov = TRUE)
      vcov <- as.matrix(prds$vcov)
      
      # Handle Ailised estimates in Weight Calculation
      vcov <- data.frame(vcov)
      NA_positions <- which(is.na(vcov$X1))
      if (length(NA_positions > 0)) {
        vcov <- vcov[-NA_positions, -NA_positions]
      }
      weights <- diag(solve(vcov))
      v <- c(weights, NA_positions)
      weights <- replace(rep(NA, length(v)), !seq_along(v) %in% NA_positions, weights)
      #Add to the dataframe
      BLUP$w <- weights
      assign(df_name, BLUP)
      BLUP_list <- append(BLUP_list, df_name)
    } else {
    print(trials[i])
  }
  rm(df_name, temp, BLUP, norm_bloc)
}

first_df <- get(BLUP_list[1])[, c("Name1", "predicted.value", "w")]
first_df$ENV <- BLUP_list[1]  

blup_df <- first_df 
for (df_name in BLUP_list[-1]) {
  df <- get(df_name)[, c("Name1", "predicted.value", "w")]  
  df$ENV <- df_name
  blup_df <- rbind(blup_df, df) 
}

yld_blup <- blup_df
rm(first_df)

#### WINSUR ####
test <- AllPhenos %>%
  group_by(Experiment.Name, Year, Location) %>% summarise(ProportionWinSur = (sum(WINSUR, na.rm = T)/n()))

WinSur <- test %>% filter(ProportionWinSur > 10) #more than 10% WINSUR

trials_to_keep <- c("2022Mead BDUP", "2023Lincoln BDUP", "2008Lincoln BVT", "2011Lincoln BVT", "2014Mead BVT", "2015Lincoln BVT", "2015Mead BVT", "2016Mead BVT", "2022Mead BVT")
pheno <- AllPhenos[AllPhenos$trial %in% trials_to_keep, ]

trials <- unique(pheno$trial)
pheno <- pheno[order(pheno$Experiment.Name, pheno$BLOC),]
factors <- c("Name1","trial", "Year", "Experiment.Name", "BLOC")
pheno[,factors] <- lapply(pheno[,factors] , factor)
pheno <- pheno %>% filter(Experiment.Name == "BVT")

#Create a list of DFs that will contain BLUPs
BLUP_list <- c()

for (i in 1:length(trials)) {
  temp <- pheno[pheno$trial == trials[i],]
  df_name <- paste(trials[i], "BLUP")
  if (length(unique(temp$BLOC)) > 1) {
    norm_bloc <- asreml(fixed = WINSUR ~ BLOC + Name1, 
                        data = temp,
                        maxiter = 25,
                        na.action = na.method(y = "include", x = "include"))
    norm_bloc <- update.asreml(norm_bloc)
    
    BLUP <- predict(norm_bloc, classify = "Name1")$pvals 
    prds <- predict(norm_bloc, classify = "Name1", vcov = TRUE)
    vcov <- as.matrix(prds$vcov)
    
    # Handle Ailised estimates in Weight Calculation
    vcov <- data.frame(vcov)
    NA_positions <- which(is.na(vcov$X1))
    if (length(NA_positions > 0)) {
      vcov <- vcov[-NA_positions, -NA_positions]
    }
    weights <- diag(solve(vcov))
    v <- c(weights, NA_positions)
    weights <- replace(rep(NA, length(v)), !seq_along(v) %in% NA_positions, weights)
    #Add to the dataframe
    BLUP$w <- weights
    assign(df_name, BLUP)
    BLUP_list <- append(BLUP_list, df_name)
  } else {
    print(trials[i])
  }
  rm(df_name, temp, BLUP, norm_bloc)
}

first_df <- get(BLUP_list[1])[, c("Name1", "predicted.value", "w")]
first_df$ENV <- BLUP_list[1]  

blup_df <- first_df 
for (df_name in BLUP_list[-1]) {
  df <- get(df_name)[, c("Name1", "predicted.value", "w")]  
  df$ENV <- df_name
  blup_df <- rbind(blup_df, df) 
}

winsur_blup <- blup_df
rm(first_df)

#### Height ####
test <- AllPhenos %>%
  group_by(Experiment.Name, Year, Location) %>%
  summarise(ProportionHeight = (sum(HEIGHT, na.rm = T)/n()))

Height <- test
Height <- Height %>% filter(ProportionHeight < 10) #less than 10% Height

Height$trial <- paste(Height$Year, Height$Location, " ", Height$Experiment.Name, sep = "")

trials_to_drop <- Height$trial
pheno <- AllPhenos[!AllPhenos$trial %in% trials_to_drop, ]

trials <- unique(pheno$trial)
pheno <- pheno[order(pheno$Experiment.Name, pheno$BLOC),]
factors <- c("Name1","trial", "Year", "Experiment.Name", "BLOC")
pheno[,factors] <- lapply(pheno[,factors] , factor)
pheno <- pheno %>% mutate(HEIGHT = replace(HEIGHT, HEIGHT == 0, NA)) #Set Zeros to NA
pheno <- pheno %>% filter(Experiment.Name == "BVT")

#Create a list of DFs that will contain BLUPs
BLUP_list <- c()

for (i in 1:length(trials)) {
  temp <- pheno[pheno$trial == trials[i],]
  df_name <- paste(trials[i], "BLUP")
  if (length(unique(temp$BLOC)) > 1) {
    norm_bloc <- asreml(fixed = HEIGHT ~ BLOC + Name1, 
                        data = temp,
                        maxiter = 25,
                        na.action = na.method(y = "include", x = "include"))
    norm_bloc <- update.asreml(norm_bloc)
    
    BLUP <- predict(norm_bloc, classify = "Name1")$pvals 
    prds <- predict(norm_bloc, classify = "Name1", vcov = TRUE)
    vcov <- as.matrix(prds$vcov)
    
    # Handle Ailised estimates in Weight Calculation
    vcov <- data.frame(vcov)
    NA_positions <- which(is.na(vcov$X1))
    if (length(NA_positions > 0)) {
      vcov <- vcov[-NA_positions, -NA_positions]
    }
    weights <- diag(solve(vcov))
    v <- c(weights, NA_positions)
    weights <- replace(rep(NA, length(v)), !seq_along(v) %in% NA_positions, weights)
    #Add to the dataframe
    BLUP$w <- weights
    assign(df_name, BLUP)
    BLUP_list <- append(BLUP_list, df_name)
  } else {
    print(trials[i])
  }
  rm(df_name, temp, BLUP, norm_bloc)
}

first_df <- get(BLUP_list[1])[, c("Name1", "predicted.value", "w")]
first_df$ENV <- BLUP_list[1]  

blup_df <- first_df 
for (df_name in BLUP_list[-1]) {
  df <- get(df_name)[, c("Name1", "predicted.value", "w")]  
  df$ENV <- df_name
  blup_df <- rbind(blup_df, df) 
}

ht_blup <- blup_df
colnames(ht_blup) <- c("Name1", "Height", "W.Height", "ENV")
colnames(yld_blup) <- c("Name1", "Yield", "W.Yield", "ENV")
colnames(winsur_blup) <- c("Name1", "WinSur", "W.WinSur", "ENV")
rm(first_df)

BLUPS <- merge(ht_blup, yld_blup, by = c("Name1", "ENV"), all = T)
BLUPS <- merge(BLUPS, winsur_blup, by = c("Name1","ENV"), all = T)

#### Merge Sample ID ####
samples <- read.csv("/Users/sydneygraham/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Barley/Barley Genotyping/2023 3K Data/UNL2023Frels-B3K_Adjusted-selected/UNL2023FrelsAdjusted-Samples.csv")
samples <- samples %>% select(c("Sample.ID", "Sample.Name"))
samples <- samples %>%
  mutate(Sample.ID = case_when(
    Sample.Name == "P-954" ~ 23098,
    Sample.Name == "NB99845" ~ 23097,
    TRUE ~ Sample.ID  
  ))

#Use same sampleID for all duplicates, just selecting 1st one
samples <- samples %>%
  group_by(Sample.Name) %>%
  mutate(Sample.ID = first(Sample.ID)) %>%
  ungroup()

samples <- samples %>% distinct()

pheno_final <- merge(BLUPS, samples, by.x = "Name1", by.y = "Sample.Name", all.x = T)

geno_list <- unique(pheno_final$Sample.ID)
pheno_final <- pheno_final %>% filter(!Sample.ID %in% c(23270, 23236)) #Removed during geno filtering

pheno_final <- pheno_final[order(pheno_final$Sample.ID),]

#### Calculate BLUPs for BS4R8 w/ spatial adjustments ####
#yield
pheno_obs$Name1 <- as.factor(pheno_obs$Name1)
norm_col <- asreml(fixed = YLDBUA ~ Name1,
                   random = ~ Plot.Range,
                   data = pheno_obs,
                   maxiter = 25,
                   na.action = na.method(y = "include", x = "include"))
norm_col <- update.asreml(norm_col)
OBS <- predict(norm_col, classify = "Name1")$pvals 

OBS$ENV <- "2024Lincoln"
OBS$Trial <- "BS4R8"
OBS$Yield <- OBS$predicted.value
OBS$std.error <- NULL
OBS$predicted.value <- NULL

Obs_yld <- OBS

modH2 <- asreml(fixed = YLDBUA ~ 1,
                   random = ~ Name1 + Plot.Range,
                   data = pheno_obs,
                   maxiter = 25,
                   na.action = na.method(y = "include", x = "include"))
modH2 <- update.asreml(modH2)

vc.g <- summary(modH2)$varcomp['Name1','component']
vdBLUP.mat <- predict(modH2, classify="Name1", only="Name1", sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
h2 <- 1 - (vdBLUP.avg / 2 / vc.g)
h2 #heritability

#height
pheno_obs$Name1 <- as.factor(pheno_obs$Name1)
norm_col <- asreml(fixed = HEIGHT ~ Name1,
                   random = ~ Plot.Range,
                   data = pheno_obs,
                   maxiter = 25,
                   na.action = na.method(y = "include", x = "include"))
norm_col <- update.asreml(norm_col)
OBS <- predict(norm_col, classify = "Name1")$pvals 

OBS$ENV <- "2024Lincoln"
OBS$Trial <- "BS4R8"
OBS$Height <- OBS$predicted.value
OBS$std.error <- NULL
OBS$predicted.value <- NULL

Obs_ht <- OBS

#WinSur
pheno_obs$Name1 <- as.factor(pheno_obs$Name1)
pheno_obs$WINSUR <- pheno_obs$WINSUR * 10 #change to percent!
norm_col <- asreml(fixed = WINSUR ~ Name1,
                   random = ~ Plot.Range,
                   data = pheno_obs,
                   maxiter = 25,
                   na.action = na.method(y = "include", x = "include"))
norm_col <- update.asreml(norm_col)
OBS <- predict(norm_col, classify = "Name1")$pvals 

OBS$ENV <- "2024Lincoln"
OBS$Trial <- "BS4R8"
OBS$WinSur <- OBS$predicted.value
OBS$std.error <- NULL
OBS$predicted.value <- NULL

Obs_WinSur <- OBS

modH2 <- asreml(fixed = WINSUR ~ 1,
                random = ~ Name1 + Plot.Range,
                data = pheno_obs,
                maxiter = 25,
                na.action = na.method(y = "include", x = "include"))
modH2 <- update.asreml(modH2)

vc.g <- summary(modH2)$varcomp['Name1','component']
vdBLUP.mat <- predict(modH2, classify="Name1", only="Name1", sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
h2 <- 1 - (vdBLUP.avg / 2 / vc.g)
h2 #heritability

S4R8 <- merge(Obs_ht, Obs_yld, by = c("Name1", "status", "ENV", "Trial"), all = T)
S4R8 <- merge(S4R8, Obs_WinSur, by = c("Name1", "status", "ENV", "Trial"), all = T)
S4R8 <- S4R8 %>% filter(Name1 != "DUMMY")

samples_obs <- read.csv("/Users/sydneygraham/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Barley/Barley Genotyping/2024 3K Data/UNL2024Frels-B3K_01-02/UNL2024Frels-B3K_01-02-Samples.csv")
samples_obs <- samples_obs %>% select(c("Sample.ID", "Sample.Name"))
samples_obs <- samples_obs %>%
  mutate(Sample.ID = case_when(
    Sample.Name == "P-954" ~ 24078,
    Sample.Name == "NB15420" ~ 24098,
    TRUE ~ Sample.ID  
  ))

samples_obs <- samples_obs %>% distinct()

pheno_final_obs <- merge(S4R8, samples_obs, by.x = "Name1", by.y = "Sample.Name", all.y = T)

pheno_final_obs$status <- NULL

pheno_final_obs$W.Height <- NA
pheno_final_obs$W.Yield <- NA
pheno_final_obs$W.WinSur <- NA
pheno_final$Trial <- "BVT"
data <- rbind(pheno_final, pheno_final_obs)

write.csv(data, "all_phenoV7.csv")

#### Import Geno Data - Yield Trials ####

#VCF input code from Practical 2.R shared from Sheryl
infileVCF <- "/Users/sydneygraham/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Barley/Barley Genotyping/2023 3K Data/UNL2023Frels-B3K_Adjusted-selected/UNL2023FrelsAdjusted.MorexV3.IDnum.vcf"
vcf <- read.vcfR(infileVCF, verbose = FALSE)
impMethod <- "naive"

### Format, manipulate and filter genotype data ----
# Converting VCF file format to numerical matrix format that can be fit in statistical models
gt <- extract.gt(vcf, element = "GT", as.numeric = F)
fix_T <- as_tibble(getFIX(vcf))
gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
colnames(gt2) <- colnames(gt)
gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))
gt2b <- gsub("0[/|]0","0",gt2a)
gt2c <- gsub("[10][/|][10]","0.5",gt2b)
gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
gt2d_num <- apply(gt2d, 2, as.numeric)
#Adding row names back in
rownames(gt2d_num) <- rownames(gt2d)
geno_num_yt <- t(gt2d_num)

#### Import Geno Data - Obs ####
infileVCF <- "/Users/sydneygraham/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Barley/Barley Genotyping/2024 3K Data/UNL2024Frels-B3K_01-02/UNL2024Frels-B3K_01-02.MorexV3.IDnum.vcf"
vcf <- read.vcfR(infileVCF, verbose = FALSE)
impMethod <- "naive"

### Format, manipulate and filter genotype data ----
# Converting VCF file format to numerical matrix format that can be fit in statistical models
gt <- extract.gt(vcf, element = "GT", as.numeric = F)
fix_T <- as_tibble(getFIX(vcf))
gt2 <- matrix(0, ncol = ncol(gt), nrow = nrow(gt))
colnames(gt2) <- colnames(gt)
gt2a <- apply(gt,2, function(x) gsub("1/1","1",x))
gt2b <- gsub("0[/|]0","0",gt2a)
gt2c <- gsub("[10][/|][10]","0.5",gt2b)
gt2d <- gsub("\\.[/|]\\.","NA",gt2c)
gt2d_num <- apply(gt2d, 2, as.numeric)
#Adding row names back in
rownames(gt2d_num) <- rownames(gt2d)
geno_num_obs <- t(gt2d_num)

#merge geno data
geno_num_obs_df <- as.data.frame(geno_num_obs)
geno_num_obs_df$Samples <- rownames(geno_num_obs_df)

geno_num_yt_df <- as.data.frame(geno_num_yt)
geno_num_yt_df$Samples <- rownames(geno_num_yt_df)

common_col_names <- intersect(colnames(geno_num_obs_df), colnames(geno_num_yt_df))
geno_num <- merge(geno_num_obs_df, geno_num_yt_df, by=common_col_names, all=TRUE)
rownames(geno_num) <- geno_num$Samples
geno_num <- as.matrix(geno_num)

# Filtering out markers on proportion of missing data (50% missing)
miss <- function(x){length(which(is.na(x)))}
mrkNA <- (apply(geno_num, MARGIN=2, FUN=miss))/dim(geno_num)[1]
ndx <- which(mrkNA > 0.5)

if (length(ndx)>0) geno_num2 <- geno_num[, -ndx] else geno_num2 <- geno_num

# Filtering out individuals on proportion of missing data (50% missing)
indNA <- (apply(geno_num2, MARGIN=1, FUN=miss))/dim(geno_num2)[2]
ndx2 <- which(indNA > 0.5)

if (length(ndx2)>0) geno_num3 <- geno_num2[-ndx2, ] else geno_num3 <- geno_num2 

# Filter markers based on MAF (0.05 < MAF < 0.95)
maf <- apply(geno_num3, MARGIN=2, FUN=mean, na.rm=T)
ndx3 <- which(maf<0.05 | maf>0.95) 

if (length(ndx3)>0) geno_num4 <- geno_num3[, -ndx3] else geno_num4 <- geno_num3
rownames <- rownames(geno_num4)
rownames2 <- gsub("^.*\\.","", rownames)
rownames(geno_num4) <- rownames2 #these rownames will match sample names

# Match phenotypic data to marker data
ndx4 <- match(rownames(geno_num4), data$Sample.ID)
ndxNA <- which(is.na(ndx4))
ndx5 <- ndx4[-ndxNA]

geno_num5 <- geno_num4[-ndxNA, ]

# Impute genotype data using either naive imputation or Markov chain implemented in the NAM package
rownames <- rownames(geno_num5)
geno_num6 <- apply(geno_num5, 2, as.numeric)
rownames(geno_num6) <- rownames
if (impMethod == "naive") geno_imp <- replaceNAwithMean(geno_num6)

geno_list <- rownames(geno_imp)
pheno2 <- data[data$Sample.ID %in% geno_list, ]

write.csv(geno_imp, "geno_allV2.csv") #376 x 3055

