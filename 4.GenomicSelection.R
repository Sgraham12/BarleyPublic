#### 4. Genomic Selection for Yield & WinSur
# S. Graham
# 2/5/2025

#### Required files: ####
# 1. all_phenoV7.csv (2.PhenoGenoCleaning.R)
# 2. W_SummerOnly.csv
# 3. geno_allV2.csv (2.PhenoGenoCleaning.R)

#### Output files: ####
# 1. BGLR_Yield_Results.csv
# 2. ML_Yield_Results.csv
# 3. BGLR_WinSur_Results.csv
# 4. ML_WinSur_Results.csv

#Load Packages
library(corrr)
library(BGLR)
library(tidyverse)
library(vcfR)
library(Metrics)
library(caret)
library(purrr)
library(randomForest)
library(xgboost)

#setwd("~/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Statistics MS/Barley")

#### Define Functions ####
mean_imputation <- function(df){
  n_cols <- ncol(df)
  for(i in 1:n_cols){
    column_i <- df[, i]
    mean_i <- mean(column_i, na.rm = TRUE)
    NAs_i <- which(is.na(column_i))
    N_NAs <- length(NAs_i)
    column_i[NAs_i] <- mean_i
    df[, i] <- column_i
  }
  return(df)
}

#### Define GS Parameters ####
burnIn <- 3000
nIter <- 30000
reps <- 1

#### Yield ####
### Import all Pheno Data ###
rep_df <- data.frame()
for (i in 1:8) {
  pheno <- read.csv("all_phenoV7.csv")
  pheno <- pheno %>%
    mutate(Expt = case_when(
      grepl("BDUP", Trial) ~ "BDUP",
      grepl("BVT", Trial) ~ "BVT",
      grepl("BS4R8", Trial) ~ "BS4R8",
      TRUE ~ NA_character_
    ))
  
  pheno <- pheno %>% filter(!is.na(Yield)) %>% filter(Expt != "BDUP") %>% filter (!is.na(Sample.ID))
  pheno$ENV <- as.factor(pheno$ENV)
  pheno <- pheno %>% select("Sample.ID", "ENV", "Yield")
  pheno$ENV <- gsub( " .*$", "", pheno$ENV)
  
  if (i == 2) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln"))
  } else if (i == 3) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 4) {
    pheno <- pheno %>% filter(str_detect(ENV, "202"))
  } else if (i == 5) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead"))
  } else if (i ==6) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 7) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln"))
  } else if (i == 8) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln")) %>% filter(str_detect(ENV, "202"))
  }
  envs <- length(unique(pheno$ENV))
  pheno$Sample.ID <- as.factor(pheno$Sample.ID)
  pheno_list <- unique(pheno$Sample.ID)
  pheno_list <- data.frame(pheno_list)
  
  #### Geno Data for GBLUP ####
  X <- read.csv("geno_allV2.csv")
  X <- X[order(X[, 1]), ]
  X_merge <- merge(pheno_list, X, by.x = "pheno_list", by.y = "X", all.x = T)
  
  geno_order <- X_merge[, 1] # Store ordered genotype sample names
  X_merge <- X_merge[, -1] # Remove the first column (sample names) from X
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X_merge)) > 0))
  X_merge <- X_merge[, !colnames(X_merge) %in% NAMkrs]
  
  # Scale the data (center and scale each column)
  X_merge <- scale(X_merge, center = TRUE, scale = TRUE)
  
  # Set row names to the ordered genotype sample names
  rownames(X_merge) <- geno_order
  X <- X_merge
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X)) > 0))
  X <- X[, !colnames(X) %in% NAMkrs]
  
  n <- nrow(X)
  p <- ncol(X)
  cat("Dimensions of X:", n, "genotypes and", p, "markers\n")
  
  G <- tcrossprod(X)/p
  dim(G)
  
  #### Geno Data for Bayes ####
  rownames(X)
  X_Bayes <- merge(pheno, X, by.x = "Sample.ID", by.y = "row.names")
  X_Bayes <- X_Bayes[,-1:-3]
  
  #### Import all Weather Data ####
  weather <- read.csv("W_SummerOnly.csv") #Tried with Weather PCs- it was worse.
  W <- merge(pheno, weather, by = "ENV", all.x = T)
  W <- W[order(W$Sample.ID),] #sort to match pheno
  W <- W %>% select(where(~mean(is.na(.)) < .01)) # Drop columns with more than 1% missing values (30 columns dropped)
  
  #### Define Training/ Testing ####
  y <- W[,"Yield"]
  yNA <- y
  test <- which(W$ENV == "2024Lincoln")
  yNA[test] <- NA
  length(test)
  obs <- W %>% filter(ENV == "2024Lincoln")
  
  #### Matrix Creation for GBLUP ####
  # incidence matrix of main effects of genotype
  Z0.L <- model.matrix(~0 + Sample.ID, data=W)
  Z.L <- as.matrix(Z0.L[1:dim(Z0.L)[1],1:dim(Z0.L)[2]]) #same info but cleaner format
  
  # incidence matrix of main effects of environment
  Z0.E <- model.matrix(~0 + ENV, data=W)
  Z.E <- as.matrix(Z0.E[1:dim(Z0.E)[1],1:dim(Z0.E)[2]]) #same info but cleaner format
  ZEZE <- tcrossprod(Z.E)
  
  # Augmented G Matrix
  G_augE <-Z.L %*% G %*% t(Z.L)
  G_augGE <- G_augE*ZEZE #Hadamard Product
  
  #Env covariate matrix
  Wmat <- W %>% select(-ENV, -Sample.ID, -X, -year, -Location, -Yield)
  Wmat <- sapply(Wmat, as.numeric)
  Wmat <- mean_imputation(Wmat)
  Wmat <- Wmat/sqrt(ncol(Wmat)) #Normalize
  Wmat2 <- tcrossprod(Wmat)/ncol(Wmat) #Omega
  
  GW <- G_augE * Wmat2
  
  # Machine Learning Data
  X_ML <- merge(pheno, X, by.x = "Sample.ID", by.y = "row.names")
  X_ML <- X_ML %>% select(Yield, everything())
  X_ML$Sample.ID <- as.factor(X_ML$Sample.ID)
  X_ML$ENV <- as.factor(X_ML$ENV)
  
  cvdf <- data.frame()
  for (j in 1:reps) {
    
    ### Model 1: G + E
    EtaM1 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(K = G_augE, model = "RKHS"))
    
    Model1 <- BGLR(y = yNA, ETA = EtaM1, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod1 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model1$yHat[test],
      model = 1,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod1)
    
    ### Model 2: G + E + GxE
    EtaM2 <- list(Line =list(X = Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(K = G_augE, model = "RKHS"),
                  GE = list(K = G_augGE, model = "RKHS"))
    
    Model2 <- BGLR(y = yNA, ETA = EtaM2, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod2 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model2$yHat[test],
      model = 2,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod2)
    
    ### Model 3: G + W
    EtaM3 <- list(Line =list(X = Z.L, model = "BRR"),
                  W = list(K = Wmat2, model = "RKHS"),
                  G = list(K = G_augE, model = "RKHS"))
    
    Model3 <- BGLR(y = yNA, ETA = EtaM3, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod3 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model3$yHat[test],
      model = 3,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod3)
    
    ### Model 4: G + W + GxW
    EtaM4 <- list(Line =list(X = Z.L, model = "BRR"),
                  W = list(K = Wmat2, model = "RKHS"),
                  G = list(K = G_augE, model = "RKHS"),
                  GW = list(K = GW, model = "RKHS"))
    
    Model4 <- BGLR(y = yNA, ETA = EtaM4, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod4 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model4$yHat[test],
      model = 4,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod4)
    
    ### Model 5: BayesA
    EtaM5 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BayesA"))
    
    Model5 <- BGLR(y = yNA, ETA = EtaM5, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod5 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model5$yHat[test],
      model = 5,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod5)
    
    ### Model 6: BayesB
    EtaM6 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BayesB"))
    
    Model6 <- BGLR(y = yNA, ETA = EtaM6, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod6 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model6$yHat[test],
      model = 6,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod6)
    
    ### Model 7: BayesC
    EtaM7 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BayesC"))
    
    Model7 <- BGLR(y = yNA, ETA = EtaM7, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod7 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model7$yHat[test],
      model = 7,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod7)
    
    ### Model 8: Bayes LASSO
    EtaM8 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BL"))
    
    Model8 <- BGLR(y = yNA, ETA = EtaM8, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod8 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model8$yHat[test],
      model = 8,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod8)
  }
  rep_df <- rbind(rep_df, cvdf)
  rm(cvdf)
  write.csv(rep_df, "BGLR_Yield_Results2.csv")
}

#### Machine Learning #####
rep_df <- data.frame()
for (i in 1:8) {
  pheno <- read.csv("all_phenoV7.csv")
  pheno <- pheno %>%
    mutate(Expt = case_when(
      grepl("BDUP", Trial) ~ "BDUP",
      grepl("BVT", Trial) ~ "BVT",
      grepl("BS4R8", Trial) ~ "BS4R8",
      TRUE ~ NA_character_
    ))
  
  pheno <- pheno %>% filter(!is.na(Yield)) %>% filter(Expt != "BDUP") %>% filter (!is.na(Sample.ID))
  pheno$ENV <- as.factor(pheno$ENV)
  pheno <- pheno %>% select("Sample.ID", "ENV", "Yield")
  pheno$ENV <- gsub( " .*$", "", pheno$ENV)
  
  if (i == 2) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln"))
  } else if (i == 3) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 4) {
    pheno <- pheno %>% filter(str_detect(ENV, "202"))
  } else if (i == 5) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead"))
  } else if (i ==6) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 7) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln"))
  } else if (i == 8) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln")) %>% filter(str_detect(ENV, "202"))
  }
  envs <- length(unique(pheno$ENV))
  pheno$Sample.ID <- as.factor(pheno$Sample.ID)
  pheno_list <- unique(pheno$Sample.ID)
  pheno_list <- data.frame(pheno_list)
  
  #### Geno Data for GBLUP ####
  X <- read.csv("geno_allV2.csv")
  X <- X[order(X[, 1]), ]
  X_merge <- merge(pheno_list, X, by.x = "pheno_list", by.y = "X", all.x = T)
  
  geno_order <- X_merge[, 1] # Store ordered genotype sample names
  X_merge <- X_merge[, -1] # Remove the first column (sample names) from X
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X_merge)) > 0))
  X_merge <- X_merge[, !colnames(X_merge) %in% NAMkrs]
  
  # Scale the data (center and scale each column)
  X_merge <- scale(X_merge, center = TRUE, scale = TRUE)
  
  # Set row names to the ordered genotype sample names
  rownames(X_merge) <- geno_order
  X <- X_merge
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X)) > 0))
  X <- X[, !colnames(X) %in% NAMkrs]
  
  n <- nrow(X)
  p <- ncol(X)
  cat("Dimensions of X:", n, "genotypes and", p, "markers\n")
  
  G <- tcrossprod(X)/p
  dim(G)
  
  #### Geno Data for Bayes ####
  rownames(X)
  X_Bayes <- merge(pheno, X, by.x = "Sample.ID", by.y = "row.names")
  X_Bayes <- X_Bayes[,-1:-3]
  
  #### Import all Weather Data ####
  weather <- read.csv("W_SummerOnly.csv") #Tried with Weather PCs- it was worse.
  W <- merge(pheno, weather, by = "ENV", all.x = T)
  W <- W[order(W$Sample.ID),] #sort to match pheno
  W <- W %>% select(where(~mean(is.na(.)) < .01)) # Drop columns with more than 1% missing values (30 columns dropped)
  
  #### Define Training/ Testing ####
  y <- W[,"Yield"]
  yNA <- y
  test <- which(W$ENV == "2024Lincoln")
  yNA[test] <- NA
  length(test)
  obs <- W %>% filter(ENV == "2024Lincoln")
  
  #### Matrix Creation for GBLUP ####
  # incidence matrix of main effects of genotype
  Z0.L <- model.matrix(~0 + Sample.ID, data=W)
  Z.L <- as.matrix(Z0.L[1:dim(Z0.L)[1],1:dim(Z0.L)[2]]) #same info but cleaner format
  
  # incidence matrix of main effects of environment
  Z0.E <- model.matrix(~0 + ENV, data=W)
  Z.E <- as.matrix(Z0.E[1:dim(Z0.E)[1],1:dim(Z0.E)[2]]) #same info but cleaner format
  ZEZE <- tcrossprod(Z.E)
  
  # Augmented G Matrix
  G_augE <-Z.L %*% G %*% t(Z.L)
  G_augGE <- G_augE*ZEZE #Hadamard Product
  
  #Env covariate matrix
  Wmat <- W %>% select(-ENV, -Sample.ID, -X, -year, -Location, -Yield)
  Wmat <- sapply(Wmat, as.numeric)
  Wmat <- mean_imputation(Wmat)
  Wmat <- Wmat/sqrt(ncol(Wmat)) #Normalize
  Wmat2 <- tcrossprod(Wmat)/ncol(Wmat) #Omega
  
  GW <- G_augE * Wmat2
  
  # Machine Learning Data
  X_ML <- merge(pheno, X, by.x = "Sample.ID", by.y = "row.names")
  X_ML <- X_ML %>% select(Yield, everything())
  X_ML$Sample.ID <- as.factor(X_ML$Sample.ID)
  X_ML$ENV <- as.factor(X_ML$ENV)
  
  cvdf <- data.frame()
  for (j in 1:reps) {
    
    # ### Model 9: Random Forest
    t <- ncol(X_ML)
    #col_list <- nearZeroVar(X_ML)
    col_list <- which(apply(X_ML[,4:t], 2, var) == 0) #Exactly 0 variance
    
    if (length(col_list)>0) {
      X_ML <- X_ML[-col_list]
    }
    
    test <- X_ML %>% filter(ENV == "2024Lincoln")
    y_ML <- test$Yield
    test$Yield <- NULL
    
    train <- X_ML %>% filter(ENV != "2024Lincoln")
    
    k_rf = c(250, 500, 1000) #mtry grid
    mod_name <- c("9A", "9B", "9C")
    ntree_list <- c(250, 500, 1000) #ntrees grid
    
    for (h in 1:3) {
      mod <- mod_name[h]
      ntree <- ntree_list[h]
      Model9 <- train(Yield ~ .,
                      data = train,
                      method = "rf",
                      tuneGrid = expand.grid(mtry = k_rf),
                      ntree = ntree)
      assign(paste0("Model_", mod), Model9) #Save models with each value of ntree
      
      pred_rf_mod9 <- data.frame(
        name = test$Sample.ID,
        env = test$ENV,
        y = y_ML,
        yhat = predict(Model9, newdata = test),
        model = 9
      )
      
      pred_mod9 <- cor(pred_rf_mod9$y, pred_rf_mod9$yhat)
      assign(paste0("Pred_", mod), pred_mod9) #Save pred_accuracy with each value of ntree
      rm(pred_rf_mod9)
    }
    
    best_mod <- c("9A", "9B", "9C")[which.max(c(Pred_9A, Pred_9B, Pred_9C))] #Get model with highest pred_accuracy
    best_mod_name <- get(paste0("Model_", best_mod)) #call the best model by name
    
    assign("Model9", best_mod_name) #rename best model to generic Model 9
    
    pred_rf_mod9 <- data.frame(
      name = test$Sample.ID,
      env = test$ENV,
      y = y_ML,
      yhat = predict(Model9, newdata = test),
      model = 9,
      rep = j,
      CV = i,
      tune = paste0(Model9$bestTune[[1]], "_", Model9$finalModel['ntree'][[1]])
    )
    
    cvdf <- rbind(cvdf, pred_rf_mod9)
    
    ### Model 10: Support Vector Machine
    k_svm = c(0.01, 0.1, 1, 10) #c grid
    
    ## Hot-one encoding for ENV and Sample ID
    dummy <- dummyVars(~ENV + Sample.ID, X_ML)
    DummyMat <- predict(dummy, X_ML)
    genos <- X_ML %>% filter(ENV == "2024Lincoln") %>% select(Sample.ID)
    locs <- X_ML %>% filter(ENV == "2024Lincoln") %>% select(ENV)
    
    SVM_dat <- cbind(DummyMat, X_ML[,-2:-3])
    SVM_dat <- SVM_dat %>% select(Yield, everything())
    
    test <- SVM_dat %>% filter(ENV.2024Lincoln == 1)
    y_ML <- test$Yield
    test$Yield <- NULL
    
    train <- SVM_dat %>% filter(ENV.2024Lincoln != 1)
    
    Model10 <- train(Yield ~ .,
                     data = train,
                     method = "svmLinear2",
                     scale = F,
                     tuneGrid = expand.grid(cost = k_svm))
    
    pred_svm_mod10 <- data.frame(
      name = unname(genos),
      env = unname(locs),
      y = y_ML,
      yhat = predict(Model10, newdata = test),
      model = 10,
      rep = j,
      CV = i,
      tune = Model10$bestTune[[1]]
    )
    
    cvdf <- rbind(cvdf, pred_svm_mod10)
  }
  rep_df <- rbind(rep_df, cvdf)
  rm(cvdf)
  write.csv(rep_df, "ML_Yield_Results2.csv")
}

#### Winter Survival ####
## Import all Pheno Data ###
rep_df <- data.frame()
for (i in c(2, 5)) {
  pheno <- read.csv("all_phenoV7.csv")
  pheno <- pheno %>%
    mutate(Expt = case_when(
      grepl("BDUP", Trial) ~ "BDUP",
      grepl("BVT", Trial) ~ "BVT",
      grepl("BS4R8", Trial) ~ "BS4R8",
      TRUE ~ NA_character_
    ))
  
  pheno <- pheno %>% filter(!is.na(WinSur)) %>% filter(Expt != "BDUP") %>% filter (!is.na(Sample.ID))
  pheno$ENV <- as.factor(pheno$ENV)
  pheno <- pheno %>% select("Sample.ID", "ENV", "WinSur")
  pheno$ENV <- gsub( " .*$", "", pheno$ENV)
  
  if (i == 2) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln"))
  } else if (i == 3) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 4) {
    pheno <- pheno %>% filter(str_detect(ENV, "202"))
  } else if (i == 5) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead"))
  } else if (i ==6) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 7) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln"))
  } else if (i == 8) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln")) %>% filter(str_detect(ENV, "202"))
  }
  envs <- length(unique(pheno$ENV))
  pheno$Sample.ID <- as.factor(pheno$Sample.ID)
  pheno_list <- unique(pheno$Sample.ID)
  pheno_list <- data.frame(pheno_list)
  
  #### Geno Data for GBLUP ####
  X <- read.csv("geno_allV2.csv")
  X <- X[order(X[, 1]), ]
  X_merge <- merge(pheno_list, X, by.x = "pheno_list", by.y = "X", all.x = T)
  
  geno_order <- X_merge[, 1] # Store ordered genotype sample names
  X_merge <- X_merge[, -1] # Remove the first column (sample names) from X
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X_merge)) > 0))
  X_merge <- X_merge[, !colnames(X_merge) %in% NAMkrs]
  
  # Scale the data (center and scale each column)
  X_merge <- scale(X_merge, center = TRUE, scale = TRUE)
  
  # Set row names to the ordered genotype sample names
  rownames(X_merge) <- geno_order
  X <- X_merge
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X)) > 0))
  X <- X[, !colnames(X) %in% NAMkrs]
  
  n <- nrow(X) 
  p <- ncol(X)
  cat("Dimensions of X:", n, "genotypes and", p, "markers\n")
  
  G <- tcrossprod(X)/p
  dim(G)
  
  #### Geno Data for Bayes ####
  rownames(X)
  X_Bayes <- merge(pheno, X, by.x = "Sample.ID", by.y = "row.names")
  X_Bayes <- X_Bayes[,-1:-3]
  
  #### Import all Weather Data ####
  weather <- read.csv("W_SummerOnly.csv") #Tried with Weather PCs- it was worse.
  W <- merge(pheno, weather, by = "ENV", all.x = T)
  W <- W[order(W$Sample.ID),] #sort to match pheno
  W <- W %>% select(where(~mean(is.na(.)) < .01)) # Drop columns with more than 1% missing values (30 columns dropped)
  
  #### Define Training/ Testing ####
  y <- W[,"WinSur"]
  yNA <- y
  test <- which(W$ENV == "2024Lincoln")
  yNA[test] <- NA
  length(test)
  obs <- W %>% filter(ENV == "2024Lincoln")
  
  #### Matrix Creation for GBLUP ####
  # incidence matrix of main effects of genotype
  Z0.L <- model.matrix(~0 + Sample.ID, data=W)
  Z.L <- as.matrix(Z0.L[1:dim(Z0.L)[1],1:dim(Z0.L)[2]]) #same info but cleaner format
  
  # incidence matrix of main effects of environment
  Z0.E <- model.matrix(~0 + ENV, data=W)
  Z.E <- as.matrix(Z0.E[1:dim(Z0.E)[1],1:dim(Z0.E)[2]]) #same info but cleaner format
  ZEZE <- tcrossprod(Z.E)
  
  # Augmented G Matrix
  G_augE <-Z.L %*% G %*% t(Z.L)
  G_augGE <- G_augE*ZEZE #Hadamard Product
  
  #Env covariate matrix
  Wmat <- W %>% select(-ENV, -Sample.ID, -X, -year, -Location, -WinSur)
  Wmat <- sapply(Wmat, as.numeric)
  Wmat <- mean_imputation(Wmat)
  Wmat <- Wmat/sqrt(ncol(Wmat)) #Normalize
  Wmat2 <- tcrossprod(Wmat)/ncol(Wmat) #Omega
  
  GW <- G_augE * Wmat2
  
  # Machine Learning Data
  X_ML <- merge(pheno, X, by.x = "Sample.ID", by.y = "row.names")
  X_ML <- X_ML %>% select(WinSur, everything())
  X_ML$Sample.ID <- as.factor(X_ML$Sample.ID)
  X_ML$ENV <- as.factor(X_ML$ENV)
  
  cvdf <- data.frame()
  for (j in 1:reps) {
    
    ### Model 1: G + E
    EtaM1 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(K = G_augE, model = "RKHS"))
    
    Model1 <- BGLR(y = yNA, ETA = EtaM1, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod1 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model1$yHat[test],
      model = 1,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod1)
    
    ### Model 2: G + E + GxE
    EtaM2 <- list(Line =list(X = Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(K = G_augE, model = "RKHS"),
                  GE = list(K = G_augGE, model = "RKHS"))
    
    Model2 <- BGLR(y = yNA, ETA = EtaM2, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod2 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model2$yHat[test],
      model = 2,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod2)
    
    ### Model 3: G + W
    EtaM3 <- list(Line =list(X = Z.L, model = "BRR"),
                  W = list(K = Wmat2, model = "RKHS"),
                  G = list(K = G_augE, model = "RKHS"))
    
    Model3 <- BGLR(y = yNA, ETA = EtaM3, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod3 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model3$yHat[test],
      model = 3,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod3)
    
    ### Model 4: G + W + GxW
    EtaM4 <- list(Line =list(X = Z.L, model = "BRR"),
                  W = list(K = Wmat2, model = "RKHS"),
                  G = list(K = G_augE, model = "RKHS"),
                  GW = list(K = GW, model = "RKHS"))
    
    Model4 <- BGLR(y = yNA, ETA = EtaM4, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod4 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model4$yHat[test],
      model = 4,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod4)
    
    ### Model 5: BayesA
    EtaM5 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BayesA"))
    
    Model5 <- BGLR(y = yNA, ETA = EtaM5, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod5 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model5$yHat[test],
      model = 5,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod5)
    
    ### Model 6: BayesB
    EtaM6 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BayesB"))
    
    Model6 <- BGLR(y = yNA, ETA = EtaM6, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod6 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model6$yHat[test],
      model = 6,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod6)
    
    ### Model 7: BayesC
    EtaM7 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BayesC"))
    
    Model7 <- BGLR(y = yNA, ETA = EtaM7, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod7 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model7$yHat[test],
      model = 7,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod7)
    
    ### Model 8: Bayes LASSO
    EtaM8 <- list(Line =list(X=Z.L, model = "BRR"),
                  ENV = list(X = Z.E, model = "FIXED"),
                  G = list(X = X_Bayes, model = "BL"))
    
    Model8 <- BGLR(y = yNA, ETA = EtaM8, nIter = nIter, burnIn = burnIn, verbose = F)
    
    pred_gblup_bglr_mod8 <- data.frame(
      name = pheno[test, "Sample.ID"],
      env = pheno[test, "ENV"],
      y = y[test],
      yhat = Model8$yHat[test],
      model = 8,
      rep = j,
      CV = i
    )
    
    cvdf <- rbind(cvdf, pred_gblup_bglr_mod8)
  }
  rep_df <- rbind(rep_df, cvdf)
  rm(cvdf)
  write.csv(rep_df, "BGLR_WinSur_Results2.csv")
}

#### Machine Learning #####
rep_df <- data.frame()
for (i in c(2,5)) {
  pheno <- read.csv("all_phenoV7.csv")
  pheno <- pheno %>%
    mutate(Expt = case_when(
      grepl("BDUP", Trial) ~ "BDUP",
      grepl("BVT", Trial) ~ "BVT",
      grepl("BS4R8", Trial) ~ "BS4R8",
      TRUE ~ NA_character_
    ))
  
  pheno <- pheno %>% filter(!is.na(WinSur)) %>% filter(Expt != "BDUP") %>% filter (!is.na(Sample.ID))
  pheno$ENV <- as.factor(pheno$ENV)
  pheno <- pheno %>% select("Sample.ID", "ENV", "WinSur")
  pheno$ENV <- gsub( " .*$", "", pheno$ENV)
  
  if (i == 2) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln"))
  } else if (i == 3) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 4) {
    pheno <- pheno %>% filter(str_detect(ENV, "202"))
  } else if (i == 5) {
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead"))
  } else if (i ==6) {
    next
    pheno <- pheno %>% filter(str_detect(ENV, "Lincoln") | str_detect(ENV, "Mead")) %>% filter(str_detect(ENV, "202"))
  } else if (i == 7) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln"))
  } else if (i == 8) {
    pheno <- pheno %>% filter(str_detect(ENV, "Sidney") | str_detect(ENV, "2024Lincoln")) %>% filter(str_detect(ENV, "202"))
  }
  envs <- length(unique(pheno$ENV))
  pheno$Sample.ID <- as.factor(pheno$Sample.ID)
  pheno_list <- unique(pheno$Sample.ID)
  pheno_list <- data.frame(pheno_list)
  
  #### Geno Data for GBLUP ####
  X <- read.csv("geno_allV2.csv")
  X <- X[order(X[, 1]), ]
  X_merge <- merge(pheno_list, X, by.x = "pheno_list", by.y = "X", all.x = T)
  
  geno_order <- X_merge[, 1] # Store ordered genotype sample names
  X_merge <- X_merge[, -1] # Remove the first column (sample names) from X
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X_merge)) > 0))
  X_merge <- X_merge[, !colnames(X_merge) %in% NAMkrs]
  
  # Scale the data (center and scale each column)
  X_merge <- scale(X_merge, center = TRUE, scale = TRUE)
  
  # Set row names to the ordered genotype sample names
  rownames(X_merge) <- geno_order
  X <- X_merge
  
  # Remove markers (columns) with NA values
  NAMkrs <- names(which(colSums(is.na(X)) > 0))
  X <- X[, !colnames(X) %in% NAMkrs]
  
  n <- nrow(X) 
  p <- ncol(X)
  cat("Dimensions of X:", n, "genotypes and", p, "markers\n")
  
  G <- tcrossprod(X)/p
  dim(G)
  
  #### Import all Weather Data ####
  weather <- read.csv("W_SummerOnly.csv") #Tried with Weather PCs- it was worse.
  W <- merge(pheno, weather, by = "ENV", all.x = T)
  W <- W[order(W$Sample.ID),] #sort to match pheno
  W <- W %>% select(where(~mean(is.na(.)) < .01)) # Drop columns with more than 1% missing values (30 columns dropped)
  
  #### Define Training/ Testing ####
  y <- W[,"WinSur"]
  yNA <- y
  test <- which(W$ENV == "2024Lincoln")
  yNA[test] <- NA
  length(test)
  obs <- W %>% filter(ENV == "2024Lincoln")
  
  # Machine Learning Data
  X_ML <- merge(pheno, X, by.x = "Sample.ID", by.y = "row.names")
  X_ML <- X_ML %>% select(WinSur, everything())
  X_ML$Sample.ID <- as.factor(X_ML$Sample.ID)
  X_ML$ENV <- as.factor(X_ML$ENV)
  
  cvdf <- data.frame()
  for (j in 1:reps) {
    
    # ### Model 9: Random Forest
    t <- ncol(X_ML)
    #col_list <- nearZeroVar(X_ML)
    col_list <- which(apply(X_ML[,4:t], 2, var) == 0) #Exactly 0 variance
    
    if (length(col_list)>0) {
      X_ML <- X_ML[-col_list]
    }
    
    test <- X_ML %>% filter(ENV == "2024Lincoln")
    y_ML <- test$WinSur
    test$WinSur <- NULL
    
    train <- X_ML %>% filter(ENV != "2024Lincoln")
    
    k_rf = c(250, 500, 1000) #mtry grid
    mod_name <- c("9A", "9B", "9C")
    ntree_list <- c(250, 500, 1000) #ntrees grid
    
    for (h in 1:3) {
      mod <- mod_name[h]
      ntree <- ntree_list[h]
      Model9 <- train(WinSur ~ .,
                      data = train,
                      method = "rf",
                      tuneGrid = expand.grid(mtry = k_rf),
                      ntree = ntree)
      assign(paste0("Model_", mod), Model9) #Save models with each value of ntree
      
      pred_rf_mod9 <- data.frame(
        name = test$Sample.ID,
        env = test$ENV,
        y = y_ML,
        yhat = predict(Model9, newdata = test),
        model = 9
      )
      
      pred_mod9 <- cor(pred_rf_mod9$y, pred_rf_mod9$yhat)
      assign(paste0("Pred_", mod), pred_mod9) #Save pred_accuracy with each value of ntree
      rm(pred_rf_mod9)
    }
    
    best_mod <- c("9A", "9B", "9C")[which.max(c(Pred_9A, Pred_9B, Pred_9C))] #Get model with highest pred_accuracy
    best_mod_name <- get(paste0("Model_", best_mod)) #call the best model by name
    
    assign("Model9", best_mod_name) #rename best model to generic Model 9
    
    pred_rf_mod9 <- data.frame(
      name = test$Sample.ID,
      env = test$ENV,
      y = y_ML,
      yhat = predict(Model9, newdata = test),
      model = 9,
      rep = j,
      CV = i,
      tune = paste0(Model9$bestTune[[1]], "_", Model9$finalModel['ntree'][[1]])
    )
    
    cvdf <- rbind(cvdf, pred_rf_mod9)
    
    ### Model 10: Support Vector Machine
    k_svm = c(0.01, 0.1, 1, 10) #c grid
    
    ## Hot-one encoding for ENV and Sample ID
    dummy <- dummyVars(~ENV + Sample.ID, X_ML)
    DummyMat <- predict(dummy, X_ML)
    genos <- X_ML %>% filter(ENV == "2024Lincoln") %>% select(Sample.ID)
    locs <- X_ML %>% filter(ENV == "2024Lincoln") %>% select(ENV)
    
    SVM_dat <- cbind(DummyMat, X_ML[,-2:-3])
    SVM_dat <- SVM_dat %>% select(WinSur, everything())
    
    test <- SVM_dat %>% filter(ENV.2024Lincoln == 1)
    y_ML <- test$WinSur
    test$WinSur <- NULL
    
    train <- SVM_dat %>% filter(ENV.2024Lincoln != 1)
    
    Model10 <- train(WinSur ~ .,
                     data = train,
                     method = "svmLinear2",
                     scale = F,
                     tuneGrid = expand.grid(cost = k_svm))
    
    pred_svm_mod10 <- data.frame(
      name = unname(genos),
      env = unname(locs),
      y = y_ML,
      yhat = predict(Model10, newdata = test),
      model = 10,
      rep = j,
      CV = i,
      tune = Model10$bestTune[[1]]
    )
    
    cvdf <- rbind(cvdf, pred_svm_mod10)
  }
  rep_df <- rbind(rep_df, cvdf)
  rm(cvdf)
  write.csv(rep_df, "ML_WinSur_Results2.csv")
}  