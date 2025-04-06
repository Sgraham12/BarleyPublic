#### 3. Analysis of Phenotypic Data
# S. Graham
# 2/5/2025

#### Required files: ####
# 1. all_phenoV8.csv (2.PhenoGenoCleaning.R)
# 2. Cleaned_YieldTrials.csv (1.DataCleaning.R)

#### Output files: ####
# 1. GenTrend.png
# 2. HeatMap.png

#Load Packages
library(asreml)
library(tidyverse)

setwd("~/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Documents/Statistics MS/Barley")

#### Test Factors for Significance ####
data <- read.csv("all_phenoV8.csv")
data$Year <- substr(data$ENV, 1, 4)
data$Loc <- gsub('[0-9]+', '', data$ENV)
data$Loc <- gsub( " .*$", "", data$Loc )
data$Name1 <- as.factor(data$Name1)
data$Year <- as.factor(data$Year)
data$Loc <- as.factor(data$Loc)
data$site <- paste0(data$Loc, " ", data$Year)
data$site <- as.factor(data$site)
BVT <- data %>% filter(Trial == "BVT")

### Yield
mod <- asreml(fixed = Yield ~ site,
              random = ~ Name1 + Name1:site,
              data = BVT, 
              weights = W.Yield,
              family = asr_gaussian(dispersion = 1.0),
              maxiter = 25,
              na.action = na.method(y = "omit", x = "include"))
mod <- update.asreml(mod)
  
modName <- asreml(fixed = Yield ~ site,
              random = ~ Name1:site,
              data = BVT, 
              weights = W.Yield,
              family = asr_gaussian(dispersion = 1.0),
              maxiter = 25,
              na.action = na.method(y = "omit", x = "include"))
modName <- update.asreml(modName)
  
modNameSite <- asreml(fixed = Yield ~ site,
              random = ~ Name1,
              data = BVT, 
              weights = W.Yield,
              family = asr_gaussian(dispersion = 1.0),
              maxiter = 25,
              na.action = na.method(y = "omit", x = "include"))
modNameSite <- update.asreml(modNameSite)

wald(mod)  
lrt(mod, modName)
lrt(mod, modNameSite)

modH2 <- asreml(fixed = Yield ~ 1,
              random = ~ site + Name1 + Name1:site,
              data = BVT, 
              weights = W.Yield,
              family = asr_gaussian(dispersion = 1.0),
              maxiter = 25,
              na.action = na.method(y = "omit", x = "include"))
modH2 <- update.asreml(modH2)

vc.g <- summary(modH2)$varcomp['Name1','component']
vdBLUP.mat <- predict(modH2, classify="Name1", only="Name1", sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
h2 <- 1 - (vdBLUP.avg / 2 / vc.g)
h2 #heritability

### Winsur
mod <- asreml(fixed = WinSur ~ site,
              random = ~ Name1 + Name1:site,
              data = BVT, 
              weights = W.WinSur,
              family = asr_gaussian(dispersion = 1.0),
              maxiter = 25,
              na.action = na.method(y = "omit", x = "include"))
mod <- update.asreml(mod)

modName <- asreml(fixed = WinSur ~ site,
                  random = ~ Name1:site,
                  data = BVT, 
                  weights = W.WinSur,
                  family = asr_gaussian(dispersion = 1.0),
                  maxiter = 25,
                  na.action = na.method(y = "omit", x = "include"))
modName <- update.asreml(modName)

modNameSite <- asreml(fixed = WinSur ~ site,
                      random = ~ Name1,
                      data = BVT, 
                      weights = W.WinSur,
                      family = asr_gaussian(dispersion = 1.0),
                      maxiter = 25,
                      na.action = na.method(y = "omit", x = "include"))
modNameSite <- update.asreml(modNameSite)

wald(mod)  
lrt(mod, modName)
lrt(mod, modNameSite)

modH2 <- asreml(fixed = WinSur ~ 1,
                random = ~ site + Name1 + Name1:site,
                data = BVT, 
                weights = W.WinSur,
                family = asr_gaussian(dispersion = 1.0),
                maxiter = 25,
                na.action = na.method(y = "omit", x = "include"))
modH2 <- update.asreml(modH2)

vc.g <- summary(modH2)$varcomp['Name1','component']
vdBLUP.mat <- predict(modH2, classify="Name1", only="Name1", sed=TRUE)$sed^2 # obtain squared s.e.d. matrix 
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
h2 <- 1 - (vdBLUP.avg / 2 / vc.g)
h2 #heritability

#### Genetic Gain ####
pheno <- read.csv("Cleaned_YieldTrials.csv")

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

## Yield
#No yield data for BDUP: 2016LN, 2017S, 2018S, 2020M, 2023L #Limited Data for BVT: 2007M, 2007S, BDUP: 2015L, 2015M, 2016S, 2019L
trials_to_remove <- c("2016Lincoln BDUP", "2017Sidney BDUP", "2018Sidney BDUP", "2020Mead BDUP", "2023Lincoln BDUP", "2007Mead BVT", "2007Sidney BVT", "2015Lincoln BDUP", "2015Mead BDUP", "2016Sidney BDUP", "2019Lincoln BDUP", "2019Sidney BDUP")
pheno <- AllPhenos[!AllPhenos$trial %in% trials_to_remove, ]

pheno$check <- ifelse(grepl("P-954", pheno$Name1), 1, 0)
pheno$design <- "RCBD"
pheno <- pheno %>% filter(Experiment.Name == "BVT")

#### Calculate Data Overlap ####
consecutive_years <- pheno %>%
  mutate(Year = as.numeric(Year)) %>%             # Ensure Year is numeric
  distinct(Name1, Year) %>%                       # Get unique Name1-Year combinations
  arrange(Year) %>%                               # Sort by Year
  rowwise() %>%
  mutate(
    ClosestPrevYear = {
      prev_years <- pheno$Year[pheno$Year < Year]
      if (length(prev_years) > 0) {
        max(prev_years, na.rm = TRUE)             # Find the closest previous year
      } else {
        NA                                       # Handle cases with no previous year
      }
    }
  ) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarize(
    LinesInYear = n(),                            # Total lines in the current year
    LinesFromPrevYear = sum(Name1 %in% pheno$Name1[pheno$Year == ClosestPrevYear], na.rm = TRUE), # Compare with closest previous year
    PercentFromPrevYear = ifelse(
      LinesInYear > 0,
      (LinesFromPrevYear / LinesInYear) * 100,   # Calculate percentage
      NA                                        # Avoid division by zero
    ),
    .groups = 'drop'
  )

DurationOfGenos <- pheno %>% group_by(Name1) %>% summarise(First = min(Year), Last = max(Year), Duration = (Last - First))

mean(consecutive_years$PercentFromPrevYear)

#### EBV METHOD  ####
#https://doi.org/10.1002/tpg2.20471
#https://github.com/Biometrics-IITA/Estimating-Realized-Genetic-Gain/blob/main/realized_genetic_gain.R
#### STEP 1: estimate and discard env. with low heritability (Cullis)
#### Linear Mixed Model (LMM) for each env.

#### Run single trial analysis: geno as random
datCleaned <- pheno |>
  mutate(across(
    c(Experiment.Name, BLOC, Year, Location, Name1, env, trial, check, design),
    as.factor
  ))

single.study.random <- datCleaned |>
  nest_by(trial, design) |>
  mutate(
    data = list(droplevels(data)),
    model.rando = ifelse(design == "Alpha",
                         list(tryCatch(
                           asreml(
                             fyld ~ rep,
                             random = ~ geno + block,
                             workspace = "3gb",
                             data = data
                           ),
                           error = function(e) {
                             NULL
                           }
                         )),
                         list(tryCatch(
                           # RCBD
                           asreml(
                             YLDBUA ~ BLOC,
                             random = ~ Name1,
                             workspace = "3gb",
                             data = data
                           ),
                           error = function(e) {
                             NULL
                           }
                         )))
  )

#### Filter null and estimate H2
single.study.random <- single.study.random |>
  filter(!is.null(model.rando)) |>
  mutate(
    model.rando = list(eval(
      parse(text = "update.asreml(model.rando)")
    )),
    model.rando = list(eval(
      parse(text = "update.asreml(model.rando)")
    )),
    model.rando = list(eval(
      parse(text = "update.asreml(model.rando)")
    )),
    vg = list(
      summary(model.rando)$varcomp |>
        rownames_to_column(var = "Effect") |>
        select(Effect, everything()) |>
        filter(grepl("Name1", Effect)) |>
        pull(component)
    ),
    vd.mat = list(
      predict(
        model.rando,
        classify = "Name1",
        only = "Name1",
        sed = TRUE
      )$sed ^ 2
    ),
    vd.avg = list(mean(vd.mat[upper.tri(vd.mat, diag = FALSE)])),
    H2 = list(1 - (vd.avg / (vg * 2)))
  )

#### STEP 2: example of a one-stage approach
#### Run a combined Linear Mixed Model (LMM) and estimate combined BLUEs and weights

#### Filter low H2 and prepare the data
list <- single.study.random %>% filter(H2 < 0.2) %>% select(trial, H2)
Remove <- droplevels(list$trial)

#### Run combined trial analysis: geno as fixed - USE STAGE 1 BLUES
data$trial <- paste0(data$Year,data$Loc, " ", data$Trial)
BVT <- data %>% filter(Trial == "BVT") %>% filter(!(trial %in% Remove)) #Filter by H2
#BVT <- data %>% filter(Trial == "BVT") #use this if we want to keep all trials regardless of H2
checks <- c("P-954", "P-845", "P-721", "TAMBAR 501", "P-713") #More than ten years
BVT <- BVT %>% filter(!(Name1 %in% checks))

model.combined <- asreml(fixed = Yield ~ Year + Name1,
                         random = ~ Loc + Year:Loc + Name1:Loc + Name1:Year + Name1:Year:Loc,
                         family = asr_gaussian(dispersion = 1.0),
                         weights = W.Yield,
                         workspace = "3gb",
                         data = BVT)

if (!is.null(model.combined)) {
  model.combined <- eval(parse(text = "update.asreml(model.combined)"))
  model.combined <- eval(parse(text = "update.asreml(model.combined)"))
  model.combined <- eval(parse(text = "update.asreml(model.combined)"))
  
  #### estimate combined BLUEs and weights
  predictions.combined <- predict.asreml(model.combined, pworkspace = "5gb", classify = "Name1", vcov = TRUE)
  blues.combined <- predictions.combined$pvals |>
    mutate(weight = 1 / (std.error) ^ 2) |>
    left_join(
      pheno |>
        select(Name1, Year) |>
        filter(Year != 0) |>
        droplevels() |>
        group_by(Name1) |>
        summarize(Year = min(Year, na.rm = TRUE))
    )
  
  blues.combined$W <- diag(solve(predictions.combined$vcov))
  blues.combined <- blues.combined %>% filter(predicted.value > 0) |>
    filter(!is.na(Year)) |>
    filter(!is.na(predicted.value)) |>
    left_join(pheno |>
                select(Name1, check) |>
                distinct(Name1, .keep_all = TRUE)) |>
    filter(check != 1) |> filter(Year != 2002) |>
    droplevels()
  
  #### STEP 3: weighted-linear regression of the combined BLUEs on the year of origin

  nY = length(unique(blues.combined$Year))
  
  if (nY >= 5) {
    fit.regression <- lm(predicted.value ~ Year,
                         weights = W,
                         data = blues.combined)
    slope <-  fit.regression$coefficients[2]
    intercept <-  fit.regression$coefficients[1]
    first.Year.geno <-  min(blues.combined$Year)
    last.Year.geno <-  max(blues.combined$Year)
    baseline <-  intercept + first.Year.geno * slope
    gg.percent <-  100 * slope / baseline
    pvalue <- summary(fit.regression)$coefficients |>
      as.data.frame() |>
      rownames_to_column("Effect") |>
      filter(Effect == "Year") |>
      pull("Pr(>|t|)")
    slope.sd <- summary(fit.regression)$coefficients |>
      as.data.frame() |>
      rownames_to_column("Effect") |>
      filter(Effect == "Year") |>
      pull("Std. Error")
    out.gg <- tibble(
      trait = "Yield",
      slope = slope,
      slope.sd = slope.sd,
      pvalue = pvalue,
      gg.percent = gg.percent
    )
  }
  out.gg
}
yield <- out.gg

#### Genetic Trend ####
trend <- ggplot(data = blues.combined, aes(x = Year, y = predicted.value)) + 
  geom_point() + 
  geom_smooth(method = 'lm') +
  ylab("Yield (bu/ac)") +
  xlab("First Year in AYT") +
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_text(size = 14),  # Axis title size
    axis.text = element_text(size = 12),   # Axis label text size
  )

ggsave(trend, filename = "GenTrend.png", units = "in", width = 7.25, height = 5, dpi = 800)

heatmap_data <- datCleaned %>%
  count(Name1, Year) %>%
  spread(key = Year, value = n, fill = 0)  # Spread Year as columns with the count as values

# Convert the data into a long format for ggplot
heatmap_data_long <- heatmap_data %>%
  gather(key = "Year", value = "Observations", -Name1)

# Plot the heatmap
heatmap <- ggplot(heatmap_data_long, aes(x = Name1, y = Year, fill = Observations)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),   # Remove labels on x-axis
        axis.title.x = element_text(size = 18),  # Remove x-axis title
        axis.title.y = element_text(size = 18), 
        axis.text.y = element_text(size = 12)) +
  labs(fill = "Number of Plots/Year") +
  xlab("Genotype")

ggsave(heatmap, filename = "HeatMap.png", units = "in", width = 7.25, height = 5, dpi = 800)
