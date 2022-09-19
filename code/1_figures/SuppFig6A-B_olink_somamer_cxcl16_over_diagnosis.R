# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----     Single cell transcriptomics reveals CD8 T cell recruitment     -----
# -----       to the cerebrospinal fluid during cognitive impairment       -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 09-06-2022
# Written by: Natalie Piehl
# Summary: Run ANCOVA on Olink and SomaLogical CXCL16 data
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("xlsx")
  library("rstatix")
  library("broom")
  library("emmeans")
  library("multcomp")
})

# Initialize paths
somalogical_dir <- "path/to/somalogic/protein/data/"
olink_path <- "path/to/olink-quanterix/protein/data/"
output_dir <- "path/to/export/results"

# Source helper functions
source("../0_preprocessing/00_helper_functions.R")

# Generate output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# OLINK
# Run ANCOVA b/w HC, MCI-AD, and AD

# Load in data
data <- read.xlsx(olink_path, sheetIndex = 2)

# Convert Gender to binary
data$Gender <- mapvalues(data$Gender,
                         from = c("F", "M"),
                         to = c(0,1))
data$Gender <- as.numeric(data$Gender)

# Define diagnoses
diagnoses <- list("HC", "MCI-AD", "AD", c("MCI-AD", "AD"))
diagnoses_labels <- c("HC", "MCI-AD", "AD", "CI")

# Run ANCOVA
data$Diagnosis <- factor(data$Diagnosis, levels = diagnoses[1:3])
res.aov <- aov(CXCL16 ~ Diagnosis + Age + Gender, data = data)

# Generate ANCOVA table
anova_table <- Anova(res.aov)
write.csv(anova_table, file = paste0(output_dir, "olink_HC_MCIAD_AD_anova_table.csv"))

# Run post-hoc test (Tukey)
postHocs <- glht(res.aov, linfct = mcp(Diagnosis = "Tukey"))
tukey_pvals <- data.frame(summary(postHocs)$test$pvalues)
tukey_pvals$Comparison <- names(summary(postHocs)$test$coefficients)
write.csv(tukey_pvals, file = paste0(output_dir, "olink_HC_MCIAD_AD_tukey_pvals.csv"))

#-------------------------------------------------------------------------------
# Run ANCOVA b/w HC and CI

# Run ANCOVA
data$Diagnosis_mod <- mapvalues(data$Diagnosis,
                                from = diagnoses[1:3],
                                to = c("HC", "CI", "CI"))
data$Diagnosis_mod <- factor(data$Diagnosis_mod, levels = c("HC", "CI"))
res.aov <- aov(CXCL16 ~ Diagnosis_mod + Age + Gender, data = data)

# Generate ANCOVA table
anova_table <- Anova(res.aov)
write.csv(anova_table, file = paste0(output_dir, "olink_HC_CI_anova_table.csv"))

# Run post-hoc test (Tukey)
postHocs <- glht(res.aov, linfct = mcp(Diagnosis_mod = "Tukey"))
tukey_pvals <- data.frame(summary(postHocs)$test$pvalues)
tukey_pvals$Comparison <- names(summary(postHocs)$test$coefficients)
write.csv(tukey_pvals, file = paste0(output_dir, "olink_HC_CI_tukey_pvals.csv"))

#-------------------------------------------------------------------------------
# Generate plot with all four groups (SuppFig 6A)

# Add duplicate rows for MCI and AD
data_ci <- data
data_ci$Diagnosis <- data_ci$Diagnosis_mod
data_ci <- data_ci[ which(data_ci$Diagnosis == "CI"),]
data <- rbind(data, data_ci)

# Generate boxplot
p_box <- ggplot(data, aes(x = Diagnosis, y = CXCL16)) +
  geom_boxplot(aes(fill = Diagnosis), fatten = 1, width = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 2) +
  geom_vline(xintercept = 3.5, color = "gray") +
  theme_Publication_blank() +
  scale_fill_manual(values = c("gray", "pink", "red", "lightblue")) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  labs(y = "Olink CXCL16",
       x = NULL) +
  theme(text = element_text(size = 25),
        legend.title=element_blank(),
        legend.position="right",
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 8))

# Export plot
set_panel_size(p_box, file = paste0(output_dir, "olink_cxcl16_boxplot.pdf"),
               width = unit(4, "in"), height = unit(4, "in"))

#-------------------------------------------------------------------------------
# SOMALOGICAL
# Run ANCOVA b/w HC, MCI-AD, and AD

# Load in data
data <- read.xlsx(paste0(somalogical_dir, "ADRC_cleaned_CSF_Wagner.xlsx"), sheetIndex = 1)
cdr <- read.xlsx(paste0(somalogical_dir, "ADRC_cleaned_CDR.xlsx"), sheetIndex = 1)
moca <- read.xlsx(paste0(somalogical_dir, "ADRC_cleaned_MOCA.xlsx"), sheetIndex = 1)

# Only complete rows
cdr <- cdr[complete.cases(cdr),]
moca <- moca[complete.cases(moca),]

# Clean up cdr and moca
cdr <- cdr[,c("ADRC_ID", "Age_at_dementia_scoring", "Clinical.Dementia.Rating..CDR..CDRSUM")]
names(cdr) <- c("ADRC_ID", "Age_CDR", "CDR")
moca <- moca[,c("ADRC_ID", "Age_at_MOCA_scoring", "MOCATOTS")]
names(moca) <- c("ADRC_ID", "Age_MOCA", "MoCA")

# Merge onto data
data <- merge(data, cdr, all.x = TRUE, by = "ADRC_ID")
data <- merge(data, moca, all.x = TRUE, by = "ADRC_ID")

# Remove columns not of interest
data <- data[,-c(1:6,10:11)]

# Rename columns
names(data) <- c("Age_GENE", "Sex", "Diagnosis",
                 "CXCL16", "ApoE", "ApoE3", "ApoE4", "ApoE2",
                 "PLTP", "APOC1", "NEFL", "GFAP", "TREM2", "UCHL1",
                 "Age_CDR", "CDR", "Age_MOCA", "MoCA")

# Convert Sex to binary
data$Sex <- mapvalues(data$Sex,
                      from = c("F", "M"),
                      to = c(0,1))
data$Sex <- as.numeric(data$Sex)

# Fix dumb "MCI " value
data[which(data$Diagnosis == "MCI "), "Diagnosis"] <- "MCI"

# Define diagnoses
diagnoses <- list("HC", "MCI", "AD", c("MCI", "AD"))
diagnoses_labels <- c("HC", "MCI", "AD", "CI")

# Subset for diagnoses of interest
data <- data[which(data$Diagnosis %in% diagnoses[1:3]),]

# Run ANCOVA
data$Diagnosis <- factor(data$Diagnosis, levels = diagnoses[1:3])
res.aov <- aov(CXCL16 ~ Diagnosis + Age_GENE + Sex, data = data)

# Generate ANCOVA table
anova_table <- Anova(res.aov)
write.csv(anova_table, file = paste0(output_dir, "somalogical_HC_MCI_AD_anova_table.csv"))

# Run post-hoc test (Tukey)
postHocs <- glht(res.aov, linfct = mcp(Diagnosis = "Tukey"))
tukey_pvals <- data.frame(summary(postHocs)$test$pvalues)
tukey_pvals$Comparison <- names(summary(postHocs)$test$coefficients)
write.csv(tukey_pvals, file = paste0(output_dir, "somalogical_HC_MCI_AD_tukey_pvals.csv"))

#-------------------------------------------------------------------------------
# Run ANCOVA b/w HC and CI

# Run ANCOVA
data$Diagnosis_mod <- mapvalues(data$Diagnosis,
                                from = diagnoses[1:3],
                                to = c("HC", "CI", "CI"))
data$Diagnosis_mod <- factor(data$Diagnosis_mod, levels = c("HC", "CI"))
res.aov <- aov(CXCL16 ~ Diagnosis_mod + Age_GENE + Sex, data = data)

# Generate ANCOVA table
anova_table <- Anova(res.aov)
write.csv(anova_table, file = paste0(output_dir, "somalogical_HC_CI_anova_table.csv"))

# Run post-hoc test (Tukey)
postHocs <- glht(res.aov, linfct = mcp(Diagnosis_mod = "Tukey"))
tukey_pvals <- data.frame(summary(postHocs)$test$pvalues)
tukey_pvals$Comparison <- names(summary(postHocs)$test$coefficients)

#-------------------------------------------------------------------------------
# Generate plot with all four groups (SuppFig 6B)

# Add duplicate rows for MCI and AD
data_ci <- data
data_ci$Diagnosis <- data_ci$Diagnosis_mod
data_ci <- data_ci[ which(data_ci$Diagnosis == "CI"),]
data <- rbind(data, data_ci)

# Generate boxplot
p_box <- ggplot(data, aes(x = Diagnosis, y = CXCL16)) +
  geom_boxplot(aes(fill = Diagnosis), fatten = 1, width = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 2) +
  geom_vline(xintercept = 3.5, color = "gray") +
  theme_Publication_blank() +
  scale_fill_manual(values = c("gray", "pink", "red", "lightblue")) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  labs(y = "SomaLogical CXCL16",
       x = NULL) +
  theme(text = element_text(size = 25),
        legend.title=element_blank(),
        legend.position="right",
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 8))

# Export plot
set_panel_size(p_box, file = paste0(output_dir, "somalogical_cxcl16_boxplot.pdf"),
               width = unit(4, "in"), height = unit(4, "in"))