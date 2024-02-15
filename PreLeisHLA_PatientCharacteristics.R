#################################################################
##                    Loading in the packages                   #
#################################################################

suppressMessages(library(here))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(pheatmap))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(gplots))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(hrbrthemes))
suppressMessages(library(extrafont))
suppressMessages(library(Cairo))

##################################################################
##                    Setting global variables                   #
##################################################################

options(ggrepel.max.overlaps = Inf)

text_font <- 'Roboto Condensed'

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#7AC5FF", "#C6FDEC", "#0348A6"))

#################################################################
##                      Loading in the data                     #
#################################################################

characteristics_data <- suppressMessages(read_csv(here("data", "meta_data", "240109_Patient_Characteristics.csv")))

characteristics_data$PreLeisHLA_ID <- as.factor(characteristics_data$PreLeisHLA_ID)
characteristics_data$Age <- as.numeric(characteristics_data$Age)
characteristics_data$Sex <- as.factor(characteristics_data$Sex)
characteristics_data$Occupation <- as.factor(characteristics_data$Occupation)
characteristics_data$BMI <- as.numeric(characteristics_data$BMI)
characteristics_data$VL_History <- as.factor(characteristics_data$VL_History)
characteristics_data$Months_From_Prev_VL <- as.numeric(characteristics_data$Months_From_Prev_VL)
characteristics_data$CD4_Baseline <- as.numeric(characteristics_data$CD4_Baseline)
characteristics_data$HIV_TimeOfVL <- as.factor(characteristics_data$HIV_TimeOfVL)
characteristics_data$ART_TimeOfVL <- as.factor(characteristics_data$ART_TimeOfVL)
characteristics_data$ART_Baseline <- as.factor(characteristics_data$ART_Baseline)
characteristics_data$ART_Regimen_Baseline <- as.factor(characteristics_data$ART_Regimen_Baseline)
characteristics_data$VL_Episodes_Baseline <- as.numeric(characteristics_data$VL_Episodes_Baseline)
characteristics_data$Concomitant_Disease_Present <- as.factor(characteristics_data$Concomitant_Disease_Present)

##################################################################
##                          Statistics                           #
##################################################################

##################################################################
##                        Total Socio-Demo                       #
##################################################################

Total_Age <- summary(characteristics_data$Age)
Total_Male <- characteristics_data %>% filter(Sex == 'Male') %>% nrow()
Total_BMI <- summary(characteristics_data$BMI)

Total_Occupation <- cbind("Count" = rowSums(table(characteristics_data$Occupation, characteristics_data$VL_History)),
                          "Proportion" = rowSums(table(characteristics_data$Occupation, characteristics_data$VL_History))/sum(table(characteristics_data$Occupation, characteristics_data$VL_History))*100)

#################################################################
##                        Group Socio-Demo                      #
#################################################################


##----------------------------------------------------------------
##                      Summary statistics                       -
##----------------------------------------------------------------

All_groups_Age <- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_age <- characteristics_data %>% filter(VL_History == x) %>% select(Age) %>% unlist(.) %>% summary(.)
  paste0(round(summary_age['Median'],1), ' (', round(summary_age['1st Qu.'],1), '-', round(summary_age['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Male <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_male <- characteristics_data %>% filter(VL_History == x) %>% filter(Sex == 'Male') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_male, ' (', round((num_male / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_BMI <- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_BMI <- characteristics_data %>% filter(VL_History == x) %>% select(BMI) %>% unlist(.) %>% summary(.)
  paste0(round(summary_BMI['Median'],1), ' (', round(summary_BMI['1st Qu.'],1), '-', round(summary_BMI['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Occupation <- table(characteristics_data$Occupation, characteristics_data$VL_History)
All_groups_Occupation <- cbind(All_groups_Occupation, (All_groups_Occupation[,1] / colSums(All_groups_Occupation)[1])*100, (All_groups_Occupation[,2] / colSums(All_groups_Occupation)[2])*100)

##---------------------------------------------------------------
##                      Statistical testing                     -
##---------------------------------------------------------------

All_groups_wilcox_Age <- wilcox.test(Age ~ VL_History, data = characteristics_data)$p.value

All_groups_fisher_Male <- fisher.test(table(characteristics_data$Sex, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()

All_groups_wilcox_BMI <- wilcox.test(BMI ~ VL_History, data = characteristics_data)$p.value

All_groups_fisher_Occupation <- fisher.test(table(characteristics_data$Occupation, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()

##################################################################
##                        Total Clinical                         #
##################################################################

Total_Months_Followed <- summary(characteristics_data$Months_Followed)

Total_VL_History <- table(characteristics_data$VL_History)
Total_Past_VL_Episodes <- summary(characteristics_data$VL_Episodes_Baseline)
Total_Months_Since_Past_VL <- summary(characteristics_data$Months_From_Prev_VL)

Total_ART <- table(characteristics_data$ART_Baseline)
Total_Concomitant_Disease_Present <- table(characteristics_data$Concomitant_Disease_Present)

Total_RK39RDT_Study <- table(characteristics_data$RK39RDT_Pos_Study)
Total_RK39ELISA_Study <- table(characteristics_data$RK39ELISA_Pos_Study)
Total_KATEX_Study <- table(characteristics_data$KATEX_Pos_Study)
Total_DAT_Study <- table(characteristics_data$DAT_Pos_Study)
Total_PCR_Study <- table(characteristics_data$PCR_Pos_Study)

Total_CD4_D0 <- summary(characteristics_data$CD4_Baseline)
Total_Lymphocytes_D0 <- summary(characteristics_data$Lymphocytes_Baseline)
Total_Platelet_D0 <- summary(characteristics_data$Platelet_Count_Baseline)
Total_HB_D0 <- summary(characteristics_data$Hemoglobin_Count_Baseline)

##----------------------------------------------------------------
##                      Summary statistics                       -
##----------------------------------------------------------------

All_groups_VL_History <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_VL_History <- characteristics_data %>% filter(VL_History == x) %>% filter(VL_History == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_VL_History, ' (', round((num_VL_History / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)


All_groups_episodes<- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_episodes <- characteristics_data %>% filter(VL_History == x) %>% select(VL_Episodes_Baseline) %>% unlist(.) %>% summary(.)
  paste0(round(summary_episodes['Median'],1), ' (', round(summary_episodes['1st Qu.'],1), '-', round(summary_episodes['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Months_Prev_VL <- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_Months_Prev_VL <- characteristics_data %>% filter(VL_History == x) %>% select(Months_From_Prev_VL) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Months_Prev_VL['Median'],1), ' (', round(summary_Months_Prev_VL['1st Qu.'],1), '-', round(summary_Months_Prev_VL['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_ART <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_ART <- characteristics_data %>% filter(VL_History == x) %>% filter(ART_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_ART, ' (', round((num_ART / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_ART_regime <- table(characteristics_data$ART_Regimen_Baseline, characteristics_data$VL_History)
All_groups_ART_regime <- cbind(All_groups_ART_regime, (All_groups_ART_regime[,1] / colSums(All_groups_ART_regime)[1])*100, (All_groups_ART_regime[,2] / colSums(All_groups_ART_regime)[2])*100)

All_groups_RK39RDT_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_RK39RDT_Pos_Baseline <- characteristics_data %>% filter(VL_History == x) %>% filter(RK39RDT_Pos_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_RK39RDT_Pos_Baseline, ' (', round((num_RK39RDT_Pos_Baseline / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_RK39ELISA_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_RK39ELISA_Pos_Baseline <- characteristics_data %>% filter(VL_History == x) %>% filter(RK39ELISA_Pos_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_RK39ELISA_Pos_Baseline, ' (', round((num_RK39ELISA_Pos_Baseline / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_KATEX_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_KATEX_Pos_Baseline <- characteristics_data %>% filter(VL_History == x) %>% filter(KATEX_Pos_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_KATEX_Pos_Baseline, ' (', round((num_KATEX_Pos_Baseline / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_DAT_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_DAT_Pos_Baseline <- characteristics_data %>% filter(VL_History == x) %>% filter(DAT_Pos_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_DAT_Pos_Baseline, ' (', round((num_DAT_Pos_Baseline / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_PCR_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_PCR_Pos_Baseline <- characteristics_data %>% filter(VL_History == x) %>% filter(PCR_Pos_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_PCR_Pos_Baseline, ' (', round((num_PCR_Pos_Baseline / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_RK39RDT_Study <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_RK39RDT_Pos_Study <- characteristics_data %>% filter(VL_History == x) %>% filter(RK39RDT_Pos_Baseline == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_RK39RDT_Pos_Study, ' (', round((num_RK39RDT_Pos_Study / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_RK39ELISA_Study <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_RK39ELISA_Pos_Study <- characteristics_data %>% filter(VL_History == x) %>% filter(RK39ELISA_Pos_Study == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_RK39ELISA_Pos_Study, ' (', round((num_RK39ELISA_Pos_Study / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_KATEX_Study <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_KATEX_Pos_Study <- characteristics_data %>% filter(VL_History == x) %>% filter(KATEX_Pos_Study == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_KATEX_Pos_Study, ' (', round((num_KATEX_Pos_Study / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_DAT_Study <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_DAT_Pos_Study <- characteristics_data %>% filter(VL_History == x) %>% filter(DAT_Pos_Study == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_DAT_Pos_Study, ' (', round((num_DAT_Pos_Study / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_PCR_Study <- sapply(levels(characteristics_data$VL_History), function(x) {
  num_PCR_Pos_Study <- characteristics_data %>% filter(VL_History == x) %>% filter(PCR_Pos_Study == 'Yes') %>% nrow()
  total_patients <- characteristics_data %>% filter(VL_History == x) %>% nrow()
  paste0(num_PCR_Pos_Study, ' (', round((num_PCR_Pos_Study / total_patients) * 100, 1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_CD4_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_CD4 <- characteristics_data %>% filter(VL_History == x) %>% select(CD4_Baseline) %>% unlist(.) %>% summary(.)
  paste0(round(summary_CD4['Median'],0), ' (', round(summary_CD4['1st Qu.'],0), '-', round(summary_CD4['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Lymphocytes_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_Lymphocytes <- characteristics_data %>% filter(VL_History == x) %>% select(Lymphocytes_Baseline) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Lymphocytes['Median'],2), ' (', round(summary_Lymphocytes['1st Qu.'],2), '-', round(summary_Lymphocytes['3rd Qu.'],2), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Platelets_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_Platelets <- characteristics_data %>% filter(VL_History == x) %>% select(Platelet_Count_Baseline) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Platelets['Median'],0), ' (', round(summary_Platelets['1st Qu.'],0), '-', round(summary_Platelets['3rd Qu.'],0), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

All_groups_Hemoglobin_Baseline <- sapply(levels(characteristics_data$VL_History), function(x) {
  summary_Hemoglobin <- characteristics_data %>% filter(VL_History == x) %>% select(Hemoglobin_Count_Baseline) %>% unlist(.) %>% summary(.)
  paste0(round(summary_Hemoglobin['Median'],1), ' (', round(summary_Hemoglobin['1st Qu.'],1), '-', round(summary_Hemoglobin['3rd Qu.'],1), ')')
}, simplify = FALSE, USE.NAMES = TRUE)

##---------------------------------------------------------------
##                      Statistical testing                     -
##---------------------------------------------------------------

All_groups_fisher_ART <- fisher.test(table(characteristics_data$ART_Baseline, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_ART_Regimen_Baseline <- fisher.test(table(characteristics_data$ART_Regimen_Baseline, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()

# All_groups_fisher_Concomitant_Diseases_Present <- fisher.test(table(characteristics_data$Concomitant_Disease_Present, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_RK39RDT_Pos_Baseline <- fisher.test(table(characteristics_data$RK39RDT_Pos_Baseline, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_RK39ELISA_Pos_Baseline <- fisher.test(table(characteristics_data$RK39ELISA_Pos_Baseline, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_KATEX_Pos_Baseline <- fisher.test(table(characteristics_data$KATEX_Pos_Baseline, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_DAT_Pos_Baseline <- fisher.test(table(characteristics_data$DAT_Pos_Baseline, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_PCR_Pos_Baseline <- fisher.test(table(characteristics_data$PCR_Pos_Baseline, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()

All_groups_fisher_RK39RDT_Pos_Study <- fisher.test(table(characteristics_data$RK39RDT_Pos_Study, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_RK39ELISA_Pos_Study <- fisher.test(table(characteristics_data$RK39ELISA_Pos_Study, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_KATEX_Pos_Study <- fisher.test(table(characteristics_data$KATEX_Pos_Study, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_DAT_Pos_Study <- fisher.test(table(characteristics_data$DAT_Pos_Study, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()
All_groups_fisher_PCR_Pos_Study <- fisher.test(table(characteristics_data$PCR_Pos_Study, characteristics_data$VL_History))['p.value'] %>% unlist() %>% unname()

All_groups_wilcox_CD4_Baseline <- wilcox.test(CD4_Baseline ~ VL_History, data = characteristics_data)$p.value
All_groups_wilcox_Lymphocytes_Baseline <- wilcox.test(Lymphocytes_Baseline ~ VL_History, data = characteristics_data)$p.value
All_groups_wilcox_Platelet_Baseline <- wilcox.test(Platelet_Count_Baseline ~ VL_History, data = characteristics_data)$p.value
All_groups_wilcox_HB_Baseline <- wilcox.test(Hemoglobin_Count_Baseline ~ VL_History, data = characteristics_data)$p.value

