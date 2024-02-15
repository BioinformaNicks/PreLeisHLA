#################################################################
##                    Loading in the packages                   #
#################################################################

suppressMessages(library(here))
suppressMessages(library(readr))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(gplots))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(hrbrthemes))
suppressMessages(library(extrafont))
suppressMessages(library(Cairo))
suppressMessages(library(midasHLA))

here::i_am("scripts/PreLeisHLA.R")

text_font <- 'Roboto Condensed'

#################################################################
##                      Loading in the data                     #
#################################################################

PreLeisHLA_Omixon <- read_excel(here("data", "raw_data", "NanoTYPE", "231201_PreLeisHLA_Final_Results.xlsx"))

colnames(PreLeisHLA_Omixon) <- str_replace(colnames(PreLeisHLA_Omixon), "HLA-", "")
PreLeisHLA_Omixon <- PreLeisHLA_Omixon %>% mutate_all(~str_replace(., "HLA-", ""))

# SpatialCL_HLA <- read_excel(here("Analyses", "omixon", "SpatialCL_HLA", "results", "230425_Intermediate_10_Sample_Results.xlsx"))
# 
# colnames(SpatialCL_HLA) <- str_replace(colnames(SpatialCL_HLA), "HLA-", "")
# SpatialCL_HLA <- SpatialCL_HLA %>% mutate_all(~str_replace(., "HLA-", ""))

PreLeisHLA_CareDx <- read_csv(here("data", "raw_data", "AlloSeq", "230421_AlloSeq_Simplified.csv"))

PreLeisHLA_CareDx <- mutate_all(PreLeisHLA_CareDx, as.character)
colnames(PreLeisHLA_CareDx) <- str_replace(colnames(PreLeisHLA_CareDx), "HLA-", "")
PreLeisHLA_CareDx$Allele <- str_replace(PreLeisHLA_CareDx$Allele, "Allele ", "")

for(x in colnames(PreLeisHLA_CareDx)[-c(1,2)]) {
  PreLeisHLA_CareDx[[x]] <- str_c(x, PreLeisHLA_CareDx[[x]], sep = "*")
}

PreLeisHLA_CareDx <- PreLeisHLA_CareDx %>% mutate(across(colnames(PreLeisHLA_CareDx)[-c(1,2)], ~str_extract(., "^[A+-Z+]+(\\*\\d+(:\\d+){0,2})")))

AFND_Amhara <- read_csv(here("data", "raw_data", "AFDB_Amhara_Population.csv"))

##################################################################
##                  Pivot tables to MiDAS format                 #
##################################################################

PreLeisHLA_Omixon <- pivot_wider(PreLeisHLA_Omixon, names_from = Allele, values_from = colnames(PreLeisHLA_Omixon)[-c(1:4)])
# SpatialCL_HLA <- pivot_wider(SpatialCL_HLA, names_from = Allele, values_from = colnames(SpatialCL_HLA)[-c(1,2)])
PreLeisHLA_CareDx <- pivot_wider(PreLeisHLA_CareDx, names_from = Allele, values_from = colnames(PreLeisHLA_CareDx)[-c(1:2)])

colnames(PreLeisHLA_Omixon)[1] <- "ID"
# colnames(SpatialCL_HLA)[1] <- "ID"
colnames(PreLeisHLA_CareDx)[1] <- "ID"

PreLeisHLA_Omixon <- PreLeisHLA_Omixon %>% mutate(across(colnames(PreLeisHLA_Omixon), ~str_remove(., "#1")))

#################################################################
##                      AlloSeq vs NanoTYPE                     #
#################################################################

Omixon_1st_field <- reduceHlaCalls(PreLeisHLA_Omixon[-c(2:3)], resolution = 2)
Omixon_2nd_field <- reduceHlaCalls(PreLeisHLA_Omixon[-c(2:3)], resolution = 4)
Omixon_3rd_field <- reduceHlaCalls(PreLeisHLA_Omixon[-c(2:3)], resolution = 6)

CareDx_1st_field <- reduceHlaCalls(PreLeisHLA_CareDx, resolution = 2)
CareDx_2nd_field <- reduceHlaCalls(PreLeisHLA_CareDx, resolution = 4)
CareDx_3rd_field <- reduceHlaCalls(PreLeisHLA_CareDx, resolution = 6)

#################################################################
##                        Get frequencies                       #
#################################################################

##----------------------------------------------------------------
##                  General population vs AFND                   -
##----------------------------------------------------------------

Omixon_2nd_field_freqs <- getHlaFrequencies(hla_calls = Omixon_2nd_field)
Omixon_3rd_field_freqs <- getHlaFrequencies(hla_calls = Omixon_3rd_field)

# SpatialCL_2nd_field <- reduceHlaCalls(SpatialCL_HLA, resolution = 4)
# SpatialCL_2nd_field_freqs <- getHlaFrequencies(hla_calls = SpatialCL_2nd_field)

AFND_versus_Omixon <- full_join(AFND_Amhara, Omixon_2nd_field_freqs, by = 'allele')
colnames(AFND_versus_Omixon) <- c("allele", "AFND_Amhara", "Counts_PreLeisH", "PreLeisH_Cohort")
# AFND_versus_Omixon <- full_join(AFND_versus_Omixon, SpatialCL_2nd_field_freqs, by = 'allele')
# colnames(AFND_versus_Omixon) <- c("allele", "AFND_Amhara", "Counts_PreLeisH", "PreLeisH_Cohort", "Counts_SCL", "SpatialCL")
# AFND_versus_Omixon <- mutate(AFND_versus_Omixon, across(c("AFND_Amhara", "Counts_PreLeisH", "PreLeisH_Cohort", "Counts_SCL", "SpatialCL"), ~replace_na(.x, 0)))
# AFND_versus_Omixon <- pivot_longer(AFND_versus_Omixon, cols = c("AFND_Amhara", "PreLeisH_Cohort", "SpatialCL"), names_to = 'Population', values_to = 'Freq')
AFND_versus_Omixon <- mutate(AFND_versus_Omixon, across(c("AFND_Amhara", "Counts_PreLeisH", "PreLeisH_Cohort"), ~replace_na(.x, 0)))
AFND_versus_Omixon <- pivot_longer(AFND_versus_Omixon, cols = c("AFND_Amhara", "PreLeisH_Cohort"), names_to = 'Population', values_to = 'Freq')
# AFND_versus_Omixon <- AFND_versus_Omixon %>% filter(!grepl("^DQA1*03", allele)) %>% filter(!grepl("^DQB1*02", allele))

AFND_versus_Omixon_HLAI <- filter(AFND_versus_Omixon, grepl("^[A|B|C]", allele))
AFND_versus_Omixon_HLAII <- filter(AFND_versus_Omixon, grepl("^(DP[AB]|DQ[AB]1|DRB1)", allele))

HLA_A_allelefreq <- filter(Omixon_2nd_field_freqs, grepl("^(A)", allele))
HLA_B_allelefreq <- filter(Omixon_2nd_field_freqs, grepl("^(B)", allele))
HLA_C_allelefreq <- filter(Omixon_2nd_field_freqs, grepl("^(C)", allele))
HLA_DPA_allelefreq <- filter(Omixon_2nd_field_freqs, grepl("^(DPA)", allele))
HLA_DPB_allelefreq <- filter(Omixon_2nd_field_freqs, grepl("^(DPB)", allele))
HLA_DQA_alellefreq <- filter(AFND_versus_Omixon, grepl("^(DQA)", allele))
HLA_DQB_alellefreq <- filter(AFND_versus_Omixon, grepl("^(DQB)", allele))
HLA_DRB1_alellefreq <- filter(AFND_versus_Omixon, grepl("^(DRB1)", allele))
HLA_DRB345_alellefreq <- filter(Omixon_2nd_field_freqs, grepl("^(DRB[3|4|5])", allele))

##----------------------------------------------------------------
##                    VL History vs no History                   -
##----------------------------------------------------------------

VL_history <- filter(PreLeisHLA_Omixon, VL_History == "Yes")
No_VL_history <- filter(PreLeisHLA_Omixon, VL_History == "No")

phenotype_dat <- select(PreLeisHLA_Omixon, c("ID", "VL_History")) %>% mutate(VL_History = case_when(VL_History == "Yes" ~ 1,
                                                                                                    VL_History == "No" ~ 0))
test <- prepareMiDAS(hla_calls = reduceHlaCalls(PreLeisHLA_Omixon[-c(2:3)], resolution = 4),
                     colData = phenotype_dat,
                     experiment = "hla_alleles")

VL_History_freqs <- getHlaFrequencies(hla_calls = reduceHlaCalls(VL_history[-c(2:3)], resolution = 4))
No_VL_History_freqs <- getHlaFrequencies(hla_calls = reduceHlaCalls(No_VL_history[-c(2:3)], resolution = 4))

VL_History_vs_No_VL_History <- full_join(VL_History_freqs, No_VL_History_freqs, by = 'allele')
colnames(VL_History_vs_No_VL_History) <- c("allele", "Counts_VLHistory", "VL_History", "Counts_NoVLHistory", "No_VL_History")
VL_History_vs_No_VL_History <- mutate(VL_History_vs_No_VL_History, across(c("Counts_VLHistory", "VL_History", "Counts_NoVLHistory", "No_VL_History"), ~replace_na(.x, 0)))
VL_History_vs_No_VL_History <- pivot_longer(VL_History_vs_No_VL_History, cols = c("VL_History", "No_VL_History"), names_to = 'Population', values_to = 'Freq')

# test_AA <- prepareMiDAS(hla_calls = PreLeisHLA_Omixon[-c(2:3)],
#                         colData = phenotype_dat,
#                         experiment = "hla_aa")

##################################################################
##                            Plotting                           #
##################################################################

##----------------------------------------------------------------
##                General population frequencies                 -
##----------------------------------------------------------------

HLA_freq_plotter <- function(HLA_gene_freqs, HLA_gene, xlims, break_value, allele_text_size) {
  
  if (HLA_gene %in% c("A", "B", "C", "DPA1", "DPB1", "DRB3/4/5")) {
  HLA_freq_plot <- ggplot(data = HLA_gene_freqs, aes(x = allele, y = Freq)) +
    geom_bar(stat = "identity", width = 0.7, colour = "black", fill = "#FF4B20") +
    ggtitle(paste0("HLA-", HLA_gene, " ", "allele frequencies")) +
    labs(x = paste0("HLA-", HLA_gene, " ", "Allele"), y = "Allele Frequency") +
    coord_flip() + 
    scale_y_continuous(labels = scales::percent, limits = c(0, xlims), breaks = seq(0, xlims, by = break_value), expand = c(0, 0)) + 
    scale_x_discrete(limits = rev) +
    theme_classic() + 
    theme(axis.text.y = element_text(size = allele_text_size, family = text_font),
          axis.title.y = element_text(size = 14, family = text_font, face = 'bold',),
          axis.text.x = element_text(color='black', size=12, family=text_font),
          axis.title.x = element_text(color='black', size=14, family=text_font, face = 'bold',),
          plot.title = element_text(size = 16, family = text_font, face = 'bold', hjust = 0.5),
          legend.text = element_text(size = 14, family = text_font),
          legend.position="bottom")
  } else {
    HLA_freq_plot <- ggplot(data = HLA_gene_freqs, aes(x = allele, y = Freq, fill = Population)) +
      geom_bar(stat = "identity", position = ggplot2::position_dodge(0.7), width = 0.7,colour = "black") +
      ggtitle(paste0("HLA-", HLA_gene, " ", "allele frequencies compared to AFND")) +
      labs(x = paste0("HLA-", HLA_gene, " ", "Allele"), y = "Allele Frequency") +
      coord_flip() + 
      scale_y_continuous(labels = scales::percent, limits = c(0, xlims), breaks = seq(0, xlims, by = break_value), expand = c(0, 0)) + 
      scale_x_discrete(limits = rev) +
      # scale_fill_manual(name = '', values = c("#C6FDEC","#FF4B20", "#7AC5FF")) +
      scale_fill_manual(name = '', values = c("#7AC5FF","#FF4B20")) +
      theme_classic() + 
      theme(axis.text.y = element_text(size = allele_text_size, family = text_font),
            axis.title.y = element_text(size = 14, family = text_font, face = 'bold',),
            axis.text.x = element_text(color='black', size=12, family=text_font),
            axis.title.x = element_text(color='black', size=14, family=text_font, face = 'bold',),
            plot.title = element_text(size = 16, family = text_font, face = 'bold', hjust = 0.5),
            legend.text = element_text(size = 14, family = text_font),
            legend.position="bottom")
  }
  return(HLA_freq_plot)
}

arranged_plot <- plot_grid(HLA_freq_plotter(HLA_A_allelefreq, "A", 0.21, 0.025, 10),
    HLA_freq_plotter(HLA_B_allelefreq, "B", 0.16, 0.025, 10),
    HLA_freq_plotter(HLA_C_allelefreq, "C", 0.265, 0.025, 10),
    HLA_freq_plotter(HLA_DPA_allelefreq, "DPA1", 0.65, 0.1, 11),
    HLA_freq_plotter(HLA_DPB_allelefreq, "DPB1", 0.325, 0.05, 11),
    HLA_freq_plotter(HLA_DRB345_alellefreq, "DRB3/4/5", 0.25, 0.05, 11),
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol = 3, nrow = 2
) + theme(plot.tag = element_text(size = 54, family=text_font, face = 'bold'))

ggsave(here("analyses", "final_output", "supplementary", "SuppFig1_PreLeisHLA_AlleleFreqs.pdf"), arranged_plot, dev = cairo_pdf, height = 2160, width = 3840, units = 'px', dpi = 200)
ggsave(here("analyses", "final_output", "supplementary", "SuppFig1_PreLeisHLA_AlleleFreqs.png"), arranged_plot, type = 'cairo-png', height = 2160, width = 3840, units = 'px', dpi = 200)

##---------------------------------------------------------------
##                        AFND_Separate                         -
##---------------------------------------------------------------

arranged_plot_AFND <- plot_grid(HLA_freq_plotter(HLA_DRB1_alellefreq, "DRB1", 0.275, 0.05, 12),
                           HLA_freq_plotter(HLA_DQA_alellefreq, "DQA1", 0.375, 0.05, 12),
                           HLA_freq_plotter(HLA_DQB_alellefreq, "DQB1", 0.375, 0.05, 12),
                           labels = c("A", "B", "C"),
                           ncol = 3, nrow = 1
) + theme(plot.tag = element_text(size = 60, family=text_font, face = 'bold'))

ggsave(here("analyses", "final_output", "Fig3_PreLeisHLA_AlleleFreqs.pdf"), arranged_plot_AFND, dev = cairo_pdf, height = 1440, width = 3840, units = 'px', dpi = 200)
ggsave(here("analyses", "final_output", "Fig3_PreLeisHLA_AlleleFreqs.png"), arranged_plot_AFND, type = 'cairo-png', height = 1440, width = 3840, units = 'px', dpi = 200)


##---------------------------------------------------------------
##                  HLA-I and HLA-II separate                   -
##---------------------------------------------------------------

plot_HLAIfreq <- ggplot(data = AFND_versus_Omixon_HLAI, aes(x = allele, y = Freq)) +
  geom_bar(stat = "identity", width = 0.7, colour = "black", fill = "#FF4B20") +
  # ggtitle("HLA class I Allele frequencies of PreLeisH versus SpatialCL cohorts") +
  ggtitle("HLA class I allele frequencies") +
  labs(x = "HLA Allele", y = "Allele Frequency") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent, limits = c(0,0.3), breaks = seq(0,0.3, by = 0.025), expand = c(0, 0)) + 
  scale_x_discrete(limits = rev) +
  # scale_fill_manual(name = '', values = c("#C6FDEC","#FF4B20", "#7AC5FF")) +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 19, family = text_font),
        axis.title.y = element_text(size = 26, family = text_font, face = 'bold',),
        axis.text.x = element_text(color='black', size=20, family=text_font),
        axis.title.x = element_text(color='black', size=26, family=text_font, face = 'bold',),
        plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 26, family = text_font, face = 'bold'),
        legend.text = element_text(size = 26, family = text_font),
        legend.position="bottom")

plot_HLAIIfreq <- ggplot(data = AFND_versus_Omixon_HLAII, aes(x = allele, y = Freq, fill = Population)) +
  geom_bar(stat = "identity", position = ggplot2::position_dodge(0.7), width = 0.7,colour = "black") +
  # ggtitle("HLA class II Allele frequencies of PreLeisH and SpatialCL cohorts versus AFND") +
  ggtitle("HLA class II allele frequencies compared to AFND") +
  labs(x = "HLA Allele", y = "Allele Frequency") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent, limits = c(0,0.7), breaks = seq(0,0.7, by = 0.05), expand = c(0, 0)) + 
  scale_x_discrete(limits = rev) +
  # scale_fill_manual(name = '', values = c("#C6FDEC","#FF4B20", "#7AC5FF")) +
  scale_fill_manual(name = '', values = c("#7AC5FF","#FF4B20")) +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 19, family = text_font),
        axis.title.y = element_text(size = 26, family = text_font, face = 'bold',),
        axis.text.x = element_text(color='black', size=20, family=text_font),
        axis.title.x = element_text(color='black', size=26, family=text_font, face = 'bold',),
        plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 26, family = text_font, face = 'bold'),
        legend.text = element_text(size = 26, family = text_font),
        legend.position="bottom")

ggsave(here("analyses", "final_output", "PreLeisHLA_HLAI.pdf"), plot_HLAIfreq, dev = cairo_pdf, height = 17.5, width = 20)
ggsave(here("analyses", "final_output", "PreLeisHLA_HLAI.png"), plot_HLAIfreq, type = 'cairo-png', height = 17.5, width = 20)

ggsave(here("analyses", "final_output", "AFND_versus_PreLeisHLA_HLAII.pdf"), plot_HLAIIfreq, dev = cairo_pdf, height = 20, width = 20)
ggsave(here("analyses", "final_output", "AFND_versus_PreLeisHLA_HLAII.png"), plot_HLAIIfreq, type = 'cairo-png', height = 20, width = 20)

##----------------------------------------------------------------
##              VL History vs No History frequencies             -
##----------------------------------------------------------------

plot_VLHistory_Vs_NoHistory_HLAI <- ggplot(data = filter(VL_History_vs_No_VL_History, grepl("^[A|B|C]", allele)), aes(x = allele, y = Freq, fill = Population)) +
  geom_bar(stat = "identity", position = ggplot2::position_dodge(0.7), width = 0.7,colour = "black") +
  # ggtitle("HLA class I Allele frequencies of PreLeisH versus SpatialCL cohorts") +
  ggtitle("HLA class I Allele frequencies of participants with versus without VL History") +
  labs(x = "HLA Allele", y = "Allele Frequency") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent, limits = c(0,0.3), breaks = seq(0,0.3, by = 0.05), expand = c(0, 0)) + 
  scale_x_discrete(limits = rev) +
  # scale_fill_manual(name = '', values = c("#C6FDEC","#FF4B20", "#7AC5FF")) +
  scale_fill_manual(name = '', values = c("#7AC5FF","#FF4B20")) +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 17, family = text_font),
        axis.title.y = element_text(size = 26, family = text_font, face = 'bold',),
        axis.text.x = element_text(color='black', size=20, family=text_font),
        axis.title.x = element_text(color='black', size=26, family=text_font, face = 'bold',),
        plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 26, family = text_font, face = 'bold'),
        legend.text = element_text(size = 26, family = text_font),
        legend.position="bottom")

ggsave(here("analyses", "final_output", "PreLeisHLA_Hist_NoHist_HLAI.png"), plot_VLHistory_Vs_NoHistory_HLAI, type = 'cairo-png', height = 17.5, width = 20)
ggsave(here("analyses", "final_output", "PreLeisHLA_Hist_NoHist_HLAI.pdf"), plot_VLHistory_Vs_NoHistory_HLAI, dev = cairo_pdf, height = 17.5, width = 20)

plot_VLHistory_Vs_NoHistory_HLAII <- ggplot(data = filter(VL_History_vs_No_VL_History, grepl("^(DP[AB]|DQ[AB]1|DRB1)", allele)), aes(x = allele, y = Freq, fill = Population)) +
  geom_bar(stat = "identity", position = ggplot2::position_dodge(0.7), width = 0.7,colour = "black") +
  # ggtitle("HLA class I Allele frequencies of PreLeisH versus SpatialCL cohorts") +
  ggtitle("HLA class I Allele frequencies of participants with versus without VL History") +
  labs(x = "HLA Allele", y = "Allele Frequency") +
  coord_flip() + 
  scale_y_continuous(labels = scales::percent, limits = c(0,0.7), breaks = seq(0,0.7, by = 0.05), expand = c(0, 0)) + 
  scale_x_discrete(limits = rev) +
  # scale_fill_manual(name = '', values = c("#C6FDEC","#FF4B20", "#7AC5FF")) +
  scale_fill_manual(name = '', values = c("#7AC5FF","#FF4B20")) +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 17, family = text_font),
        axis.title.y = element_text(size = 26, family = text_font, face = 'bold',),
        axis.text.x = element_text(color='black', size=20, family=text_font),
        axis.title.x = element_text(color='black', size=26, family=text_font, face = 'bold',),
        plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5),
        legend.title = element_text(size = 26, family = text_font, face = 'bold'),
        legend.text = element_text(size = 26, family = text_font),
        legend.position="bottom")

ggsave(here("analyses", "final_output", "PreLeisHLA_Hist_NoHist_HLAII.png"), plot_VLHistory_Vs_NoHistory_HLAII, type = 'cairo-png', height = 17.5, width = 20)
ggsave(here("analyses", "final_output", "PreLeisHLA_Hist_NoHist_HLAII.pdf"), plot_VLHistory_Vs_NoHistory_HLAII, dev = cairo_pdf, height = 17.5, width = 20)

