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
suppressMessages(library(rstatix))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(gplots))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(hrbrthemes))
suppressMessages(library(extrafont))
suppressMessages(library(Cairo))
suppressMessages(library(UpSetR))
suppressMessages(library(ggplotify))

##################################################################
##                    Setting global variables                   #
##################################################################

text_font <- 'Roboto Condensed'

##################################################################
##                     UpsetPlot VL Developers                   #
##################################################################

"
rK39RDT: Yes
1                                                        
rK39RDT: Yes PCR: Yes
2                                                        
rK39RDT: Yes DAT: Yes
5 
rK39RDT: Yes rK39ELISA: Yes 
3 
rK39ELISA: Yes DAT: Yes
1                    
rK39RDT: Yes DAT: Yes PCR: Yes
1                                                        
rK39RDT: Yes KATEX: Yes DAT: Yes 
1                                                        
rK39RDT: Yes rK39ELISA: Yes DAT: Yes
50 
rK39RDT: Yes rK39ELISA: Yes DAT: Yes PCR: Yes
8                                                        
rK39RDT: Yes rK39ELISA: Yes KATEX: Yes DAT: Yes
4                                                        
rK39RDT: Yes rK39ELISA: Yes KATEX: Yes DAT: Yes PCR: Yes 
2 
"

intersections_VLdev <- c("rK39-ELISA" = 0,
                       "DAT" = 0,
                       "KATEX" = 0,
                       "PCR" = 0,
                       "rK39-RDT" = 1,
                       "rK39-RDT&rK39-ELISA" = 3,
                       "rK39-RDT&KATEX" = 0,
                       "rK39-RDT&DAT" = 5,
                       "rK39-RDT&PCR" = 2,
                       "rK39-ELISA&KATEX" = 0,
                       "rK39-ELISA&DAT" = 1,
                       "rK39-ELISA&PCR" = 0,
                       "KATEX&DAT" = 0,
                       "KATEX&PCR" = 0,
                       "DAT&PCR" = 0,
                       "rK39-RDT&DAT&PCR" = 1,
                       "rK39-RDT&DAT&KATEX" = 1,
                       "rK39-RDT&rK39-ELISA&KATEX" = 0,
                       "rK39-RDT&rK39-ELISA&DAT" = 50,
                       "rK39-RDT&rK39-ELISA&DAT&PCR" = 8,
                       "rK39-ELISA&KATEX&PCR" = 0,
                       "rK39-RDT&rK39-ELISA&KATEX&DAT" = 4,
                       "rK39-RDT&rK39-ELISA&KATEX&DAT&PCR" = 2,
                       "rK39-RDT&rK39-ELISA&PCR" = 0,
                       "rK39-ELISA&KATEX&DAT" = 0,
                       "KATEX&DAT&PCR" = 0,
                       "rK39-RDT&rK39-ELISA&KATEX&PCR" = 0,
                       "rK39-ELISA&KATEX&DAT&PCR" = 0
)

UpsetPlot_VLdev <- as.ggplot(upset(fromExpression(intersections_VLdev),
                                 nsets = 5, 
                                 order.by = "degree",
                                 decreasing = F,
                                 mb.ratio = c(0.6, 0.4),
                                 number.angles = 0, 
                                 text.scale = 2.5, 
                                 point.size = 8, 
                                 line.size = 1.25)) + 
  ggtitle("Leishmania infection marker positivity of past VL group") +
  theme(plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5))

# ggsave(here("analyses", "final_output", "Fig2_UpsetPlot_VLdev.png"), UpsetPlot_VLdev, type = 'cairo-png', height = 7.5, width = 15)
# ggsave(here("analyses", "final_output", "Fig2_UpsetPlot_VLdev.pdf"), UpsetPlot_VLdev, dev = cairo_pdf, height = 7.5, width = 15)


# intersections_all <- c("rK39-ELISA" = 0,
#                        "DAT" = 0,
#                        "KATEX" = 0,
#                        "PCR" = 0,
#                        "rK39-RDT" = 1,
#                        "rK39-RDT&rK39-ELISA" = 8,
#                        "rK39-RDT&KATEX" = 1,
#                        "rK39-RDT&DAT" = 14,
#                        "rK39-RDT&PCR" = 2,
#                        "rK39-ELISA&KATEX" = 1,
#                        "rK39-ELISA&DAT" = 5,
#                        "rK39-ELISA&PCR" = 0,
#                        "KATEX&DAT" = 0,
#                        "KATEX&PCR" = 0,
#                        "DAT&PCR" = 0,
#                        "rK39-RDT&DAT&PCR" = 1,
#                        "rK39-RDT&DAT&KATEX" = 1,
#                        "rK39-RDT&rK39-ELISA&KATEX" = 1,
#                        "rK39-RDT&rK39-ELISA&DAT" = 69,
#                        "rK39-RDT&rK39-ELISA&DAT&PCR" = 9,
#                        "rK39-ELISA&KATEX&PCR" = 1,
#                        "rK39-RDT&rK39-ELISA&KATEX&DAT" = 6,
#                        "rK39-RDT&rK39-ELISA&KATEX&DAT&PCR" = 4,
#                        "rK39-RDT&rK39-ELISA&PCR" = 0,
#                        "rK39-ELISA&KATEX&DAT" = 0,
#                        "KATEX&DAT&PCR" = 0,
#                        "rK39-RDT&rK39-ELISA&KATEX&PCR" = 0,
#                        "rK39-ELISA&KATEX&DAT&PCR" = 0
# )
# 
# UpsetPlot_All <- as.ggplot(upset(fromExpression(intersections_all),
#                                  nsets = 5, 
#                                  order.by = "degree",
#                                  decreasing = F,
#                                  mb.ratio = c(0.6, 0.4),
#                                  number.angles = 0, 
#                                  text.scale = 2.5, 
#                                  point.size = 8, 
#                                  line.size = 1.25)) + 
#   ggtitle("Leishmania infection marker positivity of all participants") +
#   theme(plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5))
# 
# # ggsave(here("analyses", "final_output", "Fig2_UpsetPlot_All.png"), UpsetPlot_All, type = 'cairo-png', height = 7.5, width = 15)
# # ggsave(here("analyses", "final_output", "Fig2_UpsetPlot_All.pdf"), UpsetPlot_All, dev = cairo_pdf, height = 7.5, width = 15)

#################################################################
##                    UpsetPlot Asymptomatics                   #
#################################################################

"
rK39-ELISA: Yes. KATEX: Yes
1
rK39-ELISA: Yes. DAT: Yes.
4
rK39-ELISA: Yes. KATEX: Yes. PCR: Yes.   
1
rK39-RDT: Yes. KATEX: Yes
1
rK39-RDT: Yes. DAT: Yes.
9             
rK39-RDT: Yes. rK39-ELISA: Yes
5
rK39-RDT: Yes.rK39-ELISA: Yes.KATEX: Yes
1                         
rK39-RDT: Yes.rK39-ELISA: Yes.DAT: Yes.
19             
rK39-RDT: Yes. rK39-ELISA: Yes. DAT: Yes. PCR: Yes. 
1 

rK39-RDT: Yes. rK39-ELISA: Yes. DAT: Yes. KATEX: Yes
2                                                           
rK39-RDT: Yes. rK39-ELISA: Yes. DAT: Yes. KATEX: Yes. PCR: Yes. 
2 "

intersections_asymptomatics <- c("rK39-ELISA" = 0,
                                 "DAT" = 0,
                                 "KATEX" = 0,
                                 "PCR" = 0,
                                 "rK39-RDT" = 0,
                                 "rK39-RDT&rK39-ELISA" = 5,
                                 "rK39-RDT&KATEX" = 1,
                                 "rK39-RDT&DAT" = 9,
                                 "rK39-RDT&PCR" = 0,
                                 "rK39-ELISA&KATEX" = 1,
                                 "rK39-ELISA&DAT" = 4,
                                 "rK39-ELISA&PCR" = 0,
                                 "KATEX&DAT" = 0,
                                 "KATEX&PCR" = 0,
                                 "DAT&PCR" = 0,
                                 "rK39-RDT&DAT&PCR" = 0,
                                 "rK39-RDT&DAT&KATEX" = 0,
                                 "rK39-RDT&rK39-ELISA&KATEX" = 1,
                                 "rK39-RDT&rK39-ELISA&DAT" = 19,
                                 "rK39-RDT&rK39-ELISA&DAT&PCR" = 1,
                                 "rK39-ELISA&KATEX&PCR" = 1,
                                 "rK39-RDT&rK39-ELISA&KATEX&DAT" = 2,
                                 "rK39-RDT&rK39-ELISA&KATEX&DAT&PCR" = 2,
                                 "rK39-RDT&rK39-ELISA&PCR" = 0,
                                 "rK39-ELISA&KATEX&DAT" = 0,
                                 "KATEX&DAT&PCR" = 0,
                                 "rK39-RDT&rK39-ELISA&KATEX&PCR" = 0,
                                 "rK39-ELISA&KATEX&DAT&PCR" = 0
)

UpsetPlot_Asymptomatics <- as.ggplot(upset(fromExpression(intersections_asymptomatics),
                                           nsets = 5, 
                                           order.by = "degree",
                                           decreasing = F,
                                           mb.ratio = c(0.6, 0.4),
                                           number.angles = 0, 
                                           text.scale = 2.5, 
                                           point.size = 8, 
                                           line.size = 1.25)) + 
  ggtitle("Leishmania infection marker positivity of asymptomatic Leishmania controllers") +
  theme(plot.title = element_text(size = 32, family = text_font, face = 'bold', hjust = 0.5))

# ggsave(here("analyses", "final_output", "Fig3_UpsetPlot_Asymptomatics.png"), UpsetPlot_Asymptomatics, type = 'cairo-png', height = 7.5, width = 15)
# ggsave(here("analyses", "final_output", "Fig3_UpsetPlot_Asymptomatics.pdf"), UpsetPlot_Asymptomatics, dev = cairo_pdf, height = 7.5, width = 15)

#################################################################
##                        Combined plot                         #
#################################################################

layout <- "
AAAA
BBBB
"

arranged_plot <- UpsetPlot_Asymptomatics + UpsetPlot_VLdev +
  #celltype_barchart_cross_sectional + celltype_barchart_longitudinal +
  plot_layout(design = layout) & plot_annotation(tag_levels = c('A')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", "Fig2_UpsetPlot_PreLeisHLA.png"), arranged_plot, type = 'cairo-png', height = 15, width = 17.5)
ggsave(here("analyses", "final_output", "Fig2_UpsetPlot_PreLeisHLA.pdf"), arranged_plot, dev = cairo_pdf, height = 15, width = 17.5)
