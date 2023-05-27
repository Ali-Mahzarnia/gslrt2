options(java.parameters = "-Xmx2048m")  


library(xlsx)


AB_path = '/Users/ali/Desktop/april23/GSLRT/van/amyloid_mf_intrc_apoe_no_PMI.xlsx'
AB = as.data.frame(read.xlsx2(AB_path, sheetIndex = 1))
LT_path = '/Users/ali/Desktop/april23/GSLRT/van/log_tangles_mf_intrc_apoe_no_PMI.xlsx'
LT = as.data.frame(read.xlsx2(LT_path, sheetIndex = 1))
ME_path = '/Users/ali/Desktop/april23/GSLRT/van/mmse_mf_intrc_apoe_no_PMI.xlsx'
ME = as.data.frame(read.xlsx2(ME_path, sheetIndex = 1))


# # filtered :
 AB = AB [ abs(as.numeric(AB$s))>2, ] 
 AB$s =as.numeric(AB$s)
 AB$GS =gsub( "_", " ", AB$GS)
 LT = LT [ abs(as.numeric(LT$s))>1.5, ] 
 LT$s =as.numeric(LT$s)
 LT$GS =gsub( "_", " ", LT$GS)
 ME = ME [ abs(as.numeric(ME$s))>3, ] 
 ME$s =as.numeric(ME$s)
 ME$GS =gsub( "_", " ", ME$GS)
 
 



library(GOplot)

 
 AB$adjPvalue <- ifelse(AB$p <= 0.1, "significant", "non-significant")
 cols <- c("non-significant" = "grey", "significant" = "red")
 ggplot(AB, aes(reorder(GS, s), s, fill = adjPvalue)) +
   geom_col() +
   scale_fill_manual(values = cols) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
         panel.background = element_rect(fill = "white")) +
   coord_flip() +
   labs(x="Pathway", y="Normalized Enrichment Score (Test statistcs)",
        title="Hallmark pathways Enrichment Score from Amyloid beta") 

 ggsave( "AB.pdf", plot = last_plot(), device='pdf', scale=1, width=20, height=10, unit=c("in"), dpi=200)
 
LT$adjPvalue <- ifelse(LT$p <= 0.1, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(LT, aes(reorder(GS, s), s, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.background = element_rect(fill = "white")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score(Test statistcs)",
       title="Hallmark pathways Enrichment Score from Log(Tangles)")
ggsave( "Log_tangles.pdf", plot = last_plot(), device='pdf', scale=1, width=20, height=10, unit=c("in"), dpi=200)




ME$adjPvalue <- ifelse(ME$p <= 0.1, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(ME, aes(reorder(GS, s), s, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.background = element_rect(fill = "white")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score(Test statistcs)",
       title="Hallmark pathways Enrichment Score from MMSE")

ggsave( "MMSE.pdf", plot = last_plot(), device='pdf', scale=1, width=20, height=10, unit=c("in"), dpi=200)






