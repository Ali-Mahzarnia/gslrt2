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
 AB$GS =gsub( "_", " ", AB$GS)
 AB$s =as.numeric(AB$s)
 AB$p =as.numeric(AB$p)
 AB$size =as.numeric(AB$size)
 
 LT = LT [ abs(as.numeric(LT$s))>1.5, ] 
 LT$GS =gsub( "_", " ", LT$GS)
 LT$s =as.numeric(LT$s)
 LT$p =as.numeric(LT$p)
 LT$size =as.numeric(LT$size)
 
 
 ME = ME [ abs(as.numeric(ME$s))>3, ]
 ME$GS =gsub( "_", " ", ME$GS)
 ME$s =as.numeric(ME$s)
 ME$p =as.numeric(ME$p)
 ME$size =as.numeric(ME$size)
 
 



library(GOplot)

 
 AB$Pvalue <- ifelse(AB$p <= 0.1, "significant", "non-significant")
 #cols <- c("non-significant" = "grey", "significant" = "red")
 cols <- c("non-significant" = "red", "significant" = "red")
 ggplot(AB, aes(reorder(GS, -log10(p) ), -log10(p) , fill = Pvalue)) +
   geom_col() +
   scale_fill_manual(values = cols) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
         panel.background = element_rect(fill = "white")) +
   coord_flip() +
   labs(x="Pathways", y="-Log10(p-value)",
        title="Gene Ontology (GO) molecular function pathways for Amyloid beta") +
   geom_hline(yintercept = 1)+
   geom_text(aes(y=1.2, label="-log10(0.1)", x=2), colour="black", angle=0, vjust = 1)+theme(text = element_text(size = 20))
 #geom_text(aes(x=200, label="the weak cars", y=20), colour="red", angle=90, vjust = -1, text=element_text(size=11))
 

 ggsave( "AB.pdf", plot = last_plot(), device='pdf', scale=1, width=20, height=10, unit=c("in"), dpi=200)
 
LT$Pvalue <- ifelse(LT$p <= 0.1, "significant", "non-significant")
#cols <- c("non-significant" = "grey", "significant" = "red")
cols <- c("non-significant" = "red", "significant" = "red")
ggplot(LT, aes(reorder(GS, -log10(p)), -log10(p), fill = Pvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.background = element_rect(fill = "white")) +
  coord_flip() +
  labs(x="Pathways", y="-Log10(p-value)",
       title="Gene Ontology (GO) molecular function pathways for Log(Tangles)")+
  geom_hline(yintercept = 1)+
  geom_text(aes(y=1.2, label="-log10(0.1)", x=2), colour="black", angle=0, vjust = 1)+theme(text = element_text(size = 20))

ggsave( "Log_tangles.pdf", plot = last_plot(), device='pdf', scale=1, width=20, height=10, unit=c("in"), dpi=200)




ME$Pvalue <- ifelse(ME$p <= 0.1, "significant", "non-significant")
#cols <- c("non-significant" = "grey", "significant" = "red")
cols <- c("non-significant" = "red", "significant" = "red")
ggplot(ME, aes(reorder(GS, -log10(p)), -log10(p), fill = Pvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.background = element_rect(fill = "white")) +
  coord_flip() +
  labs(x="Pathways", y="-Log10(p-value)",
       title="Gene Ontology (GO) molecular function pathways for MMSE")+
  geom_hline(yintercept = 1)+
  geom_text(aes(y=1.5, label="-log10(0.1)", x=2), colour="black", angle=0, vjust = 1)+theme(text = element_text(size = 20))


ggsave( "MMSE.pdf", plot = last_plot(), device='pdf', scale=1, width=24, height=10, unit=c("in"), dpi=200)






