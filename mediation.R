require(magrittr)
require(plyr)
require(tidyverse)
library(mediation)

G=read.delim( "ROSMAP rnaseq matrix G.txt")
names=colnames(G)
G=matrix(unlist(G), ncol = dim(G)[2], nrow =dim(G)[1])
colnames(G) = names



Zab <- read.delim( "1Z_amyloid_continuous_outcome.txt" , header = F)
Zab=matrix(unlist(Zab), ncol = dim(Zab)[2], nrow =dim(Zab)[1])
Zab= as.double(Zab) 

Ztau <- read.delim( "2Z_log_tangles_continuous_outcome.txt" , header = F)
Ztau=matrix(unlist(Ztau), ncol = dim(Ztau)[2], nrow =dim(Ztau)[1])
Ztau= as.double(Ztau) 


Zmmse <- read.delim( "3Z_mmse_continous.txt" , header = F)
Zmmse=matrix(unlist(Zmmse), ncol = dim(Zmmse)[2], nrow =dim(Zmmse)[1])
Zmmse= as.double(Zmmse) 


X <- read.delim( "ROSMAP rnaseq matrix X1.txt" , header = F) %>% select(-c('V5')) #no pmi or last column to be read
X=matrix(unlist(X), ncol = dim(X)[2], nrow =dim(X)[1])
colnames(X) = c( "msex", "education", "apoe", "age_death")
####
#### recode apoe of X matrix
apoe=X[,3];
#22-> 0 , 23-> 1, 24-> 2,  33-> 1,  34-> 2,  44-> 2
apoe[apoe==22]=0; apoe[apoe==23]=1; apoe[apoe==24]=2; apoe[apoe==33]=1; apoe[apoe==34]=2; apoe[apoe==44]=2;
X[,3]=apoe
####recode apoe of X matrix





library(xlsx)





gene_sets_read <- read.delim("Bc5.go.mf.v2023.1.Hs.symbols.txt" , header = F)
VEGF_ind = which(gene_sets_read[,1] == 'GOMF_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_RECEPTOR_2_BINDING')
VEGF = gene_sets_read[VEGF_ind,1:10 ]
VEGF_genes = VEGF[2:10]

index_numerical_values = which( colnames(G) %in% unlist(VEGF_genes)   )
G2 = G[,index_numerical_values]


covars = c("msex", "education", "apoe", "age_death")
mediators = c( "CCDC88A"  ,  "PDCL3"  ,  "VEGFA"  , "DAB2IP"    , "ITGA5"   ,  "GREM1"  ,   "CDH5"  ,  "CADM4")




pvals= matrix(NA, 1, 8)
data3 = as.data.frame(na.omit(cbind(Zab, G2,Zmmse )))
summary =vector(mode = "list", length = 8)
mediation_analysis = vector(mode = "list", length = 8)
for (j in 1:8) {
  

# j=5
med = G2[,j]
data4 = as.data.frame(na.omit(cbind(Zab, med,Zmmse )))


 model_0 = lm ( Zmmse ~ Zab , data=data4 )
 summary(model_0)
 
 
 #step 2 predict medi by variable of interest
 model_m = lm( med ~ Zmmse ,  data=data4 )
 summary(model_m)
 
 #step 3 full model
 model_full = lm( Zmmse~  med +Zab ,  data=data4) 
 summary(model_full)


mediation_analysis[[j]]= mediate(model_m,model_full, treat = 'Zmmse', mediator ='med' , boot = T, sims = 1000 )
summary[[j]]  = summary(mediation_analysis[[j]]) 
 pvals[j]= summary[[j]]$d1.p
}


pval_adj = p.adjust(pvals, method = "fdr")
index_sig = which(pval_adj<0.05)
plot(mediation_analysis[[index_sig]])
summary[[index_sig]]





summary =vector(mode = "list", length = 8)
mediation_analysis = vector(mode = "list", length = 8)
for (j in index_sig) {
  med = G2[,j]
  data4 = as.data.frame(na.omit(cbind(Zab, med,Zmmse )))
  
  
  model_0 = lm ( Zmmse ~ Zab , data=data4 )
  summary(model_0)
  
  
  #step 2 predict medi by variable of interest
  model_m = lm( med ~ Zmmse ,  data=data4 )
  summary(model_m)
  
  #step 3 full model
  model_full = lm( Zmmse~  med +Zab ,  data=data4) 
  summary(model_full)

  
  mediation_analysis[[j]]= mediate(model_m,model_full, treat = 'Zmmse', mediator ='med' , boot = T, sims = 1000 )
  summary[[j]]  = summary(mediation_analysis[[j]])  
}


# library(flexplot)
# visualize(model_full)
# 
# mediate_plot(Zmmse~  med +Zab ,  data=data4)







############################### 


#j=1
# med = G2[,j]
# 
# data = cbind(Z, X, med)
# data = as.data.frame(na.omit(data))
# data2 =  as.data.frame(na.omit(cbind(Z, X, med)))

pvals2= matrix(NA, 1, 8)
data32 = as.data.frame(na.omit(cbind(Ztau, G2,Zmmse )))



summary2 =vector(mode = "list", length = 8)
mediation_analysis2 = vector(mode = "list", length = 8)
for (j in 1:8) {
  
  
  # j=5
  med = G2[,j]
  data4 = as.data.frame(na.omit(cbind(Ztau, med,Zmmse )))
  
  
  model_0 = lm ( Zmmse ~ Ztau , data=data4 )
  summary(model_0)
  
  
  #step 2 predict medi by variable of interest
  model_m = lm( med ~ Zmmse ,  data=data4 )
  summary(model_m)
  
  #step 3 full model
  model_full = lm( Zmmse~  med +Ztau ,  data=data4) 
  summary(model_full)
  
  
  mediation_analysis2[[j]]= mediate(model_m,model_full, treat = 'Zmmse', mediator ='med' , boot = T, sims = 1000 )
  summary2[[j]]  = summary(  mediation_analysis2[[j]] ) 
  pvals2[j]= summary2[[j]]$d1.p
}






pval_adj2 = p.adjust(pvals2, method = "fdr")
index_sig = which(pval_adj2<0.05)
plot(mediation_analysis2[[index_sig]])
summary2[[index_sig]]

for (j in index_sig) {
  # j=5
  med = G2[,j]
  data4 = as.data.frame(na.omit(cbind(Ztau, med,Zmmse )))
  
  
  model_0 = lm ( Zmmse ~ Ztau , data=data4 )
  summary(model_0)
  
  
  #step 2 predict medi by variable of interest
  model_m = lm( med ~ Zmmse ,  data=data4 )
  summary(model_m)
  
  #step 3 full model
  model_full = lm( Zmmse~  med +Ztau ,  data=data4) 
  summary(model_full)
  visualize(model_full)
  
  
  mediation_analysis2= mediate(model_m,model_full, treat = 'Zmmse', mediator ='med' , boot = T, sims = 1000 )
  summary2  = summary(mediation_analysis2) 
  mediate_plot( Zmmse~  med +Ztau ,  data=data4)
  
  
}


library(flexplot)
visualize(model_full)

mediate_plot(Zmmse~  med +Zab ,  data=data4)







##### 1st pc of vegf as med

summary(princomp(G2))



partial_cor = ppcor::pcor( G2 )
 min(partial_cor$p.value[lower.tri( partial_cor$p.value, diag = F ) ])
max( abs( partial_cor$estimate[lower.tri( partial_cor$estimate, diag = F ) ]))
 
# library(ggdag)
# 
# tidy_ggdag <- dagify(
#   MMSE ~ Tau_or_Ab + g,
#   g ~ Tau_or_Ab,
#   exposure = "Tau_or_Ab", 
#   outcome = "g"
# ) %>%
#   tidy_dagitty()
# ggdag(tidy_ggdag) +
#   theme_dag()



library(diagram)
library(diagram)
data <- c(0, "", 0,
          0, 0, 0, 
          "", "", 0)
M<- matrix (nrow=3, ncol=3, byrow = TRUE, data=data)
plot<- plotmat (M, pos=c(1,2), 
                name= c( expression(g^j),  expression(paste("Tau or A",beta)), "MMSE"), 
                box.type = "rect", box.size = 0.12, box.prop=0.5,  curve=0)
ggsave( "Mediation.pdf", plot = plot, device='pdf', scale=1, width=24, height=10, unit=c("in"), dpi=200)
