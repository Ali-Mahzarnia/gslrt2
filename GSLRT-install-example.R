# install.packages("https://github.com/j-g-b/gsLRT/archive/master.tar.gz", repos = NULL, type="source")


require(magrittr)
require(plyr)
require(tidyverse)

## Read R objects
# GG <- readRDS("G.rds")
G=read.delim( "ROSMAP rnaseq matrix G.txt")
names=colnames(G)
G=matrix(unlist(G), ncol = dim(G)[2], nrow =dim(G)[1])
colnames(G) = names

# ZZ <- readRDS("Z.rds")
#Z <- read.delim( "ROSMAP rnaseq matrix z1.txt" , header = F)
#Z <- read.delim( "Z_braak_continuous_outcome.txt" , header = F)
#Z <- read.delim( "Z_cerad_binary_outcome.txt" , header = F)
#Z[Z==-1]=0

Z <- read.delim( "Z_braak_continuous_outcome.txt" , header = F)
Z=matrix(unlist(Z), ncol = dim(Z)[2], nrow =dim(Z)[1])
Z= as.double(Z) 



# XX <- readRDS("X.rds")
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
#write.xlsx(c(0,1), 'test.xlsx')
#permutation
#C2
# gene_sets_EXAMPLE <- readRDS("gene_sets.rds")
gene_sets_read <- read.delim("c5.all.v2022.1.Hs.symbols.txt" , header = F)
gene_sets = vector(mode = "list", length = dim(gene_sets_read)[1])
names(gene_sets) = gene_sets_read[,1]
for (i in 1: dim(gene_sets_read)[1]) {
  temp=gene_sets_read[i,]
  temp =temp [2 : length(temp)]
  temp= unlist(temp)
  temp = temp[temp!=""]
  names(temp) =NULL
  gene_sets[[i]]=temp
}
#CONTINIOUS
#ES_df <-gsLRT2(Z[395:400], X[395:400,], G[395:400, ], interaction_term = "apoe", gene_sets = gene_sets, n_perm = 2)
ES_df_C2 <- gsLRT2(Z, X, G, interaction_term = "apoe", gene_sets = gene_sets, n_perm = 1000)
write.xlsx(ES_df_C2, 'braak_c5_intrc_apoe_no_PMI.xlsx')

#BINARY
#ES_df <-gsLRT3(Z[1:10], X[1:10,], G[1:10, ], interaction_term = "apoe", gene_sets = gene_sets, n_perm = 2)
#ES_df_C2 <- gsLRT3(Z, X, G, interaction_term = "apoe", gene_sets = gene_sets, n_perm = 1000)
#  write.xlsx(ES_df_C2, 'amyloid_C2_intrc_apoe_no_PMI.xlsx')

#Halmark
gene_sets_read <- read.delim("Gene sets Msigdb Hallmark.txt" , header = F, sep = ",")
gene_sets = vector(mode = "list", length = dim(gene_sets_read)[1])
names(gene_sets) = gene_sets_read[,1]
for (i in 1: dim(gene_sets_read)[1]) {
  temp=gene_sets_read[i,]
  temp =temp [2 : length(temp)]
  temp= unlist(temp)
  temp = temp[temp!=""]
  names(temp) =NULL
  gene_sets[[i]]=temp
}

#CONTINIOUS 
ES_df_hall <- gsLRT2(Z, X, G, interaction_term = "apoe", gene_sets = gene_sets, n_perm = 1000)
write.xlsx(ES_df_hall, 'log_tangles_Hallmark_apoe_no_PMI.xlsx')

#BINARY
#ES_df_hall <- gsLRT3(Z, X, G, interaction_term = "apoe", gene_sets = gene_sets, n_perm = 1000)
#write.xlsx(ES_df_hall, 'amyloid_C2_Hallmark_apoe_no_PMI.xlsx')











# 
# # 
# #debugging
# ES_df <- gsLRT(Z[1:10], X[1:10,], G[1:10, ], interaction_term = "sex", gene_sets = gene_sets, n_perm = 3)
# 
# j=0
# for (i in 1:length(gene_sets)) {
#   temp=gene_sets[[i]]
#   prod=prod( temp %in% colnames(G))
#   if( prod==1) j=j+1
#   { cat(i, " ")}
# }
# j
# 
# j=0
# for (i in 1:length(gene_sets)) {
#   temp=gene_sets[[i]]
#   prod=prod( temp %in% colnames(G))
#   if( prod!=1) j=j+1
#   { cat(i, " ")}
# }
# j

