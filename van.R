options(java.parameters = "-Xmx2048m")  


library(xlsx)


AB_path = '/Users/ali/Desktop/april23/GSLRT/van/amyloid_mf_intrc_apoe_no_PMI.xlsx'
AB = as.data.frame(read.xlsx2(AB_path, sheetIndex = 1))
LT_path = '/Users/ali/Desktop/april23/GSLRT/van/log_tangles_mf_intrc_apoe_no_PMI.xlsx'
LT = as.data.frame(read.xlsx2(LT_path, sheetIndex = 1))
ME_path = '/Users/ali/Desktop/april23/GSLRT/van/mmse_mf_intrc_apoe_no_PMI.xlsx'
ME = as.data.frame(read.xlsx2(ME_path, sheetIndex = 1))

# filtered :
AB = AB [ as.numeric(AB$p)<=0.1, ] 
LT = LT [ as.numeric(LT$p)<=0.1, ] 
ME = ME [ as.numeric(ME$p)<=0.1, ] 

intersect(AB$GS, LT$GS) 
intersect(LT$GS, ME$GS) 
intersect(ME$GS, AB$GS) 

intersect(LT$GS , intersect(ME$GS, AB$GS) )
