library(RMINC)
library(lme4)
library(nlme)
library(mni.cortical.statistics)
library(ggplot2)
library(stats)

## Data organization
data <- read.csv('/data/chamal/projects/atalpala/stats/Organized_data_thickness_age_with_FLASH_fBIRN_2_clust.csv')
data$X2_cluster_sol=as.factor(data$X2_cluster_sol)
data$X3_cluster_sol=as.factor(data$X3_cluster_sol)
data$ct_regions=data[c(8:87)]
data$ct_regions_l=data[c(8:46)]
data$ct_regions_r=data[c(48:86)]

data$ct_Relregions_l=data$ct_regions_l/data$Mean.L_CT
data$ct_Relregions_r=data$ct_regions_r/data$Mean.R_CT

# data$left_thickness <- paste("/data/chamal/projects/atalpala/processed_data/from_Schizconnect_T1_NC_SZ_SAPS_SANS_XC/CIVET/NUSDAST/averages/output/",data$Folder_name, "/thickness/NUSDAST_averages_",data$Folder_name,"_native_rms_rsl_tlink_28.28mm_left.txt", sep="")
# write.csv(data$left_thickness, "/data/chamal/projects/atalpala/stats/NUSDAST_thickness_locations.csv")
#data$left_thickness <- paste("/data/chamal/projects/atalpala/processed_data/from_Schizconnect_T1_NC_SZ_SAPS_SANS_XC/CIVET/NUSDAST/FLASH/output/",data$Folder_name, "/thickness/FLASH_",data$Folder_name,"_native_rms_rsl_tlink_28.28mm_left.txt", sep="")
#write.csv(data$left_thickness, "/data/chamal/projects/atalpala/stats/NUSDAST_FLASH_thickness_locations.csv")
# data$left_thickness <- paste("/data/chamal/projects/atalpala/processed_data/from_Schizconnect_T1_NC_SZ_SAPS_SANS_XC/CIVET/NMorph/3T/output/",data$Folder_name, "/thickness/NMorph_XC_",data$Folder_name,"_native_rms_rsl_tlink_28.28mm_left.txt", sep="")
# write.csv(data$left_thickness, "/data/chamal/projects/atalpala/stats/NMorph_thickness_locations.csv")
# data$left_thickness <- paste("/data/chamal/projects/atalpala/processed_data/from_Schizconnect_T1_NC_SZ_SAPS_SANS_XC/CIVET/fBIRN/T1_5/fixed/output/",data$Folder_name, "/thickness/fBIRN_T1_5_",data$Folder_name,"_native_rms_rsl_tlink_28.28mm_left.txt", sep="")
# write.csv(data$left_thickness, "/data/chamal/projects/atalpala/stats/fBIRN_15_thickness_locations.csv")
# data$left_thickness <- paste("/data/chamal/projects/atalpala/processed_data/from_Schizconnect_T1_NC_SZ_SAPS_SANS_XC/CIVET/fBIRN/T3/fixed/output/",data$Folder_name, "/thickness/fBIRN_T3_",data$Folder_name,"_native_rms_rsl_tlink_28.28mm_left.txt", sep="")
# write.csv(data$left_thickness, "/data/chamal/projects/atalpala/stats/fBIRN_3_thickness_locations.csv")
# data$left_thickness <- paste("/data/chamal/projects/atalpala/processed_data/from_Schizconnect_T1_NC_SZ_SAPS_SANS_XC/CIVET/fBIRN/T4/fixed/output/",data$Folder_name, "/thickness/fBIRN_T4_",data$Folder_name,"_native_rms_rsl_tlink_28.28mm_left.txt", sep="")
# write.csv(data$left_thickness, "/data/chamal/projects/atalpala/stats/fBIRN_4_thickness_locations.csv")

## AAL atlas analysis
# Calculates statistics and coefficients for linear model of specified anat structure
ana<-anatLm(~ X2_cluster_sol+age+sex, data, data$ct_Relregions_l)
anaFDR<-anatFDR(ana)
colnames(ana)
anaFDR
#write.csv(anaFDR, "/data/chamal/projects/atalpala/stats/AAL-analysis_CT_age.csv")
#write.csv(ana, "/data/chamal/projects/atalpala/stats/AAL-analysis_CT-ana_age.csv")

ana<-anatLm(~ X3_cluster_sol+age+sex, data, data$ct_regions)
anaFDR<-anatFDR(ana)
colnames(ana)
anaFDR

ana<-anatLm(~ X2_cluster_sol+age+sex+dataset, data, data$ct_regions)
anaFDR<-anatFDR(ana)
anaFDR
#write.csv(anaFDR, "/data/chamal/projects/atalpala/stats/AAL-analysis_CT_age.csv")
#write.csv(ana, "/data/chamal/projects/atalpala/stats/AAL-analysis_CT-ana_age.csv")

# interface to running linear mixed effects models at every vertex. Unlike standard linear models testing hypotheses in 
#linear mixed effects models is more difficult, since the denominator degrees of freedom are more difficult to determine. 
#RMINC provides estimating degrees of freedom using the anatLmerEstimateDF function.

# 2 and 3 cluster
ana<-anatLmer(~ X2_cluster_sol+age+sex+ (1|dataset), data, data$ct_regions)
ana
#anaFDR<-anatFDR(ana)
# anaFDR
#write.csv(anaFDR, "/data/chamal/projects/atalpala/stats/AAL-analysis_SA_age.csv")
#write.csv(ana, "/data/chamal/projects/atalpala/stats/AAL-analysis_SA-ana_age.csv")
ps_ana <- 2 * (1 - pnorm(abs(ana[,6])))
range(ps_ana)

data$X2_cluster_names<-relevel(data$X3_cluster_sol, ref="2")
ana3<-anatLmer(~ X2_cluster_names+age+sex+ (1|dataset), data, data$ct_regions)
ana3
colnames(ana3)
range(ana3[,7])
range(ana3[,8])
pvals_1v2<-pt2(ana3[,7],247)
range(pvals_1v2)
pvals_1v3<-pt2(ana3[,8],247)
range(pvals_1v3)
pvals_1v2[pvals_1v2<0.05]
pvals_1v3[pvals_1v3<0.05]

ps_ana3 <- 2 * (1 - pnorm(abs(ana3[,8])))
range(ps_ana3)
write.csv(ana3, "/data/chamal/projects/atalpala/stats/AAL-anatLmer_age_sex_1dataset.csv")

ana3<-anatLmer(~ X3_cluster_sol+age+sex+ (1|dataset), data, data$ct_Relregions_l)
ana3
colnames(ana3)
range(ana3[,7])
range(ana3[,8])

# # group 1 vs group 2
# data <- read.csv('/data/chamal/projects/atalpala/stats/one_v_2.csv')
# data$one_v_two=as.factor(data$one_v_two)
# data$ct_regions=data[c(7:86)]
# data$ct_regions_l=data[c(7:45)]
# data$ct_regions_r=data[c(47:85)]
# 
# data$ct_Relregions_l=data$ct_regions_l/data$Mean.L_CT
# data$ct_Relregions_r=data$ct_regions_r/data$Mean.R_CT

# c_one_v_two_l<-anatLmer(~ one_v_two + age + sex + (1|dataset), data, data$ct_Relregions_l)
# ps_c_one_v_two_l <- 2 * (1 - pnorm(abs(c_one_v_two_l[,7])))
# range(ps_c_one_v_two_l)
# c_one_v_two_r<-anatLmer(~ one_v_two+age+sex+ (1|dataset), data, data$ct_Relregions_r)
# ps_c_one_v_two_r <- 2 * (1 - pnorm(abs(c_one_v_two_r[,6])))
# range(ps_c_one_v_two_r)
# c_one_v_two<-anatLmer(~ one_v_two+age+sex+ (1|dataset), data, data$ct_regions)
# ps_c_one_v_two <- 2 * (1 - pnorm(abs(c_one_v_two[,6])))
# range(ps_c_one_v_two)


# # group 1 vs group 3
# data <- read.csv('/data/chamal/projects/atalpala/stats/one_v_3.csv')
# data$one_v_three=as.factor(data$one_v_three)
# data$ct_regions=data[c(7:86)]
# data$ct_regions_l=data[c(7:45)]
# data$ct_regions_r=data[c(47:85)]
# 
# data$ct_Relregions_l=data$ct_regions_l/data$Mean.L_CT
# data$ct_Relregions_r=data$ct_regions_r/data$Mean.R_CT
# 
# c_one_v_three_l<-anatLmer(~ one_v_three+age+sex+ (1|dataset), data, data$ct_Relregions_l)
# ps_c_one_v_three_l <- 2 * (1 - pnorm(abs(c_one_v_three_l[,6])))
# range(ps_c_one_v_three_l)
# c_one_v_three_r<-anatLmer(~ one_v_three+age+sex+ (1|dataset), data, data$ct_Relregions_r)
# ps_c_one_v_three_r <- 2 * (1 - pnorm(abs(c_one_v_three_r[,6])))
# range(ps_c_one_v_three_r)
# c_one_v_three<-anatLmer(~ one_v_three+age+sex+ (1|dataset), data, data$ct_regions)
# ps_c_one_v_three <- 2 * (1 - pnorm(abs(c_one_v_three[,6])))
# range(ps_c_one_v_three)
# 
# # group 2 vs group 3
# data <- read.csv('/data/chamal/projects/atalpala/stats/two_v_3.csv')
# data$two_v_three=as.factor(data$two_v_three)
# data$ct_regions=data[c(7:86)]
# data$ct_regions_l=data[c(7:45)]
# data$ct_regions_r=data[c(47:85)]
# 
# data$ct_Relregions_l=data$ct_regions_l/data$Mean.L_CT
# data$ct_Relregions_r=data$ct_regions_r/data$Mean.R_CT
# 
# c_two_v_three_l<-anatLmer(~ two_v_three+age+sex+ (1|dataset), data, data$ct_Relregions_l)
# ps_c_two_v_three_l <- 2 * (1 - pnorm(abs(c_two_v_three_l[,6])))
# range(ps_c_two_v_three_l)
# c_two_v_three_r<-anatLmer(~ two_v_three+age+sex+ (1|dataset), data, data$ct_Relregions_r)
# ps_c_two_v_three_r <- 2 * (1 - pnorm(abs(c_two_v_three_r[,6])))
# range(ps_c_two_v_three_r)
# c_two_v_three<-anatLmer(~ two_v_three+age+sex+ (1|dataset), data, data$ct_regions)
# ps_c_two_v_three <- 2 * (1 - pnorm(abs(c_two_v_three[,6])))
# range(ps_c_two_v_three)


# other analyses
ana<-anatLmer(~ X2_cluster_sol+age+sex+(1|dataset), data, data$ct_Relregions_l)
ana[,6]
ana<-anatLmer(~ X2_cluster_sol+age+sex+(1|dataset), data, data$ct_Relregions_r)
ana[,6]
ana<-anatLmer(~ X2_cluster_sol+age+sex+(1|dataset), data, data$Mean.L_CT)
ana[,6]
ana<-anatLmer(~ X2_cluster_sol+age+sex+(1|dataset), data, data$Mean.R_CT)
ana[,6]


summary(lm(data$Mean.L_CT~data$X2_cluster_sol))
summary(lm(data$Mean.R_CT~data$X2_cluster_sol))

## Vertex LM
data <- read.csv('/data/chamal/projects/atalpala/stats/Organized_data_thickness_age_with_FLASH_fBIRN_2_clust.csv')
data$X2_cluster_sol=as.factor(data$X2_cluster_sol)
data$X3_cluster_sol=as.factor(data$X3_cluster_sol)

names(data)
data$X2_cluster_sol

data$X2_cluster_names<-relevel(data$X2_cluster_sol, ref="1")
vs <- vertexLm(left_thickness ~ X2_cluster_names + sex + age, data)
range(vs[,10])
vertexFDR(vs)

data$X3_cluster_names<-relevel(data$X3_cluster_sol, ref="1")
vs <- vertexLm(left_thickness ~ X3_cluster_names + sex + age, data)
colnames(vs)
range(vs[,9])
range(vs[,10])
vertexFDR(vs)
write.csv(vs, "/data/chamal/projects/atalpala/stats/vertexLm_sex_age.csv")
write.table(x=vs, col.names=FALSE, row.names=FALSE, file="vertexLm_sex_age.txt")

# vs <- vertexLm(left_thickness ~ X2_cluster_sol + Mean.L_CT + dataset + sex + age, data)
# vertexFDR(vs)
# vs <- vertexLm(left_thickness ~ X2_cluster_sol + Mean.L_CT + sex + age, data)
# range(vs[,10])
# vertexFDR(vs)

#write.table(vs, row.names = FALSE, col.names = FALSE, "/data/chamal/projects/atalpala/stats/groupxsex-age.txt")

## vertex LMER
data <- read.csv('/data/chamal/projects/atalpala/stats/Organized_data_thickness_age_with_FLASH_fBIRN_2_clust.csv')
data$left_thickness=as.character(data$left_thickness)
data$X3_cluster_sol = as.factor(data$X3_cluster_sol)

data$clust_1<-relevel(data$X3_cluster_sol, ref="2")
l_thick = vertexTable(data$left_thickness)

#ana<-vertexLmer(l_thick ~ sex+age+(1|dataset), data)

#ana<-vertexLmer(data,'y~ clust_1 + age + sex + (1|dataset)',left_thickness)

data_left_thickness = mni.build.data.table(data, "left_thickness")
vs <-mni.vertex.mixed.model(data,'y ~ clust_1 + age + sex' , '~ 1 | dataset', data_left_thickness)


pvals_2v1<-pt2(vs$t.value[,"clust_11"],247)
range(pvals_2v1)
write.table(x=pvals_2v1, col.names=FALSE, row.names=FALSE, file="pvals_2v1.txt")
pvals_2v3<-pt2(vs$t.value[,"clust_13"],247)
write.table(x=pvals_2v3, col.names=FALSE, row.names=FALSE, file="pvals_2v3.txt")
range(pvals_2v3)
qvals<- p.adjust(pvals, "fdr")

ana<-vertexLmer('left_thickness ~ X3_cluster_sol+age+sex+(1|dataset)', data)#, l_thick)
# things after + are covariates (sex) and | expressions are random effects to be taken out
ana
ps_ana <- 2 * (1 - pnorm(abs(ana[,7])))
range(ps_ana)

# # group 1 vs group 3
# data <- read.csv('/data/chamal/projects/atalpala/stats/one_v_3.csv')
# data$one_v_three=as.factor(data$one_v_three)
# vs <- vertexLm(left_thickness ~ one_v_three + age + sex + dataset, data)
# range(vs[,10])
# vertexFDR(vs)

data$DX_bl<-relevel(data$one_v_three, ref="1")
l_thick = vertexTable(data$left_thickness)
ana<-vertexLmer(filtered,'y ~ X3_cluster_sol+sex+age+(1|dataset)', l_thick)

L.thal.sa.long = vertexLmer(filtered,'y~ DX_bl*Month_bl + APOE4_STATUS + PTGENDER + Current_age + (1|PTID)',lThal_sa)
lThal_sa = vertexTable(filtered$leftthalamus_sa)

## Plots
qplot(age,SFGdor.L_CT , color=X2_cluster_sol, facets=.~sex, data=subset(data), geom=c("smooth","point"))

qplot(age,ANG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ORBsup.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,PreCG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,SFGdor.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,MFG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ORBmid.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ORBinf.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ROL.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,SMA.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,OLF.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,SFGmed.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ORBsupmed.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,REC.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,INS.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ACG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,DCG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,PCG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,CUN.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,LING.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,SOG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,IOG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,FFG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,PoCG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,SPG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,IPL.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,SMG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ANG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,PCUN.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,PCL.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,HES.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,STG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,TPOsup.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,MTG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,TPOmid.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,ITG.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))
qplot(age,Mean.L_CT , color=X2_cluster_sol, data=subset(data), geom=c("smooth","point"))


qplot(age,SMG.L_CT , color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))
qplot(age,SFGdor.R_CT , color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))
qplot(age,SMG.L_CT/Mean.L_CT , color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))
qplot(age,STG.L_CT/Mean.L_CT , color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))
qplot(age,MTG.L_CT/Mean.L_CT , color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))

qplot(age,MTG.L_CT/Mean.L_CT , color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))
hist(data$age)

qplot(age^2, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm)
qplot(age^3, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm)
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 2))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 3))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 4))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 5))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 1))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 2))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST'), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 3))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST' $ data$age<20), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 1))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$dataset=='NUSDAST' & data$age<20), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 1))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$age<20), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 1))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$age<40), geom=c("smooth","point"))+geom_smooth(method=lm, formula =  y ~ poly(x, 1))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$age<20), geom=c("smooth","point"))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$age<40), geom=c("smooth","point"))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$age>40), geom=c("smooth","point"))
qplot(age, Mean.L_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$age<30), geom=c("smooth","point"))
qplot(age, Mean.R_CT, color=X2_cluster_sol, facets=.~sex, data=subset(data, data$age<30), geom=c("smooth","point"))
