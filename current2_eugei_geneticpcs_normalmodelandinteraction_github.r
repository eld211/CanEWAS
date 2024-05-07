#cannabis_current2
 cd /mnt/data1/Emma/Cannabis/EUGI/

setwd("/mnt/data1/Emma/Cannabis/EUGI/")
load("eugi_phenobeta_newcan_snpshybremoved.RData")
load("pheno_020420.RData")
rownames(pheno)<-pheno[,1]
pheno2<-subset(pheno, pheno$cannabis_current2>=0)

betas2<-betas[,match((pheno2$Basename),colnames(betas))]



##read in genetic pcs
eug<-pheno2

pcs.eugei<-read.table("/mnt/data1/EuGEI/Genotypes/MergedWithHapMap/Attempt2_LDprune_PCA.eigenvec",stringsAsFactors = FALSE)
map<-read.csv("/mnt/data1/EuGEI/Genotypes/EUGEI_SNP_genotypes_matchedIDs.csv",stringsAsFactors = FALSE)
map<-map[match(pcs.eugei$V1, map$Geno.Plate.ID),]
eug<-cbind(eug, pcs.eugei[match(eug$Geno.CHIP.Location, map$Geno.Plate.ID), c("V3", "V4")])

colnames(eug)[(ncol(eug)-1):ncol(eug)]<-c("PC1", "PC2")


mat<-matrix(data = NA, ncol =5, nrow = nrow(betas2))
rownames(mat)<-rownames(betas2)
colnames(mat)<-c("Eugi_Beta", "Eugi_SE", "Eugi_P","Eugi_Beta_smoke","Eugi_smokeP")
for(i in 1:nrow(betas2)){
	### standard linear regression model
	

	model<-lm(betas2[i,] ~ as.factor(eug$cannabis_current2)+as.factor(eug$case_control)+eug$SmokingScore+as.factor(eug$Cohort)+eug$Age.x+eug$Bcell+eug$CD8T+eug$CD4T+eug$NK+eug$Mono+eug$Gran+as.factor(eug$Sex)+as.factor(eug$MethPlate)+ eug$PC1 +eug$PC2)
	mat[i,1:3]<-summary(model)$coefficients["as.factor(eug$cannabis_current2)1", c(1,2,4)]
	mat[i,4:5]<-summary(model)$coefficients["eug$SmokingScore", c(1,4)]
}


epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7)
epicManifest[,c("Name","CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]->epi
epi<-epi[match(rownames(mat), epi$Name),]

write.csv(mat,"can_current2_eug_geneticpc_dec23.csv")




##interaction model



mat<-matrix(data = NA, ncol =9, nrow = nrow(betas2))
rownames(mat)<-rownames(betas2)
colnames(mat)<-c("Eugi_Beta", "Eugi_SE", "Eugi_P","Eugi_Beta_smoke","Eugi_smokeP")
for(i in 1:nrow(betas2)){
	### standard linear regression model
	

	model<-lm(betas2[i,] ~ as.factor(eug$cannabis_current2)+as.factor(eug$cannabis_current2)* as.factor(eug$case_control) + as.factor(eug$case_control)+eug$SmokingScore+as.factor(eug$Cohort)+eug$Age.x+eug$Bcell+eug$CD8T+eug$CD4T+eug$NK+eug$Mono+eug$Gran+as.factor(eug$Sex)+as.factor(eug$MethPlate)+ eug$PC1 +eug$PC2)
	mat[i,1:3]<-summary(model)$coefficients["as.factor(eug$cannabis_current2)1", c(1,2,4)]
	mat[i,4:6]<-summary(model)$coefficients["eug$SmokingScore", c(1,2,4)]
	mat[i,7:9]<-summary(model)$coefficients["as.factor(eug$cannabis_current2)1:as.factor(eug$case_control)1", c(1,2,4)]
}


epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7)
epicManifest[,c("Name","CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]->epi
epi<-epi[match(rownames(mat), epi$Name),]

write.csv(mat,"can_current2_eug_geneticpc_interaction_dec23.csv")







