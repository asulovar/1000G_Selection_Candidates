#Author: Arvis Sulovari, PhD
#GWAS-based enrichment for highly-differentiated SNVs

setwd("C:/Users/arvis/[...]")
#Read in all files post QC

all <- read.delim("simple_BEST_LATEST_FINAL_afterQC_GWAScatalogue_Nov13_2015.txt",header=F)

afr_afr <- read.delim("afr_afr_FINAL_matched.gwascat",header=F)
sas_sas <- read.delim("sas_sas_FINAL_matched.gwascat",header=F)
eur_eur <- read.delim("eur_eur_FINAL_matched.gwascat",header=F)
eas_eas <- read.delim("eas_eas_FINAL_matched.gwascat",header=F)
amr_amr <- read.delim("amr_amr_FINAL_matched.gwascat",header=F)


#Start enrichment analysis using hypergeometric model


gwas_snps_nr <- length(unique(all[,2]))

#AFRICANS (AFR)

l=length(unique(afr_afr$V1))
afr_snps_nr <- length(unique(afr_afr$V2))
afr_profile <- array(NA, dim=c(l,7))


for (i in 1:l) {
  #m - high Fst SNPs for disease i
  m <- nrow(afr_afr[afr_afr$V1==toString(unique((afr_afr$V1))[i]),])
  
  #n - number of GWAS catalogue SNPs for disease i
  n <- nrow(all[all$V1==toString(unique((afr_afr$V1))[i]),])
  
  obs <- m
  exp <- afr_snps_nr*(n/gwas_snps_nr)
  
  afr_profile[i,1] <- toString(unique((afr_afr$V1))[i])
  afr_profile[i,2] <- n
  afr_profile[i,3] <- m
  afr_profile[i,4] <- gwas_snps_nr
  afr_profile[i,5] <- (obs/exp)
  

  afr_profile[i,6] <- (phyper(m-1,
                              n,
                              gwas_snps_nr-n,
                              afr_snps_nr,
                              lower.tail = F,log.p = F))
  
  
}

# 
# afr_profile[,6] <- (1-phyper(as.numeric(afr_profile[,3]),
#                              as.numeric(afr_profile[,2]),
#                              gwas_snps_nr,
#                              sum(as.numeric(afr_profile[,3]))))
# 
# 

#######CORRECTED PHYPER CODE############

afr_profile[,7] <- p.adjust(as.numeric(afr_profile[,6]),method = "BH")

colnames(afr_profile) <- c("Disease_Trait","GWAS_SNPs","AFR_SNPs","Total_GWAS_SNPs","R","AFR_UnAdj_P","AFR_adj_P")

write.table(afr_profile,"AFR_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",sep="\t")

######################################

#Set variables to zero before moving to the next POP
exp <- 0
obs <- 0
m <- 0
n <- 0


# 
# #afr_profile[,8] <- p.adjust(afr_profile[,6],method = "fdr")
# 
# afr_profile_ordered <- afr_profile[order(as.numeric(as.character(afr_profile[,6])), as.numeric(as.character(afr_profile[,5])),decreasing = F),]
# colnames(afr_profile_ordered) <- c("Disease/Trait","GWAS_SNPs","AFR_SNPs","Total GWAS_SNPs","R","AFR_unadj_P","AFR_adj_P")
# 
# #View(afr_profile[order(as.numeric(as.character(afr_profile[,6])), as.numeric(as.character(afr_profile[,5])),decreasing = F),])
# 
# write.table(afr_profile_ordered,"AFR_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",sep="\t")



#EUROPEANS (EUR)

l=length(unique(eur_eur$V1))
eur_snps_nr <- length(unique(eur_eur$V2))
eur_profile <- array(NA, dim=c(l,7))


for (i in 1:l) {
  #m - high Fst SNPs for disease i
  m <- nrow(eur_eur[eur_eur$V1==toString(unique((eur_eur$V1))[i]),])
  
  #n - number of GWAS catalogue SNPs for disease i
  n <- nrow(all[all$V1==toString(unique((eur_eur$V1))[i]),])
  
  obs <- m
  exp <- eur_snps_nr*(n/gwas_snps_nr)
  
  #Populate output array 
  eur_profile[i,1] <- toString(unique((eur_eur$V1))[i])
  eur_profile[i,2] <- n
  eur_profile[i,3] <- m
  eur_profile[i,4] <- gwas_snps_nr
  eur_profile[i,5] <- (obs/exp)
  
  eur_profile[i,6] <- (phyper(m-1,
                              n,
                              gwas_snps_nr-n,
                              eur_snps_nr,
                              lower.tail = F,log.p = F))
  
  
  
  
}


#######CORRECTED PHYPER CODE############
# 
# eur_profile[,6] <- -1*(log10(exp(phyper(as.numeric(eur_profile[,3])-1,
#                                         as.numeric(eur_profile[,2]),
#                                         gwas_snps_nr-as.numeric(eur_profile[,2]),
#                                         eur_total_hits,lower.tail = F,log.p = T))))


eur_profile[,7] <- p.adjust(as.numeric(eur_profile[,6]),method = "BH")

colnames(eur_profile) <- c("Disease_Trait","GWAS_SNPs","EUR_SNPs","Total_GWAS_SNPs","R","EUR_UnAdj_P","EUR_adj_P")

write.table(eur_profile,"EUR_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",sep="\t")

######################################


#Set variables to zero before moving to the next POP
exp <- 0
obs <- 0
m <- 0
n <- 0
# 
# eur_profile[,6] <- (1-phyper(as.numeric(eur_profile[,3]),as.numeric(eur_profile[,2]),gwas_snps_nr,sum(as.numeric(eur_profile[,3]))))
# eur_profile[,7] <- p.adjust(as.numeric(eur_profile[,6]),method = "BH")
# 
# #eur_profile[,8] <- p.adjust(eur_profile[,6],method = "fdr")
# 
# eur_profile_ordered <- eur_profile[order(as.numeric(as.character(eur_profile[,6])), as.numeric(as.character(eur_profile[,5])),decreasing = F),]
# colnames(eur_profile_ordered) <- c("Disease/Trait","GWAS_SNPs","EUR_SNPs","Total GWAS_SNPs","R","EUR_unadj_P","EUR_adj_P")
# 
# #View(eur_profile[order(as.numeric(as.character(eur_profile[,6])), as.numeric(as.character(eur_profile[,5])),decreasing = F),])
# 
# write.table(eur_profile_ordered,"EUR_DISEASE_PROFILE_HYPERGEO_FINAL.txt",sep="\t")




#EAST ASIANS (EAS)

l=length(unique(eas_eas$V1))
eas_snps_nr <- length(unique(eas_eas$V2))
eas_profile <- array(NA, dim=c(l,7))

for (i in 1:l) {
  #m - high Fst SNPs for disease i
  m <- nrow(eas_eas[eas_eas$V1==toString(unique((eas_eas$V1))[i]),])
  
  #n - number of GWAS catalogue SNPs for disease i
  n <- nrow(all[all$V1==toString(unique((eas_eas$V1))[i]),])
  
  obs <- m
  exp <- eas_snps_nr*(n/gwas_snps_nr)
  
  eas_profile[i,1] <- toString(unique((eas_eas$V1))[i])
  eas_profile[i,2] <- n
  eas_profile[i,3] <- m
  eas_profile[i,4] <- gwas_snps_nr
  eas_profile[i,5] <- (obs/exp)
  
  eas_profile[i,6] <- (phyper(m-1,
                              n,
                              gwas_snps_nr-n,
                              eas_snps_nr,
                              lower.tail = F,log.p = F))
  
  
}



#######CORRECTED PHYPER CODE############

# eas_profile[,6] <- -1*(log10(exp(phyper(as.numeric(eas_profile[,3])-1,
#                                         as.numeric(eas_profile[,2]),
#                                         gwas_snps_nr-as.numeric(eas_profile[,2]),
#                                         eas_total_hits,lower.tail = F,log.p = T))))


eas_profile[,7] <- p.adjust((eas_profile[,6]),method = "BH")

colnames(eas_profile) <- c("Disease_Trait","GWAS_SNPs","EAS_SNPs","Total_GWAS_SNPs","R","EAS_UnAdj_P","EAS_adj_P")

write.table(eas_profile,"EAS_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",sep="\t")

######################################


#Set variables to zero before moving to the next POP
exp <- 0
obs <- 0
m <- 0
n <- 0

# 
# eas_profile[,6] <- (1-phyper(as.numeric(eas_profile[,3]),as.numeric(eas_profile[,2]),gwas_snps_nr,sum(as.numeric(eas_profile[,3]))))
# eas_profile[,7] <- p.adjust(as.numeric(eas_profile[,6]),method = "BH")
# 
# #eas_profile[,8] <- p.adjust(eas_profile[,6],method = "fdr")
# 
# eas_profile_ordered <- eas_profile[order(as.numeric(as.character(eas_profile[,6])), as.numeric(as.character(eas_profile[,5])),decreasing = F),]
# colnames(eas_profile_ordered) <- c("Disease/Trait","GWAS_SNPs","eas_SNPs","Total GWAS_SNPs","R","EAS_unadj_P","EAS_adj_P")
# 
# #View(eas_profile[order(as.numeric(as.character(eas_profile[,6])), as.numeric(as.character(eas_profile[,5])),decreasing = F),])
# 
# write.table(eas_profile_ordered,"eas_DISEASE_PROFILE_HYPERGEO_FINAL.txt",sep="\t")
# 

#SOUTH ASIANS (SAS)
l=length(unique(sas_sas$V1))
sas_snps_nr <- length(unique(sas_sas$V2))
sas_profile <- array(NA, dim=c(l,7))

i=1
for (i in 1:l) {
  #m - high Fst SNPs for dissase i
  m <- nrow(sas_sas[sas_sas$V1==toString(unique((sas_sas$V1))[i]),])
  
  #n - number of GWAS catalogue SNPs for dissase i
  n <- nrow(all[all$V1==toString(unique((sas_sas$V1))[i]),])
  
  obs <- m
  exp <- sas_snps_nr*(n/gwas_snps_nr)
  
  
  sas_profile[i,1] <- toString(unique((sas_sas$V1))[i])
  sas_profile[i,2] <- n
  sas_profile[i,3] <- m
  sas_profile[i,4] <- gwas_snps_nr
  sas_profile[i,5] <- (obs/exp)
  sas_profile[i,6] <- (phyper(m-1,
                              n,
                              gwas_snps_nr-n,
                              sas_snps_nr,
                              lower.tail = F,log.p = F))
  
  
}



#######CORRECTED PHYPER CODE############


# sas_profile[,6] <- -1*(log10(exp(phyper(as.numeric(sas_profile[,3])-1,
#                                         as.numeric(sas_profile[,2]),
#                                         gwas_snps_nr-as.numeric(sas_profile[,2]),
#                                         sas_total_hits,lower.tail = F,log.p = T))))
# 

sas_profile[,7] <- p.adjust((sas_profile[,6]),method = "BH")

colnames(sas_profile) <- c("Disease_Trait","GWAS_SNPs","SAS_SNPs","Total_GWAS_SNPs","R","SAS_UnAdj_P","SAS_adj_P")

write.table(sas_profile,"SAS_DISSASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",sep="\t")


######################################


#Set variables to zero before moving to the next POP
exp <- 0
obs <- 0
m <- 0
n <- 0


# sas_profile[,6] <- (1-phyper(as.numeric(sas_profile[,3]),as.numeric(sas_profile[,2]),gwas_snps_nr,sum(as.numeric(sas_profile[,3]))))
# sas_profile[,7] <- p.adjust(as.numeric(sas_profile[,6]),method = "BH")
# 
# #sas_profile[,8] <- p.adjust(sas_profile[,6],method = "fdr")
# 
# sas_profile_ordered <- sas_profile[order(as.numeric(as.character(sas_profile[,6])), as.numeric(as.character(sas_profile[,5])),decreasing = F),]
# colnames(sas_profile_ordered) <- c("Disease/Trait","GWAS_SNPs","SAS_SNPs","Total GWAS_SNPs","R","SAS_unadj_P","SAS_adj_P")
# 
# #View(sas_profile[order(as.numeric(as.character(sas_profile[,6])), as.numeric(as.character(sas_profile[,5])),decreasing = F),])
# 
# write.table(sas_profile_ordered,"SAS_DISEASE_PROFILE_HYPERGEO_FINAL.txt",sep="\t")
# 

#AMERICANS (AMR)

l=length(unique(amr_amr$V1))
amr_snps_nr <- length(unique(amr_amr$V2))
amr_profile <- array(NA, dim=c(l,7))

i=1
for (i in 1:l) {
  #m - high Fst SNPs for disease i
  m <- nrow(amr_amr[amr_amr$V1==toString(unique((amr_amr$V1))[i]),])
  
  #n - number of GWAS catalogue SNPs for disease i
  n <- nrow(all[all$V1==toString(unique((amr_amr$V1))[i]),])
  
  obs <- m
  exp <- amr_snps_nr*(n/gwas_snps_nr)
  
  #Probability of picking m/n from the whole GWAS catalogue
  # p <- (choose(n,m) * choose(length(unique(all$V1))-n,nrow(amr_amr)-m))/choose(length(unique(all$V1)),nrow(amr_amr))
  
  #Binomial distribution for probability of hitting m SNPs from n SNPs of disease i (in 14,612 total SNPs).
  #p <- choose(n,m)*(n/nrow(all))^m*(1-n/nrow(all))^(n-m)
  #test <- matrix(c(m,n,nrow(amr_amr)-m,8641-n),ncol=2)
  
  #p <- choose(nrow(amr_amr),m)*(n/8641)^m*(1-n/8641)^(nrow(amr_amr)-m)
  
  amr_profile[i,1] <- toString(unique((amr_amr$V1))[i])
  amr_profile[i,2] <- n
  amr_profile[i,3] <- m
  amr_profile[i,4] <- gwas_snps_nr
  amr_profile[i,5] <- (obs/exp)
  amr_profile[i,6] <- (phyper(m-1,
                              n,
                              gwas_snps_nr-n,
                              amr_snps_nr,
                              lower.tail = F,log.p = F))
  
  
}


#######CORRECTED PHYPER CODE############

# amr_profile[,6] <- -1*(log10(exp(phyper(as.numeric(amr_profile[,3])-1,
#                                         as.numeric(amr_profile[,2]),
#                                         gwas_snps_nr-as.numeric(amr_profile[,2]),
#                                         amr_total_hits,lower.tail = F,log.p = T))))
# 

amr_profile[,7] <- p.adjust((amr_profile[,6]),method = "BH")

colnames(amr_profile) <- c("Disease_Trait","GWAS_SNPs","AMR_SNPs","Total_GWAS_SNPs","R","AMR_UnAdj_P","AMR_adj_P")

write.table(amr_profile,"AMR_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",sep="\t")

######################################


#Set variables to zero before moving to the next POP
exp <- 0
obs <- 0
m <- 0
n <- 0

# 
# 
# amr_profile[,6] <- (1-phyper(as.numeric(amr_profile[,3]),as.numeric(amr_profile[,2]),gwas_snps_nr,sum(as.numeric(amr_profile[,3]))))
# amr_profile[,7] <- p.adjust(as.numeric(amr_profile[,6]),method = "BH")
# 
# #amr_profile[,8] <- p.adjust(amr_profile[,6],method = "fdr")
# 
# amr_profile_ordered <- amr_profile[order(as.numeric(as.character(amr_profile[,6])), as.numeric(as.character(amr_profile[,5])),decreasing = F),]
# colnames(amr_profile_ordered) <- c("Disease/Trait","GWAS_SNPs","amr_SNPs","Total GWAS_SNPs","R","AMR_unadj_P","AMR_adj_P")
# 
# #View(amr_profile[order(as.numeric(as.character(amr_profile[,6])), as.numeric(as.character(amr_profile[,5])),decreasing = F),])
# 
# write.table(amr_profile_ordered,"amr_DISEASE_PROFILE_HYPERGEO_FINAL.txt",sep="\t")
# 


###MERGE RESULTS FILE

##GWASCAT disease/trait merging

# afr_gwas <- read.delim("AFR_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",header=T)
# eur_gwas <- read.delim("EUR_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",header=T)
# eas_gwas <- read.delim("EAS_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",header=T)
# sas_gwas <- read.delim("SAS_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",header=T)
# amr_gwas <- read.delim("AMR_DISEASE_PROFILE_HYPERGEO_FINAL_fixedBug.txt",header=T)

temp <- merge(afr_profile,eur_profile,by="Disease_Trait",all=T)
temp2 <- merge(temp,eas_profile,by="Disease_Trait",all=T)
temp3 <- merge(temp2,sas_profile,by="Disease_Trait",all=T)
temp4 <- merge(temp3,amr_profile,by="Disease_Trait",all=T)

write.csv(temp4, "Merged_gwasCatalogue_phenotypes_Dec1.csv")


############################COMPARISON TO FISHER&Chi square TEST############################


gwas_snps_nr <- length(unique(all[,2]))

#AFRICANS (AFR)

l=length(unique(afr_afr$V1))
afr_snps_nr <- length(unique(afr_afr$V2))
afr_profile <- array(NA, dim=c(l,9))


for (i in 1:l) {
  #m - high Fst SNPs for disease i
  m <- nrow(afr_afr[afr_afr$V1==toString(unique((afr_afr$V1))[i]),])
  
  #n - number of GWAS catalogue SNPs for disease i
  n <- nrow(all[all$V1==toString(unique((afr_afr$V1))[i]),])
  
  obs <- m
  exp <- afr_snps_nr*(n/gwas_snps_nr)
  
  afr_profile[i,1] <- toString(unique((afr_afr$V1))[i])
  afr_profile[i,2] <- n
  afr_profile[i,3] <- m
  afr_profile[i,4] <- gwas_snps_nr
  afr_profile[i,5] <- (obs/exp)
  
  
  afr_profile[i,6] <- (phyper(m-1,
                              n,
                              gwas_snps_nr-n,
                              afr_snps_nr,
                              lower.tail = F,log.p = F))
  
  #Fisher exact test, 2x2 table
  
  ft <- fisher.test(matrix(c(m,afr_snps_nr,n,gwas_snps_nr),nrow=2))
  afr_profile[i,8] <- ft$p.value
  
  #Chi square
  chi <- chisq.test(matrix(c(m,afr_snps_nr,n,gwas_snps_nr),nrow=2))
  
  afr_profile[i,9] <- chi$p.value
  
}


afr_profile[,7] <- p.adjust(as.numeric(afr_profile[,6]),method = "BH")




#