#Author: Arvis Sulovari, PhD
#Gene-based enrichment for highly-differentiated SNVs


setwd("C:/Users/arvis/OneDrive/Documents/WorkComp/Papers_v1/Under preparation/1000Genomes_Fst/F_st/Manuscript/Annotation/HighFst_genome_wide_annotation/Gene_based_enrichment")


#Read in the input file

gene_counts <- read.delim("HighFST_hitsPerGene.txt",header=T)

##HACKING THE GWAS BASED ENRICHMENT

#Start enrichment analysis using hypergeometric model

all_gene_variants <- 47240409

#AFRICANS (AFR)

l=29062
#l=length(na.omit(gene_counts$AFR_AFR))
#gene_list <- toString((gene_counts$Gene))
#afr_genelist <- (na.omit(cbind(toString(gene_counts$Gene),gene_counts$AFR_AFR)))
afr_GeneSNPs_nr <- sum(na.omit(gene_counts$AFR_AFR))
afr_geneEnrichment_profile <- array(NA, dim=c(l,7))

#i=7449
#1:l
for (i in 1:l) {
  #m - high Fst variants for gene i
  m <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$AFR_AFR)
  
  #n - number of 1000 Genomes  variants for gene i
  n <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$SNP_total)
  
  obs <- m
  exp <- afr_GeneSNPs_nr*(n/all_gene_variants)
  
  #Probability of picking m/n from the whole GWAS catalogue
  # p <- (choose(n,m) * choose(length(unique(all$V1))-n,nrow(afr_afr)-m))/choose(length(unique(all$V1)),nrow(afr_afr))
  
  #Binomial distribution for probability of hitting m SNPs from n SNPs of disease i (in 14,612 total SNPs).
  #p <- choose(n,m)*(n/nrow(all))^m*(1-n/nrow(all))^(n-m)
  #test <- matrix(c(m,n,nrow(afr_afr)-m,8641-n),ncex=2,col=2)
  
  #p <- choose(nrow(afr_afr),m)*(n/8641)^m*(1-n/8641)^(nrow(afr_afr)-m)
  
  afr_geneEnrichment_profile[i,1] <- toString(gene_counts[i,1])
  afr_geneEnrichment_profile[i,2] <- as.numeric(n)
  afr_geneEnrichment_profile[i,3] <- as.numeric(m)
  afr_geneEnrichment_profile[i,4] <- exp
  afr_geneEnrichment_profile[i,5] <- (obs/exp)
  afr_geneEnrichment_profile[i,6] <- -1*(log10(exp(phyper(m-1,
                                                          n,
                                                          all_gene_variants-n,
                                                          afr_GeneSNPs_nr,
                                                          lower.tail = F,log.p = T))))
}

#Reset these variables to 0
obs <- 0
exp <- 0
m <- 0
n <- 0

#afr_geneEnrichment_profile[,6] <- (1-phyper(as.numeric(afr_geneEnrichment_profile[,3]),as.numeric(afr_geneEnrichment_profile[,2]),all_gene_variants,sum(as.numeric(afr_geneEnrichment_profile[,3]))))
#afr_geneEnrichment_profile[,7] <- p.adjust(as.numeric(afr_geneEnrichment_profile[,6]),method = "BH")

#afr_profile[,8] <- p.adjust(afr_profile[,6],method = "fdr")

#afr_geneEnrichment_profile_ordered <- afr_geneEnrichment_profile[order(as.numeric(as.character(afr_geneEnrichment_profile[,6])), as.numeric(as.character(afr_geneEnrichment_profile[,5])),decreasing = F),]
colnames(afr_geneEnrichment_profile) <- c("GENE","1KG_Variants","AFR_Variants","Expected","R","-Log_10(AFR_unadj_P)","XXXAFR_adj_PXXX")

#View(afr_profile[order(as.numeric(as.character(afr_profile[,6])), as.numeric(as.character(afr_profile[,5])),decreasing = F),])

write.table(afr_geneEnrichment_profile,"AFR_GENE_ENRICHMENT_PROFILE_HYPERGEO_FINAL.txt",sep="\t")



#EUROPEANS
all_gene_variants <- 47240409


l=29062
eur_GeneSNPs_nr <- sum(na.omit(gene_counts$EUR_EUR))
eur_geneEnrichment_profile <- array(NA, dim=c(l,7))

for (i in 1:l) {
  m <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$EUR_EUR)
  
  n <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$SNP_total)
  
  obs <- m
  exp <- eur_GeneSNPs_nr*(n/all_gene_variants)
  
  eur_geneEnrichment_profile[i,1] <- toString(gene_counts[i,1])
  eur_geneEnrichment_profile[i,2] <- as.numeric(n)
  eur_geneEnrichment_profile[i,3] <- as.numeric(m)
  eur_geneEnrichment_profile[i,4] <- exp
  eur_geneEnrichment_profile[i,5] <- (obs/exp)
  eur_geneEnrichment_profile[i,6] <- -1*(log10(exp(phyper(m-1,
                                                          n,
                                                          all_gene_variants-n,
                                                          eur_GeneSNPs_nr,
                                                          lower.tail = F,log.p = T))))
  
}

#Reset these variables to 0
obs <- 0
exp <- 0
m <- 0
n <- 0

colnames(eur_geneEnrichment_profile) <- c("GENE","1KG_Variants","EUR_Variants","Expected","R","-Log_10(EUR_unadj_P)","XXXEUR_adj_PXXX")

write.table(eur_geneEnrichment_profile,"EUR_GENE_ENRICHMENT_PROFILE_HYPERGEO_FINAL.txt",sep="\t")



#EAST ASIANS
all_gene_variants <- 47240409


l=29062
eas_GeneSNPs_nr <- sum(na.omit(gene_counts$EAS_EAS))
eas_geneEnrichment_profile <- array(NA, dim=c(l,7))

for (i in 1:l) {
  m <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$EAS_EAS)
  
  n <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$SNP_total)
  
  obs <- m
  exp <- eas_GeneSNPs_nr*(n/all_gene_variants)
  
  eas_geneEnrichment_profile[i,1] <- toString(gene_counts[i,1])
  eas_geneEnrichment_profile[i,2] <- as.numeric(n)
  eas_geneEnrichment_profile[i,3] <- as.numeric(m)
  eas_geneEnrichment_profile[i,4] <- exp
  eas_geneEnrichment_profile[i,5] <- (obs/exp)
  eas_geneEnrichment_profile[i,6] <- -1*(log10(exp(phyper(m-1,n,all_gene_variants-n,eas_GeneSNPs_nr,lower.tail = F,log.p = T))))
  
}


#Reset these variables to 0
obs <- 0
exp <- 0
m <- 0
n <- 0


colnames(eas_geneEnrichment_profile) <- c("GENE","1KG_Variants","EAS_Variants","Expected","R","-Log_10(EAS_unadj_P)","XXXEAS_adj_PXXX")

write.table(eas_geneEnrichment_profile,"EAS_GENE_ENRICHMENT_PROFILE_HYPERGEO_FINAL.txt",sep="\t")



#SOUTH ASIANS
all_gene_variants <- 47240409


l=29062
sas_GeneSNPs_nr <- sum(na.omit(gene_counts$SAS_SAS))
sas_geneEnrichment_profile <- array(NA, dim=c(l,7))

for (i in 1:l) {
  m <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$SAS_SAS)
  
  n <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$SNP_total)
  
  obs <- m
  exp <- sas_GeneSNPs_nr*(n/all_gene_variants)
  
  sas_geneEnrichment_profile[i,1] <- toString(gene_counts[i,1])
  sas_geneEnrichment_profile[i,2] <- as.numeric(n)
  sas_geneEnrichment_profile[i,3] <- as.numeric(m)
  sas_geneEnrichment_profile[i,4] <- exp
  sas_geneEnrichment_profile[i,5] <- (obs/exp)
  sas_geneEnrichment_profile[i,6] <- -1*(log10(exp(phyper(m-1,n,all_gene_variants-n,sas_GeneSNPs_nr,lower.tail = F,log.p = T))))
  
}


colnames(sas_geneEnrichment_profile) <- c("GENE","1KG_Variants","SAS_Variants","Expected","R","-Log_10(SAS_unadj_P)","XXXSAS_adj_PXXX")

write.table(sas_geneEnrichment_profile,"SAS_GENE_ENRICHMENT_PROFILE_HYPERGEO_FINAL.txt",sep="\t")


#AMERICANS
all_gene_variants <- 47240409


l=29062
amr_GeneSNPs_nr <- sum(na.omit(gene_counts$AMR_AMR))
amr_geneEnrichment_profile <- array(NA, dim=c(l,7))

for (i in 1:l) {
  m <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$AMR_AMR)
  
  n <- as.numeric(gene_counts[gene_counts$Gene==toString(gene_counts[i,1]),]$SNP_total)
  
  obs <- m
  exp <- amr_GeneSNPs_nr*(n/all_gene_variants)
  
  amr_geneEnrichment_profile[i,1] <- toString(gene_counts[i,1])
  amr_geneEnrichment_profile[i,2] <- as.numeric(n)
  amr_geneEnrichment_profile[i,3] <- as.numeric(m)
  amr_geneEnrichment_profile[i,4] <- exp
  amr_geneEnrichment_profile[i,5] <- (obs/exp)
  amr_geneEnrichment_profile[i,6] <- -1*(log10(exp(phyper(m-1,n,all_gene_variants-n,amr_GeneSNPs_nr,lower.tail = F,log.p = T))))
  
}


colnames(amr_geneEnrichment_profile) <- c("GENE","1KG_Variants","AMR_Variants","Expected","R","-Log_10(AMR_unadj_P)","XXXAMR_adj_PXXX")

write.table(amr_geneEnrichment_profile,"AMR_GENE_ENRICHMENT_PROFILE_HYPERGEO_FINAL.txt",sep="\t")



#####PLOTTING AREA#############

#PLOT enrichment P-values against density

bonf_thresh <- -log(0.05/29062,10)
dens_thresh <- log(0.01,10)


par(mfrow=c(2,3))

#AFRICANS
plot((log(as.numeric(afr_geneEnrichment_profile[,3])/as.numeric(afr_geneEnrichment_profile[,2]),base=10)),
     afr_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="AFR",
     col=rgb(0, 0, 0,0.5),pch=15)

abline(h = bonf_thresh,v = dens_thresh,col="red")

text("HLA-DRB1",y=298,x=log10(260/1650))
text("HLA-DQB1",y=262.8,x=log10(242/1766))
text("SLC38A9",y=207.2,x=log10(249/3242))
text("SMC2",y=186.5,x=log10(195/1911))
text("HLA-T",y=175.5,x=log10(122/428))
text("HLA-DRB9",y=161.3,x=log10(149/1104))
text("HLA-DRA",y=141.5,x=log10(108/494))

#AC104389.28
text("1",y=130.6,x=log10(238/6196))

text("C6orf10",y=130.47,x=log10(171/2637))
text("PLD1",y=127.3,x=log10(233/6108))
text("HLA-H",y=125.2,x=log10(99/499))
text("DEPTOR",y=122.2,x=log10(218/5515))
text("ACER3",y=115.7071977,x=log10(0.043875686))
text("CD36",y=115.0948581,x=log10(0.033304746))
text("KIF20B",y=109.4082328,x=log10(0.062580369))
text("FRRS1",y=109.2074003,x=log10(0.075170843))

#RP11-624M8.1
text("2",y=103.87857,x=log10(0.037791848))

text("OR12D1",y=101.7743371,x=log10(0.278431373))
text("GOLGA8A",y=94.75016732,x=log10(0.049538777))


#EUROPEANS
plot((log(as.numeric(eur_geneEnrichment_profile[,3])/as.numeric(eur_geneEnrichment_profile[,2]),base=10)),
     eur_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="EUR",
     col=rgb(0, 0, 0,0.5),pch=16)

#abline(h = bonf_thresh,v = dens_thresh,col="red")

text("MICU1",y=259.5,x=log10(0.0386))
text("HLA-DPA1",y=254.1,x=log10(0.161))
text("HLA-DRB5",y=242.7,x=log10(0.11))
text("PLD1",y=236.9,x=log10(0.039))
text("HLA-DQA1",y=235.5,x=log10(0.069332355))
text("GSTCD",y=225.7883456,x=log10(0.054470426))

#RP11-94P11.4
text("3",y=205.133229,x=log10(0.025313913))

#HLA-DRB6
text("4",y=195.8376343,x=log10(0.123106061))

text("NCAM1",y=185.3597047,x=log10(0.02500546))
text("HLA-DRB1",y=174.9454439,x=log10(0.080606061))
text("MAP3K19",y=166.8831687,x=log10(0.059821807))
text("PBRM1",y=161.8879348,x=log10(0.045563549))
text("TMEM131",y=161.558489,x=log10(0.028752512))
text("MKL2",y=155.0322105,x=log10(0.032251908))
text("IRGM",y=153.0337845,x=log10(0.077124183))
text("TLR6",y=139.5583017,x=log10(0.080424886))
text("HERC2",y=120.2903388,x=log10(0.022918823))


#EAST ASIAN
plot((log(as.numeric(eas_geneEnrichment_profile[,3])/as.numeric(eas_geneEnrichment_profile[,2]),base=10)),
     eas_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="EAS",
     col=rgb(0, 0, 0,0.5),pch=17)


text("CNTLN",y=275.1980861,x=log10(0.024252891))
text("HLA-A",y=262.6173326,x=log10(0.139294927))
text("ABI1",y=242.5992067,x=log10(0.066037736))

#HCG4B
text("5",y=240.1904905,x=log10(0.199701937))

text("CLEC16A",y=209.6543317,x=log10(0.027615063))
text("LARP1B",y=207.3548072,x=log10(0.043609371))
text("OR5H5P",y=196.3319325,x=log10(0.177394035))
text("DOCK9",y=194.3661501,x=log10(0.025805698))
text("FRMD4B",y=175.9537992,x=log10(0.019106191))
text("MUC4",y=171.0447224,x=log10(0.046561825))
text("MSH4",y=167.9487578,x=log10(0.051217765))

#HLA-W
text("6",y=166.9676562,x=log10(0.138069705))

#snoU13
text("7",y=166.9099521,x=log10(0.117263844))

text("ANKRD36B",y=164.7744489,x=log10(0.066149871))
text("C7orf50",y=158.9973764,x=log10(0.028623731))
text("POC1B",y=155.6551548,x=log10(0.047752809))
text("STRN",y=140.8944602,x=log10(0.032956903))
text("WFS1",y=140.3146029,x=log10(0.056622354))
text("MEF2C",y=137.598827,x=log10(0.03330139))
text("NBAS",y=130.7644909,x=log10(0.016148224))
text("KRT6C",y=129.7747198,x=log10(0.131707317))
text("FADS2",y=124.4138731,x=log10(0.04))




#abline(h = bonf_thresh,v = dens_thresh,col="red")




#SOUTH ASIANS
plot((log(as.numeric(sas_geneEnrichment_profile[,3])/as.numeric(sas_geneEnrichment_profile[,2]),base=10)),
     sas_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="SAS",
     col=rgb(0, 0, 0,0.5),pch=8)

#HLA-T
text("8",y=264.71,x=log10(0.322429907))

text("MAST2",y=252.21,x=log10(0.042212042))
text("EGLN1",y=235.32,x=log10(0.085354478))
text("DDX11",y=230.20,x=log10(0.115056818))

#RP11-202K23.1
text("9",y=206.94,x=log10(0.035071668))

#HLA-H
text("10",y=196.85,x=log10(0.226452906))

#snoU13
text("11",y=181.51,x=log10(0.132464712))

#RP11-385F5.4
text("12",y=160.08,x=log10(0.162640902))

text("ENTPD1",y=144.85,x=log10(0.037329505))

#Y_RNA
text("13",y=143.49,x=log10(0.03418259))

#LINC00535
text("14",y=139.42,x=log10(0.018041447))

text("FCGR2B",y=125.20,x=log10(0.064966606))
text("CDK8",y=122.93,x=log10(0.03163017))

#HLA-L
text("15",y=120.68,x=log10(0.109137056))

text("MRPS6",y=118.94,x=log10(0.046586345))
text("LRFN5",y=118.76,x=log10(0.019839501))
text("FAM184B",y=115.45,x=log10(0.027512419))
text("MICU1",y=112.30,x=log10(0.022814248))
text("FCGR3B",y=112.06,x=log10(0.117117117))
text("CENPK",y=110.40,x=log10(0.063587684))
text("FOXJ3",y=109.14,x=log10(0.029177719))
text("GRM5",y=108.09,x=log10(0.012744413))


#abline(h = bonf_thresh,v = dens_thresh,col="red")




#AMERICANS
plot((log(as.numeric(amr_geneEnrichment_profile[,3])/as.numeric(amr_geneEnrichment_profile[,2]),base=10)),
     amr_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="AMR",
     col=rgb(0, 0, 0,0.5),pch=20)

text("TBC1D5",y=308.92,x=log10(0.013980365))

#DPY19L2P2
text("16",y=282.10,x=log10(0.069524131))

text("ATF6",y=235.53,x=log10(0.040917177))
text("HLA-DQA1",y=217.84,x=log10(0.060895084))
text("PLXNA4",y=209.70,x=log10(0.017043326))
text("DCC",y=192.43,x=log10(0.009519264))
text("TRIM40",y=184.58,x=log10(0.136138614))

#LINC00243
text("17",y=182.43,x=log10(0.08665269))

#AC073869.20
text("18",y=176.03,x=log10(0.196687371))

text("ROBO1",y=137.48,x=log10(0.008555192))

#RP11-398M15.1
text("19",y=137.42,x=log10(0.063384615))
text("LMTK2",y=117.75,x=log10(0.033810539))

#LINC00383
text("20",y=116.44,x=log10(0.041804636))

#RP11-142C4.6
text("21",y=113.94,x=log10(0.057704918))

text("RFX8",y=108.20,x=log10(0.037290715))

#HLA-DRB6
text("22",y=107.78,x=log10(0.072916667))

#RP5-991G20.1
text("23",y=106.78,x=log10(0.039030403))

text("DEPTOR",y=106.19,x=log10(0.021940163))

#RP11-98J9.1
text("24",y=106.06,x=log10(0.040522876))

text("FAM13A",y=97.49,x=log10(0.013322308))
text("PTPRD",y=96.22,x=log10(0.004254612))
text("IQGAP1",y=94.10,x=log10(0.028537098))


#abline(h = bonf_thresh,v = dens_thresh,col="red")











#######ALL IN ONE PANEL#################

par(mfrow=c(1,1))

#AFRICANS
plot((log(as.numeric(afr_geneEnrichment_profile[,3])/as.numeric(afr_geneEnrichment_profile[,2]),base=10)),
     afr_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="",
     cex=2,col=rgb(0, 0, 0,0.5),pch=15,xlim=c(-3,0),ylim=c(0,320))


text("HLA-DRB1",y=298.873515,x=log10(260/1650),cex=1.2,font=3)
points(y=298.87,x=log10(260/1650),cex=2,col=rgb(215/255,25/255,28/255,0.5),pch=15)

text("HLA-DQB1",y=262.8,x=log10(242/1766),cex=1.2,font=3)
points(y=262.8,x=log10(242/1766),cex=2,col=rgb(215/255,25/255,28/255,0.5),pch=15)

text("SLC38A9",y=207.2,x=log10(249/3242),cex=1.2,font=3)
points(y=207.2,x=log10(249/3242),cex=2,col=rgb(215/255,25/255,28/255,0.5),pch=15)

text("SMC2",y=186.5,x=log10(195/1911),cex=1.2, font=3)
points(y=186.5,x=log10(195/1911),cex=2,col=rgb(215/255,25/255,28/255,0.5),pch=15)

#HLA-T
text("HLA-T",y=175.5,x=log10(122/428),cex=1.2, font=3)
points(y=175.5,x=log10(122/428),cex=2,col=rgb(215/255,25/255,28/255,0.5),pch=15)

text("HLA-DRB9",y=161.3,x=log10(149/1104),cex=1.2, font=3)
points(y=161.3,x=log10(149/1104),cex=2,col=rgb(215/255,25/255,28/255,0.5),pch=15)

text("HLA-DRA",y=141.5,x=log10(108/494),cex=1.2, font=3)
points(y=141.5,x=log10(108/494),cex=2,col=rgb(215/255,25/255,28/255,0.5),pch=15)

#AC104389.28
#text("2",y=130.6,x=log10(238/6196),cex=1.2, font=3)

#text("C6orf10",y=130.47,x=log10(171/2637),cex=1.2, font=3)
#text("PLD1",y=127.3,x=log10(233/6108),cex=1.2, font=3)

#HLA-H
#text("3",y=125.2,x=log10(99/499),cex=1.2, font=3)

#text("DEPTOR",y=122.2,x=log10(218/5515),cex=1.2, font=3)
#text("ACER3",y=115.7071977,x=log10(0.043875686),cex=1.2, font=3)
#text("CD36",y=115.0948581,x=log10(0.033304746),cex=1.2, font=3)
#text("KIF20B",y=109.4082328,x=log10(0.062580369),cex=1.2, font=3)
#text("FRRS1",y=109.2074003,x=log10(0.075170843),cex=1.2, font=3)

#RP11-624M8.1
#text("4",y=103.87857,x=log10(0.037791848),cex=1.2, font=3)

#text("OR12D1",y=101.7743371,x=log10(0.278431373),cex=1.2, font=3)
#text("GOLGA8A",y=94.75016732,x=log10(0.049538777),cex=1.2, font=3)

par(new=T)


#EUROPEANS
plot((log(as.numeric(eur_geneEnrichment_profile[,3])/as.numeric(eur_geneEnrichment_profile[,2]),base=10)),
     eur_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="",
     cex=2,col=rgb(0, 0, 0,0.5),pch=16,xlim=c(-3,0),ylim=c(0,320))

#abline(h = bonf_thresh,v = dens_thresh,cex=2,col="red")

text("MICU1",y=259.5,x=log10(0.0386),cex=1.2, font=3)
points(y=259.5,x=log10(0.0386),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)


text("HLA-DPA1",y=254.1,x=log10(0.161),cex=1.2, font=3)
points(y=254.1,x=log10(0.161),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("HLA-DRB5",y=242.7,x=log10(0.11),cex=1.2, font=3)
points(y=242.7,x=log10(0.11),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("PLD1",y=236.9,x=log10(0.039),cex=1.2, font=3)
points(y=236.9,x=log10(0.039),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("HLA-DQA1",y=235.5,x=log10(0.069332355),cex=1.2, font=3)
points(y=235.5,x=log10(0.069332355),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("GSTCD",y=225.7883456,x=log10(0.054470426),cex=1.2, font=3)
points(y=225.7883456,x=log10(0.054470426),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

#RP11-94P11.4
text("RP11-94P11.4",y=205.133229,x=log10(0.025313913),cex=1.2, font=3)
points(y=205.133229,x=log10(0.025313913),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

#HLA-DRB6
text("HLA-DRB6",y=195.8376343,x=log10(0.123106061),cex=1.2, font=3)
points(y=195.8376343,x=log10(0.123106061),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("NCAM1",y=185.3597047,x=log10(0.02500546),cex=1.2, font=3)
points(y=185.3597047,x=log10(0.02500546),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("HLA-DRB1",y=174.9454439,x=log10(0.080606061),cex=1.2, font=3)
points(y=174.9454439,x=log10(0.080606061),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("MAP3K19",y=166.8831687,x=log10(0.059821807),cex=1.2, font=3)
points(y=166.8831687,x=log10(0.059821807),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("PBRM1",y=161.8879348,x=log10(0.045563549),cex=1.2, font=3)
points(y=161.8879348,x=log10(0.045563549),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("TMEM131",y=161.558489,x=log10(0.028752512),cex=1.2, font=3)
points(y=161.558489,x=log10(0.028752512),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("MKL2",y=155.0322105,x=log10(0.032251908),cex=1.2, font=3)
points(y=155.0322105,x=log10(0.032251908),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("IRGM",y=153.0337845,x=log10(0.077124183),cex=1.2, font=3)
points(y=153.0337845,x=log10(0.077124183),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("ACAP2",y=150.62,x=log10(0.03318992),cex=1.2, font=3)
points(y=150.62,x=log10(0.03318992),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("Y_RNA",y=144.94,x=log10(0.033121019),cex=1.2, font=3)
points(y=144.94,x=log10(0.033121019),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("TLR6",y=139.5583017,x=log10(0.080424886),cex=1.2, font=3)
points(y=139.5583017,x=log10(0.080424886),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

text("HERC2",y=120.2903388,x=log10(0.022918823),cex=1.2, font=3)
#points(y=225.7883456,x=log10(0.054470426),cex=2,col=rgb(253/255,174/255,97/255,0.7),pch=16)

par(new=T)

#EAST ASIAN
plot((log(as.numeric(eas_geneEnrichment_profile[,3])/as.numeric(eas_geneEnrichment_profile[,2]),base=10)),
     eas_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="",
     cex=2,col=rgb(0, 0, 0,0.5),pch=17,xlim=c(-3,0),ylim=c(0,320))


text("CNTLN",y=275.1980861,x=log10(0.024252891),cex=1.2, font=3)
points(y=275.1980861,x=log10(0.024252891),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)


text("HLA-A",y=262.6173326,x=log10(0.139294927),cex=1.2, font=3)
points(y=262.6173326,x=log10(0.139294927),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)


text("ABI1",y=242.5992067,x=log10(0.066037736),cex=1.2, font=3)
points(y=242.5992067,x=log10(0.066037736),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

#HCG4B
text("HCG4B",y=240.1904905,x=log10(0.199701937),cex=1.2, font=3)
points(y=240.1904905,x=log10(0.199701937),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("CLEC16A",y=209.6543317,x=log10(0.027615063),cex=1.2, font=3)
points(y=209.6543317,x=log10(0.027615063),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("LARP1B",y=207.3548072,x=log10(0.043609371),cex=1.2, font=3)
points(y=207.3548072,x=log10(0.043609371),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("OR5H5P",y=196.3319325,x=log10(0.177394035),cex=1.2, font=3)
points(y=196.3319325,x=log10(0.177394035),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("DOCK9",y=194.3661501,x=log10(0.025805698),cex=1.2, font=3)
points(y=194.3661501,x=log10(0.025805698),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("FRMD4B",y=175.9537992,x=log10(0.019106191),cex=1.2, font=3)
points(y=175.9537992,x=log10(0.019106191),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("MUC4",y=171.0447224,x=log10(0.046561825),cex=1.2, font=3)
points(y=171.0447224,x=log10(0.046561825),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("MSH4",y=167.9487578,x=log10(0.051217765),cex=1.2, font=3)
points(y=167.9487578,x=log10(0.051217765),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)
#HLA-W
text("HLA-W",y=166.9676562,x=log10(0.138069705),cex=1.2, font=3)
points(y=166.9676562,x=log10(0.138069705),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)
#snoU13

text("snoU13",y=166.9099521,x=log10(0.117263844),cex=1.2, font=3)
points(y=166.9099521,x=log10(0.117263844),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("ANKRD36B",y=164.7744489,x=log10(0.066149871),cex=1.2, font=3)
points(y=164.7744489,x=log10(0.066149871),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("C7orf50",y=158.9973764,x=log10(0.028623731),cex=1.2, font=3)
points(y=158.9973764,x=log10(0.028623731),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("POC1B",y=155.6551548,x=log10(0.047752809),cex=1.2, font=3)
points(y=155.6551548,x=log10(0.047752809),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("STRN",y=140.8944602,x=log10(0.032956903),cex=1.2, font=3)
points(y=140.8944602,x=log10(0.032956903),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

text("WFS1",y=140.3146029,x=log10(0.056622354),cex=1.2, font=3)
points(y=140.3146029,x=log10(0.056622354),cex=2,col=rgb(255/255,255/255,191/255,0.7),pch=17)

#text("MEF2C",y=137.598827,x=log10(0.03330139),cex=1.2, font=3)
#text("NBAS",y=130.7644909,x=log10(0.016148224),cex=1.2, font=3)
text("KRT6C",y=129.7747198,x=log10(0.131707317),cex=1.2, font=3)
text("FADS2",y=124.4138731,x=log10(0.04),cex=1.2, font=3)


par(new=T)

#SOUTH ASIANS
plot((log(as.numeric(sas_geneEnrichment_profile[,3])/as.numeric(sas_geneEnrichment_profile[,2]),base=10)),
     sas_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="",
     cex=2,col=rgb(0, 0, 0,0.5),pch=25,bg=rgb(0, 0, 0,0.5),xlim=c(-3,0),ylim=c(0,320))

#HLA-T
text("HLA-T",y=264.71,x=log10(0.322429907),cex=1.2, font=3)
points(x=log10(0.322429907),y=264.71,cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))

text("MAST2",y=252.21,x=log10(0.042212042),cex=1.2, font=3)
points(y=252.21,x=log10(0.042212042),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))

text("EGLN1",y=235.32,x=log10(0.085354478),cex=1.2, font=3)
points(y=235.32,x=log10(0.085354478),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))

text("DDX11",y=230.20,x=log10(0.115056818),cex=1.2, font=3)
points(y=230.20,x=log10(0.115056818),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))


#RP11-202K23.1
text("RP11-202K23.1",y=206.94,x=log10(0.035071668),cex=1.2, font=3)
points(y=206.94,x=log10(0.035071668),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))


#HLA-H
text("HLA-H",y=196.85,x=log10(0.226452906),cex=1.2, font=3)
points(y=196.85,x=log10(0.226452906),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))

text("RP11-202K23.1",y=206.94,x=log10(0.035071668),cex=1.2, font=3)
points(y=206.94,x=log10(0.035071668),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))


#snoU13
text("snoU13",y=181.51,x=log10(0.132464712),cex=1.2, font=3)
points(y=181.51,x=log10(0.132464712),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))


#RP11-385F5.4
text("RP11-385F5.4",y=160.08,x=log10(0.162640902),cex=1.2, font=3)
points(y=160.08,x=log10(0.162640902),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))


text("ENTPD1",y=144.85,x=log10(0.037329505),cex=1.2, font=3)
points(y=144.85,x=log10(0.037329505),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))

#Y_RNA
text("Y_RNA",y=143.49,x=log10(0.03418259),cex=1.2, font=3)
points(y=143.49,x=log10(0.03418259),cex=2,col=rgb(171/255,221/255,164/255,0.5),
       pch=25,bg=rgb(171/255,221/255,164/255,0.5))


#LINC00535
#text("LINC00535",y=139.42,x=log10(0.018041447),cex=1.2, font=3)

#text("FCGR2B",y=125.20,x=log10(0.064966606),cex=1.2, font=3)
#text("CDK8",y=122.93,x=log10(0.03163017),cex=1.2, font=3)

#HLA-L
#text("17",y=120.68,x=log10(0.109137056),cex=1.2, font=3)

# text("MRPS6",y=118.94,x=log10(0.046586345),cex=1.2, font=3)
# text("LRFN5",y=118.76,x=log10(0.019839501),cex=1.2, font=3)
# text("FAM184B",y=115.45,x=log10(0.027512419),cex=1.2, font=3)
# text("MICU1",y=112.30,x=log10(0.022814248),cex=1.2, font=3)
# text("FCGR3B",y=112.06,x=log10(0.117117117),cex=1.2, font=3)
# text("CENPK",y=110.40,x=log10(0.063587684),cex=1.2, font=3)
# text("FOXJ3",y=109.14,x=log10(0.029177719),cex=1.2, font=3)
text("GRM5",y=108.09,x=log10(0.012744413),cex=1.2, font=3)


#AMERICANS

par(new=T)


#AMERICANS
plot((log(as.numeric(amr_geneEnrichment_profile[,3])/as.numeric(amr_geneEnrichment_profile[,2]),base=10)),
     amr_geneEnrichment_profile[,6], xlab= "Log_10(highFST_density)" ,ylab= "-Log_10(P)", main="",
     cex=2,col=rgb(0, 0, 0,0.5),pch=8,xlim=c(-3,0),ylim=c(0,320))

text("TBC1D5",y=308.92,x=log10(0.013980365),cex=1.2, font=3)
points(y=308.92,x=log10(0.013980365),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)

#DPY19L2P2
text("DPY19L2P2",y=282.10,x=log10(0.069524131),cex=1.2, font=3)
points(y=282.10,x=log10(0.069524131),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


text("ATF6",y=235.53,x=log10(0.040917177),cex=1.2, font=3)
points(y=235.53,x=log10(0.040917177),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


text("HLA-DQA1",y=217.84,x=log10(0.060895084),cex=1.2, font=3)
points(y=217.84,x=log10(0.060895084),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


text("PLXNA4",y=209.70,x=log10(0.017043326),cex=1.2, font=3)
points(y=209.70,x=log10(0.017043326),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


text("DCC",y=192.43,x=log10(0.009519264),cex=1.2, font=3)
points(y=192.43,x=log10(0.009519264),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


text("TRIM40",y=184.58,x=log10(0.136138614),cex=1.2, font=3)
points(y=184.58,x=log10(0.136138614),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


#LINC00243
text("LINC00243",y=182.43,x=log10(0.08665269),cex=1.2, font=3)
points(y=182.43,x=log10(0.08665269),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


#AC073869.20
text("AC073869.20",y=176.03,x=log10(0.196687371),cex=1.2, font=3)
points(y=176.03,x=log10(0.196687371),cex=2,col=rgb(43/255,131/255,186/255,0.5),pch=8)


#text("ROBO1",y=137.48,x=log10(0.008555192),cex=1.2, font=3)

#RP11-398M15.1
#text("21",y=137.42,x=log10(0.063384615),cex=1.2, font=3)
#text("LMTK2",y=117.75,x=log10(0.033810539),cex=1.2, font=3)

#LINC00383
#text("22",y=116.44,x=log10(0.041804636),cex=1.2, font=3)
# 
# #RP11-142C4.6
# text("23",y=113.94,x=log10(0.057704918),cex=1.2, font=3)
# 
# text("RFX8",y=108.20,x=log10(0.037290715),cex=1.2, font=3)
# 
# #HLA-DRB6
# text("24",y=107.78,x=log10(0.072916667),cex=1.2, font=3)
# 
# #RP5-991G20.1
# text("25",y=106.78,x=log10(0.039030403),cex=1.2, font=3)
# 
# text("DEPTOR",y=106.19,x=log10(0.021940163),cex=1.2, font=3)
# 
# #RP11-98J9.1
# text("26",y=106.06,x=log10(0.040522876),cex=1.2, font=3)
# 
# text("FAM13A",y=97.49,x=log10(0.013322308),cex=1.2, font=3)
# text("PTPRD",y=96.22,x=log10(0.004254612),cex=1.2, font=3)
# text("IQGAP1",y=94.10,x=log10(0.028537098),cex=1.2, font=3)
# 


abline(h = bonf_thresh,v = dens_thresh,cex=2,col="azure3")

legend(-2.7,250,c("AFR","EUR","EAS","SAS","AMR"),pch=c(15,16,17,25,8),cex=2,col=c(rgb(215/255,25/255,28/255,0.5),rgb(253/255,174/255,97/255,0.7),
                                                  rgb(255/255,255/255,191/255,0.7),rgb(171/255,221/255,164/255,0.5),
                                                  rgb(43/255,131/255,186/255,0.5)))


#Additional Genes
text("EDAR",y=21.53,x=log10(0.010906825),cex=1.2, font=3)
text("LCT",y=49.89,x=log10(0.05134189),cex=1.2, font=3)
text("TLR1",y=49.80351113,x=log10(0.077777778),cex=1.2, font=3)



