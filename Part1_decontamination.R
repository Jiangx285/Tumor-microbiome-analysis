---
title: "Microbiome decontamination pipeline"
author: "Xuan Jiang"
---
  
  
#de1_quantile test
#de2_decontam_freq_0.5
#merge1_quan_fre
#de3_correlation
#merge2_all
##1.separate ABFV file
##2.merge Bacteria file to RData



library(tidyverse)
library(reshape2)



#overall
setwd("D:/Projects/tumor/ESCC/data")
meta=read.csv("table/metadata.csv")
abfv=read.delim("table/bracken_ABFV.txt") %>% 
  filter(rank=="S") %>% column_to_rownames("tax_name") %>% select(!c("rank","taxid")) %>% filter(rowSums(.)>0)
abfvrel=apply(abfv, 2, function(x){x/sum(x)}) %>% as.data.frame() 
  

#------------de1_quantile test---------
temp1 <- abfvrel %>% select(starts_with("NC")) %>% filter(rowSums(.)>0)
temp2 <- abfvrel %>% select(!starts_with("NC")) %>% filter(rowSums(.)>0)
genus_test=intersect(rownames(temp1),rownames(temp2))

two_sample_quantail_test=function(x){
  s1=temp1[x,] %>% as.numeric() %>% na.omit()
  s2=temp2[x,] %>% as.numeric() %>% na.omit()
  test=EnvStats::quantileTest(x = s1,y = s2,
                              target.quantile = 0.95,
                              exact.p = T,
                              alternative = "greater")
  return(test$p.value)
}
p_value=lapply(genus_test,two_sample_quantail_test)
genus_contami=data.frame(genus=genus_test,p_value=unlist(p_value))
genus_contami$p.adj=p.adjust(genus_contami$p_value,method = "fdr")
taxid_con_qua=genus_contami[genus_contami$p.adj <0.05,] %>% .$genus
temp=genus_contami[genus_contami$p.adj <0.05,] %>% rename(taxa=genus) %>% left_join(.,abundance)

#abfv
temp1 <- abfv %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfv %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfv3= temp2 %>% filter(!taxa %in% taxid_con_qua) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfv3[is.na(abfv3)]=0

#abfvrel
temp1 <- abfvrel %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfvrel %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfvrel3= temp2 %>% filter(!taxa %in% taxid_con_qua) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfvrel3[is.na(abfvrel3)]=0


#------------de2_decontam_freq_0.5---------
library(decontam)
library(scales)
#
meta2=meta %>% mutate(check=if_else(sample_type=="NC",TRUE,FALSE))
abfv44=abfv
abfv44[abfv44==0]=10
abfv44=abfv44[,meta2$sample] %>%  filter(rowSums(.)>0)
temp <- isContaminant(t(abfv44), method="frequency",neg=meta2$check,
                      conc=meta2$Amount_Bacteria_ng_50ul,threshold = 0.5)
taxid_con_fre=temp %>% filter(contaminant=="TRUE") %>% rownames(.)

#abfv
temp1 <- abfv %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfv %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfv4= temp2 %>% filter(!taxa %in% taxid_con_fre) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfv4[is.na(abfv4)]=0

#abfvrel
temp1 <- abfvrel %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfvrel %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfvrel4= temp2 %>% filter(!taxa %in% taxid_con_fre) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfvrel4[is.na(abfvrel4)]=0




#------------merge1_quan_fre---------
#
taxid_con_mer=unique(c(taxid_con_fre,taxid_con_qua))
#abfv
temp1 <- abfv %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfv %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfv5= temp2 %>% filter(!taxa %in% taxid_con_mer) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfv5[is.na(abfv5)]=0

#abfvrel
temp1 <- abfvrel %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfvrel %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfvrel5= temp2 %>% filter(!taxa %in% taxid_con_mer) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfvrel5[is.na(abfvrel5)]=0




#------------de3_correlation---------
##1.caculate
#row:sample, column:taxa
fast_spearman_fdr <- function(otu_niche){
  r <- cor(otu_niche,method="spearman")
  n = nrow(otu_niche)
  t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
  p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
  return(list(r=r,p=p))
}

abfvrel2=abfvrel
genus_corr <- fast_spearman_fdr(t(abfvrel2))
genus_corr$p.adj <- p.adjust(genus_corr$p,method = "fdr")
#
tmp1=matrix(genus_corr$p.adj,nr=nrow(abfvrel2))
colnames(tmp1)=colnames(genus_corr$p)
rownames(tmp1)=rownames(genus_corr$p)

plotdata1 <- genus_corr$p  #p value 
plotdata1[lower.tri(plotdata1,diag = T)]=NA
plotdata1= melt(as.matrix(plotdata1),varnames=c("taxa1","taxa2")) %>% filter(!is.na(.$value)) %>% rename(p_value=value)
#
plotdata2 <- tmp1 %>% as.matrix() #this p value is fdr
plotdata2[lower.tri(plotdata2,diag = T)]=NA
plotdata2= melt(as.matrix(plotdata2),varnames=c("taxa1","taxa2")) %>% filter(!is.na(.$value)) %>% rename(p_adjust=value)
#
plotdata <- genus_corr$r
plotdata[lower.tri(plotdata,diag = T)]=NA
plotdata= melt(as.matrix(plotdata),varnames=c("taxa1","taxa2")) %>% filter(!is.na(.$value)) %>% rename(r_value=value) %>% 
  left_join(.,plotdata1) %>% left_join(.,plotdata2)
#
rm(genus_corr,plotdata1,plotdata2)


taxid_con_cor=plotdata %>% filter(taxa1 %in% taxid_con_mer) %>% filter(!taxa2 %in% taxid_con_mer) %>% 
  filter(p_adjust>0) %>% filter(p_adjust<0.05) %>% 
  filter(r_value>0.7) %>% .$taxa2 %>% unique() %>% as.character()



#------------merge2_all---------
taxid_con_all=unique(c(taxid_con_fre,taxid_con_qua,taxid_con_cor))
#abfv
temp1 <- abfv %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfv %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfv6= temp2 %>% filter(!taxa %in% taxid_con_all) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfv6[is.na(abfv6)]=0

#abfvrel
temp1 <- abfvrel %>% select(starts_with("NC")) %>% rownames_to_column("taxa")
temp2 <- abfvrel %>% select(!starts_with("NC")) %>% rownames_to_column("taxa")
abfvrel6= temp2 %>% filter(!taxa %in% taxid_con_all) %>% full_join(.,temp1) %>% column_to_rownames("taxa")
abfvrel6[is.na(abfvrel6)]=0

#save
#save(meta,abfv,abfvrel,abfv3,abfvrel3,abfv4,abfvrel4,abfv5,abfvrel5,abfv6,abfvrel6,
#     taxid_con_qua,taxid_con_fre,taxid_con_mer,taxid_con_all,taxid_con_cor,
#     file = "1_decontam_S/data_abfvabfvrel_bracken.RData")


###1.separate ABFV file-----------
#abfvid.label.RData is the taxonomy.tbl from ncbi
load("abfvid.label.RData")
taxid_de=data.frame(taxa=setdiff(rownames(abfvrel),taxid_con_all)) %>% left_join(taxid[,c("taxa","label")])
#
temp1=taxid_de %>% filter(label == "Bacteria")
temp2=taxid_de %>% filter(label == "Virus")
temp3=taxid_de %>% filter(label == "Archaea")
temp4=taxid_de %>% filter(label == "Fungi")

abfv11=abfv[temp1$taxa,] %>% select(!starts_with("NC")) %>% filter(rowSums(.)>0)
abfv22=abfv[temp2$taxa,] 
abfv33=abfv[temp3$taxa,]
abfv44=abfv[temp4$taxa,]

#write.table(abfv11,"1_decontam_S/bracken.S.254.bacteria.txt",sep="\t",row.names = T,quote = F)
#write.table(abfv22,"1_decontam_S/bracken.S.254.virus.txt",sep="\t",row.names = T,quote = F)
#write.table(abfv33,"1_decontam_S/bracken.S.254.arhcaea.txt",sep="\t",row.names = T,quote = F)
#write.table(abfv44,"1_decontam_S/bracken.S.254.fungi.txt",sep="\t",row.names = T,quote = F)


###2.merge Bacteria file to RData---------------------
#threshold Bacteria > 50000
meta=read.csv("table/metadata.csv",row.names = 1) %>% filter(!sample_type=="NC")
abfv=read.delim(paste("1_decontam_S/bracken.S.254.bacteria.txt",sep="")) %>% 
  select(meta$sample) %>% filter(rowSums(.)>0)
id1=t(abfv) %>% as.data.frame() %>% transmute(de_reads=rowSums(.)) %>% rownames_to_column("sample") %>% 
  filter(de_reads>50000)
meta=meta[id1$sample,]
abfv=abfv[,id1$sample] %>% filter(rowSums(.)>0)
abfvrel=apply(abfv, 2, function(x){x/sum(x)}) %>% as.data.frame() 
abfvrel=abfvrel[rowSums(abfvrel>0.001)>0,]
abfv=abfv[rownames(abfvrel),]
abfvrel=apply(abfv, 2, function(x){x/sum(x)}) %>% as.data.frame() 
#
tmp1=apply(abfvrel,1,median) %>% as.data.frame()
colnames(tmp1)="median"
tmp2=apply(abfvrel,1,mean) %>% as.data.frame()
colnames(tmp2)="mean"
tmp3=apply(abfvrel,1,function(x) sum(x > 0)) %>% as.data.frame()
colnames(tmp3)="prevalence"
abundance=cbind(tmp1,tmp2,tmp3) %>% rownames_to_column("taxa") 
#save(meta,abfv,abfvrel,abundance,file = "table/use_taxa_bracken_data.RData")
#2.paired sample
meta=meta %>% group_by(ID) %>% mutate(check=n()) %>% filter(check==2) %>% select(!"check") %>% 
  mutate(tmp=sample) %>% column_to_rownames("tmp")
abfv=abfv[,meta$sample] %>% filter(rowSums(.)>0)
abfvrel=apply(abfv, 2, function(x){x/sum(x)}) %>% as.data.frame() 
#save(meta,abfv,abfvrel,abundance,file = "table/use_taxa_bracken_paired_data.RData")


