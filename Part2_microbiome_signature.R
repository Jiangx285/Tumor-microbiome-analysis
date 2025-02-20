---
title: "Microbiome signature analysis"
author: "Xuan Jiang"
note: "Figure1 & Supp.Figure1"
---



library(readxl)
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggpubr)
library(Hmisc)




setwd("D:/Projects/Tumor/ESCC/data")
#1_NC & sample-------------
##S1A tSNE-----------
library(Rtsne)
library(vegan)
library(phyloseq)
library(ade4)
library(picante)
conflicted::conflicts_prefer(phyloseq::distance)
set.seed(123)
#1.before
load("table/use_taxa_bracken_data.RData")
abfvrel=read.delim("table/abfv_sort/bracken_ABFV.txt") %>% 
  filter(rank=="S") %>% column_to_rownames("tax_name") %>% select(!colnames(.)[1:2]) %>% 
  filter(rowSums(.)>0)
abfv_jsd <- vegdist(t(abfvrel), method = "bray")
tsne_out <- Rtsne(abfv_jsd, dims = 2) 
meta=read.csv("table/metadata.csv",row.names = 1) 
meta=meta[colnames(abfvrel),] %>% mutate(sample_type=ifelse(sample_type=="NC","NC","samples"))
plot1=tsne_out$Y %>% as.data.frame() %>% mutate(sample=colnames(abfvrel)) %>%
  left_join(.,meta[,c("sample","sample_type")])
#2.after
load("table/use_taxa_bracken_data.RData")
meta=read.csv("table/metadata.csv",row.names = 1) 
tmp1=read.delim("1_decontam_S/bracken_ABFV.txt") %>% 
  select(starts_with("NC")) %>% filter(rowSums(.)>0) 
tmp1=apply(tmp1,2,function(x){x/sum(x)}) %>% as.data.frame() %>% rownames_to_column("taxa")
abfvrel=abfvrel %>% rownames_to_column("taxa") %>% full_join(.,tmp1) %>% column_to_rownames("taxa")
abfvrel[is.na(abfvrel)]=0
meta=meta[colnames(abfvrel),] %>% mutate(sample_type=ifelse(sample_type=="NC","NC","samples"))
#
abfv_jsd <- vegdist(t(abfvrel), method = "bray")
tsne_out <- Rtsne(abfv_jsd, dims = 2) 
plot2=tsne_out$Y %>% as.data.frame() %>% mutate(sample=colnames(abfvrel)) %>%
  left_join(.,meta[,c("sample","sample_type")])
#
ggplot(plot1, aes(V1, V2,group=sample_type,color=sample_type)) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.2) + 
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.2) +
  geom_point(size = 0.8)+
  stat_ellipse(level = 0.95,linewidth = 0.1)+
  scale_color_manual(values=c("#FB8072","#80B1D3"))+
  theme_classic()+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, linewidth = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        axis.title.x=element_text(size=12,color = 'black'),
        axis.title.y=element_text(size=12,color = 'black'),
        axis.text.x=element_text(size=12,color = 'black'),
        axis.text.y=element_text(size=12,color = 'black'),
        legend.title=element_blank(),
        legend.position = "none")+
  labs(#title = paste(" PERMANOVA: R2=",round(temp1$R2[1],2),", p=",temp1$`Pr(>F)`[1],sep=""),
    title="Before decontamination",
    x = "t-SNE 1",
    y = "t-SNE 2")+
ggplot(plot2, aes(V1, V2,group=sample_type,color=sample_type)) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.2) + 
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.2) +
  geom_point(size = 0.8)+
  stat_ellipse(level = 0.95,linewidth = 0.1)+
  scale_color_manual(values=c("#FB8072","#80B1D3"))+
  theme_classic()+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, linewidth = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        axis.title.x=element_text(size=12,color = 'black'),
        axis.title.y=element_text(size=12,color = 'black'),
        axis.text.x=element_text(size=12,color = 'black'),
        axis.text.y=element_text(size=12,color = 'black'),
        legend.title=element_blank(),
        legend.text=element_text(size=12,color = 'black'),
        legend.position = "right")+
  labs(#title = paste(" PERMANOVA: R2=",round(temp1$R2[1],2),", p=",temp1$`Pr(>F)`[1],sep=""),
    title="After decontamination",
       x = "t-SNE 1",
       y = "t-SNE 2")




#2_16S & meta----------------------------------------------------------
##S1B Procruster---------
library(vegan)
library(ade4)
load("table/ampli/abfv_16S_taxa_decontam.RData")
abfvrel_16S=abfvrel
load("table/use_taxa_bracken_data.RData")
id1=intersect(colnames(abfvrel_16S),colnames(abfvrel))
abfvrel_16S=abfvrel_16S[,id1] %>% filter(rowSums(.)>0)
abfvrel_meta=abfvrel[,id1] %>% filter(rowSums(.)>0)
meta=meta[id1,]
abfvrel_16S2=abfvrel_16S %>% setNames(paste("16S",colnames(.),sep="__")) %>% rownames_to_column("taxa")
abfvrel_meta2=abfvrel_meta %>% setNames(paste("meta",colnames(.),sep="__")) %>% rownames_to_column("taxa")
abfvrel2=full_join(abfvrel_16S2,abfvrel_meta2) %>% column_to_rownames("taxa")
abfvrel2[is.na(abfvrel2)]=0
##
##
abfvrel_16S3=abfvrel_16S2 %>% column_to_rownames("taxa") %>% t()
abfvrel_meta3=abfvrel_meta2 %>% column_to_rownames("taxa") %>% t()
##
env.bray <- vegdist(abfvrel_16S3, method = "bray")
add <-  !(is.euclid(env.bray))
pcoa.env <- cmdscale(env.bray, k = nrow(abfvrel_16S3)-1, eig = TRUE, add = add)
ordiplot(pcoa.env, type = "text", main = "PCoA for 16S Data, Bray Distances")
spe.bray <- vegdist(abfvrel_meta3, method = "bray")
add <-  !(is.euclid(spe.bray))
pcoa.spe <- cmdscale(spe.bray, k = nrow(abfvrel_meta3)-1, eig = TRUE, add = add)
ordiplot(pcoa.spe, type = "text", main = "PCoA for meta Data, Bray Distances")
#
pro.1 <- procrustes(pcoa.env, pcoa.spe, symmetric = FALSE, scores = "sites",  
                    choices = c(1,2))
prot.1 <- protest(X = pcoa.env, Y = pcoa.spe, permutations = how(nperm = 999))
prot.1$signif  #p 值
prot.1$ss  #偏差平方和 M2 统计量

#> prot.1$signif  #p 值
#[1] 0.001
#> prot.1$ss  #偏差平方和 M2 统计量
#[1] 0.372828


#plot
Y <- cbind(data.frame(pro.1$Yrot), data.frame(pro.1$X))
X <- data.frame(pro.1$rotation)
group <- data.frame(samples=colnames(abfvrel2)) %>% 
  separate(samples,c("type","ID"),sep="__",remove = F) 
Y$samples <- rownames(Y)
plot1=Y %>% mutate(type1="16S",type2="meta")

#ggplot2
ggplot(plot1) +
  geom_segment(aes(x = X1, y = X2, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.1, 'cm')),
               color = 'grey', size = 0.3) +
  geom_point(aes(X1, X2), size = 1.5, shape = 16,color="#80B1D3") +
  geom_point(aes(Dim1, Dim2), size = 1.5, shape = 16,color="#FB8072") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'PCoA1', y = 'PCoA2') +
  annotate('text', label = sprintf('Protest: M^2 == 0.373'),
           x = -0.4, y = 0.6, size = 3, parse = TRUE) +
  annotate('text', label = paste0("P = ",prot.1$signif),
           x = -0.4, y = 0.53, size = 3)




##S1C correlation---------
load("table/ampli/abfv_16S_taxa_decontam.RData")
abfvrel_16S=abfvrel
load("table/use_taxa_bracken_data.RData")
abfvrel_meta=abfvrel %>% rownames_to_column("taxa") %>% mutate(genus=word(taxa,start = 1)) %>%   
  group_by(genus) %>%  
  summarise(across(grep("^Y|^L", colnames(.), value = TRUE), sum)) %>% 
  column_to_rownames("genus")
id1=intersect(colnames(abfvrel_16S),colnames(abfvrel))
abfvrel_16S=abfvrel_16S[,id1] %>% filter(rowSums(.)>0)
abfvrel_meta=abfvrel[,id1] %>% filter(rowSums(.)>0)
meta=meta[id1,]
tmp=t(abfv) %>% as.data.frame() %>% transmute(sum=rowSums(.)) %>% rownames_to_column("sample")
##
id1=intersect(rownames(abfvrel_16S),rownames(abfvrel_meta))
#
#correlation
abfvrel_16S2=abfvrel_16S[id1,] %>% t() %>% as.data.frame() %>% setNames(paste("16S",colnames(.),sep="_")) 
abfvrel_meta2=abfvrel_meta[id1,] %>% t() %>% as.data.frame() %>% setNames(paste("meta",colnames(.),sep="_"))
tmp=cbind(abfvrel_16S2,abfvrel_meta2)
result1=rcorr(as.matrix(tmp),type = "spearman")
result1$P.adj=matrix(p.adjust(result1$P, method = 'fdr'),nr=ncol(tmp))
colnames(result1$P.adj)=colnames(result1$P)
rownames(result1$P.adj)=rownames(result1$P)
heatdata=result1$r[colnames(abfvrel_16S2),colnames(abfvrel_meta2)]
p_sum=result1$P.adj[colnames(abfvrel_16S2),colnames(abfvrel_meta2)]
#
tmp1=apply(abfvrel_meta,1,median) %>% as.data.frame() %>% setNames("median")
tmp2=apply(abfvrel_meta,1,mean) %>% as.data.frame() %>% setNames("mean")
tmp3=apply(abfvrel_meta,1,function(x) sum(x > 0)) %>% as.data.frame() %>% setNames("prevalence")
abundance=cbind(tmp1,tmp2,tmp3) %>% rownames_to_column("taxa") 
#
tmp1=apply(abfvrel_16S,1,median) %>% as.data.frame() %>% setNames("median")
tmp2=apply(abfvrel_16S,1,mean) %>% as.data.frame() %>% setNames("mean")
tmp3=apply(abfvrel_16S,1,function(x) sum(x > 0)) %>% as.data.frame() %>% setNames("prevalence")
abundance2=cbind(tmp1,tmp2,tmp3) %>% rownames_to_column("taxa") 
tmp=melt(p_sum,varnames = c("id1","id2")) %>% 
  separate(id1, c("type1", "taxa1"), "_") %>% 
  separate(id2, c("type2", "taxa2"), "_") %>% filter(taxa1==taxa2) %>% rename(taxa=taxa1) %>% rename(p=value)
plot1=melt(heatdata,varnames = c("id1","id2")) %>% 
  separate(id1, c("type1", "taxa1"), "_") %>% 
  separate(id2, c("type2", "taxa2"), "_") %>% filter(taxa1==taxa2) %>% 
  rename(taxa=taxa1) %>% left_join(.,abundance) %>% #select(c("taxa","value","mean")) %>% 
  rename(r=value) %>% 
  left_join(tmp[,c("taxa","p")]) %>% mutate(label=ifelse(p <= 0.0001,"****",
                                                         ifelse(p > 0.0001& p <= 0.001,"***",
                                                         ifelse(p > 0.001& p <= 0.01,"**",
                                                                ifelse(p > 0.01 & p < 0.05,"*",""))))) %>% 
  .[order(-.$r),] %>% mutate(taxa=factor(taxa,levels=(taxa)))

tmp=cor.test(plot1$mean,plot1$r,method = "spearman")
ggplot(data=plot1,aes(x=mean, y=r))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(title = paste0("r = ",tmp$estimate,", p = ",tmp$p.value),x="Mean abundance in metagenome", y= "Spearman'r")+
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.title.y=element_text(size=14,color = "black"),
        strip.text = element_text(size = 14),
        legend.position = "right")


p1=ggplot(data=plot1,aes(x=taxa, y=r,fill=log10(mean)))+
  geom_bar(stat="identity")+
  scale_fill_gradient(low = "#D7EDF9", high = "#4390BA")+
  labs(x="", y= "Spearman r value")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=11,angle=90,vjust = 0.5,hjust = 0.98,color = "black"),
        axis.text.y=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.title.y=element_text(size=14,color = "black"),
        strip.text = element_text(size = 14),
        legend.position = "right")
p2=ggplot(plot1) +
  geom_text(aes(x = taxa, y = 0, label = label),hjust = 0,size=6,angle=90) +
  theme_void()
p=p1/p2+
  plot_layout(heights = c(5, 2))
p



#3_PERMANOVA-----------------
##1B PERMANOVA--------
load("table/use_taxa_bracken_paired_data.RData")
ID=meta$ID
abfvrel2=abfvrel[,rownames(meta)]
tmp1=adonis2(as.data.frame(t(abfvrel2)) ~ sample_type, 
             data = meta, permutations = 999, 
             strata = ID,
             method="bray",by="margin")
load("table/use_taxa_bracken_data.RData")
meta22=meta %>% filter(sample_type=="T")
abfvrel2=abfvrel[,rownames(meta22)]
tmp2=adonis2(as.data.frame(t(abfvrel2)) ~ sex, 
             data = meta22, permutations = 999, 
             method="bray",by="margin")
meta22=meta %>% filter(sample_type=="T")
abfvrel2=abfvrel[,rownames(meta22)]
tmp3=adonis2(as.data.frame(t(abfvrel2)) ~ age, 
             data = meta22, permutations = 999, 
             method="bray",by="margin")
tmp4=adonis2(as.data.frame(t(abfvrel2)) ~ stagegp, 
             data = meta22, permutations = 999, 
             method="bray",by="margin")
meta22=meta %>% filter(!is.na(survivalgp)) %>% filter(sample_type=="T")
abfvrel2=abfvrel[,rownames(meta22)]
tmp5=adonis2(as.data.frame(t(abfvrel2)) ~ survivalgp, 
             data = meta22, permutations = 999, 
             method="bray",by="margin")
tmp=read.csv("table/metadata_smoke_drink.csv") 
meta2=meta %>% filter(sample_type=="T") %>% left_join(.,tmp)
rownames(meta2)=meta2$sample
meta22=meta2 %>% filter(!is.na(Smoking.status))
abfvrel2=abfvrel[,rownames(meta22)]
tmp6=adonis2(as.data.frame(t(abfvrel2)) ~ Smoking.status, 
             data = meta22, permutations = 999, 
             method="bray",by="margin")
tmp7=adonis2(as.data.frame(t(abfvrel2)) ~ Drinking.status, 
             data = meta22, permutations = 999, 
             method="bray",by="margin")
tmp11=data.frame(type="Sample type",R2=round(tmp1$R2[1],3),p=tmp1$`Pr(>F)`[1])
tmp22=data.frame(type="Gender",R2=round(tmp2$R2[1],3),p=tmp2$`Pr(>F)`[1])
tmp33=data.frame(type="Age",R2=round(tmp3$R2[1],3),p=tmp3$`Pr(>F)`[1])
tmp44=data.frame(type="Stage group",R2=round(tmp4$R2[1],3),p=tmp4$`Pr(>F)`[1])
tmp55=data.frame(type="Survival group",R2=round(tmp5$R2[1],3),p=tmp5$`Pr(>F)`[1])
tmp66=data.frame(type="Smoking",R2=round(tmp6$R2[1],3),p=tmp6$`Pr(>F)`[1])
tmp77=data.frame(type="Alcohol consumption",R2=round(tmp7$R2[1],3),p=tmp7$`Pr(>F)`[1])
plot1=rbind(tmp11,tmp22,tmp33,tmp44,tmp55,tmp66,tmp77) %>% 
  mutate(label=ifelse(p <= 0.001,"***",
                      ifelse(p > 0.001& p <= 0.01,"**",
                             ifelse(p > 0.01 & p < 0.05,"*","")))) %>% 
  .[order(-.$p),] %>% mutate(type=factor(type,levels=type))

ggplot(plot1,aes(x=R2,y=(type),fill=-log10(p)))+
  geom_bar(stat = "identity",width = 0.5)+
  labs(title="PERMANOVA",y="")+
  theme_bw()+
  theme(
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    axis.text.x=element_text(size=12,angle=0,color="black"),
    axis.text.y=element_text(size=12,angle=0,color="black")) 





#4_DNA load qPCR-------------------
##1E qPCR------------
#Raw_bacteria normalized to DNA qubit.
color1=c("#80B1D3", "#FB8072")
plot1 <- read_excel("table/qPCR_all.xlsx",sheet = "Raw_stdCT")
#
plot2=plot1 %>% mutate(sample_type=ifelse(sample_type=="N","Normal","Tumor"))
ggplot(plot2,aes(sample_type,`log10((1/2^Ct)/qubit*10^14)`,fill=sample_type))+
  #geom_point(size=0.5)+
  geom_line(aes(group=ID),color="grey",linewidth=0.2)+
  geom_boxplot(width=0.8,linewidth=0.3,outlier.size = 0.5)+
  #scale_y_log10()+
  labs(y="Bacteria /total DNA ratio(log10)",x="")+
  scale_color_manual(values=color1)+
  scale_fill_manual(values=color1)+
  stat_compare_means(paired = T,label="p.signif")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=14,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.title.y=element_text(size=14,color = "black"),
        strip.text = element_text(size = 14),
        legend.position = "none")



#5_alpha diversity-----------
##S1D alpha------
library(vegan)
library(picante)
library(ggpubr)
library(psych)
load("table/use_taxa_bracken_data.RData")
tmp=read.csv("table/metadata.csv") 
meta=meta %>% left_join(.,tmp)
rownames(meta)=meta$sample

abfv_rare <-  as.data.frame(t(rrarefy(t(abfv), min(colSums(abfv)))))
tmp1 <- diversity(t(abfv_rare), index = 'shannon') %>% as.data.frame() %>% rownames_to_column("sample") %>% 
  mutate(type="Shannon") %>% rename(value=".")
tmp2 <- rowSums(t(abfv_rare) >0) %>% as.data.frame() %>% rownames_to_column("sample") %>% 
  mutate(type="Richness") %>% rename(value=".")
plotdata1 <-  rbind(tmp1,tmp2) %>% left_join(.,meta) %>% 
  mutate(sample_type=ifelse(sample_type=="N","Normal","Tumor"))
plotdata2 <- plotdata1 %>% group_by(ID) %>% mutate(check=n()) %>% filter(check==4)
plotdata3 <- plotdata1 %>% filter(!is.na(Survival.time))

color1=c("#80B1D3", "#FB8072")
ggplot(plotdata2,aes(sample_type,value,fill=sample_type))+
  geom_violin(width=0.5,linewidth=0.3)+
  geom_boxplot(width=0.2,fill="white",linewidth=0.3)+
  #geom_jitter(size=0.9)+
  stat_compare_means(paired = T,label = "p.signif")+
  scale_color_manual(values=color1)+
  scale_fill_manual(values=color1)+
  facet_wrap(.~type,scales = "free")+
  labs(x="")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=14,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.title.y=element_text(size=14,color = "black"),
        strip.text = element_text(size = 14),
        legend.position = "none")


#6_Top composition -----------
##S1E compo-----------
##1 species
source("R_script/store/1_composition.R")
source("R_script/store/00_color.R")
color=c("#B7B7B7",color_21,color_21,color_20)
#
load("table/use_taxa_bracken_paired_data.RData")
#N
id1=meta %>% filter(sample_type=="N") %>% .$sample
abfvrel1=abfvrel[,id1] %>% filter(rowSums(.)>0)
tmp1=abfvrel1 %>% transmute(value=rowSums(.)/ncol(.)) %>% .[order(-.$value),,drop=F] %>% rename(Normal=value)
t1=rownames(tmp1)[1:15] 
#T
id1=meta %>% filter(sample_type=="T") %>% .$sample
abfvrel1=abfvrel[,id1] %>% filter(rowSums(.)>0)
tmp2=abfvrel1 %>% transmute(value=rowSums(.)/ncol(.)) %>% .[order(-.$value),,drop=F] %>% rename(Tumor=value)
t2=rownames(tmp2)[1:15] 
id_taxa=unique(c(t1,t2))
tmp3=cbind(tmp1[id_taxa,,drop=F],tmp2[id_taxa,,drop=F])
tmp3["Others",]=apply(tmp3, 2, function(x){1-sum(x)})
plot1 = melt(as.matrix(tmp3),varnames=c("taxa","sample")) %>% 
  mutate(taxa=factor(.$taxa,levels = rev(rownames(tmp3)))) %>% 
  left_join(.,meta[,c("sample","sample_type")])

##2 genus
source("R_script/store/00_color.R")
color=c("#B7B7B7",color_21,color_21,color_20)
#
load("table/use_taxa_bracken_paired_data.RData")
abfvrel=abfvrel %>% rownames_to_column("taxa") %>% mutate(genus=word(taxa,start = 1)) %>%   
  group_by(genus) %>%  
  summarise(across(grep("^Y|^L", colnames(.), value = TRUE), sum)) %>% 
  column_to_rownames("genus")
meta=meta %>% group_by(ID) %>% mutate(n=n()) %>% filter(n==2)
abfvrel=abfvrel[,meta$sample]
#N
id1=meta %>% filter(sample_type=="N") %>% .$sample
abfvrel1=abfvrel[,id1] %>% filter(rowSums(.)>0)
tmp1=abfvrel1 %>% transmute(value=rowSums(.)/ncol(.)) %>% .[order(-.$value),,drop=F] %>% rename(Normal=value)
t1=rownames(tmp1)[1:15] 
#T
id1=meta %>% filter(sample_type=="T") %>% .$sample
abfvrel1=abfvrel[,id1] %>% filter(rowSums(.)>0)
tmp2=abfvrel1 %>% transmute(value=rowSums(.)/ncol(.)) %>% .[order(-.$value),,drop=F] %>% rename(Tumor=value)
t2=rownames(tmp2)[1:15] 
id_taxa=unique(c(t1,t2))
tmp3=cbind(tmp1[id_taxa,,drop=F],tmp2[id_taxa,,drop=F])
tmp3["Others",]=apply(tmp3, 2, function(x){1-sum(x)})
plot2 = melt(as.matrix(tmp3),varnames=c("taxa","sample")) %>% 
  mutate(taxa=factor(.$taxa,levels = rev(rownames(tmp3)))) %>% 
  left_join(.,meta[,c("sample","sample_type")])


##3 phylum
load("table/use_taxa_bracken_paired_data.RData")
load("abfvid.label.RData")
species_df=taxid %>% filter(taxa %in% rownames(abfvrel))
phylum_df=taxid %>% filter(rank=="phylum")
species_df2=c()
for (i in 1:nrow(phylum_df)) {
  tmp=species_df %>%
    rowwise() %>% #Ensures that the mutate() operation is applied row by row.
    mutate(phylum_taxa = ifelse(any(str_detect(path, phylum_df$path[i])), phylum_df$taxa[i], NA_character_)) %>% 
    filter(!is.na(phylum_taxa))
  species_df2=rbind(species_df2,tmp)
}
#
abfvrel0=abfvrel %>% rownames_to_column("taxa") %>% 
  left_join(.,species_df2[,c("taxa","phylum_taxa")]) %>% 
  select("phylum_taxa",everything()) %>% filter(!is.na(phylum_taxa)) %>% 
  select(!"taxa") %>% group_by(phylum_taxa) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup()  %>% column_to_rownames("phylum_taxa")
#N
id1=meta %>% filter(sample_type=="N") %>% .$sample
abfvrel1=abfvrel0[,id1] %>% filter(rowSums(.)>0)
tmp1=abfvrel1 %>% transmute(value=rowSums(.)/ncol(.)) %>% .[order(-.$value),,drop=F] %>% rename(Normal=value)
t1=rownames(tmp1)
#T
id1=meta %>% filter(sample_type=="T") %>% .$sample
abfvrel1=abfvrel0[,id1] %>% filter(rowSums(.)>0)
tmp2=abfvrel1 %>% transmute(value=rowSums(.)/ncol(.)) %>% .[order(-.$value),,drop=F] %>% rename(Tumor=value)
t2=rownames(tmp2)
id_taxa=unique(c(t1,t2))
tmp3=cbind(tmp1[id_taxa,,drop=F],tmp2[id_taxa,,drop=F])
tmp3["Others",]=apply(tmp3, 2, function(x){1-sum(x)})
plot3 = melt(as.matrix(tmp3),varnames=c("taxa","sample")) %>% 
  mutate(taxa=factor(.$taxa,levels = rev(rownames(tmp3)))) %>% 
  left_join(.,meta[,c("sample","sample_type")])
#
source("R_script/store/00_color.R")
color=c("#B7B7B7",color_21,color_21,color_20)
#
ggplot(data=plot1,aes(x=sample, y=value, fill=taxa))+
  geom_bar(stat='identity')+
  labs(title="",y="Relative abundance",x="")+
  scale_fill_manual(values=c(color))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black",size=12),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.title.y=element_text(size=12,color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "right")+
ggplot(data=plot2,aes(x=sample, y=value, fill=taxa))+
  geom_bar(stat='identity')+
  labs(title="",y="Relative abundance",x="")+
  scale_fill_manual(values=c(color))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black",size=12),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.title.y=element_text(size=12,color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "right")+
ggplot(data=plot3,aes(x=sample, y=value, fill=taxa))+
  geom_bar(stat='identity',width = 0.5)+
  labs(title="",y="Relative abundance",x="")+
  scale_fill_manual(values=c(color))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black",size=12),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.title.y=element_text(size=12,color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "right")




#7_aero-------------------
##1.process NT diff using ALDEx2--------
library(ALDEx2)
cacu_DE_NT=function(meta,abfv,abfvrel){
  #need count data
  meta2 <- meta[,c("sample","sample_type")] %>% .[order(.$sample_type),] 
  abfv2 <- abfv[,meta2$sample] %>% round()
  aldex.sample <- aldex.clr(abfv2, meta2$sample_type, mc.samples=128, verbose=TRUE)
  plotdata_p <- aldex.ttest(aldex.sample, paired.test=T, hist.plot = F) %>% rownames_to_column("taxa") %>% 
    mutate(p_value=wi.ep,p_adjust=wi.eBH,type="ALDEx2") %>% dplyr::select(c("taxa","p_value","p_adjust","type")) 
  plotdata_p1=plotdata_p
  plotdata_p1.effect <- aldex.effect(aldex.sample, CI=T, verbose=F, include.sample.summary=F, 
                                     paired.test=T) %>% rownames_to_column("taxa")
  res_list=list(plotdata_p1=plotdata_p1,plotdata_p1.effect=plotdata_p1.effect)
  return(res_list)
}
#
load("table/use_taxa_bracken_paired_data.RData")
meta2=meta
tmp=abundance %>% mutate(prev=prevalence/ncol(abfvrel)) %>% filter(prev>0.3)
abfv2=round(abfv[tmp$taxa,meta2$sample]) 
abfvrel2=abfvrel[rownames(abfv2),colnames(abfv2)]
plot_NT=cacu_DE_NT(meta=meta2,abfv = abfv2,abfvrel=abfvrel2)
res1=plot_NT$plotdata_p1
res2=plot_NT$plotdata_p1.effect
plot1=res2 %>% left_join(res1,.) %>% left_join(.,abundance) %>% 
  mutate(type2=ifelse(p_value > 0.05,"Unchanged",ifelse(diff.btw > 0,"ESCC > NOR","ESCC < NOR"))) %>% 
  mutate(type2=gsub("ESCC > NOR","Tumor > Normal",type2)) %>%
  mutate(type2=gsub("ESCC < NOR","Tumor < Normal",type2)) 
save(res1,res2,"res.DE.RData")


##1F aerophilicity----------
#1.oxygen
library(bugphyzz)
dat_bugphyzz=bugphyzz::importBugphyzz()
tmp1=dat_bugphyzz$aerophilicity %>% mutate(taxa=Taxon_name,genus=Taxon_name) %>% filter(Score > 0.5)
tmp2=plot1[,"taxa",drop=F] %>% left_join(.,tmp1[,c("taxa", "Attribute","Attribute_value")]) %>%
  filter(!is.na(Attribute)) %>% 
  group_by(taxa) %>% mutate(n=n()) %>% 
  summarise(type_phy = paste(Attribute_value, collapse = ",")) %>% ungroup() %>% 
  left_join(.,tmp1[,c("taxa", "Score")]) 
#function of summarise is make every taxa only one row, make two feature to one row
tmp3=plot1[,"taxa",drop=F] %>% filter(!taxa %in% tmp2$taxa) %>% separate(taxa,c("genus"),sep=" ",remove = F) %>% 
  left_join(.,tmp1[,c("genus", "Attribute","Attribute_value")]) %>%
  filter(!is.na(Attribute)) %>% 
  group_by(taxa) %>% mutate(n=n()) %>% 
  summarise(type_phy = paste(Attribute_value, collapse = ",")) %>% ungroup() %>% 
  separate(taxa,c("genus"),sep=" ",remove = F) %>% left_join(.,tmp1[,c("genus", "Score")]) %>%  
  select(colnames(tmp2)) %>% rbind(.,tmp2) 
plot2=plot1 %>% left_join(.,tmp3) 
plot2$type_phy[is.na(plot2$type_phy)]="unknown"
plot22=plot2[,c("taxa","type2","type_phy")] %>% left_join(.,tmp1[,c("taxa", "Attribute","Attribute_value","Score")]) 
plot3=plot2 %>% group_by(type2,type_phy) %>% transmute(n=n()) %>% unique() %>% ungroup() %>% 
  group_by(type2) %>% mutate(sum=sum(n)) %>% ungroup() %>% mutate(ratio=n/sum) %>% mutate(type="Aerophilicity") %>% 
  mutate(type_phy=gsub("facultatively anaerobic","Facultative anaerobe",type_phy)) %>% 
  mutate(type_phy=gsub("anaerobic","Anaerobe",type_phy)) %>% 
  mutate(type_phy=gsub("aerobic","Aerobe",type_phy))

ggplot(plot3,aes(x=type2,y=ratio,fill=type_phy))+
  geom_bar(stat = "identity",width = 0.8)+
  labs(title="",y="",x="")+
  facet_wrap(.~type)+
  scale_fill_manual(values = c("#80B1D3","#f4a39d", "#f96b6b","#B3B3B3"))+
  theme_bw(base_line_size = 0.2)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        strip.text.x.top = element_text(size=12,color = 'black'),
        axis.title.x=element_text(size=12,color = 'black'),
        axis.title.y=element_text(size=12,color = 'black'),
        axis.text.x=element_text(size=12,angle = 30,vjust = 0.98,hjust = 0.98,color = 'black'),
        axis.text.y=element_text(size=12,color = 'black'),
        legend.title=element_blank(),
        legend.text=element_text(size=12,color = 'black'),
        legend.position = "right")




##Table S1------------
#164 DE sig taxa oxygen table
plot33=plot2 %>% filter(!type2=="Unchanged") %>% 
  mutate(type_phy=gsub("facultatively anaerobic","Facultative anaerobe",type_phy)) %>% 
  mutate(type_phy=gsub("anaerobic","Anaerobe",type_phy)) %>% 
  mutate(type_phy=gsub("aerobic","Aerobe",type_phy)) %>% 
  select(c("taxa","type2","type_phy","Score")) %>% filter(!type2=="Unchanged") %>% 
  rename(Species=taxa,Group=type2,Aerophilicity=type_phy)
#write.csv(plot33,"Table S1.csv",row.names = F)












