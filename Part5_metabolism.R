---
title: "Metabolism analysis"
author: "Xuan Jiang"
note: "Figure6"
---
  
#metabolism
##6B PCA
##6C Bacteria derived Top metabolism
##6D Top metabolism
  
library(readxl)
library(tidyverse)
library(reshape2)
library(patchwork)
library(conflicted)
library(ggpubr)
library(Hmisc)


setwd("D:/Projects/Tumor/ESCC/data")



#metabolism------------
#138 metabolism,Be carful about the [*]
t1=read_excel("5_Metabolome/MWY-24-3452-a_2024-08-27-11-51-28/2.Basic_Analysis/Difference_analysis/Pm_vs_PBS/Pm_vs_PBS_info.xlsx") %>% 
  filter(`P-value`<0.05) %>% filter(!grepl("[*]",Index)) %>% filter(!Level == "3")


##6B PCA---------------
t1=read_excel("5_Metabolome/MWY-24-3452-a_2024-08-27-11-51-28/2.Basic_Analysis/Difference_analysis/Pm_vs_PBS/Pm_vs_PBS_info.xlsx") %>% 
  filter(`P-value`<0.05) %>% filter(!grepl("[*]",Index)) %>% filter(!Level == "3")
dat1=t1[,c("Index","Pm-1","Pm-2","Pm-3","PBS-1","PBS-2","PBS-3")] %>% 
  setNames(gsub("Pm","P. m.",colnames(.))) %>% column_to_rownames("Index")
meta1=data.frame(sample=colnames(dat1)) %>% 
  separate(sample,c("group","ID"),sep="-",remove = F)
rownames(meta1)=meta1$sample
pca_result <- prcomp(as.data.frame(t(dat1)), scale. = TRUE)
pca_scores1 <- as.data.frame(pca_result$x) %>% mutate(sample=rownames(.)) %>% left_join(.,meta1)
explained_variance1 <- summary(pca_result)$importance[2, ]

ggplot(pca_scores1, aes(x = PC1, y = PC2, label = sample)) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.4) +
  geom_point(size=3,aes(color=group)) +
  #geom_text_repel() + 
  scale_color_manual(values = c("#FF8080", "#A0A0A4","#FDB462","#377EB8"))+
  labs(title = "Tumor interstitial fluid",
       x = paste0("PC1 (", round(explained_variance1[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(explained_variance1[2] * 100, 2), "%)")) +
  theme_bw(base_line_size = 0.2)+
  theme(
    legend.text = element_text(color = 'black',size = 12, face = 'plain'),
    axis.text = element_text(color = 'black',size = 12,  face = 'plain'),
    axis.title = element_text(color = 'black',size = 12,  face = 'plain'),
    legend.position = "right") 



##6C Bacteria derived Top metabolism----------
tmp=read_excel("5_Metabolome/MWY-24-3452-a_2024-08-27-11-51-28/2.Basic_Analysis/Difference_analysis/Pm_vs_PBS/trace/Pm_vs_PBS_trace.xlsx") %>% 
  select(c("Index","Origin")) 
t1=read_excel("5_Metabolome/MWY-24-3452-a_2024-08-27-11-51-28/2.Basic_Analysis/Difference_analysis/Pm_vs_PBS/Pm_vs_PBS_info.xlsx") %>% 
  filter(`P-value`<0.05) %>% filter(!grepl("[*]",Index)) %>% filter(!Level == "3") %>% left_join(.,tmp)

plot1=t1 %>% filter(grepl("Bacteria",Origin)) %>%  filter(!grepl("Homo",Origin)) %>% 
  mutate(Log2FC_abs=abs(Log2FC)) %>% mutate(t=`物质`)
  top_n(10, Log2FC_abs) %>% .[order(.$Log2FC),] %>% 
  mutate(Compounds=factor(Compounds,levels=.$Compounds)) %>% 
  

p1=ggplot(plot1,aes(Log2FC,Compounds,color=Type)) + 
  geom_segment(aes(x = 0, xend = Log2FC, y = Compounds, yend = Compounds), 
               color = "grey") + 
  geom_point(aes(fill=Type,size=VIP)) +
  labs(x="",y="",title="Culture medium")+
  #labs(x="",y="",title="Tumor interstitial fluid")+
  scale_color_manual(values=c("#80B1D3", "#FB8072"))+
  theme_bw( base_line_size = 0.2)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,color = "black"),
        axis.title.y=element_text(size=12,color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "right")

p2=ggplot()+
  geom_tile(data = plot1, aes(x = 1,y=Compounds, fill = `Class I`), width = 1) +
  scale_fill_manual(values = c("Amino acid and Its metabolites"="#377EB8","Nucleotide and Its metabolites"="#66C2A5","Others"="#B3B3B3",
                               "Carbohydrates and Its metabolites"="#FDB462", "Organic acid and Its derivatives"="#FF7F00",
                               "Heterocyclic compounds"="#E78AC3","Benzene and substituted derivatives"="#4DAF4A",
                               "Aldehyde,Ketones,Esters"="#FFED6F"))+
  labs(x="",y="")+
  theme_bw( base_line_size = 0.2)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.98,hjust = 0.98,color = "black"),
        axis.text.y=element_blank(),
        legend.position = "right",
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"))

p3=p1+p2+plot_layout(widths = c(5,1))





##6D Top metabolism----------
t1=read_excel("5_Metabolome/MWY-24-3452-a_2024-08-27-11-51-28/2.Basic_Analysis/Difference_analysis/Pm_vs_PBS/Pm_vs_PBS_info.xlsx") %>% 
  filter(`P-value`<0.05) %>% filter(!grepl("[*]",Index)) %>% filter(!Level == "3")
plot1=t1 %>% mutate(Log2FC_abs=abs(Log2FC)) %>% 
  group_by(Type) %>% top_n(10, Log2FC_abs) %>% ungroup() %>% .[order(.$Log2FC),] %>% 
  mutate(Compounds=factor(Compounds,levels=.$Compounds))
  

p1=ggplot(plot1,aes(Log2FC,Compounds,color=Type)) + 
  geom_segment(aes(x = 0, xend = Log2FC, y = Compounds, yend = Compounds), 
               color = "grey") + 
  geom_point(aes(fill=Type,size=VIP)) +
  labs(x="",y="",title="Culture medium")+
  #labs(x="",y="",title="Tumor interstitial fluid")+
  scale_color_manual(values=c("#80B1D3", "#FB8072"))+
  theme_bw( base_line_size = 0.2)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.5,hjust = 0.5,color = "black"),
        axis.text.y=element_text(size=12,angle=0,color = "black"),
        axis.title.y=element_text(size=12,color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "right")

p2=ggplot()+
  geom_tile(data = plot1, aes(x = 1,y=Compounds, fill = `Class I`), width = 1) +
  scale_fill_manual(values = c("Amino acid and Its metabolites"="#377EB8","Nucleotide and Its metabolites"="#66C2A5","Others"="#B3B3B3",
                               "Carbohydrates and Its metabolites"="#FDB462", "Organic acid and Its derivatives"="#FF7F00",
                               "Heterocyclic compounds"="#E78AC3","Benzene and substituted derivatives"="#4DAF4A",
                               "Aldehyde,Ketones,Esters"="#FFED6F"))+
  labs(x="",y="")+
  theme_bw( base_line_size = 0.2)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=12,angle=0,vjust = 0.98,hjust = 0.98,color = "black"),
        axis.text.y=element_blank(),
        legend.position = "right",
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"))

p3=p1+p2+plot_layout(widths = c(5,1))















