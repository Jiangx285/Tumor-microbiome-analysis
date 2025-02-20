---
title: "taxa correlation with single cell data"
author: "Xuan Jiang"
note: "Figure3"
---
  
#1_Single Cell signature
##3A NT
###process
###plot
##3B survival
#2_correlation
##process
##3C corr  
  
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggpubr)
library(Hmisc)


setwd("D:/Projects/Tumor/ESCC/data")




#1_Single Cell signature-----------
##3A NT-------------
###process-----------
sing=read.csv("cell.score.csv")
meta=read.csv("meta_cell.csv") %>% filter(sample %in% colnames(sing))
plot1=melt(as.matrix(sing),varnames = c("sing","sample")) %>% as.data.frame() %>% left_join(.,meta[,c("sample","sample_type")])
tmp1=plot1 %>% group_by(sing,sample_type) %>% mutate(median=median(value)) %>% 
  select(c("sing","sample_type","median")) %>% unique() %>% 
  spread(.,key=sample_type,value=median) %>% mutate(T_N=T-N)
#
plot2 <- c()
id1=meta %>% filter(sample_type=="T") %>% .$sample
id2=meta %>% filter(sample_type=="N") %>% .$sample
for (i in rownames(sing)) {
  group_1 = plot1 %>% filter(sing==i) %>% filter(sample %in% id1) %>%  .$value
  group_2 = plot1 %>% filter(sing==i) %>% filter(sample %in% id2) %>%  .$value
  result_wilcoxon = wilcox.test(group_1, group_2, paried=F)
  p <- data.frame(
    sing=i,
    p_value=result_wilcoxon$p.value)
  plot2 <- rbind(plot2,p)
}
plot2 <- plot2 %>% mutate(p_adjust = p.adjust(p_value, method = "fdr")) 
#write.csv(plot2,"cell_NT_diff.csv",row.names = F)
#write.csv(plot1,"cell_score_merge.csv",row.names = F)


###plot---------------
meta=read.csv("meta_cell.csv") 
dat1=read.csv("cell_score_merge.csv")
plot2=c()
for (i in unique(dat1$sing)) {
  tmp1=dat1 %>% filter(sing==i,sample_type=="T")
  tmp2=dat1 %>% filter(sing==i,sample_type=="N")
  tmp=data.frame(sing=i,fd=mean(tmp1$value)/mean(tmp2$value))
  plot2=rbind(plot2,tmp)
}
plot1=read.csv("cell_NT_diff.csv") %>% 
  filter(p_value<0.05) %>% .[order(.$p_value),] %>% left_join(plot2) %>% 
  mutate(type=ifelse(fd > 1,"up","down")) %>% 
  #.[order(.$fd),] %>% 
  .[order(.$type,-.$p_value),] %>% 
  mutate(sing=factor(sing,levels=sing)) %>% 
  mutate(label=ifelse(p_value <= 0.001,"***",
                    ifelse(p_value > 0.001& p_value <= 0.01,"**",
                           ifelse(p_value > 0.01 & p_value < 0.05,"*","n.s.")))) 

color1=c("#80B1D3", "#FB8072")
ggplot(plot1,aes(log2(fd),sing)) +
  geom_bar(stat = "identity",width = 0.05)+
  geom_point(data=plot1,aes(log2(fd),sing,color=type,size=-log10(p_value)))+
  labs(x="log2(Fold change)",y="")+
  scale_color_manual(values = color1)+
  theme_bw()+
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



##3B survival-------------
library(survival)
library(survminer)
library(gridExtra)
#data1.taxa
fit2 <- survfit(Surv(Survival.time, type_survival) ~ taxa, data = sing_taxa)
p1=ggsurvplot(fit2,title=k,
                 pval = TRUE, #conf.int = TRUE,
                 #risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 #linetype = "strata", # Change line type by groups
                 #surv.median.line = "hv", # Specify median survival
                 xlab = "Time (Months)",
                 ylab = "",
                 legend.title = "",
                 legend.labs = labels,
                 ggtheme = theme_survminer(), # Change ggplot2 theme
                 palette = c("#00468B","#ED0000"))
write.csv("cell_survival.csv")



#2_correlation-----------
##process---------------
load("table/use_taxa_bracken_data.RData")
sing=read.csv("cell.score.csv") %>% select(ends_with("_T"))
id=intersect(colnames(abfvrel),colnames(sing))
abfvrel2=abfvrel[,id] %>% filter(rowSums(.)>0)
tmp=cbind(t(abfvrel2),t(sing[,id]))
result1=rcorr(as.matrix(tmp),type = "spearman")
result1$P.adj=matrix(p.adjust(result1$P, method = 'fdr'),nr=ncol(tmp))
colnames(result1$P.adj)=colnames(result1$P)
rownames(result1$P.adj)=rownames(result1$P)
p_sum=result1$P[rownames(abfvrel2),rownames(sing)] 
heatdata=result1$r[rownames(p_sum),colnames(p_sum)]
tmp1=melt(as.matrix(p_sum),varnames = c("taxa","sing")) %>% as.data.frame() %>% rename(p=value)
tmp2=melt(as.matrix(heatdata),varnames = c("taxa","sing")) %>% as.data.frame() %>% rename(r=value) %>% 
  left_join(.,tmp1) %>% mutate(type=id1)
plot1=rbind(plot1,tmp2)



##3C corr----------------
#select single-cell p.NT < 0.5
t1=read.csv("cell_survival.csv") %>% filter(p_KM<0.05) %>% 
tmp1=read.csv("cell_NT_diff.csv") %>% 
  filter(p_value<0.05) %>% .[order(.$p_value),] %>% left_join(.,t1)%>% 
  select(!"type") %>% left_join(.,unique(t0[,c("sing","type")]))
tmp2=tmp1 %>% filter(T_N < 0) %>% .[order(.$type),] %>% .[order(.$sing),]
t2=tmp1 %>% filter(T_N > 0)  %>% .[order(.$type),] %>% .[order(.$sing),] %>% rbind(.,tmp2) %>% 
  mutate(NT=ifelse(T_N>0,"ESCC > NOR","ESCC < NOR")) %>% 
  mutate(sur=ifelse(is.na(HR),"n.s.",ifelse(HR>1,"Risk","Protect")))
t20=data.frame(type=unique(t2$type),
               Cell_type=c("#4E79A7","#999999","#4DAF4A","#BC80BD","#FF7F00","#d40012"))
t21=data.frame(NT=unique(t2$NT),
               Enrichment=c("#FB8072","#80B1D3"))
t22=data.frame(sur=unique(t2$sur),
               Survival=c("#FFBE7D","grey","#86BCB6"))
t2=t2 %>% left_join(t20) %>% left_join(t21) %>% left_join(t22)
#select 17 TPPMs and 8 NFPMs
load("Micro_overlap.Rdata")
tmp1=plot3 %>% filter(change=="Tumor") %>% mutate(median=median_T) %>% .[order(-.$median),] 
tmp2=plot3 %>% filter(change=="Normal") %>% mutate(median=median_N) %>% .[order(-.$median),]
tmp3=tmp2 %>% rbind(.,tmp1)%>% .[order(-.$median),] #%>% .[1:10,]
tmp11= tmp1 %>% filter(taxa %in% tmp3$taxa)
t3=tmp2 %>% filter(taxa %in% tmp3$taxa) %>% rbind(tmp11,.) %>% 
  mutate(NT=ifelse(change=="Tumor","ESCC > NOR","ESCC < NOR")) %>% 
  mutate(sur=ifelse(is.na(HR),"n.s.",ifelse(HR>1,"Risk","Protect"))) %>% 
  left_join(t21) %>% left_join(t22)
#select cell & taxa
plot1=plot1
  filter(sing %in% t2$sing) %>% filter(taxa %in% t3$taxa) %>% 
  mutate(taxa=factor(taxa,levels=as.character(t3$taxa))) %>% 
  mutate(sing=factor(sing,levels=t2$sing))
plot11=plot1 %>% filter(p<0.05)
plot2=t2 %>% mutate(sing=factor(sing,levels=.$sing))
plot3=t3 %>% mutate(taxa=factor(taxa,levels=.$taxa))
#
p1=ggplot() + 
  theme_classic()+
  geom_point(data=plot1,aes(x = (sing),y=fct_rev(taxa),size=-log10(p),fill=r),shape=21, color=rgb(0,0,0,0)) +
  geom_point(data=plot11,aes(x = (sing),y=fct_rev(taxa),size=-log10(p), fill=r),shape = 21, color="black",stroke = 0.2) +
  scale_fill_gradientn(colors = c("#0A94DE","white","#D01515"))+
  labs(x="",y="")+
  theme_bw( base_line_size = 0.2)+
  theme(#panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        plot.caption = element_text(color = "black"),
        axis.text.x=element_text(size=12,angle=45,vjust = 0.98,hjust = 0.98,color = "black"),
        axis.text.y=element_text(size=12,angle=0,color = "black"),
        axis.title.y=element_text(size=12,color = "black"),
        strip.text = element_text(size = 12),
        legend.position = "right")

p2=ggplot()+
  geom_tile(data = plot2, aes(x = factor(sing),y=3, fill = Enrichment), width = 1) +
  geom_tile(data = plot2, aes(x = factor(sing),y=2, fill = Survival), width = 1) +
  geom_tile(data = plot2, aes(x = factor(sing),y=1, fill = Cell_type), width = 1) +
  scale_fill_identity() +
  labs(x="",y="")+
  theme_bw( base_line_size = 0.2)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    axis.text.x=element_text(size=12,angle=45,vjust = 0.98,hjust = 0.98,color = "black"),
    legend.position = "none")

p3=ggplot()+
  geom_tile(data = plot3, aes(x = 2,y = fct_rev(taxa), fill = Enrichment), width = 1) +
  geom_tile(data = plot3, aes(x = 1,y = fct_rev(taxa), fill = Survival), width = 1) +
  scale_fill_identity() +
  labs(x="",y="")+
  theme_bw( base_line_size = 0.2)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none")
p4=p3+p1+plot_layout(widths = c(1,10))
p5=p4/p2+plot_layout(heights = c(5,1))



