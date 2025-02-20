---
title: "NT and survival analysis"
author: "Xuan Jiang"
note: "Figure2 & Supp.Figure2"
---

#1_taxa with NT & survival significance
##caculate Survival
##2A venn
##2A order
##2B plot
#2_P. m. correlation in saliva & tumor
##S2A corr
#3_P. m. correlation in NT
##S2B corr
#4_P. m. in spatial transcriptomic data
##S2C spatial
##S2D HE
#5_P. m. in CRC & GC 
##1.process data
##2F& 2G
  
library(tidyverse)
library(reshape2)
library(patchwork)


setwd("D:/Projects/Tumor/ESCC/data")




#1_taxa with NT & survival significance--------
##caculate Survival----------
library(survival)
library(survminer)
cacu_survival = function(abfvrel2,meta2){
  name1 = data.frame(taxa=rownames(abfvrel2),id=paste("Taxa_",seq(1:nrow(abfvrel2)),sep=""))
  rownames(name1) = name1$taxa
  name1 = name1[rownames(abfvrel2),]
  rownames(abfvrel2)=name1$id
  data2=t(abfvrel2) %>% as.data.frame() %>% rownames_to_column("sample") %>% 
    left_join(.,meta2[,c("sample","age","sex","Survival.time","type_survival")]) %>% column_to_rownames("sample")
  cutoff=surv_cutpoint(data2, 
                       time="Survival.time",
                       event="type_survival",
                       minprop = 0.25,
                       variables=rownames(abfvrel2))
  data_group=surv_categorize(cutoff) %>% rownames_to_column("sample") %>% 
    left_join(meta2[,c("sample","age","sex","stage2","Smoking.status","Drinking.status")],.) %>% 
    dplyr::select(c("Survival.time","type_survival","sample","age","sex","stage2","Smoking.status","Drinking.status"),everything()) %>% 
    column_to_rownames("sample")
  num_col=7

  for (i in (num_col+1):ncol(data_group)) {
    data_group[,i]=factor(data_group[,i],levels = c("low","high"))
  }
  #multiple Cox
  p_sum=data.frame(id="",p_value="",HR="",HR_low="",HR_high="",count_high="",count_low="")
  p_sum=p_sum[-1,]
  for (i in 1:(ncol(data_group)-num_col)) {
    j=colnames(data_group)[i+num_col]
    res.cox = coxph(Surv(Survival.time, type_survival) ~ age + sex + stage2 + Smoking.status + Drinking.status+as.factor(get(j)), data =  data_group)
    temp=broom::tidy(res.cox) 
    mul_cox1 = summary(res.cox)
    temp2=data_group[,j,drop=F] %>% group_by(get(j)) %>% mutate(check=n()) %>% unique()
    p_sum[i,1]=j
    p_sum[i,2]= temp[num_col-1,5] 
    p_sum[i,3]= mul_cox1$conf.int[num_col-1,1] 
    p_sum[i,4:5]= mul_cox1$conf.int[num_col-1,3:4] 
    p_sum[i,6]= temp2[1,3]
    p_sum[i,7]= temp2[2,3]
  }
  p_sum[,2:7]=apply(p_sum[,2:7],2,as.numeric)
  result_multi_Cox = p_sum %>% mutate(p_adj=p.adjust(p_value,method="fdr")) %>% left_join(.,name1) %>% 
    .[order(.$p_value),] %>% rename(p_multip=p_value) 
  
  id1=c("data_group","result_multi_Cox")
  res_list=list(data_group,result_multi_Cox)
  for (i in 1:length(res_list)) {names(res_list)[i]=id1[i]}
  return(res_list)
}


#data1.1 taxa
load("table/use_taxa_bracken_data.RData")
tmp=read.csv("table/metadata.csv") 
meta2 <- meta %>% filter(!is.na(Survival.time))  %>% left_join(.,tmp) %>% filter(sample_type=="T")
abfv2 <- abfv[,meta2$sample] %>% filter(rowSums(.)>0) 
abfvrel2 <- abfvrel[rownames(abfv2),meta2$sample] %>% filter(rowSums(.>0)>(0.3*ncol(.))) 
res_list <- cacu_survival(abfvrel2 = abfvrel2,meta2 = meta2) 
save(res_list,file = "res.survival.RData")



##2A venn------
library(ggvenn)
load("res.DE.RData")
#NT
plot1=res2 %>% left_join(.,res1) %>% 
  mutate(type2=ifelse(p_value > 0.05,"Unchanged",ifelse(diff.btw > 0,"ESCC > NOR","ESCC < NOR"))) 
t1=plot1 %>% filter(type2=="ESCC > NOR") %>% mutate(prev=prev_T) 
t2=plot1 %>% filter(type2=="ESCC < NOR") %>% mutate(prev=prev_N) 
#survival
load("res.survival.RData")
plot2=res_list$result_multi_Cox %>% rename(p_value=p_multip) %>% left_join(tmp2,.) %>% na.omit() %>% 
  mutate(type2=ifelse(p_value > 0.05,"n.s.",ifelse(HR > 1,"HR > 1","HR < 1"))) 
t3=plot2 %>% filter(type2=="HR > 1") 
t4=plot2 %>% filter(type2=="HR < 1") 
#prev
t5=plot1 %>% filter(prev_T>0.6)
t6=plot1 %>% filter(prev_N>0.6)


#data1-bad
col1=t3$taxa
col2=t1$taxa
col3=t5$taxa
#data2-good
col1=t4$taxa
col2=t2$taxa
col3=t6$taxa
#

plot3 <- tibble(  
  survival = c(col1, rep(NA, max(length(col1), length(col2), length(col3)) - length(col1))),  
  NT = c(col2, rep(NA, max(length(col1), length(col2), length(col3)) - length(col2))),  
  prev = c(col3, rep(NA, max(length(col1), length(col2), length(col3)) - length(col3))))  
venn_list <- as.list(plot3)              
venn_list <- purrr::map(venn_list, na.omit)
venn_list <- purrr::map(venn_list, function(x){x[x!=""]})

ggvenn(
  data = venn_list,        
  columns = NULL,          
  show_elements = F,       
  label_sep = "\n",        
  show_percentage = F,     
  digits = 1,              
  fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"), 
  fill_alpha = 0.5,         
  stroke_color = "white",   
  stroke_alpha = 0.5,       
  stroke_size = 0.5,        
  stroke_linetype = "solid",
  set_name_color = "black", 
  set_name_size = 6,        
  text_color = "black",     
  text_size = 4            
)
#


##2A order--------
load("res.NT.RData")
plot1=res2 %>% left_join(.,res1) %>% 
  mutate(type2=ifelse(p_value > 0.05,"Unchanged",ifelse(diff.btw > 0,"ESCC > NOR","ESCC < NOR"))) 
t1=plot1 %>% filter(type2=="ESCC > NOR") 
t2=plot1 %>% filter(type2=="ESCC < NOR") 
load("res.survival.RData")
plot2=result_multi_Cox %>% rename(p_value=p_multip) %>% left_join(tmp2,.) %>% na.omit() %>% 
  mutate(type2=ifelse(p_value > 0.05,"n.s.",ifelse(HR > 1,"HR > 1","HR < 1")))
t3=plot2 %>% filter(type2=="HR > 1") %>% filter(taxa %in% t1$taxa) %>% 
  filter(prev_T >0.6) %>% mutate(median=median_T)
t4=plot2 %>% filter(type2=="HR < 1") %>% filter(taxa %in% t2$taxa) %>%  
  filter(prev_N >0.6) %>% mutate(median=median_N)
plot3=rbind(t3,t4) %>% .[order(-.$median),] %>% 
  mutate(taxa=factor(taxa,levels=taxa))
save(plot3,"Micro_overlap.Rdata")
color1=c("#80B1D3", "#FB8072")

ggplot(plot3,aes(100*median,fct_rev(taxa),fill=type2)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=color1)+
  labs(x="Relative abundance(%)",y="Species")+
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
        legend.position = "none")


##2B plot------------
load("res.NT.RData")
plot1=res2 %>% left_join(.,res1) %>% 
  mutate(type2=ifelse(p_value > 0.05,"Unchanged",ifelse(diff.btw > 0,"ESCC > NOR","ESCC < NOR"))) 
t1=plot1 %>% filter(type2=="ESCC > NOR") 
t2=plot1 %>% filter(type2=="ESCC < NOR") 
load("res.survival.RData")
plot2=result_multi_Cox %>% rename(p_value=p_multip) %>% left_join(tmp2,.) %>% na.omit() %>% 
  mutate(type2=ifelse(p_value > 0.05,"n.s.",ifelse(HR > 1,"HR > 1","HR < 1"))) %>% rename(p_value_HR=p_value)
t3=plot2 %>% filter(type2=="HR > 1") %>% filter(taxa %in% t1$taxa) %>% 
  filter(prev_T>0.6) %>% mutate(median=median_T)
t4=plot2 %>% filter(type2=="HR < 1") %>% filter(taxa %in% t2$taxa) %>%  
  filter(prev_N>0.6) %>% mutate(median=median_N)
plot3=rbind(t3,t4) %>% left_join(.,plot1[,c("taxa","diff.btw")]) %>% 
  .[order(-.$median),] %>% mutate(taxa=factor(taxa,levels=taxa))
#
plot4=plot2 %>% select(all_of(c("taxa","HR","p_value_HR"))) %>%
  left_join(.,plot1) %>% na.omit()
t5=plot3 %>% filter(type2=="HR > 1") %>% .[order(-.$median_T),] %>% .[1:5,] %>% 
  mutate(taxa=factor(taxa,levels=taxa))
t6=plot3 %>% filter(type2=="HR < 1") %>% .[order(-.$median_N),] %>% .[1:5,] %>% 
  mutate(taxa=factor(taxa,levels=taxa))
tmp=rbind(t5,t6)
#
#oxygen
library(bugphyzz)
dat_bugphyzz=bugphyzz::importBugphyzz()
tmp1=dat_bugphyzz$aerophilicity %>% mutate(taxa=Taxon_name,genus=Taxon_name) %>% filter(Score > 0.5)
tmp2=plot1[,"taxa",drop=F] %>% left_join(.,tmp1[,c("taxa", "Attribute","Attribute_value")]) %>%
  filter(!is.na(Attribute)) %>% group_by(taxa) %>% mutate(n=n()) %>% 
  summarise(type_phy = paste(Attribute_value, collapse = ",")) %>% ungroup()
tmp3=plot1[,"taxa",drop=F] %>% filter(!taxa %in% tmp2$taxa) %>% separate(taxa,c("genus"),sep=" ",remove = F) %>% 
  left_join(.,tmp1[,c("genus", "Attribute","Attribute_value")]) %>%
  filter(!is.na(Attribute)) %>% group_by(taxa) %>% mutate(n=n()) %>% 
  summarise(type_phy = paste(Attribute_value, collapse = ",")) %>% ungroup() %>% rbind(.,tmp2) %>% 
  mutate(type_phy=ifelse(type_phy=="aerobic","Aerobic",
                         ifelse(type_phy=="unknown","Unknown",
                                ifelse(type_phy=="anaerobic","Anaerobic","Facultatively anaerobic"))))
plot5=plot4 %>% left_join(.,tmp3) 
plot5$type_phy[is.na(plot5$type_phy)]="unknown"
plot3=plot3 %>% left_join(.,tmp3)
#
library(ggrepel)
color1=c("#80B1D3", "#FB8072")
plot6=plot5 %>% filter(p_value<0.05)
ggplot(plot6,aes(log2(HR),diff.btw)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.3,linetype = 'dashed') + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.3,linetype = 'dashed') + 
  geom_point(aes(size=median_T,shape=type_phy),color = "grey")+
  geom_point(data=plot3,aes(log2(HR),diff.btw,size=median_T,color = type2,shape=type_phy),alpha=0.5)+
  geom_text_repel(data=tmp,aes(log2(HR),diff.btw,label=taxa))+
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



#2_P. m. correlation in saliva & tumor-----------
##S2A corr------------
load("Micro_overlap.Rdata")
plot1=readRDS("2_caculate_result/Saliva_cor/result.RDS") %>% 
  rename(taxa=spec) %>% left_join(plot3[,c("taxa","change")]) %>% 
  .[order(-(.$correlation)),] %>% 
  mutate(taxa=factor(taxa,levels=.$taxa))
tmp=plot1 %>% filter(!significance=="ns")

color1=c("#80B1D3", "#FB8072")
ggplot(plot1,aes(correlation,fct_rev(taxa))) +
  geom_bar(stat = "identity",width = 0.05)+
  geom_point(data=plot1,aes(correlation,fct_rev(taxa),size=-log10(p_value),color=change))+
  labs(x="Spearman correlation",y="")+
  scale_color_manual(values = color1)+
  theme_bw()+
  geom_text(data=tmp,aes(x=correlation,y=taxa,label=significance),size=3)+
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



#3_P. m. correlation in NT----------------
##S2B corr-----------
#Spearman'r between tumor and normal
load("Micro_overlap.Rdata")
t1=plot3 %>% filter(!grepl("sp[.]",taxa)) %>% 
  filter(!taxa %in% c("Streptococcus equi","Leptotrichia trevisanii","Prevotella scopos")) %>% 
  .[order(-.$median_T),]
load("table/use_taxa_bracken_rmMus_5wan_210samples_data.RData")
id1=meta %>% filter(sample_type=="N") %>%  .[order(.$ID),] %>% .$sample
id2=meta %>% filter(sample_type=="T") %>%  .[order(.$ID),] %>% .$sample
plot1=c()
for (i in 1:nrow(abfvrel)) {
  tmp1=abfvrel[i,id1] %>% as.numeric()
  tmp2=abfvrel[i,id2] %>% as.numeric()
  tmp3=cor.test(tmp1,tmp2,method = "spearman",adjust = "fdr")
  tmp4=data.frame(taxa=rownames(abfvrel)[i],r=as.numeric(tmp3$estimate),p=as.numeric(tmp3$p.value))
  plot1=rbind(plot1,tmp4)
}
plot4=res2 %>% left_join(.,res1) %>% 
  mutate(type2=ifelse(p_value > 0.05,"Unchanged",ifelse(diff.btw > 0,"ESCC > NOR","ESCC < NOR"))) %>% 
  filter(type2=="ESCC > NOR") %>% mutate(prev=prev_T) 
#
plot2=plot1 %>% na.omit() %>% mutate(p.adj = p.adjust(p, method = "fdr")) %>% 
  filter(p.adj<0.05) %>% 
  mutate(Type=ifelse(taxa %in% plot3$taxa,"Pro-tumor",
                     ifelse(taxa %in% micro_NG_0.5$taxa,"Anti-tumor","None"))) %>% 
  mutate(size1=ifelse(taxa %in% plot3$taxa,1.5,
                      ifelse(taxa %in% micro_NG_0.5$taxa,1.5,1))) %>% 
  left_join(.,plot4) %>% filter(type=="ALDEx2") %>% 
  mutate(Type1=ifelse(diff.btw>0,"Tumor enriched","Normal enriched"))
plot3=plot2 %>% filter(taxa %in% plot3$taxa) %>% .[order(-.$r),] %>% .[1:3,]
plot4=plot2 %>% filter(taxa %in% micro_NG_TB_0.5$taxa) %>% .[order(-.$r),] %>% 
  mutate(taxa=factor(taxa,levels = taxa))
ggplot(plot4,aes(r,fct_rev(taxa))) +
  geom_bar(stat = "identity",width = 0.05)+
  geom_point(data=plot4,aes(r,fct_rev(taxa),size=-log10(p.adj),color=Type))+
  labs(x="Spearman correlation",y="")+
  scale_color_manual(values = c("#80B1D3", "#FB8072"))+
  theme_bw()+
  #geom_text(data=plot4,aes(x=r,y=taxa,label=significance),size=3)+
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




#4_P. m. in spatial transcriptomic data------------
##S2C spatial----------
library(Seurat)
library(SeuratData)
library(data.table)
id="P32"
#
tissue_sample=readRDS(paste("ST_ESCC/",id,".filtered_SCTransform.rds",sep=""))
umi_table00=read.csv(paste("ST_ESCC/umi_table_",id,".csv",sep="")) %>% 
  filter(rank=="species") %>% mutate(genus=taxa) #species!
umi_table2=umi_table00 %>% 
  group_by(barcode,UMI,genus) %>% transmute(count=n()) %>% ungroup() %>% unique() %>% 
  group_by(barcode,genus) %>% transmute(count=n()) %>% ungroup() %>% unique() %>% 
  spread(.,key=genus,value=count) %>% column_to_rownames("barcode") 
umi_table2[is.na(umi_table2)]=0
umi_table=umi_table2
umi_table$Total <- rowSums(umi_table)
umi_table[umi_table==0] <- NA
tissue_sample2 <- AddMetaData(tissue_sample, umi_table)
#
SpatialFeaturePlot(tissue_sample2,features = c("Parvimonas micra"),pt.size=3) +
  ggtitle("Parvimonas micra nUMI") + 
  theme()+
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10))


##S2D HE---------
id="P32"
plot1=read.csv(paste0("ST_ESCC/HE_type_UMI_1NOR_",id,".csv")) %>% 
  filter(genus=="Parvimonas micra")
ggplot(plot1, aes(x = HE_type,y=count,fill = genus, 
                  stratum = genus, alluvium = genus)) +
  geom_col(width = 0.4,color=NA)+
  geom_flow(width = 0.4,alpha = 0.4,knot.pos = 0)   +
  scale_fill_manual(values="#EEC373")+
  labs(title="Parvimonas micra",y="UMI/spots",x="")+
  theme_bw()+
  theme(
    panel.border = element_blank(),  
    axis.line = element_line(colour = "black", linewidth = 0.5), 
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    panel.background = element_blank(),        
    plot.background = element_blank(),         
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(size = 11),  axis.title.y = element_text(size = 12),  
    axis.text.x=element_text(size=12,angle=0,color="black"),
    legend.position = "none")



#5_P. m. in CRC & GC -----------
##1.process data------------
#The raw sequences (PRJNA280026, PRJNA861885, and PRJNA383606 from CRC; PRJNA1032279 from GC) were downloaded and processed using the qiime2 pipeline. 
#ASVs classified under the genus Parvimonas were extracted and re-aligned to the NCBI 16S rRNA database using MEGABLAST.
#Only ASVs identified as P. micra with the highest match scores were retained for further analysis.

#Samples of data process
type="CRC"
id="PRJNA830432"
meta=read.csv(paste0("table/16S/",type,"/",id,"/SraRunTable.csv")) %>% 
  filter(isolation_source=="gut") %>% 
  mutate(Sample.Name=gsub("XEL","",Sample.Name)) %>% 
  mutate(sample_type1 = sapply(str_extract_all(Sample.Name, "[a-zA-Z]+"), paste, collapse = "")) %>% 
  mutate(sample_type = gsub("LA","",sample_type1)) %>% mutate(sample_type = gsub("RA","",sample_type)) %>% 
  mutate(sample_type2=substr(sample_type1,1,1)) %>% 
  mutate(ID = sapply(str_extract_all(Sample.Name, "\\d+"), paste, collapse = "")) %>% 
  unique()
abfv0=fread(paste0("table/16S/",type,"/",id,"/table.merge.txt")) %>% 
  select("taxonomy",everything()) %>% 
  filter(grepl("d__Bacteria;",taxonomy)) %>% filter(grepl("g__",taxonomy)) 
abfv1=abfv0 %>% filter(!`#OTU ID` %in% otu_id_g_remove$OTU_ID) %>% 
  group_by(taxonomy) %>%  summarise(across(starts_with("SRR"), sum)) %>% 
  column_to_rownames("taxonomy") %>% select(meta$Run)
abfv2=abfv0 %>% filter(`#OTU ID` %in% otu_id_g_remove$OTU_ID) %>% 
  group_by(taxonomy) %>%  summarise(across(starts_with("SRR"), sum)) %>% 
  column_to_rownames("taxonomy") %>% select(meta$Run)
rownames(abfv2)="g__Parvimonas_notPm"
abfv=rbind(abfv1,abfv2)
tmp=apply(abfv,1,function(x) sum(x > 0)) %>% as.data.frame() %>% setNames("prevalence") %>% rownames_to_column("taxonomy") %>% 
  mutate(prev=prevalence/ncol(abfv)) %>% filter(prev>0.1) %>% 
  mutate(taxa=sub(".*g__", "g__", taxonomy))
abfv=abfv[tmp$taxonomy,] %>% rownames_to_column("taxa")
abfv$taxa=tmp$taxa
abfv=abfv %>% group_by(taxa) %>%  summarise(across(starts_with("SRR"), sum)) %>% column_to_rownames("taxa")
abfvrel=apply(abfv, 2, function(x){x/sum(x)}) %>% as.data.frame() 
plot1=melt(as.matrix(abfvrel)) %>% rename(taxa=Var1,Run=Var2) %>% 
  filter(taxa=="g__Parvimonas") %>% 
  left_join(.,meta[,c("Run","sample_type","ID","sample_type2")]) %>% 
  unique() %>% mutate(taxa=as.character(taxa)) 

tmp=abfv %>% transmute(value=rowSums(.)/ncol(.)) %>% .[order(-.$value),,drop=F] %>% rownames_to_column("taxa")
tmp %>% filter(grepl("g__Parvimonas",taxa)) %>% .$taxa
unique(plot1$sample_type)
unique(plot1$taxa)
##p.value
rownames(plot1)=plot1$Run
plot1=plot1 %>% mutate(ID2=paste0(sample_type2,ID))
meta=meta %>% mutate(ID2=paste0(sample_type2,ID))
sort(unique(plot1$sample_type))
tmp1=meta %>% filter(sample_type==sort(unique(plot1$sample_type))[1]) %>% .[order(.$ID2),]
tmp2=meta %>% filter(sample_type==sort(unique(plot1$sample_type))[2]) %>% .[order(.$ID2),] 
wilcox.test(plot1[tmp1$Run,"value"],plot1[tmp2$Run,"value"],paired = T)


##2F& 2G-----
ggplot(data=plot1,aes(x=sample_type, y=(value)))+
  geom_point(shape=16)+
  geom_boxplot(alpha=0)+
  labs(title=id,y="Relative abundance",x="")+
  stat_compare_means(label = "p.signif", paired = T,
                     comparisons=list(c("N", "P"),
                                      c("N", "T"),
                                      c("P", "T")))+
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










