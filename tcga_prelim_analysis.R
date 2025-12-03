setwd("C:/Users/sweta/OneDrive - Emory University/NSF REU UCR/Project")
library(dplyr)
library(tidyr)

#Import gene expression data
gene_exp <- read.delim("C:/Users/sweta/OneDrive - Emory University/NSF REU UCR/Project/TCGA.STAD.sampleMap_HiSeqV2_PANCAN.tsv", header=FALSE)
gene_exp<-t(gene_exp)
colnames(gene_exp)<-gene_exp[1,]
sample <- rownames(gene_exp)
gene_exp<-gene_exp[-1,]
rownames(gene_exp)<-NULL
gene_exp<-as.data.frame(gene_exp)

#Import phenotype data
library(readr)
phen <-read.delim("C:/Users/sweta/OneDrive - Emory University/NSF REU UCR/Project/TCGA.STAD")
colnames(phen)[1]<-"sample"
phen<-as.data.frame(phen)

#Import survival data
surv<-read.csv("C:/Users/sweta/OneDrive - Emory University/NSF REU UCR/Project/Survival Data - Stomach Cancer.csv")

#Clean datasets

#find samples present in all datasets
common<-Reduce(intersect,list(surv[,1],phen[,1],gene_exp[,1]))

gene_exp<-
  gene_exp%>%
  filter(sample %in% common)

names(phen) <- make.unique(names(phen))
phen<-
  phen%>%
  filter(sample %in% common)

surv<-
  surv%>%
  filter(sample %in% common)

#combine datasets
surv_phen<-phen%>%
  merge(surv,by="sample",all=T)

data<-surv_phen%>%
  merge(gene_exp,by="sample",all=T)

#Exploratory Data Analysis
library(ggplot2)

#Samples in the dataset
length(data$sample) #450 samples

#Patients in the dataset
length(unique(data$patient_id)) #418 patients, 32 individuals with more than 1 sample taken 

#Number of samples per sample type
table(data$sample_type)
ggplot(data,aes(x=sample_type))+
  geom_bar(fill=c("deeppink2","skyblue3"))+
  labs(x="Sample Type",y="Number of Samples",title="Number of Samples per Sample Type")
ggsave("Stomach - Number of Samples per Sample Type.png",width=5,height=7)

#Age Distribution of all Individuals by Sample Type
data%>%
  ggplot(aes(x=sample_type,y=age_at_initial_pathologic_diagnosis,fill=sample_type))+
  geom_boxplot(fill=c("deeppink2","skyblue3"))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Age Distribution of Individuals by Sample Type") +
  labs(x="Sample Type",y="Age at Initial Pathologic Diagnosis")
ggsave("Stomach - Age Distribution by Sample Type.png",width=5,height=7)

#DSS by age - by pathologic stage and sample type - ALL STAGES - info unavailable for 9
data2<-filter(data,pathologic_stage!="[Discrepancy]")

ggplot(data=subset(data2,!is.na(pathologic_stage)),aes(y=DSS.time,x=age_at_initial_pathologic_diagnosis,color=pathologic_stage, shape = sample_type))+
  geom_point(size=6)+
  scale_color_brewer(palette="Set3")+
  labs(x="Age at Pathologic Diagnosis (years)",y="Disease Specific Survival Time (days)",title="Age at Diagnosis vs. Time until Death Caused by Disease",color="Pathologic Stage",shape="Sample Type")
ggsave("Stomach - DSS by age - by pathologic stage and sample type - ALL STAGES.png",width=9,height=7)

#DSS by age - by pathologic stage and sample type - STAGES COMBINED - info unavailable for 9
library(forcats)
data2$pathologic_stage<-fct_collapse(data2$pathologic_stage, "Stage I"=c("Stage I","Stage IA","Stage IB"),"Stage II"=c("Stage II","Stage IIA","Stage IIB"),"Stage III"=c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC"))

ggplot(data=subset(data2,!is.na(pathologic_stage)),aes(y=DSS.time,x=age_at_initial_pathologic_diagnosis,color=pathologic_stage, shape = sample_type))+
  geom_point(size=6)+
  scale_color_brewer(palette="Dark2")+
  labs(x="Age at Pathologic Diagnosis (years)",y="Disease Specific Survival Time (days)",title="Age at Diagnosis vs. Time until Death Caused by Disease",color="Pathologic Stage",shape="Sample Type")+
  scale_x_continuous(breaks=seq(30,90,by=10))+
  theme(legend.position = "bottom")
ggsave("Stomach - DSS by age - by pathologic stage and sample type - STAGES COMBINED.png",width=9,height=7)

#Disease Specific Survival by Pathologic Stage
data2%>%
  ggplot(aes(x=pathologic_stage,y=DSS.time,fill=pathologic_stage))+
  geom_boxplot()+
  scale_fill_brewer(palette="Set2")+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Distribution of Disease Specific Survival by Pathologic Stage") +
  labs(x="Pathologic Stage",y="Disease Specific Survival (days)")


#DFI by age - by pathologic stage
ggplot(filter(data2,pathologic_stage!="Stage IV"),aes(y=DFI.time,x=age_at_initial_pathologic_diagnosis,color=pathologic_stage))+
  geom_point(size=6)+
  scale_color_brewer(palette="Dark2")+
  labs(x="Age at Pathologic Diagnosis (years)",y="Disease Free Interval (days)",title="Age at Diagnosis vs. Disease Free Interval",color="Pathologic Stage")

#OS by age - by pathologic stage - by alive/deceased status
ggplot(data2,aes(y=OS.time,x=age_at_initial_pathologic_diagnosis,color=pathologic_stage,shape=vital_status))+
  geom_point(size=6)+
  scale_color_brewer(palette="Dark2")+
  scale_shape_discrete(labels=c("Alive","Deceased"))+
  labs(x="Age at Pathologic Diagnosis (years)",y="Days of Survival",title="Age at Diagnosis vs. Overall Survival",color="Pathologic Stage", shape="Alive/Deceased")+
  theme(legend.key.size=unit(1.0,'cm'))
ggsave("Stomach - OS by age - by pathologic stage and vital status",width=9,height=7)

## DSS by Age - by Pathologic Stage and Survival Status
ggplot(data2,aes(y=DSS.time,x=age_at_initial_pathologic_diagnosis,color=pathologic_stage,shape=as.character(OS)))+
  geom_point(size=6)+
  scale_color_brewer(palette="Dark2")+
  scale_shape_discrete(labels=c("Alive","Deceased"))+
  labs(x="Age at Pathologic Diagnosis (years)",y="Disease Specific Survival Time (days)",title="Age at Diagnosis vs. Overall Survival",color="Pathologic Stage", shape="Survival Status")+
  theme(legend.key.size=unit(1.0,'cm'))

## DSS by Histological Type
data3<-data%>%
  filter(histological_type!="[Discrepancy]")%>%
  drop_na(histological_type)

data3$histological_type<-as.factor(data3$histological_type)

levels(data3$histological_type)<-c("Adenocarcinoma, Signet Ring Type","Adenocarcinoma, Diffuse Type","Adenocarcinoma, Not Otherwise Specified (NOS)","Intestinal Adenocarcinoma, Mucinous Type","Intestinal Adenocarcinoma, Not Otherwise Specified (NOS)","Intestinal Adenocarcinoma, Papillary Type","Intestinal Adenocarcinoma, Tubular Type")

data3%>%
  ggplot(aes(x=histological_type,y=DSS.time,fill=histological_type))+
  geom_boxplot()+
  scale_fill_brewer(palette="Set1")+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 17,vjust=0.5, color="black"),
    axis.text.y=element_text(color="black")
  ) +
  ggtitle("Stomach Cancer Distribution of Disease Specific Survival by Histological Stage") +
  labs(x="Histological Type",y="Disease Specific Survival (days)")

ggsave("Stomach - DSS by Histological Type.png",width=15,height=9)


#Survival Analysis
library(survival) #for survival analysis
library(survminer) #for plot
library(ggfortify)
library(ranger)

#remove duplicate patients for survival analysis (keep only 1 entry)
surv<-data%>%
  slice(-c(which(duplicated(data$patient_id))))

#Add DSS_month and DSS_year as a variables from DSS.time
surv<-mutate(surv,DSS_month=DSS.time/30,DSS_year=DSS.time/365)

#Add DFI_month and DFI_year as a variables from DFI.time
surv<-mutate(surv,DFI_month=DFI.time/30,DFI_year=DFI.time/365)

##HIF3A gene

#Add variable to separate observations based on gene expression
surv$HIF3A_f<-ifelse(surv$HIF3A>=median(as.numeric(surv$HIF3A)),"High HIF3A Expression","Low HIF3A Expression")

#Plot Kaplan-Meier Curve
plot1<-survfit(formula=Surv(DSS_year,DSS)~HIF3A_f,data=surv)
ggsurvplot(plot1,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High HIF3A Expression","Low HIF3A Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by HIF3A Gene Expression",pval=TRUE,risk.table=TRUE,conf.int=TRUE) #kaplan-meier plot

##HDAC5 gene - overexpressed in stomach cancer, poor prognosis

#Add variable to separate observations based on gene expression
surv$HDAC5_f<-ifelse(surv$HDAC5>=median(as.numeric(surv$HDAC5)),"High HDAC5 Expression","Low HDAC5 Expression")

#Plot Kaplan-Meier Curve
plot2<-survfit(formula=Surv(DSS_year,DSS)~HDAC5_f,data=surv)
ggsurvplot(plot2,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High HDAC5 Expression","Low HDAC5 Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by HDAC5 Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot

#ADH1B gene - downregulated in gastric cancer 

#Add variable to separate observations based on gene expression
surv$ADH1B_f<-ifelse(surv$ADH1B>=median(as.numeric(surv$ADH1B)),"High ADH1B Expression","Low ADH1B Expression")

#Plot Kaplan-Meier Curve
plot3<-survfit(formula=Surv(DSS_year,DSS)~ADH1B_f,data=surv)
ggsurvplot(plot3,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High ADH1B Expression","Low ADH1B Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by ADH1B Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot

#ADAM12 gene - upregulated in gastric cancer 

#Add variable to separate observations based on gene expression
surv$ADAM12_f<-ifelse(surv$ADAM12>=median(as.numeric(surv$ADAM12)),"High ADAM12 Expression","Low ADAM12 Expression")

#Plot Kaplan-Meier Curve
plot4<-survfit(formula=Surv(DSS_year,DSS)~ADAM12_f,data=surv)
ggsurvplot(plot4,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High ADAM12 Expression","Low ADAM12 Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by ADAM12 Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot

#ADAM12 - add DFI

#Plot Kaplan-Meier Curve
plot5<-survfit(formula=Surv(DFI_year,DFI)~ADAM12_f,data=surv)
ggsurvplot(plot5,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High ADAM12 Expression","Low ADAM12 Expression"),xlab="Years",ylab="Survival Probability",title="Disease Free Interval of Stomach Cancer by ADAM12 Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot

#PLAU gene

#Add variable to separate observations based on gene expression
surv$PLAU_f<-ifelse(surv$PLAU>=median(as.numeric(surv$PLAU)),"High PLAU Expression","Low PLAU Expression")

#Plot Kaplan-Meier Curve
plot6<-survfit(formula=Surv(DSS_year,DSS)~PLAU_f,data=surv)
ggsurvplot(plot6,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High PLAU Expression","Low PLAU Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by PLAU Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot

#MAP2K1 gene

#Add variable to separate observations based on gene expression
surv$MAP2K1_f<-ifelse(surv$MAP2K1>=median(as.numeric(surv$MAP2K1)),"High MAP2K1 Expression","Low MAP2K1 Expression")

#Plot Kaplan-Meier Curve
plot7<-survfit(formula=Surv(DSS_year,DSS)~MAP2K1_f,data=surv)
ggsurvplot(plot7,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High MAP2K1 Expression","Low MAP2K1 Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by MAP2K1 Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot

#FN1 gene

#Add variable to separate observations based on gene expression
surv$FN1_f<-ifelse(surv$FN1>=median(as.numeric(surv$FN1)),"High FN1 Expression","Low FN1 Expression")

#Plot Kaplan-Meier Curve
plot8<-survfit(formula=Surv(DSS_year,DSS)~FN1_f,data=surv)
ggsurvplot(plot8,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High FN1 Expression","Low FN1 Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by FN1 Gene Expression",pval=TRUE,risk.table=TRUE,conf.int=TRUE) #kaplan-meier plot


#CEP55 gene

#Add variable to separate observations based on gene expression
surv$CEP55_f<-ifelse(surv$CEP55>=median(as.numeric(surv$CEP55)),"High CEP55 Expression","Low CEP55 Expression")

#Plot Kaplan-Meier Curve
plot9<-survfit(formula=Surv(DSS_year,DSS)~CEP55_f,data=surv)
ggsurvplot(plot9,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High CEP55 Expression","Low CEP55 Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by CEP55 Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot


#PRTG gene
surv$PRTG_f<-ifelse(surv$PRTG >=median(as.numeric(surv$PRTG )),"High PRTG Expression","Low PRTG Expression")
table(surv$PRTG_f)

plot9<-survfit(formula=Surv(DSS_year,DSS)~PRTG_f,data=surv)
ggsurvplot(plot9,data=surv,surv.median.line="hv",legend.title="",legend.labs=c("High PRTG Expression","Low PRTG Expression"),xlab="Years",ylab="Survival Probability",title="Disease Specific Survival of Stomach Cancer by PRTG Gene Expression",pval=TRUE,risk.table=TRUE) #kaplan-meier plot
ggsave("PRTG Gene Kaplan Meier Curve.png",width=5,height=7)


#HISTOLOGICAL TYPE

#Plot Kaplan-Meier Curve
plot6<-survfit(formula=Surv(DSS_year,DSS)~histological_type,data=(filter(surv,histological_type!="[Discrepancy]"))) #formula
ggsurvplot(plot6,data=(filter(surv,histological_type!="[Discrepancy]")),,surv.median.line="hv",title="Disease Specific Survival of Stomach Cancer by Histological Type",legend.title="",xlab="Years",ylab="Disease Specific Survival Probability",legend.labs=c("Adenocarcinoma, Signet Ring Type","Adenocarcinoma, Diffuse Type","Adenocarcinoma, Not Otherwise Specified (NOS)","Intestinal Adenocarcinoma, Mucinous Type","Intestinal Adenocarcinoma, Not Otherwise Specified (NOS)","Intestinal Adenocarcinoma, Papillary Type","Intestinal Adenocarcinoma, Tubular Type"),pval=TRUE,risk.table=TRUE) #kaplan-meier plot

#MULTIPLE TESTING CORRECTION
library(dplyr)
library(tidyr)
library(survival)

#gene expression columns
gene_columns <- 119:20648
gene_data <- sapply(data[, gene_columns], as.numeric)

#calculate median expression
median <- apply(gene_data, 2, median)

#recode gene expressions based on median threshold
recode <- as.data.frame(t(apply(gene_data, 1, function(x) ifelse(x > median, 1, 0))))
recode$DSS.time <- data$DSS.time
recode$DSS <- data$DSS

#remove genes with NA values
recode <- recode %>% select(where(~ !any(is.na(.))) | one_of("DSS.time", "DSS"))

#remove genes with low variation (not enough 0s or)
valid_genes <- names(recode)[1:(ncol(recode) - 2)]
unique_counts <- sapply(valid_genes, function(gene) length(unique(recode[[gene]])))

#filter genes
valid_genes <- valid_genes[unique_counts > 1]

#remove gene C20orf185 (not enough variation)
valid_genes <- valid_genes[!valid_genes %in% "C20orf185"]

#survival analysis function
survival_analysis <- function(gene, dss_time, dss) {
  if (length(unique(gene)) < 2) {
    return(NA)
  }
  fit <- survfit(Surv(dss_time, dss) ~ gene)
  log_rank <- survdiff(Surv(dss_time, dss) ~ gene)
  p_value <- 1 - pchisq(log_rank$chisq, length(log_rank$n) - 1)
  return(p_value)
}

#calculate p-values for each valid gene
p_values <- sapply(valid_genes, function(gene) {
  tryCatch({
    survival_analysis(recode[[gene]], recode$DSS.time, recode$DSS)
  }, error = function(e) {
    cat("Error for gene:", gene, ":", e$message, "\n")
    return(NA)
  })
})

#filter
valid_p_values <- p_values[!is.na(p_values)]

#add adjusted p-value, using Benjamini-Hochberg and Bonferroni methods
bh_p_values <- p.adjust(valid_p_values, method = "BH")

#combine p-values into a dataframe
result<-data.frame(Gene = valid_genes[!is.na(p_values)], 
                      P_Value = valid_p_values, 
                      BH_P_Value = bh_p_values)
#arrange dataframe from most to least significant based on BH estimate
arrange(result,bh_p_values,descending=FALSE)

## COX PROPORTIONAL HAZARDS MODEL (allows inclusion of multiple predictors)

#PRTG gene
summary(coxph(formula=Surv(DSS_year,DSS)~PRTG_f,data=surv))

