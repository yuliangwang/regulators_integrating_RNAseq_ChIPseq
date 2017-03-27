#Integrate histone marks from diffReps with DESeq results 
#Submitted to Github on March 27th 2017
setwd("~/Dropbox/research_projects/nathan_palpant/corrected_analysis_2017/")

###Read files that identified differential H3K4me3- and H3K27me3-bound regions between each pair of conditions
K4me3_HA_CVP<-read.table('K4me3_HA_CVP.diff.txt.annotated',header=TRUE,sep="\t")
K27me3_HA_CVP<-read.table('K27me3_HA_CVP.diff.txt.annotated',header=TRUE,sep="\t")
K4me3_HA_HP<-read.table('K4me3_HA_HP.diff.txt.annotated',header=TRUE,sep="\t")
K27me3_HA_HP<-read.table('K27me3_HA_HP.diff.txt.annotated',header=TRUE,sep="\t")
K4me3_HP_CVP<-read.table('K4me3_HP_CVP.diff.txt.annotated',header=TRUE,sep="\t")
K27me3_HP_CVP<-read.table('K27me3_HP_CVP.diff.txt.annotated',header=TRUE,sep="\t")

###Each gene may be associated with multiple differential H3K4me3- and H3K27me3-bound regions. This function assigns the region with smallest p-value to each gene.
###Each gene will have at most one differential region. 
library(dplyr)
uniq_consistent_diffReps<-function(data){
  data<-data[-grep("^$",data$GName),]
  data<-arrange(data,GName,padj,pval)
  data<-data[!duplicated(data$GName),]
}

#Assign each gene a unique differential bound region and select relevant columns from diffReps results. 
K4me3_HA_CVP<-uniq_consistent_diffReps(K4me3_HA_CVP)
K4me3_HA_CVP<-K4me3_HA_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]

K27me3_HA_CVP<-uniq_consistent_diffReps(K27me3_HA_CVP)
K27me3_HA_CVP<-K27me3_HA_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K4me3_HA_HP<-uniq_consistent_diffReps(K4me3_HA_HP)
K4me3_HA_HP<-K4me3_HA_HP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K27me3_HA_HP<-uniq_consistent_diffReps(K27me3_HA_HP)
K27me3_HA_HP<-K27me3_HA_HP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K4me3_HP_CVP<-uniq_consistent_diffReps(K4me3_HP_CVP)
K4me3_HP_CVP<-K4me3_HP_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K27me3_HP_CVP<-uniq_consistent_diffReps(K27me3_HP_CVP)
K27me3_HP_CVP<-K27me3_HP_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]

#"T5_deseq.RData" stores DESeq results compair each pair of cell lines. 
load("T5_deseq.RData")
CVP_HA<-CVP_HA[complete.cases(CVP_HA),]
CVP_HP<-CVP_HP[complete.cases(CVP_HP),]
HA_HP<-HA_HP[complete.cases(HA_HP),]

#################  First, find candidate lineage regulators in each pairwise cell line comparison#############
###HA_CVP
#An integrated candidate regulator first needs to be differentially expressed (2 fold change, FDR<0.05)
deg<-subset(CVP_HA,subset=padj<0.05&abs(log2FoldChange)>1,select=gene_name,drop=TRUE)
#Second, it needs to have significant changes in either H3K4me3 or H3K27me3
chipseq<-union(K27me3_HA_CVP$GName,K4me3_HA_CVP$GName)
#That is, significant changes in both RNA-seq and ChIP-seq
common<-intersect(deg,chipseq)
num_common<-length(common)
#Compile RNA-seq and ChIP-seq change information for these regulators
initial<-matrix(NA,num_common,1)
ind<-match(common,CVP_HA$gene_name)
regulators_HA_CVP<-data.frame(gene_name=common,
                              meanHA=CVP_HA$BaseMeanHA[ind],
                              meanCVP=CVP_HA$BaseMeanCVP[ind],
                              log2FoldChange=CVP_HA$log2FoldChange[ind],
                              deseq_padj=CVP_HA$padj[ind],
                              K27me3_meanHA=initial,
                              K27me3_meanCVP=initial,
                              K27me3_log2FoldChange=initial,
                              K27me3_padj=initial,
                              K4me3_meanHA=initial,
                              K4me3_meanCVP=initial,
                              K4me3_log2FoldChange=initial,
                              K4me3_padj=initial
)

ind<-match(common,K27me3_HA_CVP$GName)
temp<-K27me3_HA_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HA_CVP[!is.na(ind),6:9]<-temp

ind<-match(common,K4me3_HA_CVP$GName)
temp<-K4me3_HA_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HA_CVP[!is.na(ind),10:13]<-temp
ind<-c(5,9,13)
pvals<-regulators_HA_CVP[,ind]
pvals[is.na(pvals)]<-1
for (i in 1:3){
  ind<-pvals[,i]==0
  pvals[ind,i]<-min(pvals[!ind,i])
}
ind<-c(4,8,12)
fc<-regulators_HA_CVP[,ind]
fc[is.na(fc)]<-0
for (i in 1:3){
  ind<-abs(fc[,i])==Inf
  bound<-max(abs(fc[!ind,i]))
  fc[fc[,i]==Inf,i]<-bound
  fc[fc[,i]==-Inf,i]<-(-bound)
}

regulators_HA_CVP$score<-sign(fc[,1])*abs(log10(pvals[,1]))-sign(fc[,2])*abs(log10(pvals[,2]))+sign(fc[,3])*abs(log10(pvals[,3]))
#regulators_HA_CVP$score<-sign(fc[,1])*abs(log10(pvals[,1]))# When we only consider RNAseq data, and ignore ChIP-seq data, use this line
regulators_HA_CVP<-regulators_HA_CVP[order(regulators_HA_CVP$score),]


###Same analysis for HP vs. CVP
deg<-subset(CVP_HP,subset=padj<0.05&abs(log2FoldChange)>1,select=gene_name,drop=TRUE)
chipseq<-union(K27me3_HP_CVP$GName,K4me3_HP_CVP$GName)
common<-intersect(deg,chipseq)
num_common<-length(common)
initial<-matrix(NA,num_common,1)
ind<-match(common,CVP_HP$gene_name)
regulators_HP_CVP<-data.frame(gene_name=common,
                              meanHP=CVP_HP$BaseMeanHP[ind],
                              meanCVP=CVP_HP$BaseMeanCVP[ind],
                              log2FoldChange=CVP_HP$log2FoldChange[ind],
                              deseq_padj=CVP_HP$padj[ind],
                              K27me3_meanHP=initial,
                              K27me3_meanCVP=initial,
                              K27me3_log2FoldChange=initial,
                              K27me3_padj=initial,
                              K4me3_meanHP=initial,
                              K4me3_meanCVP=initial,
                              K4me3_log2FoldChange=initial,
                              K4me3_padj=initial
)


ind<-match(common,K27me3_HP_CVP$GName)
temp<-K27me3_HP_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HP_CVP[!is.na(ind),6:9]<-temp

ind<-match(common,K4me3_HP_CVP$GName)
temp<-K4me3_HP_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HP_CVP[!is.na(ind),10:13]<-temp

ind<-c(5,9,13)
pvals<-regulators_HP_CVP[,ind]
pvals[is.na(pvals)]<-1
for (i in 1:3){
  ind<-pvals[,i]==0
  pvals[ind,i]<-min(pvals[!ind,i])
}
ind<-c(4,8,12)
fc<-regulators_HP_CVP[,ind]
fc[is.na(fc)]<-0
for (i in 1:3){
  ind<-abs(fc[,i])==Inf
  bound<-max(abs(fc[!ind,i]))
  fc[fc[,i]==Inf,i]<-bound
  fc[fc[,i]==-Inf,i]<-(-bound)
}
regulators_HP_CVP$score<-sign(fc[,1])*abs(log10(pvals[,1]))-sign(fc[,2])*abs(log10(pvals[,2]))+sign(fc[,3])*abs(log10(pvals[,3]))
#regulators_HP_CVP$score<-sign(fc[,1])*abs(log10(pvals[,1]))
regulators_HP_CVP<-regulators_HP_CVP[order(regulators_HP_CVP$score),]

###Same analysis for HA vs. HP
deg<-subset(HA_HP,subset=padj<0.05&abs(log2FoldChange)>1,select=gene_name,drop=TRUE)
chipseq<-union(K27me3_HA_HP$GName,K4me3_HA_HP$GName)
common<-intersect(deg,chipseq)
num_common<-length(common)
initial<-matrix(NA,num_common,1)
ind<-match(common,HA_HP$gene_name)
regulators_HP_HA<-data.frame(gene_name=common,
                             meanHP=HA_HP$BaseMeanHP[ind],
                             meanHA=HA_HP$BaseMeanHA[ind],
                             log2FoldChange=HA_HP$log2FoldChange[ind],
                             deseq_padj=HA_HP$padj[ind],
                             K27me3_meanHP=initial,
                             K27me3_meanHA=initial,
                             K27me3_log2FoldChange=initial,
                             K27me3_padj=initial,
                             K4me3_meanHP=initial,
                             K4me3_meanHA=initial,
                             K4me3_log2FoldChange=initial,
                             K4me3_padj=initial
)


ind<-match(common,K27me3_HA_HP$GName)
temp<-K27me3_HA_HP[ind[!is.na(ind)],]
temp$log2FC<-(-temp$log2FC)
temp<-subset(temp,select=c(5,4,3,2))
regulators_HP_HA[!is.na(ind),6:9]<-temp

ind<-match(common,K4me3_HA_HP$GName)
temp<-K4me3_HA_HP[ind[!is.na(ind)],]
temp$log2FC<-(-temp$log2FC)
temp<-subset(temp,select=c(5,4,3,2))
regulators_HP_HA[!is.na(ind),10:13]<-temp

ind<-c(5,9,13)
pvals<-regulators_HP_HA[,ind]
pvals[is.na(pvals)]<-1
for (i in 1:3){
  ind<-pvals[,i]==0
  pvals[ind,i]<-min(pvals[!ind,i])
}
ind<-c(4,8,12)
fc<-regulators_HP_HA[,ind]
fc[is.na(fc)]<-0
for (i in 1:3){
  ind<-abs(fc[,i])==Inf
  bound<-max(abs(fc[!ind,i]))
  fc[fc[,i]==Inf,i]<-bound
  fc[fc[,i]==-Inf,i]<-(-bound)
}

regulators_HP_HA$score<-sign(fc[,1])*abs(log10(pvals[,1]))-sign(fc[,2])*abs(log10(pvals[,2]))+sign(fc[,3])*abs(log10(pvals[,3]))
#regulators_HP_HA$score<-sign(fc[,1])*abs(log10(pvals[,1]))
regulators_HP_HA<-regulators_HP_HA[order(regulators_HP_HA$score),]

###Identify candidate regulators for each lineage###
#For example, candidate CVP regulators are genes that scores higher for CVP in both HA vs. CVP and HP vs. CVP pairwise comparisons. 
common_CVP<-intersect(subset(regulators_HA_CVP,score<0,select=gene_name,drop=TRUE),subset(regulators_HP_CVP,score<0,select=gene_name,drop=TRUE))
common_HA<-intersect(subset(regulators_HP_HA,score<0,select=gene_name,drop=TRUE),subset(regulators_HA_CVP,score>0,select=gene_name,drop=TRUE))
common_HP<-intersect(subset(regulators_HP_HA,score>0,select=gene_name,drop=TRUE),subset(regulators_HP_CVP,score>0,select=gene_name,drop=TRUE))

#Lineage regulators common to HP and HA are genes scores higher in HA and HP in HA vs. CVP and HP vs. CVP.
common_HAHP<-intersect(subset(regulators_HP_CVP,score>0,select=gene_name,drop=TRUE),subset(regulators_HA_CVP,score>0,select=gene_name,drop=TRUE))
#Lineage regulators common to CVP and HA are genes scores higher in CVP and HA in CVP vs. HP and HA vs. HP. 
common_HACVP<-intersect(subset(regulators_HP_CVP,score<0,select=gene_name,drop=TRUE),subset(regulators_HP_HA,score<0,select=gene_name,drop=TRUE))


#Compile RNA-seq and ChIP-seq data for these regulators 
ind1<-match(common_CVP,regulators_HP_CVP$gene_name)
ind2<-match(common_CVP,regulators_HA_CVP$gene_name)
regulators_CVP<-data.frame(gene_name=common_CVP,
                           meanCVP=regulators_HP_CVP$meanCVP[ind1],
                           meanHA=regulators_HA_CVP$meanHA[ind2],
                           meanHP=regulators_HP_CVP$meanHP[ind1],
                           log2FoldChangeCVP_HA=(-regulators_HA_CVP$log2FoldChange[ind2]),
                           log2FoldChangeCVP_HP=(-regulators_HP_CVP$log2FoldChange[ind1]),           
                           deseq_padj_CVP_HA=regulators_HA_CVP$deseq_padj[ind2],
                           deseq_padj_CVP_HP=regulators_HP_CVP$deseq_padj[ind1],
                           K27me3_log2FoldChange_CVP_HA=(-regulators_HA_CVP$K27me3_log2FoldChange[ind2]),
                           K27me3_log2FoldChange_CVP_HP=(-regulators_HP_CVP$K27me3_log2FoldChange[ind1]),
                           K27me3_padj_CVP_HA=regulators_HA_CVP$K27me3_padj[ind2],
                           K27me3_padj_CVP_HP=regulators_HP_CVP$K27me3_padj[ind1],
                           K4me3_log2FoldChange_CVP_HA=(-regulators_HA_CVP$K4me3_log2FoldChange[ind2]),
                           K4me3_log2FoldChange_CVP_HP=(-regulators_HP_CVP$K4me3_log2FoldChange[ind1]),
                           K4me3_padj_CVP_HA=regulators_HA_CVP$K4me3_padj[ind2],
                           K4me3_padj_CVP_HP=regulators_HP_CVP$K4me3_padj[ind1],
                           score_CVP_HA=-regulators_HA_CVP$score[ind2],
                           score_CVP_HP=-regulators_HP_CVP$score[ind1],
                           score_aggregate=-regulators_HA_CVP$score[ind2]-regulators_HP_CVP$score[ind1],
                           Function=regulators_HA_CVP$Function[ind2]
)

regulators_CVP<-regulators_CVP[order(-regulators_CVP$score_aggregate),]

ind1<-match(common_HP,regulators_HP_CVP$gene_name)
ind2<-match(common_HP,regulators_HP_HA$gene_name)
regulators_HP<-data.frame(gene_name=common_HP,
                          meanHP=regulators_HP_CVP$meanHP[ind1],
                          meanHA=regulators_HP_HA$meanHA[ind2],
                          meanCVP=regulators_HP_CVP$meanCVP[ind1],
                          log2FoldChangeHP_HA=regulators_HP_HA$log2FoldChange[ind2],
                          log2FoldChangeHP_CVP=regulators_HP_CVP$log2FoldChange[ind1],           
                          deseq_padj_HP_HA=regulators_HP_HA$deseq_padj[ind2],
                          deseq_padj_HP_CVP=regulators_HP_CVP$deseq_padj[ind1],
                          K27me3_log2FoldChange_HP_HA=regulators_HP_HA$K27me3_log2FoldChange[ind2],                         
                          K27me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K27me3_log2FoldChange[ind1],
                          K27me3_padj_HP_HA=regulators_HP_HA$K27me3_padj[ind2],
                          K27me3_padj_HP_CVP=regulators_HP_CVP$K27me3_padj[ind1],
                          K4me3_log2FoldChange_HP_HA=regulators_HP_HA$K4me3_log2FoldChange[ind2],                         
                          K4me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K4me3_log2FoldChange[ind1],
                          K4me3_padj_HP_HA=regulators_HP_HA$K4me3_padj[ind2],
                          K4me3_padj_HP_CVP=regulators_HP_CVP$K4me3_padj[ind1],
                          score_HP_HA=regulators_HP_HA$score[ind2],
                          score_HP_CVP=regulators_HP_CVP$score[ind1],
                          score_aggregate=regulators_HP_HA$score[ind2]+regulators_HP_CVP$score[ind1],
                          Function=regulators_HP_CVP$Function[ind1]
)

regulators_HP<-regulators_HP[order(-regulators_HP$score_aggregate),]

ind1<-match(common_HA,regulators_HA_CVP$gene_name)
ind2<-match(common_HA,regulators_HP_HA$gene_name)
regulators_HA<-data.frame(gene_name=common_HA,
                          meanHA=regulators_HP_HA$meanHA[ind2],                         
                          meanHP=regulators_HP_HA$meanHP[ind2],
                          meanCVP=regulators_HA_CVP$meanCVP[ind1],
                          log2FoldChangeHA_HP=-regulators_HP_HA$log2FoldChange[ind2],
                          log2FoldChangeHA_CVP=regulators_HA_CVP$log2FoldChange[ind1],           
                          deseq_padj_HA_HP=regulators_HP_HA$deseq_padj[ind2],
                          deseq_padj_HA_CVP=regulators_HA_CVP$deseq_padj[ind1],
                          K27me3_log2FoldChange_HA_HP=-regulators_HP_HA$K27me3_log2FoldChange[ind2],                         
                          K27me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K27me3_log2FoldChange[ind1],
                          K27me3_padj_HA_HP=regulators_HP_HA$K27me3_padj[ind2],
                          K27me3_padj_HA_CVP=regulators_HA_CVP$K27me3_padj[ind1],
                          K4me3_log2FoldChange_HA_HP=-regulators_HP_HA$K4me3_log2FoldChange[ind2],                         
                          K4me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K4me3_log2FoldChange[ind1],
                          K4me3_padj_HA_HP=regulators_HP_HA$K4me3_padj[ind2],
                          K4me3_padj_HA_CVP=regulators_HA_CVP$K4me3_padj[ind1],
                          score_HA_HP=-regulators_HP_HA$score[ind2],
                          score_HA_CVP=regulators_HA_CVP$score[ind1],
                          score_aggregate=-regulators_HP_HA$score[ind2]+regulators_HA_CVP$score[ind1],
                          Function=regulators_HA_CVP$Function[ind1]
)

regulators_HA<-regulators_HA[order(-regulators_HA$score_aggregate),]

ind1<-match(common_HAHP,regulators_HA_CVP$gene_name)
ind2<-match(common_HAHP,regulators_HP_CVP$gene_name)
regulators_HAHP<-data.frame(gene_name=common_HAHP,
                            meanHA=regulators_HA_CVP$meanHA[ind1],                         
                            meanHP=regulators_HP_CVP$meanHP[ind2],
                            meanCVP=regulators_HA_CVP$meanCVP[ind1],
                            log2FoldChangeHA_CVP=regulators_HA_CVP$log2FoldChange[ind1],
                            log2FoldChangeHP_CVP=regulators_HP_CVP$log2FoldChange[ind2],           
                            deseq_padj_HA_CVP=regulators_HA_CVP$deseq_padj[ind1],
                            deseq_padj_HP_CVP=regulators_HP_CVP$deseq_padj[ind2],
                            K27me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K27me3_log2FoldChange[ind1],                         
                            K27me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K27me3_log2FoldChange[ind2],
                            K27me3_padj_HA_CVP=regulators_HA_CVP$K27me3_padj[ind1],
                            K27me3_padj_HP_CVP=regulators_HP_CVP$K27me3_padj[ind2],
                            K4me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K4me3_log2FoldChange[ind1],                         
                            K4me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K4me3_log2FoldChange[ind2],
                            K4me3_padj_HA_CVP=regulators_HA_CVP$K4me3_padj[ind1],
                            K4me3_padj_HP_CVP=regulators_HP_CVP$K4me3_padj[ind2],
                            score_HA_CVP=regulators_HA_CVP$score[ind1],
                            score_HP_CVP=regulators_HP_CVP$score[ind2],
                            score_aggregate=regulators_HA_CVP$score[ind1]+regulators_HP_CVP$score[ind2],
                            Function=regulators_HA_CVP$Function[ind1]
)

regulators_HAHP<-regulators_HAHP[order(-regulators_HAHP$score_aggregate),]

ind1<-match(common_HACVP,regulators_HP_HA$gene_name)
ind2<-match(common_HACVP,regulators_HP_CVP$gene_name)
regulators_HACVP<-data.frame(gene_name=common_HACVP,
                             meanHA=regulators_HP_HA$meanHA[ind1], 
                             meanCVP=regulators_HP_CVP$meanCVP[ind2],                           
                             meanHP=regulators_HP_HA$meanHP[ind1],
                             log2FoldChangeHA_HP=-regulators_HP_HA$log2FoldChange[ind1],
                             log2FoldChangeCVP_HP=-regulators_HP_CVP$log2FoldChange[ind2],           
                             deseq_padj_HA_HP=regulators_HP_HA$deseq_padj[ind1],
                             deseq_padj_CVP_HP=regulators_HP_CVP$deseq_padj[ind2],
                             K27me3_log2FoldChange_HA_HP=-regulators_HP_HA$K27me3_log2FoldChange[ind1],                         
                             K27me3_log2FoldChange_CVP_HP=-regulators_HP_CVP$K27me3_log2FoldChange[ind2],
                             K27me3_padj_HA_HP=regulators_HP_HA$K27me3_padj[ind1],
                             K27me3_padj_CVP_HP=regulators_HP_CVP$K27me3_padj[ind2],
                             K4me3_log2FoldChange_HA_HP=-regulators_HP_HA$K4me3_log2FoldChange[ind1],                         
                             K4me3_log2FoldChange_CVP_HP=-regulators_HP_CVP$K4me3_log2FoldChange[ind2],
                             K4me3_padj_HA_HP=regulators_HA_CVP$K4me3_padj[ind1],
                             K4me3_padj_CVP_HP=regulators_HP_CVP$K4me3_padj[ind2],
                             score_HA_HP=-regulators_HP_HA$score[ind1],
                             score_CVP_HP=-regulators_HP_CVP$score[ind2],
                             score_aggregate=-regulators_HP_HA$score[ind1]-regulators_HP_CVP$score[ind2],
                             Function=regulators_HP_CVP$Function[ind2]
)

regulators_HACVP<-regulators_HACVP[order(-regulators_HACVP$score_aggregate),]

#For regulators shared in HA/CVP, or HA/HP, remove those already identified as regulators in each individual lineage
regulators_HACVP_exclude_indivi<-regulators_HACVP[!(regulators_HACVP$gene_name %in% regulators_HA$gene_name |regulators_HACVP$gene_name %in% regulators_CVP$gene_name),]
regulators_HAHP_exclude_indivi<-regulators_HAHP[!(regulators_HAHP$gene_name %in% regulators_HA$gene_name |regulators_HAHP$gene_name %in% regulators_HP$gene_name),]


###The integrated score are derived from RNA-seq, H3K4me3 and H3K27me3. Calculate the contribution of each data source to the final score
#HA breakdown
ind1<-match(regulators_HA$gene_name,regulators_HA_CVP$gene_name)
ind2<-match(regulators_HA$gene_name,regulators_HP_HA$gene_name)
breakdown_HA<-data.frame(gene_name=regulators_HA$gene_name,exprs.HAvsCVP=-log10(regulators_HA_CVP$deseq_padj[ind1]),
                         exprs.HAvsHP=-log10(regulators_HP_HA$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HA_CVP$deseq_padj[ind1])-log10(regulators_HP_HA$deseq_padj[ind2]))/2,
                         K4me3.HAvsCVP=-log10(regulators_HA_CVP$K4me3_padj[ind1]),K4me3.HAvsHP=-log10(regulators_HP_HA$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HA_CVP$K4me3_padj[ind1])-log10(regulators_HP_HA$K4me3_padj[ind2]))/2,
                         K27me3.HAvsCVP=-log10(regulators_HA_CVP$K27me3_padj[ind1]),K27me3.HAvsHP=-log10(regulators_HP_HA$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HA_CVP$K27me3_padj[ind1])-log10(regulators_HP_HA$K27me3_padj[ind2]))/2,
                         final_score=regulators_HA$score_aggregate)
#CVP breakdown
ind1<-match(regulators_CVP$gene_name,regulators_HA_CVP$gene_name)
ind2<-match(regulators_CVP$gene_name,regulators_HP_CVP$gene_name)
breakdown_CVP<-data.frame(gene_name=regulators_CVP$gene_name,exprs.CVPvsHA=-log10(regulators_HA_CVP$deseq_padj[ind1]),
                          exprs.CVPvsHP=-log10(regulators_HP_CVP$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HA_CVP$deseq_padj[ind1])-log10(regulators_HP_CVP$deseq_padj[ind2]))/2,
                          K4me3.CVPvsHA=-log10(regulators_HA_CVP$K4me3_padj[ind1]),K4me3.CVPvsHP=-log10(regulators_HP_CVP$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HA_CVP$K4me3_padj[ind1])-log10(regulators_HP_CVP$K4me3_padj[ind2]))/2,
                          K27me3.CVPvsHA=-log10(regulators_HA_CVP$K27me3_padj[ind1]),K27me3.CVPvsHP=-log10(regulators_HP_CVP$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HA_CVP$K27me3_padj[ind1])-log10(regulators_HP_CVP$K27me3_padj[ind2]))/2,
                          final_score=regulators_CVP$score_aggregate)
#HP breakdown
ind2<-match(regulators_HP$gene_name,regulators_HP_HA$gene_name)
ind1<-match(regulators_HP$gene_name,regulators_HP_CVP$gene_name)
breakdown_HP<-data.frame(gene_name=regulators_HP$gene_name,exprs.HPvsCVP=-log10(regulators_HP_CVP$deseq_padj[ind1]),
                         exprs.HPvsHA=-log10(regulators_HP_HA$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HP_CVP$deseq_padj[ind1])-log10(regulators_HP_HA$deseq_padj[ind2]))/2,
                         K4me3.HPvsCVP=-log10(regulators_HP_CVP$K4me3_padj[ind1]),K4me3.HPvsHA=-log10(regulators_HP_HA$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HP_CVP$K4me3_padj[ind1])-log10(regulators_HP_HA$K4me3_padj[ind2]))/2,
                         K27me3.HPvsCVP=-log10(regulators_HP_CVP$K27me3_padj[ind1]),K27me3.HPvsHA=-log10(regulators_HP_HA$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HP_CVP$K27me3_padj[ind1])-log10(regulators_HP_HA$K27me3_padj[ind2]))/2,
                         final_score=regulators_HP$score_aggregate)

#HAHP breakdown
ind1<-match(regulators_HAHP_exclude_indivi$gene_name,regulators_HA_CVP$gene_name)
ind2<-match(regulators_HAHP_exclude_indivi$gene_name,regulators_HP_CVP$gene_name)
breakdown_HAHP<-data.frame(gene_name=regulators_HAHP_exclude_indivi$gene_name,exprs.CVPvsHA=-log10(regulators_HA_CVP$deseq_padj[ind1]),
                           exprs.CVPvsHP=-log10(regulators_HP_CVP$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HA_CVP$deseq_padj[ind1])-log10(regulators_HP_CVP$deseq_padj[ind2]))/2,
                           K4me3.CVPvsHA=-log10(regulators_HA_CVP$K4me3_padj[ind1]),K4me3.CVPvsHP=-log10(regulators_HP_CVP$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HA_CVP$K4me3_padj[ind1])-log10(regulators_HP_CVP$K4me3_padj[ind2]))/2,
                           K27me3.CVPvsHA=-log10(regulators_HA_CVP$K27me3_padj[ind1]),K27me3.CVPvsHP=-log10(regulators_HP_CVP$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HA_CVP$K27me3_padj[ind1])-log10(regulators_HP_CVP$K27me3_padj[ind2]))/2,
                           final_score=regulators_HAHP_exclude_indivi$score_aggregate)

#HACVP breakdown
ind2<-match(regulators_HACVP_exclude_indivi$gene_name,regulators_HP_HA$gene_name)
ind1<-match(regulators_HACVP_exclude_indivi$gene_name,regulators_HP_CVP$gene_name)
breakdown_HACVP<-data.frame(gene_name=regulators_HACVP_exclude_indivi$gene_name,exprs.HPvsCVP=-log10(regulators_HP_CVP$deseq_padj[ind1]),
                            exprs.HPvsHA=-log10(regulators_HP_HA$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HP_CVP$deseq_padj[ind1])-log10(regulators_HP_HA$deseq_padj[ind2]))/2,
                            K4me3.HPvsCVP=-log10(regulators_HP_CVP$K4me3_padj[ind1]),K4me3.HPvsHA=-log10(regulators_HP_HA$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HP_CVP$K4me3_padj[ind1])-log10(regulators_HP_HA$K4me3_padj[ind2]))/2,
                            K27me3.HPvsCVP=-log10(regulators_HP_CVP$K27me3_padj[ind1]),K27me3.HPvsHA=-log10(regulators_HP_HA$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HP_CVP$K27me3_padj[ind1])-log10(regulators_HP_HA$K27me3_padj[ind2]))/2,
                            final_score=regulators_HACVP_exclude_indivi$score_aggregate)




##########################################################################
#### Only use ChIP-seq data in ranking genes#####
K4me3_HA_CVP<-read.table('K4me3_HA_CVP.diff.txt.annotated',header=TRUE,sep="\t")
K27me3_HA_CVP<-read.table('K27me3_HA_CVP.diff.txt.annotated',header=TRUE,sep="\t")
K4me3_HA_HP<-read.table('K4me3_HA_HP.diff.txt.annotated',header=TRUE,sep="\t")
K27me3_HA_HP<-read.table('K27me3_HA_HP.diff.txt.annotated',header=TRUE,sep="\t")
K4me3_HP_CVP<-read.table('K4me3_HP_CVP.diff.txt.annotated',header=TRUE,sep="\t")
K27me3_HP_CVP<-read.table('K27me3_HP_CVP.diff.txt.annotated',header=TRUE,sep="\t")

library(dplyr)
uniq_consistent_diffReps<-function(data){
  data<-data[-grep("^$",data$GName),]
  data<-arrange(data,GName,padj,pval)
  data<-data[!duplicated(data$GName),]
}

K4me3_HA_CVP<-uniq_consistent_diffReps(K4me3_HA_CVP)
K4me3_HA_CVP<-K4me3_HA_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]

K27me3_HA_CVP<-uniq_consistent_diffReps(K27me3_HA_CVP)
K27me3_HA_CVP<-K27me3_HA_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K4me3_HA_HP<-uniq_consistent_diffReps(K4me3_HA_HP)
K4me3_HA_HP<-K4me3_HA_HP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K27me3_HA_HP<-uniq_consistent_diffReps(K27me3_HA_HP)
K27me3_HA_HP<-K27me3_HA_HP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K4me3_HP_CVP<-uniq_consistent_diffReps(K4me3_HP_CVP)
K4me3_HP_CVP<-K4me3_HP_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]
K27me3_HP_CVP<-uniq_consistent_diffReps(K27me3_HP_CVP)
K27me3_HP_CVP<-K27me3_HP_CVP[,c("GName","padj","log2FC","Treatment.avg","Control.avg","Event","Feature","D2TSS")]

#######Same analysis proceedure as in the integrated regulator appraoch, but does not require a regulator to be differentially expressed at the RNA level.
###HA_CVP
chipseq<-union(K27me3_HA_CVP$GName,K4me3_HA_CVP$GName)
common<-chipseq
num_common<-length(common)
initial<-matrix(NA,num_common,1)
ind<-match(common,CVP_HA$gene_name)
regulators_HA_CVP<-data.frame(gene_name=common,
                              meanHA=CVP_HA$BaseMeanHA[ind],
                              meanCVP=CVP_HA$BaseMeanCVP[ind],
                              log2FoldChange=CVP_HA$log2FoldChange[ind],
                              deseq_padj=CVP_HA$padj[ind],
                              K27me3_meanHA=initial,
                              K27me3_meanCVP=initial,
                              K27me3_log2FoldChange=initial,
                              K27me3_padj=initial,
                              K4me3_meanHA=initial,
                              K4me3_meanCVP=initial,
                              K4me3_log2FoldChange=initial,
                              K4me3_padj=initial
)

ind<-match(common,K27me3_HA_CVP$GName)
temp<-K27me3_HA_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HA_CVP[!is.na(ind),6:9]<-temp

ind<-match(common,K4me3_HA_CVP$GName)
temp<-K4me3_HA_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HA_CVP[!is.na(ind),10:13]<-temp
ind<-c(5,9,13)
pvals<-regulators_HA_CVP[,ind]
pvals[is.na(pvals)]<-1
for (i in 1:3){
  ind<-pvals[,i]==0
  pvals[ind,i]<-min(pvals[!ind,i])
}
ind<-c(4,8,12)
fc<-regulators_HA_CVP[,ind]
fc[is.na(fc)]<-0
for (i in 1:3){
  ind<-abs(fc[,i])==Inf
  bound<-max(abs(fc[!ind,i]))
  fc[fc[,i]==Inf,i]<-bound
  fc[fc[,i]==-Inf,i]<-(-bound)
}

regulators_HA_CVP$score<- -sign(fc[,2])*abs(log10(pvals[,2]))+sign(fc[,3])*abs(log10(pvals[,3]))
regulators_HA_CVP<-regulators_HA_CVP[order(regulators_HA_CVP$score),]
############HP_CVP
chipseq<-union(K27me3_HP_CVP$GName,K4me3_HP_CVP$GName)
common<-chipseq
num_common<-length(common)
initial<-matrix(NA,num_common,1)
ind<-match(common,CVP_HP$gene_name)
regulators_HP_CVP<-data.frame(gene_name=common,
                              meanHP=CVP_HP$BaseMeanHP[ind],
                              meanCVP=CVP_HP$BaseMeanCVP[ind],
                              log2FoldChange=CVP_HP$log2FoldChange[ind],
                              deseq_padj=CVP_HP$padj[ind],
                              K27me3_meanHP=initial,
                              K27me3_meanCVP=initial,
                              K27me3_log2FoldChange=initial,
                              K27me3_padj=initial,
                              K4me3_meanHP=initial,
                              K4me3_meanCVP=initial,
                              K4me3_log2FoldChange=initial,
                              K4me3_padj=initial
)


ind<-match(common,K27me3_HP_CVP$GName)
temp<-K27me3_HP_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HP_CVP[!is.na(ind),6:9]<-temp

ind<-match(common,K4me3_HP_CVP$GName)
temp<-K4me3_HP_CVP[ind[!is.na(ind)],]
temp<-subset(temp,select=c(4,5,3,2))
regulators_HP_CVP[!is.na(ind),10:13]<-temp

ind<-c(5,9,13)
pvals<-regulators_HP_CVP[,ind]
pvals[is.na(pvals)]<-1
for (i in 1:3){
  ind<-pvals[,i]==0
  pvals[ind,i]<-min(pvals[!ind,i])
}
ind<-c(4,8,12)
fc<-regulators_HP_CVP[,ind]
fc[is.na(fc)]<-0
for (i in 1:3){
  ind<-abs(fc[,i])==Inf
  bound<-max(abs(fc[!ind,i]))
  fc[fc[,i]==Inf,i]<-bound
  fc[fc[,i]==-Inf,i]<-(-bound)
}
regulators_HP_CVP$score<- -sign(fc[,2])*abs(log10(pvals[,2]))+sign(fc[,3])*abs(log10(pvals[,3]))
regulators_HP_CVP<-regulators_HP_CVP[order(regulators_HP_CVP$score),]
###################HP_HA
chipseq<-union(K27me3_HA_HP$GName,K4me3_HA_HP$GName)
common<-chipseq
num_common<-length(common)
initial<-matrix(NA,num_common,1)
ind<-match(common,HA_HP$gene_name)
regulators_HP_HA<-data.frame(gene_name=common,
                             meanHP=HA_HP$BaseMeanHP[ind],
                             meanHA=HA_HP$BaseMeanHA[ind],
                             log2FoldChange=HA_HP$log2FoldChange[ind],
                             deseq_padj=HA_HP$padj[ind],
                             K27me3_meanHP=initial,
                             K27me3_meanHA=initial,
                             K27me3_log2FoldChange=initial,
                             K27me3_padj=initial,
                             K4me3_meanHP=initial,
                             K4me3_meanHA=initial,
                             K4me3_log2FoldChange=initial,
                             K4me3_padj=initial
)


ind<-match(common,K27me3_HA_HP$GName)
temp<-K27me3_HA_HP[ind[!is.na(ind)],]
temp$log2FC<-(-temp$log2FC)
temp<-subset(temp,select=c(5,4,3,2))
regulators_HP_HA[!is.na(ind),6:9]<-temp

ind<-match(common,K4me3_HA_HP$GName)
temp<-K4me3_HA_HP[ind[!is.na(ind)],]
temp$log2FC<-(-temp$log2FC)
temp<-subset(temp,select=c(5,4,3,2))
regulators_HP_HA[!is.na(ind),10:13]<-temp

ind<-c(5,9,13)
pvals<-regulators_HP_HA[,ind]
pvals[is.na(pvals)]<-1
for (i in 1:3){
  ind<-pvals[,i]==0
  pvals[ind,i]<-min(pvals[!ind,i])
}
ind<-c(4,8,12)
fc<-regulators_HP_HA[,ind]
fc[is.na(fc)]<-0
for (i in 1:3){
  ind<-abs(fc[,i])==Inf
  bound<-max(abs(fc[!ind,i]))
  fc[fc[,i]==Inf,i]<-bound
  fc[fc[,i]==-Inf,i]<-(-bound)
}

regulators_HP_HA$score<- -sign(fc[,2])*abs(log10(pvals[,2]))+sign(fc[,3])*abs(log10(pvals[,3]))
regulators_HP_HA<-regulators_HP_HA[order(regulators_HP_HA$score),]


common_CVP<-intersect(subset(regulators_HA_CVP,score<0,select=gene_name,drop=TRUE),subset(regulators_HP_CVP,score<0,select=gene_name,drop=TRUE))
common_HA<-intersect(subset(regulators_HP_HA,score<0,select=gene_name,drop=TRUE),subset(regulators_HA_CVP,score>0,select=gene_name,drop=TRUE))
common_HP<-intersect(subset(regulators_HP_HA,score>0,select=gene_name,drop=TRUE),subset(regulators_HP_CVP,score>0,select=gene_name,drop=TRUE))
common_HAHP<-intersect(subset(regulators_HP_CVP,score>0,select=gene_name,drop=TRUE),subset(regulators_HA_CVP,score>0,select=gene_name,drop=TRUE))
common_HACVP<-intersect(subset(regulators_HP_CVP,score<0,select=gene_name,drop=TRUE),subset(regulators_HP_HA,score<0,select=gene_name,drop=TRUE))

ind1<-match(common_CVP,regulators_HP_CVP$gene_name)
ind2<-match(common_CVP,regulators_HA_CVP$gene_name)
regulators_CVP<-data.frame(gene_name=common_CVP,
                           meanCVP=regulators_HP_CVP$meanCVP[ind1],
                           meanHA=regulators_HA_CVP$meanHA[ind2],
                           meanHP=regulators_HP_CVP$meanHP[ind1],
                           log2FoldChangeCVP_HA=(-regulators_HA_CVP$log2FoldChange[ind2]),
                           log2FoldChangeCVP_HP=(-regulators_HP_CVP$log2FoldChange[ind1]),           
                           deseq_padj_CVP_HA=regulators_HA_CVP$deseq_padj[ind2],
                           deseq_padj_CVP_HP=regulators_HP_CVP$deseq_padj[ind1],
                           K27me3_log2FoldChange_CVP_HA=(-regulators_HA_CVP$K27me3_log2FoldChange[ind2]),
                           K27me3_log2FoldChange_CVP_HP=(-regulators_HP_CVP$K27me3_log2FoldChange[ind1]),
                           K27me3_padj_CVP_HA=regulators_HA_CVP$K27me3_padj[ind2],
                           K27me3_padj_CVP_HP=regulators_HP_CVP$K27me3_padj[ind1],
                           K4me3_log2FoldChange_CVP_HA=(-regulators_HA_CVP$K4me3_log2FoldChange[ind2]),
                           K4me3_log2FoldChange_CVP_HP=(-regulators_HP_CVP$K4me3_log2FoldChange[ind1]),
                           K4me3_padj_CVP_HA=regulators_HA_CVP$K4me3_padj[ind2],
                           K4me3_padj_CVP_HP=regulators_HP_CVP$K4me3_padj[ind1],
                           score_CVP_HA=-regulators_HA_CVP$score[ind2],
                           score_CVP_HP=-regulators_HP_CVP$score[ind1],
                           score_aggregate=-regulators_HA_CVP$score[ind2]-regulators_HP_CVP$score[ind1],
                           Function=regulators_HA_CVP$Function[ind2]
)

regulators_CVP<-regulators_CVP[order(-regulators_CVP$score_aggregate),]

ind1<-match(common_HP,regulators_HP_CVP$gene_name)
ind2<-match(common_HP,regulators_HP_HA$gene_name)
regulators_HP<-data.frame(gene_name=common_HP,
                          meanHP=regulators_HP_CVP$meanHP[ind1],
                          meanHA=regulators_HP_HA$meanHA[ind2],
                          meanCVP=regulators_HP_CVP$meanCVP[ind1],
                          log2FoldChangeHP_HA=regulators_HP_HA$log2FoldChange[ind2],
                          log2FoldChangeHP_CVP=regulators_HP_CVP$log2FoldChange[ind1],           
                          deseq_padj_HP_HA=regulators_HP_HA$deseq_padj[ind2],
                          deseq_padj_HP_CVP=regulators_HP_CVP$deseq_padj[ind1],
                          K27me3_log2FoldChange_HP_HA=regulators_HP_HA$K27me3_log2FoldChange[ind2],                         
                          K27me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K27me3_log2FoldChange[ind1],
                          K27me3_padj_HP_HA=regulators_HP_HA$K27me3_padj[ind2],
                          K27me3_padj_HP_CVP=regulators_HP_CVP$K27me3_padj[ind1],
                          K4me3_log2FoldChange_HP_HA=regulators_HP_HA$K4me3_log2FoldChange[ind2],                         
                          K4me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K4me3_log2FoldChange[ind1],
                          K4me3_padj_HP_HA=regulators_HP_HA$K4me3_padj[ind2],
                          K4me3_padj_HP_CVP=regulators_HP_CVP$K4me3_padj[ind1],
                          score_HP_HA=regulators_HP_HA$score[ind2],
                          score_HP_CVP=regulators_HP_CVP$score[ind1],
                          score_aggregate=regulators_HP_HA$score[ind2]+regulators_HP_CVP$score[ind1],
                          Function=regulators_HP_CVP$Function[ind1]
)

regulators_HP<-regulators_HP[order(-regulators_HP$score_aggregate),]

ind1<-match(common_HA,regulators_HA_CVP$gene_name)
ind2<-match(common_HA,regulators_HP_HA$gene_name)
regulators_HA<-data.frame(gene_name=common_HA,
                          meanHA=regulators_HP_HA$meanHA[ind2],                         
                          meanHP=regulators_HP_HA$meanHP[ind2],
                          meanCVP=regulators_HA_CVP$meanCVP[ind1],
                          log2FoldChangeHA_HP=-regulators_HP_HA$log2FoldChange[ind2],
                          log2FoldChangeHA_CVP=regulators_HA_CVP$log2FoldChange[ind1],           
                          deseq_padj_HA_HP=regulators_HP_HA$deseq_padj[ind2],
                          deseq_padj_HA_CVP=regulators_HA_CVP$deseq_padj[ind1],
                          K27me3_log2FoldChange_HA_HP=-regulators_HP_HA$K27me3_log2FoldChange[ind2],                         
                          K27me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K27me3_log2FoldChange[ind1],
                          K27me3_padj_HA_HP=regulators_HP_HA$K27me3_padj[ind2],
                          K27me3_padj_HA_CVP=regulators_HA_CVP$K27me3_padj[ind1],
                          K4me3_log2FoldChange_HA_HP=-regulators_HP_HA$K4me3_log2FoldChange[ind2],                         
                          K4me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K4me3_log2FoldChange[ind1],
                          K4me3_padj_HA_HP=regulators_HP_HA$K4me3_padj[ind2],
                          K4me3_padj_HA_CVP=regulators_HA_CVP$K4me3_padj[ind1],
                          score_HA_HP=-regulators_HP_HA$score[ind2],
                          score_HA_CVP=regulators_HA_CVP$score[ind1],
                          score_aggregate=-regulators_HP_HA$score[ind2]+regulators_HA_CVP$score[ind1],
                          Function=regulators_HA_CVP$Function[ind1]
)

regulators_HA<-regulators_HA[order(-regulators_HA$score_aggregate),]

ind1<-match(common_HAHP,regulators_HA_CVP$gene_name)
ind2<-match(common_HAHP,regulators_HP_CVP$gene_name)
regulators_HAHP<-data.frame(gene_name=common_HAHP,
                            meanHA=regulators_HA_CVP$meanHA[ind1],                         
                            meanHP=regulators_HP_CVP$meanHP[ind2],
                            meanCVP=regulators_HA_CVP$meanCVP[ind1],
                            log2FoldChangeHA_CVP=regulators_HA_CVP$log2FoldChange[ind1],
                            log2FoldChangeHP_CVP=regulators_HP_CVP$log2FoldChange[ind2],           
                            deseq_padj_HA_CVP=regulators_HA_CVP$deseq_padj[ind1],
                            deseq_padj_HP_CVP=regulators_HP_CVP$deseq_padj[ind2],
                            K27me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K27me3_log2FoldChange[ind1],                         
                            K27me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K27me3_log2FoldChange[ind2],
                            K27me3_padj_HA_CVP=regulators_HA_CVP$K27me3_padj[ind1],
                            K27me3_padj_HP_CVP=regulators_HP_CVP$K27me3_padj[ind2],
                            K4me3_log2FoldChange_HA_CVP=regulators_HA_CVP$K4me3_log2FoldChange[ind1],                         
                            K4me3_log2FoldChange_HP_CVP=regulators_HP_CVP$K4me3_log2FoldChange[ind2],
                            K4me3_padj_HA_CVP=regulators_HA_CVP$K4me3_padj[ind1],
                            K4me3_padj_HP_CVP=regulators_HP_CVP$K4me3_padj[ind2],
                            score_HA_CVP=regulators_HA_CVP$score[ind1],
                            score_HP_CVP=regulators_HP_CVP$score[ind2],
                            score_aggregate=regulators_HA_CVP$score[ind1]+regulators_HP_CVP$score[ind2],
                            Function=regulators_HA_CVP$Function[ind1]
)

regulators_HAHP<-regulators_HAHP[order(-regulators_HAHP$score_aggregate),]

ind1<-match(common_HACVP,regulators_HP_HA$gene_name)
ind2<-match(common_HACVP,regulators_HP_CVP$gene_name)
regulators_HACVP<-data.frame(gene_name=common_HACVP,
                             meanHA=regulators_HP_HA$meanHA[ind1], 
                             meanCVP=regulators_HP_CVP$meanCVP[ind2],                           
                             meanHP=regulators_HP_HA$meanHP[ind1],
                             log2FoldChangeHA_HP=-regulators_HP_HA$log2FoldChange[ind1],
                             log2FoldChangeCVP_HP=-regulators_HP_CVP$log2FoldChange[ind2],           
                             deseq_padj_HA_HP=regulators_HP_HA$deseq_padj[ind1],
                             deseq_padj_CVP_HP=regulators_HP_CVP$deseq_padj[ind2],
                             K27me3_log2FoldChange_HA_HP=-regulators_HP_HA$K27me3_log2FoldChange[ind1],                         
                             K27me3_log2FoldChange_CVP_HP=-regulators_HP_CVP$K27me3_log2FoldChange[ind2],
                             K27me3_padj_HA_HP=regulators_HP_HA$K27me3_padj[ind1],
                             K27me3_padj_CVP_HP=regulators_HP_CVP$K27me3_padj[ind2],
                             K4me3_log2FoldChange_HA_HP=-regulators_HP_HA$K4me3_log2FoldChange[ind1],                         
                             K4me3_log2FoldChange_CVP_HP=-regulators_HP_CVP$K4me3_log2FoldChange[ind2],
                             K4me3_padj_HA_HP=regulators_HA_CVP$K4me3_padj[ind1],
                             K4me3_padj_CVP_HP=regulators_HP_CVP$K4me3_padj[ind2],
                             score_HA_HP=-regulators_HP_HA$score[ind1],
                             score_CVP_HP=-regulators_HP_CVP$score[ind2],
                             score_aggregate=-regulators_HP_HA$score[ind1]-regulators_HP_CVP$score[ind2],
                             Function=regulators_HP_CVP$Function[ind2]
)

regulators_HACVP<-regulators_HACVP[order(-regulators_HACVP$score_aggregate),]

regulators_HACVP_exclude_indivi<-regulators_HACVP[!(regulators_HACVP$gene_name %in% regulators_HA$gene_name |regulators_HACVP$gene_name %in% regulators_CVP$gene_name),]
regulators_HAHP_exclude_indivi<-regulators_HAHP[!(regulators_HAHP$gene_name %in% regulators_HA$gene_name |regulators_HAHP$gene_name %in% regulators_HP$gene_name),]

#HA breakdown
ind1<-match(regulators_HA$gene_name,regulators_HA_CVP$gene_name)
ind2<-match(regulators_HA$gene_name,regulators_HP_HA$gene_name)
breakdown_HA<-data.frame(gene_name=regulators_HA$gene_name,exprs.HAvsCVP=-log10(regulators_HA_CVP$deseq_padj[ind1]),
                         exprs.HAvsHP=-log10(regulators_HP_HA$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HA_CVP$deseq_padj[ind1])-log10(regulators_HP_HA$deseq_padj[ind2]))/2,
                         K4me3.HAvsCVP=-log10(regulators_HA_CVP$K4me3_padj[ind1]),K4me3.HAvsHP=-log10(regulators_HP_HA$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HA_CVP$K4me3_padj[ind1])-log10(regulators_HP_HA$K4me3_padj[ind2]))/2,
                         K27me3.HAvsCVP=-log10(regulators_HA_CVP$K27me3_padj[ind1]),K27me3.HAvsHP=-log10(regulators_HP_HA$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HA_CVP$K27me3_padj[ind1])-log10(regulators_HP_HA$K27me3_padj[ind2]))/2,
                         final_score=regulators_HA$score_aggregate)
#CVP breakdown
ind1<-match(regulators_CVP$gene_name,regulators_HA_CVP$gene_name)
ind2<-match(regulators_CVP$gene_name,regulators_HP_CVP$gene_name)
breakdown_CVP<-data.frame(gene_name=regulators_CVP$gene_name,exprs.CVPvsHA=-log10(regulators_HA_CVP$deseq_padj[ind1]),
                          exprs.CVPvsHP=-log10(regulators_HP_CVP$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HA_CVP$deseq_padj[ind1])-log10(regulators_HP_CVP$deseq_padj[ind2]))/2,
                          K4me3.CVPvsHA=-log10(regulators_HA_CVP$K4me3_padj[ind1]),K4me3.CVPvsHP=-log10(regulators_HP_CVP$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HA_CVP$K4me3_padj[ind1])-log10(regulators_HP_CVP$K4me3_padj[ind2]))/2,
                          K27me3.CVPvsHA=-log10(regulators_HA_CVP$K27me3_padj[ind1]),K27me3.CVPvsHP=-log10(regulators_HP_CVP$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HA_CVP$K27me3_padj[ind1])-log10(regulators_HP_CVP$K27me3_padj[ind2]))/2,
                          final_score=regulators_CVP$score_aggregate)
#HP breakdown
ind2<-match(regulators_HP$gene_name,regulators_HP_HA$gene_name)
ind1<-match(regulators_HP$gene_name,regulators_HP_CVP$gene_name)
breakdown_HP<-data.frame(gene_name=regulators_HP$gene_name,exprs.HPvsCVP=-log10(regulators_HP_CVP$deseq_padj[ind1]),
                         exprs.HPvsHA=-log10(regulators_HP_HA$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HP_CVP$deseq_padj[ind1])-log10(regulators_HP_HA$deseq_padj[ind2]))/2,
                         K4me3.HPvsCVP=-log10(regulators_HP_CVP$K4me3_padj[ind1]),K4me3.HPvsHA=-log10(regulators_HP_HA$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HP_CVP$K4me3_padj[ind1])-log10(regulators_HP_HA$K4me3_padj[ind2]))/2,
                         K27me3.HPvsCVP=-log10(regulators_HP_CVP$K27me3_padj[ind1]),K27me3.HPvsHA=-log10(regulators_HP_HA$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HP_CVP$K27me3_padj[ind1])-log10(regulators_HP_HA$K27me3_padj[ind2]))/2,
                         final_score=regulators_HP$score_aggregate)

#HAHP breakdown
ind1<-match(regulators_HAHP_exclude_indivi$gene_name,regulators_HA_CVP$gene_name)
ind2<-match(regulators_HAHP_exclude_indivi$gene_name,regulators_HP_CVP$gene_name)
breakdown_HAHP<-data.frame(gene_name=regulators_HAHP_exclude_indivi$gene_name,exprs.CVPvsHA=-log10(regulators_HA_CVP$deseq_padj[ind1]),
                           exprs.CVPvsHP=-log10(regulators_HP_CVP$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HA_CVP$deseq_padj[ind1])-log10(regulators_HP_CVP$deseq_padj[ind2]))/2,
                           K4me3.CVPvsHA=-log10(regulators_HA_CVP$K4me3_padj[ind1]),K4me3.CVPvsHP=-log10(regulators_HP_CVP$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HA_CVP$K4me3_padj[ind1])-log10(regulators_HP_CVP$K4me3_padj[ind2]))/2,
                           K27me3.CVPvsHA=-log10(regulators_HA_CVP$K27me3_padj[ind1]),K27me3.CVPvsHP=-log10(regulators_HP_CVP$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HA_CVP$K27me3_padj[ind1])-log10(regulators_HP_CVP$K27me3_padj[ind2]))/2,
                           final_score=regulators_HAHP_exclude_indivi$score_aggregate)

#HACVP breakdown
ind2<-match(regulators_HACVP_exclude_indivi$gene_name,regulators_HP_HA$gene_name)
ind1<-match(regulators_HACVP_exclude_indivi$gene_name,regulators_HP_CVP$gene_name)
breakdown_HACVP<-data.frame(gene_name=regulators_HACVP_exclude_indivi$gene_name,exprs.HPvsCVP=-log10(regulators_HP_CVP$deseq_padj[ind1]),
                            exprs.HPvsHA=-log10(regulators_HP_HA$deseq_padj[ind2]),exprs.mean=(-log10(regulators_HP_CVP$deseq_padj[ind1])-log10(regulators_HP_HA$deseq_padj[ind2]))/2,
                            K4me3.HPvsCVP=-log10(regulators_HP_CVP$K4me3_padj[ind1]),K4me3.HPvsHA=-log10(regulators_HP_HA$K4me3_padj[ind2]),K4me3.mean=(-log10(regulators_HP_CVP$K4me3_padj[ind1])-log10(regulators_HP_HA$K4me3_padj[ind2]))/2,
                            K27me3.HPvsCVP=-log10(regulators_HP_CVP$K27me3_padj[ind1]),K27me3.HPvsHA=-log10(regulators_HP_HA$K27me3_padj[ind2]),K27me3.mean=(-log10(regulators_HP_CVP$K27me3_padj[ind1])-log10(regulators_HP_HA$K27me3_padj[ind2]))/2,
                            final_score=regulators_HACVP_exclude_indivi$score_aggregate)

