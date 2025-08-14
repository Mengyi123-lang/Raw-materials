
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
  dir.create("00_pre_data")
}

library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)

library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group'){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot()+
    scale_fill_manual(values =group_cols)+   
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) 
  return(p)
}
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  

  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) 
  return(geneList)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  # set.seed(2024)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}

genecode=read.delim('origin_datas/TCGA/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)


tcga.pancancer.cli=read.xlsx('origin_datas/TCGA/TCGA_pancancer_cli_PMID_29625055.xlsx')
head(tcga.pancancer.cli)
table(tcga.pancancer.cli$type)
tcga.cli=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='OV'),]
head(tcga.cli)

tcga.cli<-read.delim('origin_datas/TCGA/Merge_OV_clinical.txt',sep='\t',header = T)
head(tcga.cli)
tcga.cli=data.frame(Samples=paste0(tcga.cli$A0_Samples,'-01'),
                    Age=tcga.cli$A17_Age, Grade=tcga.cli$A7_Grade,Stage=tcga.cli$A6_Stage,
                    Status=tcga.cli$A2_Event,OS.time=tcga.cli$A1_OS)
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
tcga.cli$OS.time
tcga.cli=tcga.cli %>% drop_na(OS.time)
tcga.cli=tcga.cli[tcga.cli$OS.time>0,]
tcga.cli$OS=ifelse(tcga.cli$Status=='Alive',0,1)

dim(tcga.cli)
table(tcga.cli$Stage)
tcga.cli$Stage[tcga.cli$Stage=='']=NA
tcga.cli$Stage=gsub('Stage ','',tcga.cli$Stage)
tcga.cli$Stage=gsub('[ABC]','',tcga.cli$Stage)

table(tcga.cli$Grade)
tcga.cli$Grade[tcga.cli$Grade%in%c('GB','GX','Not Available')]=NA
dim(tcga.cli)


tcga.data=read.delim('origin_datas/TCGA/TCGA_OV_TPM.txt',row.names = 1,check.names = F)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))
dim(tcga.data)

sample_T=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==1)]#肿瘤样本
sample_T=intersect(sample_T,tcga.cli$Samples)

tcga.exp=tcga.data[rownames(tcga.data) %in% mrna_genecode$SYMBOL,sample_T]
range(tcga.exp)
tcga.exp=log2(tcga.exp+1)
dim(tcga.exp)

tcga.cli=tcga.cli[sample_T,]
head(tcga.cli)
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>59,'>59','<=59')

save(tcga.exp,file = '00_pre_data/tcga.exp.RData')

#GSE32062#####
GSE32062=getGEOExpData('GSE32062')
save(GSE32062,file = 'origin_datas/GEO/GSE32062.RData')
load('origin_datas/GEO/GSE32062.RData')
GSE32062.cli=GSE32062$Sample
head(GSE32062.cli)
table(GSE32062.cli$tissue)
GSE32062.cli=data.frame(Samples=GSE32062.cli$Acc,Grade=paste0('G',GSE32062.cli$grading),Stage=GSE32062.cli$Stage,
                        OS=GSE32062.cli$`death (1)`,OS.time=GSE32062.cli$`os (m)`)
GSE32062.cli=GSE32062.cli%>%drop_na(OS)
rownames(GSE32062.cli)=GSE32062.cli$Samples
head(GSE32062.cli)
dim(GSE32062.cli)
# 
GSE32062.exp=GSE32062$Exp$GPL6480_41093_Data_col1
GSE32062.exp[1:5,1:5]
dim(GSE32062.exp)
GSE32062.exp=exp_probe2symbol_v2(datExpr = GSE32062.exp,GPL = 'GPL6480')
range(GSE32062.exp)
#GSE32062.exp=log2(2^GSE32062.exp+1)
#01.TRP 
dir.create('01_TRP_gene')
TRP.geneset=read.csv("01_TRP_gene/trp_gene_ PMID37986367.csv",header = F)
TRP.geneset <- TRP.geneset$V1
TRP.geneset <- TRP.geneset[TRP.geneset!=""]
length(TRP.geneset)
#33


TRP.cox=cox_batch(dat = tcga.exp[TRP.geneset,tcga.cli$Samples],
                   time = tcga.cli$OS.time,event = tcga.cli$OS)
TRP.cox=na.omit(TRP.cox)
head(TRP.cox)
rownames(TRP.cox)=gsub('-','_',rownames(TRP.cox))
p_cutoff=0.05
table(TRP.cox$p.value<p_cutoff)
# FALSE  TRUE 
# 116    16 
TRP.cox.fit=TRP.cox[TRP.cox$p.value<p_cutoff,]
dim(TRP.cox.fit)
write.csv(TRP.cox.fit,'01_TRP_gene/TRP.cox.fit.csv')

bioForest=function(rt=null,col){

  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  

  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}

# pdf('01_TRP_gene/TRP_forest.pdf',height = 10,width = 8,onefile = F)
# bioForest(rt = TRP.cox.fit[order(TRP.cox.fit$HR),],col=c('blue','#F05A7E'))
# dev.off()

tcga.maf=getTCGAMAFByCode('OV')
pdf('01_TRP_gene/TRP_MAF.pdf',height = 10,width = 12,onefile = F)
oncoplot(maf = tcga.maf,genes = TRP.geneset)
dev.off()
##TRP.ssgsea####
library(survival)
library(survminer)
TRP.score <- t(ssGSEAScore_by_genes(gene.exp = tcga.exp,genes =rownames(TRP.cox )))

dir.create('02.WGCNA')

library(WGCNA)
allowWGCNAThreads(nThreads = 36)
enableWGCNAThreads(nThreads = 36)

my_mad <- function(x){mad(x,na.rm = TRUE)} 
wgcna_exp=t(tcga.exp)
m.mad <- apply(wgcna_exp,2,my_mad)
dim(tcga.exp)
# 19503   373

tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]


range(tpm_T2)

pdf('02.WGCNA/1.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()


tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=3,
                                 mergeCutHeight=0.3,
                                 minModuleSize=60)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('02.WGCNA/2.pdf',height = 5,width = 12)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
MODULE.gene.num <- as.data.frame(table(tpm_T2.module$Modules[,2]))
write.table(MODULE.gene.num,file = "02.WGCNA/MODULE.gene.num.txt",sep = "\t",quote = F,row.names =F) 
writeMatrix(tpm_T2.module$Modules,outpath = '02.WGCNA/tcga.wgcna.module.genes.txt')
pdf('02.WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()


# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Risktype module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('02.WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Risktypeing of module eigengenes",xlab = "", sub = "")
dev.off()

tcga_cli_use <-data.frame(TRP.score=TRP.score[tcga.cli$Samples,])
head(tcga_cli_use)
tcga_cli_use.part=tcga_cli_use
str(tcga_cli_use.part)

tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))
spms=tcga_cli_use.part

MEs_col<-tpm_T2.module$MEs
dim(MEs_col)

modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('02.WGCNA/5.pdf',width =8,height = 6)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,xLabelsAngle = 0,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))

geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "brown"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
tcga.wgcna.genes=names(which(moduleGenes))
length(tcga.wgcna.genes)
#1050
write.table(tcga.wgcna.genes,"04_model/tcga.wgcna.genes.txt",sep = "\t",col.names = F,row.names = F)
pdf('02.WGCNA/6.pdf',height = 8,width = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, "TRP.score"]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TRP.score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
dev.off()


dir.create('03_wgcna_module_gsea')

cox.gene.enrichment=mg_clusterProfiler(genes =tcga.wgcna.genes )
fig3a=list()
fig3a[[1]]=enrichplot::dotplot(cox.gene.enrichment$KEGG)+ggtitle('KEGG')
fig3a[[2]]=enrichplot::dotplot(cox.gene.enrichment$GO_BP)+ggtitle('Biological Process')
fig3a[[3]]=enrichplot::dotplot(cox.gene.enrichment$GO_CC)+ggtitle('Cellular Component')
fig3a[[4]]=enrichplot::dotplot(cox.gene.enrichment$GO_MF)+ggtitle('Molecular Function')
# 
# fig3=mg_merge_plot(fig3a[[1]],fig3a[[2]],fig3a[[4]],labels = c('A','B','C'),ncol=3,
# nrow=1)
fig3=mg_merge_plot(mg_merge_plot(fig3a[[1]],fig3a[[2]],labels = c('A','B'),ncol=2),
                   mg_merge_plot(fig3a[[3]],fig3a[[4]],ncol=2,nrow=1,labels = LETTERS[3:4],heights = c(1,1),widths = c(1,1)),
                   nrow=2,heights = c(1,1))
savePDF('03_wgcna_module_gsea/Fig3.pdf',fig3,height = 12,width = 18)


dir.create('04_model')
###

tcga.wgcna.genes<- gsub('-', '_', tcga.wgcna.genes)
TRP.geneset
intersect(TRP.geneset,rownames(tcga.exp))
comm_gene <- Reduce(intersect,list(
  
  tcga.wgcna.genes,
  rownames(GSE32062.exp),
  #rownames(GSE51088_exp)
  #rownames(GSE26712_exp)
  #rownames(GSE49997_exp)
 rownames(GSE14764_exp) ,
 #rownames(GSE17260_exp)
 #rownames(GSE26193_exp)
 rownames(GSE30161_exp)
  ))
length(comm_gene)
##COR####
set.seed(123)
tcga_trp_cor <- Hmisc::rcorr(t(tcga.exp[intersect(TRP.geneset,rownames(tcga.exp)), tcga.cli$Samples]),
                             t(tcga.exp[comm_gene, tcga.cli$Samples]))

tcga_trp_cor_r <- reshape2::melt(tcga_trp_cor$r[intersect(TRP.geneset,rownames(tcga.exp)), comm_gene])
tcga_trp_cor_r <- as.data.frame(tcga_trp_cor_r)
colnames(tcga_trp_cor_r) <- c('TRP', 'WGCNA', 'Corr')

tcga_trp_cor_p <- reshape2::melt(tcga_trp_cor$P[intersect(TRP.geneset,rownames(tcga.exp)), comm_gene])
tcga_trp_cor_p <- as.data.frame(tcga_trp_cor_p)
colnames(tcga_trp_cor_p) <- c('TRP', 'WGCNA', 'P.value')

tcga_trp_cor <- tcga_trp_cor_r
tcga_trp_cor$P.value <- tcga_trp_cor_p$P.value

tcga_trp_cor_filtered <- tcga_trp_cor[tcga_trp_cor$P.value < 0.01 & tcga_trp_cor$Corr > 0.5, ]
dim(tcga_trp_cor_filtered)
table(tcga_trp_cor_filtered$Corr > 0)
tcga_positive_genes <- as.character(tcga_trp_cor_filtered$PCG[tcga_trp_cor_filtered$Corr > 0])
length(unique(tcga_trp_cor_filtered $WGCNA))
write.csv(tcga_trp_cor_filtered[tcga_trp_cor_filtered$Corr > 0, ],
          file = '04_model/tcga_trp_cor_filtered.csv')


sig.gene.cox=cox_batch(dat = tcga.exp[intersect(unique(tcga_trp_cor_filtered $WGCNA),rownames(tcga.exp)),tcga.cli$Samples],
                       time = tcga.cli$OS.time,event = tcga.cli$OS)
sig.gene.cox
table(sig.gene.cox$p.value<0.05)
# FALSE  TRUE 
#  419    51 
pre.genes=rownames(sig.gene.cox[sig.gene.cox$p.value<0.05,])
length(pre.genes)#51
pre.genes <- gsub('-', '_', pre.genes)
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
write.csv(sig.gene.cox,'04_model/sig.cox.csv')

tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[,-c(1:2)],
                               os = tcga_model_data$OS,
                               os.time = tcga_model_data$OS.time)
length(tcga.lasso$lasso.gene)#18
tcga.lasso$plot

fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(tcga.lasso$lasso.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
"-0.391*CD40LG+0.233*HGF+0.165*NAAA+-0.262*S1PR4+-0.275*STAT4+0.187*STAT5A+-0.195*TAP1+0.191*TRPM2+0.195*VSIG4"
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)


rownames(tcga.exp) <- gsub('-', '_', rownames(tcga.exp))
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[names(lan), tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
fig_ggforest <- survminer::ggforest(cox,data=tcga_model_data)
ggsave("04_model/ggforest.pdf",fig_ggforest ,height =4,width =7 )


lan.dataframe <- as.data.frame(lan)
lan.dataframe$gene <- rownames(lan.dataframe) 
lan.dataframe$gene <- factor(lan.dataframe$gene,levels = rownames(lan.dataframe)[order(lan.dataframe$lan)])

lan.dataframe$color_group <- ifelse(lan.dataframe$lan > 0, "Positive", "Negative")
library(ggplot2)

p <- ggplot(lan.dataframe, aes(x=gene, y=lan,fill=color_group)) +
  geom_bar(stat="identity") +
  xlab("Gene Name") +
  ylab("Coefficient") +
  ggtitle("Gene Coefficients") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#FCF596", "Negative" = "#795458")) +
  theme_bw()+
  guides(fill=FALSE)
p1 <- p+geom_text(aes(label=sprintf("%.3f", lan)), hjust=-0.2, size=3, color="black")
ggsave('04_model/gene_ Coefficients.pdf',p1,height = 4,width = 9)

risktype.col=c('#DE8F5F',"#7ED4AD")

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)

# tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(tcga.risktype.cli$Riskscore),'High','Low')
tcga.cutoff<-surv_cutpoint(tcga.risktype.cli,time="OS.time",event="OS",variables='Riskscore')
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff$cutpoint$cutpoint,'High','Low')

#####tcga.km.OS#########
tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
tcga.roc

pdf("04_model/tcga.ROC.pdf",height = 6,width = 6)
tcga.roc=mg_surv_pROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
dev.off()
tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='TCGA-OV',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,
                      ylab='Overall Survival(OS)',
                      legend=c(0.85,0.8),
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

tcga.km.OS

tcga.risktype.cli$Status=ifelse(tcga.risktype.cli$OS==0,'Alive','Dead')
tcga.model.p=my_riskplot(cli_dat = tcga.risktype.cli,cols =risktype.col,xlab = 'sample',
                         a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = tcga.cutoff$cutpoint$cutpoint,labs = '')

# tcga.gene.expr <- mg_PlotMutiBoxplot(data = t(tcga.exp[names(lan),tcga.risktype.cli$Samples]),group = tcga.risktype.cli$Risktype,group_cols = risktype.col,add = "box",legend.pos = "top")
# tcga.gene.expr <- tcga.gene.expr+ggtitle(label = "TCGA")

tcga.risktype.cli$Status <- ifelse(tcga.risktype.cli$OS==0,"alive","dead")
tcga.barplot=my_mutibarplot(df=table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')

####GSE32062 ok##########

names(lan)
intersect(names(lan),rownames(GSE32062.exp))
GSE32062_model_data <- cbind(GSE32062.cli[, c("OS.time", "OS")],
                             t(GSE32062.exp[names(lan), GSE32062.cli$Samples]))
colnames(GSE32062_model_data) <- gsub('-', '_', colnames(GSE32062_model_data))

risk.GSE32062=as.numeric(lan%*%as.matrix(t(GSE32062_model_data[GSE32062.cli$Samples,names(lan)])))
GSE32062.risktype.cli=data.frame(GSE32062.cli,Riskscore=risk.GSE32062)
GSE32062.cutoff<-surv_cutpoint(GSE32062.risktype.cli,time="OS.time",event="OS",variables='Riskscore')
GSE32062.risktype.cli$Risktype=ifelse(GSE32062.risktype.cli$Riskscore>GSE32062.cutoff$cutpoint$cutpoint,'High','Low')
# GSE32062.risktype.cli$Risktype=ifelse(GSE32062.risktype.cli$Riskscore>median(risk.GSE32062),'High','Low')

# GSE32062.roc=ggplotTimeROC(GSE32062.risktype.cli$OS.time,
#                            GSE32062.risktype.cli$OS,
#                            GSE32062.risktype.cli$Riskscore,mks = c(1:5))
# GSE32062.roc
pdf("04_model/GSE32062.ROC.pdf",height = 6,width = 6)
GSE32062.roc=mg_surv_pROC(GSE32062.risktype.cli$OS.time,
                          GSE32062.risktype.cli$OS,
                          GSE32062.risktype.cli$Riskscore,mks = c(1:5))

dev.off()
GSE32062.km.OS=ggsurvplot(fit=survfit( Surv(OS.time, OS) ~ Risktype,
                                       data = GSE32062.risktype.cli),
                          data=GSE32062.risktype.cli,
                          conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                          surv.median.line = 'hv',title='GSE32062',
                          linetype = c("solid", "dashed","strata")[1],
                          palette = risktype.col,ggtheme = custom_theme(),
                          legend = c(0.8,0.85), 
                          legend.title = "Risktype",legend.labs=c('High','Low'))
GSE32062.km.OS=mg_merge_plot(GSE32062.km.OS$plot,GSE32062.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE32062.km.OS

my_mutibarplot=function(df,xlab='group',leg.title='',cols=pal_d3()(10)[5:6]){
  prop.pval=round(chisq.test(df)$p.value,2)#round(-log10(chisq.test(df)$p.value),2)
  if( prop.pval<0.001)
    prop.pval='<0.001'
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('pvalue  ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
  p
  return(p)
}
GSE32062.risktype.cli$Status <- ifelse(GSE32062.risktype.cli$OS==0,"alive","dead")
GSE32062.barplot=my_mutibarplot(df=table(GSE32062.risktype.cli$Status,GSE32062.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')



fig4=mg_merge_plot(mg_merge_plot(tcga.lasso$plot,p1,widths = c(2,1),labels = c('','C')),
                   mg_merge_plot(tcga.km.OS,tcga.roc,tcga.barplot,widths = c(1,1,1),labels = LETTERS[4:6],ncol = 3),
                   
                   mg_merge_plot(GSE32062.km.OS,GSE32062.roc,GSE32062.barplot,ncol=3,labels = LETTERS[7:9]),
                   nrow=3)

ggsave('04_model/Fig4.pdf',fig4,height = 18,width = 18)



#5.GSEA###########
dir.create('05_GSEA')
tcga.geneList=getGeneFC(gene.exp=tcga.exp,group=tcga.risktype.cli$Risktype
                        ,ulab='High',dlab = 'Low')

library(clusterProfiler)
tcga.gsea.KEGG=gseKEGG(
  tcga.geneList,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea")
saveRDS(tcga.gsea.KEGG,file = "05_GSEA/tcga.gsea.KEGG.RDS")
tcga.gsea.KEGG <- readRDS("05_GSEA/tcga.gsea.KEGG.RDS")

tcga.gsea.KEGG.res=tcga.gsea.KEGG@result
write.csv(tcga.gsea.KEGG.res,'05_GSEA/GSEA_res.csv',row.names = F)
table(tcga.gsea.KEGG.res$p.adjust<0.05 & tcga.gsea.KEGG.res$NES<0)
table(tcga.gsea.KEGG.res$p.adjust<0.05 & tcga.gsea.KEGG.res$NES>0)
library(dplyr)
ind1=tcga.gsea.KEGG.res %>% slice_max(n =5, order_by = NES)
ind2=tcga.gsea.KEGG.res %>% slice_min(n =5, order_by = NES)  


p5a=enrichplot::gseaplot2(tcga.gsea.KEGG,ind1$ID, pvalue_table = F,title ='KEGG enrichment in High group')
p5b=enrichplot::gseaplot2(tcga.gsea.KEGG,ind2$ID, pvalue_table = F,title ='KEGG enrichment in Low group')
gsea.p=mg_merge_plot(p5a,p5b,nrow = 2,labels = c('A','B'))
ggsave('05_GSEA/risktype_GSEA.pdf',gsea.p,height = 15,width = 12)



dir.create('06_nomogram')
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)
table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'
table(tcga_cox_datas$Grade)
tcga_cox_datas$Grade[which(tcga_cox_datas$Grade=='G1' |tcga_cox_datas$Grade=='G2')]='G1+G2'
tcga_cox_datas$Grade[which(tcga_cox_datas$Grade=='G3' |tcga_cox_datas$Grade=='G4')]='G3+G4'
#Age
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age1,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],4))
Age_sig_cox_dat


table(tcga_cox_datas$Stage)
#Stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],4))
Stage_sig_cox_dat

#Grade
Grade_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Grade,
                               data=tcga_cox_datas))
Grade_sig_cox_dat <- data.frame(Names=rownames(Grade_sig_cox[[8]]),
                                HR = round(Grade_sig_cox[[7]][,2],3),
                                lower.95 = round(Grade_sig_cox[[8]][,3],3),
                                upper.95 = round(Grade_sig_cox[[8]][,4],3),
                                p.value=round(Grade_sig_cox[[7]][,5],4))
Grade_sig_cox_dat

#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],4))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     # T.stage_sig_cox_dat,
                     # N.stage_sig_cox_dat,
                     # M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        # "T.stage",
                        # "N.stage",
                        # "M.stage",
                        "Stage",
                        "Grade",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value <- ifelse(data.sig$p.value<0.0001,"<0.0001",data.sig$p.value)
pdf('06_nomogram/Fig6a.pdf', width = 6, height = 5,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =10,lineheight = 9,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='#009966',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 7,graph.pos = 2)
dev.off()



#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~ Age1+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c('AGE',
                         # 'M.stage',
                         # 'T.stage',
                         # "Stage",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti
data.muti$p.value <- ifelse(data.muti$p.value<0.001,"<0.001",data.muti$p.value)
pdf('06_nomogram/Fig6b.pdf', width = 6, height = 5,onefile = F)
#mg_forestplot_v2(data.muti,xlog = T,colgap = 8,lineheight = 10)
mg_forestplot_v2(data.muti,xlog = T,colgap =10,lineheight = 9,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1,xlim = c(-1,10)
                 ,box_col='#009966',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 6,graph.pos = 2)

dev.off()


############nomogran

pdf('06_nomogram/nomogran.pdf', width = 12, height = 5)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5))
dev.off()
#mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))



################AUC
library("timeROC")
tcga_cox_auc=tcga_cox_datas
tcga_cox_auc$Riskscore=as.numeric(tcga_cox_auc$Riskscore)
# tcga_cox_auc$T.stage=as.numeric(as.factor(tcga_cox_auc$T.stage))
# tcga_cox_auc$N.stage=as.numeric(as.factor(tcga_cox_auc$N.stage))
# tcga_cox_auc$M.stage=as.numeric(as.factor(tcga_cox_auc$M.stage))
tcga_cox_auc$Stage=as.numeric(as.factor(tcga_cox_auc$Stage))
tcga_cox_auc$Grade=as.numeric(as.factor(tcga_cox_auc$Grade))
tcga_cox_auc$Age=as.numeric(as.factor(tcga_cox_auc$Age))

head(tcga_cox_auc)
ROC.DSST.Age=timeROC(T=tcga_cox_auc$OS.time/365,
                     delta=tcga_cox_auc$OS,
                     marker=tcga_cox_auc$Age,
                     cause=1,weighting="marginal",
                     times=c(1,2,3,4,5),
                     iid=TRUE)
ROC.DSST.Age$AUC

ROC.DSST.stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Stage,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)
ROC.DSST.stage$AUC
ROC.DSST.Grade=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Grade,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)
ROC.DSST.Grade$AUC

ROC.DSST.Risk=timeROC(T=tcga_cox_auc$OS.time/365,
                      delta=tcga_cox_auc$OS,
                      marker=tcga_cox_auc$Riskscore,
                      cause=1,weighting="marginal",
                      times=c(1,2,3,4,5),
                      iid=TRUE)
ROC.DSST.Risk$AUC
ROC.DSST.Nomo<-timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Riskscore,
                       other_markers=as.matrix(tcga_cox_auc[,c("Stage")]),
                       cause=1,
                       weighting="cox",
                       times=c(1,2,3,4,5),
                       iid=F)
ROC.DSST.Nomo$AUC


library(ggsci)
pdf('06_nomogram/Fig5de.pdf',height = 5,width = 12)
par(mfrow=c(1,3))
plot(ROC.DSST.Age,time=1, col=pal_nejm()(8)[1], lwd=2, title = "",)   
# plot(ROC.DSST.T,time=1, col=pal_nejm()(8)[2], add=TRUE, lwd=2) 
# plot(ROC.DSST.N,time=1, col=pal_nejm()(8)[3], add=TRUE, lwd=2) 
# plot(ROC.DSST.M,time=1, col=pal_nejm()(8)[4], add=TRUE, lwd=2) 
plot(ROC.DSST.stage,time=1, col=pal_nejm()(8)[5], add=TRUE, lwd=2) 
plot(ROC.DSST.Grade,time=1, col=pal_nejm()(8)[6], add=TRUE, lwd=2) 
plot(ROC.DSST.Risk,time=1, col=pal_nejm()(8)[7], add=TRUE, lwd=2) 
plot(ROC.DSST.Nomo,time=1, col=pal_nejm()(8)[8], add=TRUE, lwd=2) 
legend("bottomright",
       c(paste0("Age AUC at 1 year: ",round(ROC.DSST.Age[["AUC"]][1],2)), 
         # paste0("T AUC at 1 year: ",round(ROC.DSST.T[["AUC"]][1],2)), 
         # paste0("N AUC at 1 year: ",round(ROC.DSST.N[["AUC"]][1],2)), 
         # paste0("M AUC at 1 year: ",round(ROC.DSST.M[["AUC"]][1],2)), 
         paste0("Stage AUC at 1 year: ",round(ROC.DSST.stage[["AUC"]][1],2)), 
         paste0("Grade AUC at 1 year: ",round(ROC.DSST.Grade[["AUC"]][1],2)), 
         paste0("RiskScore AUC at 1 year: ",round(ROC.DSST.Risk[["AUC"]][1],2)), 
         paste0("Nomogram AUC at 1 year: ",round(ROC.DSST.Nomo[["AUC"]][1],2))),
       col=pal_nejm()(8),
       lty=1, lwd=2,bty = "n",cex = 0.7)   


#pdf('results/08.nomogram/Fig8f.pdf',height = 5,width = 5)
plot(ROC.DSST.Age,time=3, col=pal_nejm()(8)[1], lwd=2, title = "",)   
# plot(ROC.DSST.T,time=3, col=pal_nejm()(8)[2], add=TRUE, lwd=2) 
# plot(ROC.DSST.N,time=3, col=pal_nejm()(8)[3], add=TRUE, lwd=2) 
# plot(ROC.DSST.M,time=3, col=pal_nejm()(8)[4], add=TRUE, lwd=2) 
plot(ROC.DSST.stage,time=3, col=pal_nejm()(8)[5], add=TRUE, lwd=2) 
plot(ROC.DSST.Grade,time=3, col=pal_nejm()(8)[6], add=TRUE, lwd=2) 
plot(ROC.DSST.Risk,time=3, col=pal_nejm()(8)[7], add=TRUE, lwd=2) 
plot(ROC.DSST.Nomo,time=3, col=pal_nejm()(8)[8], add=TRUE, lwd=2) 
legend("bottomright",
       c(paste0("Age AUC at 3 year: ",round(ROC.DSST.Age[["AUC"]][3],2)), 
         # paste0("T AUC at 3 year: ",round(ROC.DSST.T[["AUC"]][3],2)), 
         # paste0("N AUC at 3 year: ",round(ROC.DSST.N[["AUC"]][3],2)), 
         # paste0("M AUC at 3 year: ",round(ROC.DSST.M[["AUC"]][3],2)), 
         paste0("Stage AUC at 3 year: ",round(ROC.DSST.stage[["AUC"]][3],2)), 
         paste0("Grade AUC at 3 year: ",round(ROC.DSST.Grade[["AUC"]][3],2)), 
         paste0("RiskScore AUC at 3 year: ",round(ROC.DSST.Risk[["AUC"]][3],2)), 
         paste0("Nomogram AUC at 3 year: ",round(ROC.DSST.Nomo[["AUC"]][3],2))),
       col=pal_nejm()(8),
       lty=1, lwd=2,bty = "n",cex = 0.7)   


#pdf('results/08.nomogram/Fig8g.pdf',height = 5,width = 5)
plot(ROC.DSST.Age,time=5, col=pal_nejm()(8)[1], lwd=2, title = "",)   
# plot(ROC.DSST.T,time=5, col=pal_nejm()(8)[2], add=TRUE, lwd=2) 
# plot(ROC.DSST.N,time=5, col=pal_nejm()(8)[3], add=TRUE, lwd=2) 
# plot(ROC.DSST.M,time=5, col=pal_nejm()(8)[4], add=TRUE, lwd=2) 
plot(ROC.DSST.stage,time=5, col=pal_nejm()(8)[5], add=TRUE, lwd=2) 
plot(ROC.DSST.Grade,time=5, col=pal_nejm()(8)[6], add=TRUE, lwd=2) 
plot(ROC.DSST.Risk,time=5, col=pal_nejm()(8)[7], add=TRUE, lwd=2) 
plot(ROC.DSST.Nomo,time=5, col=pal_nejm()(8)[8], add=TRUE, lwd=2) 
legend("bottomright",
       c(paste0("Age AUC at 5 year: ",round(ROC.DSST.Age[["AUC"]][5],2)), 
         # paste0("T AUC at 5 year: ",round(ROC.DSST.T[["AUC"]][5],2)), 
         # paste0("N AUC at 5 year: ",round(ROC.DSST.N[["AUC"]][5],2)), 
         # paste0("M AUC at 5 year: ",round(ROC.DSST.M[["AUC"]][5],2)), 
         paste0("Stage AUC at 5 year: ",round(ROC.DSST.stage[["AUC"]][5],2)), 
         paste0("Grade AUC at 5 year: ",round(ROC.DSST.Grade[["AUC"]][5],2)), 
         paste0("RiskScore AUC at 5 year: ",round(ROC.DSST.Risk[["AUC"]][5],2)), 
         paste0("Nomogram AUC at 5 year: ",round(ROC.DSST.Nomo[["AUC"]][5],2))),
       col=pal_nejm()(8),
       lty=1, lwd=2,bty = "n",cex = 0.7)   
dev.off()

dir.create("07_immunity")
mycolor <- risktype.col

# TCGA
tcga.exp.icg=immu_ICGs(tcga.exp)
dim(tcga.exp.icg)
tcga_icg_plot <- mg_PlotMutiBoxplot(tcga.exp.icg[rownames(tcga.risktype.cli), ]
                                    , group = tcga.risktype.cli$Risktype
                                    , legend.pos = 'top'
                                    , add = 'boxplot'
                                    , xangle = 60
                                    , ylab = 'EXP'
                                    , group_cols = risktype.col
                                    #, test_method = 'kruskal.test'
)
tcga_icg_plot
ggsave("07_immunity/tcga_icg_plot.pdf",tcga_icg_plot,height = 6,width = 14)

library(ComplexHeatmap)
library(dplyr)

colnames(tcga.risktype.cli)
tcga.risktype.cli <- arrange(tcga.risktype.cli, Risktype)

tcga_ssgsea_icg <- t(scale(tcga.exp.icg[rownames(tcga.risktype.cli), ]))
#pdf('05_immunity/tcga_ssgsea_icg_plot.pdf', width = 8, height = 8)
tcga_ssgsea_icg_plot <- Heatmap(tcga_ssgsea_icg
                                , name = "score"
                                , col = circlize::colorRamp2(c(-2, 0, 2), c(mycolor[2], 'white', mycolor[1]))
                                , border = T
                                , show_column_names = F
                                , show_column_dend = F
                                , show_row_dend = F
                                , cluster_columns = T
                                , cluster_rows = T
                                , column_split = factor(tcga.risktype.cli$Risktype)
                                , clustering_distance_rows  ='pearson'
                                , clustering_method_rows = 'ward.D2'
                                # , column_km =2
                                , row_names_gp = gpar(fontsize = 10)
                                , top_annotation = HeatmapAnnotation(Risktype = tcga.risktype.cli$Risktype
                                                                     , col=list(Risktype=c('High'=risktype.col[1],
                                                                                           'Low' = risktype.col[2])
                                                                     )
                                                                     , annotation_width = unit(c(1,2), 'cm')
                                                                     , annotation_height = unit(0.2, "cm")
                                                                     , gap = unit(1, 'mm')))
tcga_ssgsea_icg_plot
dev.off()
# pdf('05_immunity/tcga_ssgsea_icg_pheatmat_plot.pdf', width = 8, height = 8)
# tcga_ssgsea_icg_plot
# dev.off()
##mcp####
tcga.mcp <- immu_MCPcounter(exp = tcga.exp)
Fig7b <- my_mutiboxplot(dat = tcga.mcp[tcga.risktype.cli$Samples,],
                        group = tcga.risktype.cli$Risktype,
                        ylab = 'Score',group_cols = risktype.col,labs(title = "MCP"))
ggsave("07_immunity/Fig7b.pdf",Fig7b,height = 6,width = 8)
head(tcga.mcp)
mcp_cor_RS=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,
                            tcga.mcp[tcga.risktype.cli$Samples,])
mcp_cor_RS=Hmisc::rcorr(as.matrix(mcp_cor_RS),type = 'spearman')
mcp_cor_RS$P[is.na(mcp_cor_RS$P)] <- 0

pdf('07_immunity/Fig7b_2.pdf',height = 6,width = 6)
corrplot(
  mcp_cor_RS$r,
  method = c('ellipse'),
  type = c('upper'), col = rev(red_blue()),
  outline = 'grey', order = c('original'), diag = TRUE,
  tl.cex = 1, tl.col = 'black', tl.pos = 'lt',
  p.mat = mcp_cor_RS$P,
  sig.level = c(.001, .01, .05),
  insig = "label_sig", ："pch", "p-value", "blank", "n", "label_sig"
  pch.cex = 1.2,
  pch.col = 'black' 
)

corrplot(
  mcp_cor_RS$r,
  add = TRUE,
  method = c('number'), type = c('lower'), col = rev(red_blue()),
  order = c('original'), diag = FALSE, number.cex = 0.8,
  tl.pos = 'n', cl.pos = 'n',
  p.mat = mcp_cor_RS$P,
  insig = "pch"
)
dev.off()
##estimate####
tcga_estimate <- immu_estimate(exp = tcga.exp)
head(tcga_estimate)
p7c <- my_mutiboxplot(dat = tcga_estimate[tcga.risktype.cli$Samples,],group = tcga.risktype.cli$Risktype,group_cols = risktype.col)
ggsave("07_immunity/figure7c.pdf",p7c,height = 4,width = 6)
df=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,tcga_estimate[tcga.risktype.cli$Samples,-4])
estimate_cor_res <- Hmisc::rcorr(as.matrix(df),type = 'spearman')
estimate_cor_res$P[is.na(estimate_cor_res$P)] <- 0
estimate_cor_res2 <- cor(df)

pdf('07_immunity/figure7c_cor_2.pdf',height = 5,width = 5,onefile = F)
corrplot(estimate_cor_res2, p.mat = estimate_cor_res$P, 
         col=colorRampPalette(c('#F5004F', 'white','#FFAF00'))(100),
         diag = FALSE, type = 'upper',   tl.col = 'black',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', 
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[1],
         pch.col = 'grey20',order = c("original", "AOE", "FPC", "hclust", "alphabet")[1])
dev.off()

##ssgsea####
tcga.ssgsea=immu_ssgsea(exp = tcga.exp)
saveRDS(taga.ssgsea,file ='07_immunity/tcga.immu.ssgsea.RDS')
tcga.immu.ssgsea <- readRDS("07_immunity/tcga.immu.ssgsea.RDS")
pdf('07_immunity/Fig7d.pdf',height = 6,width = 14)
mg_PlotMutiBoxplot(data = tcga.immu.ssgsea[tcga.risktype.cli$Samples,],
                   group = tcga.risktype.cli$Risktype,
                   test_method = 'wilcox.test',
                   add = 'boxplot',ylab='Score',legend.pos = 'top',group_cols = risktype.col)
dev.off()

save.image('20241013_OV.Rdata')
