###########################################################################
## Make graphs showing correlation and pleiotropy findings
###########################################################################
##  data is stored:
load('/Users/wenzhang/Desktop/RData/MetaData.RData')
save(metaXcanData,file='/Users/wenzhang/Desktop/RData/MetaData.RData')
## following is done for previous MetaXcan results without moving MHC SNPs.
## for the new MHC results, no need to do removing duplicates stuff.
#***************************************
## add SLE results
#library(stringr)
#setwd('/Users/wenzhang/Desktop/SLE/')
#filenames = list.files(full.names=TRUE, pattern="SLE")
#for (i in 1:length(filenames)){
#  a<-read.csv(filenames[i],header=TRUE,stringsAsFactors=FALSE)
#  b=strsplit(gsub("_with.csv","",filenames[i]),"./")[[1]][2]
#  b1=strsplit(b,"_")[[1]][1]
#  b2=str_replace(string = b, pattern = b1, replacement = "")
#  b3=sub(pattern ='_' ,replacement = "",x = b2)
#  #a$tissue<-strsplit(gsub("./","",filenames[i]),"_")[[1]][1]
#  #a$trait<-strsplit(gsub("./","",filenames[i]),"_")[[1]][2]
#  a$tissue<-b1
#  a$trait<-b3
#  write.csv(a,file=filenames[i],quote=F,row.names=F)
#}

#rawfileSLE= do.call("rbind", lapply(filenames, read.csv, header = TRUE, stringsAsFactors=FALSE))
#rawfileSLE[rawfileSLE$tissue=='gtexAor',13]='GTEx_AOR'
#rawfileSLE[rawfileSLE$tissue=='gtexBLD',13]='GTEx_BLD'
# rawfileSLE[rawfileSLE$tissue=='gtexLiv',13]='GTEx_Liver'
# rawfileSLE[rawfileSLE$tissue=='gtexSF',13]='GTEx_SF'
# rawfileSLE[rawfileSLE$tissue=='gtexVAF',13]='GTEx_VAF'
# rawfileSLE[rawfileSLE$tissue=='gtexSKLM',13]='GTEx_SKLM'
# rawfileSLE[rawfileSLE$tissue=='StarAor',13]='STARNET_AOR'
# rawfileSLE[rawfileSLE$tissue=='StarBLD',13]='STARNET_BLD'
# rawfileSLE[rawfileSLE$tissue=='StarLiv',13]='STARNET_Liver'
# rawfileSLE[rawfileSLE$tissue=='StarSF',13]='STARNET_SF'
# rawfileSLE[rawfileSLE$tissue=='StarVAF',13]='STARNET_VAF'
# rawfileSLE[rawfileSLE$tissue=='StarSKLM',13]='STARNET_SKLM'
# rawfileSLE[rawfileSLE$tissue=='StarMam',13]='STARNET_MAM'
# dat1=dat1[,-15]
# dat1=dat1[,-15]
# dat1=rbind(dat1,rawfileSLE)

# # remove several traits:
# dat=dat1[dat1$trait!='cIMT',]
# dat=dat[dat$trait!='cIMT_PLQ',]
# dat=dat[dat$trait!='iPSYCH_ADHD',]
# dat=dat[dat$trait!='PGC2_BD',]
# dat=dat[dat$trait!='PGC2_MDD',]
# dat=dat[dat$trait!='CAD_REC',]
# dat=dat[dat$trait!='NEUROTICISM',]
# # remove the duplicated Eduyear(s)
# dat=dat[dat$trait!='EduYear',]
# dat=dat[dat$trait!='EduYears',]

# # get unique EduYear
# Edu1=dat1[dat1$trait=='EduYears',]
# Edu=dat1[dat1$trait=='EduYear',]
# # STARNET_Liver is duplicated in EduYears and EduYear: remove it
# Edu1=Edu1[Edu1$tissue!='STARNET_Liver',]
# Edu1[,14]='EduYear'
# Edu=rbind(Edu,Edu1)
# # combine the unique Edu
# dat=rbind(dat,Edu)
# # has replaced ASD with new one: iPSYCH_ASD
# dat=dat[dat$trait!='ASD',]
# # remove one duplicate SCZ: the same to PGC2_SCZ
# dat=dat[dat$trait!='SCZ_PGC2',]

#**********************************************************

# load trait names and category:
Traitname=read.csv('/Users/wenzhang/Desktop/Trait_AGene/TraitNameCategory.csv',header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Traitname)<-Traitname$V1
# order the trait according to category
a=order(Traitname[unique(dat$trait),3])
#display the 58 traits
Traitname[unique(dat$trait)[a],2]

# data that reduce SNPs in MHC region
load('/Users/wenzhang/Desktop/Trait_AGene/MetaXcan-MHC.RData')
datMHC=metaXcanData

suppressMessages(library(qvalue))
library(RColorBrewer)
myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")
my_palette=colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),space="Lab")
library(WGCNA)
mpdf=function(x,width=10,height=7)eval.parent(substitute({ pdf(paste0(outDir,"/plot_",gsub("(\\(|\\))","",gsub(" ","_",x)),".pdf"),width=width,height=height) }))
outDir='/Users/wenzhang/Desktop/'

gwasTissueSorter=function(x)
  x[  tissue[tissue %in% rownames(x)],  trait[trait %in% colnames(x)]  ,drop=F]


## all MetaXcan results
#dat1<-read.csv('/Users/wenzhang/Desktop/MetaXresult.csv',header=TRUE,stringsAsFactors = FALSE)
## filter out null pavlues
#match=is.na(dat1$pvalue)
#dat1=dat1[!match,]
## get adjusted p values -- qvalues
## using FDR<=0.01 to filter significant genes
#qobjwo <- qvalue(dat1$pvalue, fdr.level = 0.01)
#dat1$pva.qval=qobjwo$qvalues
#dat1$isSignificant=qobjwo$significant
# count number of significant genes
## without removing SNPs in MHC region:
# match=is.na(dat$pvalue)
# dat=dat[!match,]
# qobjwo <- qvalue(dat$pvalue, fdr.level = 0.01)
# dat$pva.qval=qobjwo$qvalues
# dat$isSignificant=qobjwo$significant
# dat=dat[dat$pred_perf_qval<=0.01,]

# for removing MHC results
match=is.na(datMHC$pvalue)
datMHC=datMHC[!match,]
qobjwo <- qvalue(datMHC$pvalue, fdr.level = 0.01)
datMHC$pva.qval=qobjwo$qvalues
datMHC$isSignificant=qobjwo$significant
datMHC=datMHC[datMHC$pred_perf_qval<=0.01,]

#datTab=table(datSig[datSig$isSignificant,c("tissue","trait")])

datTabMHC=table(datMHC[datMHC$isSignificant,c("tissue","trait")])


# to get the order that same to the category order
# indcategory=rep(0,58)
# for (i in 1:length(Traitname[,3])){
#   for (j in 1:length(colnames(datTab))){
#     if (colnames(datTab)[j]==Traitname[i,1])
#       indcategory[i]=j
#   }
# }

indcategory=rep(0,58)
for (i in 1:length(Traitname[,3])){
  for (j in 1:length(colnames(datTabMHC))){
    if (colnames(datTabMHC)[j]==Traitname[i,1])
      indcategory[i]=j
  }
}

Tissuename=read.csv('/Users/wenzhang/Desktop/Trait_AGene/TissueName.csv',header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Tissuename)<-Tissuename$V1


# # all of the numbers
# allTab=table(datSig[,c("tissue","trait")])
# # normalized count numbers
# NORMALIZED_SIGNIFICANT_COUNTS = (datTab/(1+allTab[rownames(datTab),colnames(datTab)]))

# all of the numbers
datSigMHC=datMHC[datMHC$isSignificant,]
#datSigMHC1=datSigMHC[abs(datSigMHC$zscore)>min(abs(datSigMHC$zscore)),]

#allTabMHC1=table(datSigMHC1[,c("tissue","trait")])
sigTabMHC=table(datSigMHC[,c("tissue","trait")])
allTabMHC=table(datMHC[,c("tissue","trait")])
# normalized count numbers
NORMALIZED_SIGNIFICANT_COUNTS_MHC = (sigTabMHC/(allTabMHC[rownames(allTabMHC),colnames(allTabMHC)]))

# # scale the normalized count
# matrix_norm=NORMALIZED_SIGNIFICANT_COUNTS
# for(i in 1:14){
#   for (j in 1:58){
#     matrix_norm[i,j]=(NORMALIZED_SIGNIFICANT_COUNTS[i,j]-mean(NORMALIZED_SIGNIFICANT_COUNTS[,j]))/sd(NORMALIZED_SIGNIFICANT_COUNTS[,j])
#   }
# }
# matrix_norm=matrix_norm[,indcategory]

# scale the normalized count
matrix_normMHC=NORMALIZED_SIGNIFICANT_COUNTS_MHC
for(i in 1:14){
  for (j in 1:58){
    matrix_normMHC[i,j]=(NORMALIZED_SIGNIFICANT_COUNTS_MHC[i,j]-mean(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j]))/sd(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j])
  }
}
# matrix_normMHC=matrix_normMHC[,indcategory]

# matrix_norm1=matrix(0,14,58)
# for(i in 1:14){
#   for (j in 1:58){
#     matrix_norm1[i,j]=datTab[i,j]
#   }
# }

# datSig=dat[dat$isSignificant,]

# datSig1=datSig[abs(datSig$zscore)>=mean(abs(datSig$zscore)),]

#save(dat1,dat,datSig,datTab,allTab,file='/Users/wenzhang/Desktop/RData/MetaData.RData')


TRAIT=unique(datSigMHC$trait)
TISSUE=unique(datSigMHC$tissue)
GENE=unique(datSigMHC$gene)

# traitGenecount=matrix(0,length(TRAIT),1)
# tissueGenecount=matrix(0,length(TISSUE),1)
# geneTraitcount=matrix(0,length(GENE),1)

traitGenecountMHC=matrix(0,length(TRAIT),1)
tissueGenecountMHC=matrix(0,length(TISSUE),1)
geneTraitcountMHC=matrix(0,length(GENE),1)


# rownames(traitGenecount)<-TRAIT
# rownames(tissueGenecount)<-TISSUE
# rownames(geneTraitcount)<-GENE

rownames(traitGenecountMHC)<-TRAIT
rownames(tissueGenecountMHC)<-TISSUE
rownames(geneTraitcountMHC)<-GENE

# for (i in 1:length(TRAIT)){
#   dattmp=subset(datSig,datSig$trait==TRAIT[i])
#   traitGenecount[i]=length(unique(dattmp$gene))
# }
# 
# for (i in 1:length(TISSUE)){
#   dattmp=subset(datSig,datSig$tissue==TISSUE[i])
#   tissueGenecount[i]=length(unique(dattmp$gene))
# }
# 
for (i in 1:length(GENE)){
  dattmp=subset(datSigMHC,datSigMHC$gene==GENE[i])
  geneTraitcountMHC[i]=length(unique(dattmp$trait))
}

for (i in 1:length(TRAIT)){
  dattmp=subset(datSigMHC,datSigMHC$trait==TRAIT[i])
  traitGenecountMHC[i]=length(unique(dattmp$gene))
}

for (i in 1:length(TISSUE)){
  dattmp=subset(datSigMHC,datSigMHC$tissue==TISSUE[i])
  tissueGenecountMHC[i]=length(unique(dattmp$gene))
}



geneTraitcoun=geneTraitcountMHC
colnames(geneTraitcoun)<-'Traitnumber'
geneTraitcoun=geneTraitcoun[order(geneTraitcoun,decreasing = TRUE)]


# colNam<-colnames(matrix_norm)
# colNam1=Traitname[colNam,2]
# colNam1<-paste(colNam1," (",traitGenecount[colNam,],") ",sep='')
# rowNam<-rownames(matrix_norm)
# rowNam1=Tissuename[rowNam,2]
# rowNam1<-paste(rowNam1," (",tissueGenecount[rowNam,],") ",sep='')
# 
# mpdf("SIGNIFICANT_COUNTS.NORMALIZED_RATIO_countNumber1_SigAll", width=3+length(unique(dat$trait))*0.30, height=4.5+length(unique(dat$tissue))*0.30)
# mar.default = c(6,5,5,3) + 0.1
# par(mar = mar.default + c(11, 13, 0, 0)) #c(bottom, left, top, right)
# labeledHeatmap(Matrix = apply(matrix_norm, 2, scale),
#                xLabels = colNam1,
#                yLabels = rowNam1,
#                colorLabels = F,
#                colors = myPalette2way(1000),
#                textMatrix = datTab,
#                setStdMargins = FALSE,
#                cex.text = 0.65, 
#                zlim = c(-3.5,3.5),
#                naColor = "white",
#                main = "")
# dev.off()


colNamMHC<-colnames(matrix_normMHC)
colNam1MHC=Traitname[colNamMHC,2]
colNam1MHC<-paste(colNam1MHC," (",traitGenecountMHC[colNamMHC,],") ",sep='')
rowNamMHC<-rownames(matrix_normMHC)
rowNam1MHC=Tissuename[rowNamMHC,2]
rowNam1MHC<-paste(rowNam1MHC," (",tissueGenecountMHC[rowNamMHC,],") ",sep='')

matrix_normMHC=matrix_normMHC[,indcategory]
datTabMHC=datTabMHC[,indcategory]
mpdf("SIGNIFICANT_COUNTS.NORMALIZED_RATIO_countNumber1-MHC", width=3+length(unique(datMHC$trait))*0.30, height=4.5+length(unique(datMHC$tissue))*0.30)
mar.default = c(6,5,5,3) + 0.1
par(mar = mar.default + c(11, 13, 0, 0)) #c(bottom, left, top, right)
labeledHeatmap(Matrix = apply(matrix_normMHC, 2, scale),
               xLabels = colNam1MHC,
               yLabels = rowNam1MHC,
               colorLabels = F,
               colors = myPalette2way(1000),
               textMatrix = datTabMHC,
               setStdMargins = FALSE,
               cex.text = 0.65, 
               zlim = c(-3.5,3.5),
               naColor = "white",
               main = "")
dev.off()
## To get gene numbers and gene names per tissue
STARNET=datSig[grep('STARNET', datSig$tissue),]
CMC=datSig[grep('CMC', datSig$tissue),]
GTEx=datSig[grep('GTEx', datSig$tissue),]

starnetgene=unique(STARNET$gene)
gtexgene=unique(GTEx$gene)
cmcgene=unique(CMC$gene)
library(VennDiagram)
library(limma)

# to make venn diagram of three gene sets
universe <- sort(unique(c(starnetgene, gtexgene, cmcgene)))
Counts <- matrix(0, nrow=length(universe), ncol=3)
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% starnetgene
  Counts[i,2] <- universe[i] %in% gtexgene
  Counts[i,3] <- universe[i] %in% cmcgene
}

# Name the columns with the tissue names
colnames(Counts) <- c("starnetgene","gtexgene","cmcgene")

# Specify the colors for the sets and make Figure 2A
cols<-c("Red", "Green", "Blue")
pdf('/Users/wenzhang/Desktop/venn')
vennDiagram(vennCounts(Counts), circle.col=cols)
dev.off()

# for removing MHC

STARNET=datSigMHC[grep('STARNET', datSigMHC$tissue),]
CMC=datSigMHC[grep('CMC', datSigMHC$tissue),]
GTEx=datSigMHC[grep('GTEx', datSigMHC$tissue),]

starnetgene=unique(STARNET$gene)
gtexgene=unique(GTEx$gene)
cmcgene=unique(CMC$gene)
library(VennDiagram)
library(limma)

# to make venn diagram of three gene sets
universe <- sort(unique(c(starnetgene, gtexgene, cmcgene)))
Counts <- matrix(0, nrow=length(universe), ncol=3)
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% starnetgene
  Counts[i,2] <- universe[i] %in% gtexgene
  Counts[i,3] <- universe[i] %in% cmcgene
}

# Name the columns with the tissue names
colnames(Counts) <- c("starnetgene","gtexgene","cmcgene")

# Specify the colors for the sets and make Figure 2A
cols<-c("Red", "Green", "Blue")
pdf('/Users/wenzhang/Desktop/venn')
vennDiagram(vennCounts(Counts), circle.col=cols)
dev.off()


mpdf("SIGNIFICANT_COUNTS.NORMALIZED_RATIO", width=2+length(unique(dat1$trait))*0.30, height=4+length(unique(dat1$tissue))*0.30)
mar.default = c(6,5,5,3) + 0.1
par(mar = mar.default + c(10, 12, 0, 0)) #c(bottom, left, top, right)
labeledHeatmap(Matrix = apply(matrix_norm, 2, scale),
               xLabels = (colnames(matrix_norm)),
               yLabels = (rownames(matrix_norm)),
               colorLabels = F,
               colors = myPalette2way(1000),
               textMatrix = datTab,
               setStdMargins = FALSE,
               cex.text = 0.65, 
               zlim = c(-3.5,3.5),
               naColor = "white",
               main = "")
dev.off()

# to count number of trait-associated genes per trait/tissue

mpdf("Density_TraitPerGene", width=3+length(unique(dat1$trait))*0.30, height=4.5+length(unique(dat1$tissue))*0.30)
hist(geneTraitcount)
dev.off()

mpdf("Density_GenePerTrait", width=3+length(unique(dat1$trait))*0.30, height=4.5+length(unique(dat1$tissue))*0.30)
hist(traitGenecount)
dev.off()

# to get traits that will be plotted with sharing associations
plotTrait=subset(traitGenecount,traitGenecount[,1]>=50)
write.table(plotTrait,file = '/Users/wenzhang/Desktop/Trait_AGene/plotTraitNode.txt',quote=F,row.names = TRUE,col.names = FALSE)
plotTrait=read.table('/Users/wenzhang/Desktop/Trait_AGene/plotTraitNode.txt',header = FALSE,stringsAsFactors = FALSE)

# for MHC
plotTrait=subset(traitGenecountMHC,traitGenecountMHC[,1]>=50)
write.table(plotTrait,file = '/Users/wenzhang/Desktop/Trait_AGene/plotTraitNodeMHC.txt',quote=F,row.names = TRUE,col.names = FALSE)
plotTrait=read.table('/Users/wenzhang/Desktop/Trait_AGene/plotTraitNodeMHC.txt',header = FALSE,stringsAsFactors = FALSE)


# to get associated genes with respect to each trait
for (i in 1:length(TRAIT)){
  dattmp=subset(datSigMHC,datSigMHC$trait==TRAIT[i])
  dattmp=dattmp[,1:2]
  write.table(dattmp,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',TRAIT[i],'_AssoGeneMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
  
}

plotTrait[,3]=Traitname[plotTrait[,1],3]
plotTrait=plotTrait[order(plotTrait[,3]),]
plotTrait[,4]=Traitname[plotTrait[,1],2]

suppressMessages(library(data.table)) 
#DT <- data.table(Edu)
#y<-unique(DT, by = "gene")
## keep the EduYear_AssoGene.txt as the final
#write.table(y,file='/Users/wenzhang/Desktop/Trait_AGene/EduYear_AssoGene.txt',quote=F,row.names=F,sep='\t')
## remove duplicated
#plotTrait=plotTrait[-33,]

sharednumber=matrix(0,dim(plotTrait)[1]*(dim(plotTrait)[1]-1)/2,1)
sharednumber=data.table(sharednumber)
write.table(sharednumber,file='/Users/wenzhang/Desktop/Trait_AGene/Trait_trait_sharedMHC.txt',quote=F,row.names=F,col.names = F,sep='\t')
sharednumber=read.table('/Users/wenzhang/Desktop/Trait_AGene/Trait_trait_sharedMHC.txt',header=FALSE,stringsAsFactors = FALSE)

sharedindex=1
for (i in 1:(dim(plotTrait)[1]-1)){
  for (j in (i+1):dim(plotTrait)[1]){
    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',plotTrait[i,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',plotTrait[j,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    a=intersect(gene1$gene,gene2$gene)
    sharednumber[sharedindex,1]=length(a)
    rownames(sharednumber)[sharedindex]=paste(plotTrait[i,1],plotTrait[j,1],sep = '+')
    sharedindex=sharedindex+1
  }
}

write.table(sharednumber,file='/Users/wenzhang/Desktop/Trait_AGene/Trait_trait_shared_edgeMHC.txt',row.names = TRUE,col.names = FALSE,quote = FALSE)

Category<-read.table('/Users/wenzhang/Desktop/Trait_AGene/TraitCategory.txt',header=TRUE,stringsAsFactors = FALSE,row.names = 'Trait')

sharedAsso=read.table('/Users/wenzhang/Desktop/Trait_AGene/Trait_trait_shared_edgeMHC_cleaned.txt',header=FALSE,stringsAsFactors = FALSE)



sharedAsso$edge=2
for (i in 1:dim(sharedAsso)[1]){
  if (Traitname[sharedAsso[i,1],3]==Traitname[sharedAsso[i,2],3]) { 
    sharedAsso[i,4]=1}
}

rownames(plotTrait)<-plotTrait$V1 

sharedAsso$node=1
sharedAsso$node2=1

for (i in 1:dim(sharedAsso)[1]){
  # set node size
  sharedAsso[i,5]=plotTrait[sharedAsso[i,1],2]
  sharedAsso[i,6]=plotTrait[sharedAsso[i,2],2]
  # change node names to full name
  sharedAsso[i,1]=plotTrait[sharedAsso[i,1],4]
  sharedAsso[i,2]=plotTrait[sharedAsso[i,2],4]
    
}
#sharedAsso[,3]=sharedAsso[,3]*100 # edge width can be adjusted in Cytoscape

write.table(sharedAsso,'/Users/wenzhang/Desktop/Trait_AGene/CategoryPlotMHC.txt',row.names = FALSE,quote = FALSE,col.names = TRUE,sep='\t')


#### Make autoimmune disorder associated gene plots
## get subsets for each autoimmune disorder

datSigMHC1=datSigMHC[abs(datSigMHC$zscore)>=mean(abs(datSigMHC$zscore)),]

ActDerm=subset(datSigMHC,datSigMHC$trait=='ActDerm')
SLE=subset(datSigMHC,datSigMHC$trait=='SLE')
CD=subset(datSigMHC,datSigMHC$trait=='CD')
UC=subset(datSigMHC,datSigMHC$trait=='UC')
RA=subset(datSigMHC,datSigMHC$trait=='RA')

ActDerm=subset(ActDerm,ActDerm$pred_perf_qval<=0.001)
SLE=subset(SLE,SLE$pred_perf_qval<=0.001)
CD=subset(CD,CD$pred_perf_qval<=0.001)
UC=subset(UC,UC$pred_perf_qval<=0.001)
RA=subset(RA,RA$pred_perf_qval<=0.001)

ActDerm=subset(ActDerm,ActDerm$pva.qval<=0.001)
SLE=subset(SLE,SLE$pva.qval<=0.001)
CD=subset(CD,CD$pva.qval<=0.001)
UC=subset(UC,UC$pva.qval<=0.001)
RA=subset(RA,RA$pva.qval<=0.001)

# number of genes whose expression was associated with multiple traits
autoimmune=matrix(0,5,1)
autoimmune[1,1]=length(unique(ActDerm$gene))
autoimmune[2,1]=length(unique(SLE$gene))
autoimmune[3,1]=length(unique(CD$gene))
autoimmune[4,1]=length(unique(UC$gene))
autoimmune[5,1]=length(unique(RA$gene))
# make node size of autoimmune disorder
rownames(autoimmune)<-c('ActDerm','SLE','CD','UC','RA')
write.table(autoimmune,'/Users/wenzhang/Desktop/Trait_AGene/autoimmuneNodeSizeMHC.txt',row.names = TRUE,col.names = FALSE,quote = FALSE)
autoimmune=read.table('/Users/wenzhang/Desktop/Trait_AGene/autoimmuneNodeSizeMHC.txt',header=FALSE,stringsAsFactors = FALSE)
rownames(autoimmune)<-autoimmune$V1

autoshared=matrix(0,dim(autoimmune)[1]*(dim(autoimmune)[1]-1)/2,1)
autoshared=data.table(autoshared)
write.table(autoshared,file='/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_sharedMHC.txt',quote=F,row.names=F,col.names = F,sep='\t')
autoshared=read.table('/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_sharedMHC.txt',header=FALSE,stringsAsFactors = FALSE)

autosharedindex=1
for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[j,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    a=intersect(gene1$gene,gene2$gene)
    autoshared[autosharedindex,1]=length(a)
    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
    autosharedindex=autosharedindex+1
  }
}

write.table(autoshared,file='/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_shared_edgenumberMHC.txt',row.names = TRUE,col.names = FALSE,quote = FALSE)


for (i in 1:dim(autoimmune)[1]){
  dattmp=subset(datSigMHC,datSigMHC$trait==autoimmune[i,1])
  # for five autoimmune disorders, filtering the genes with low pred.perf.qvalue 
  dattmp=subset(dattmp,dattmp$pred_perf_qval<=0.001)
  dattmp=subset(dattmp,dattmp$pva.qval<=0.001)
  dattmp=dattmp[,1:2]
  write.table(dattmp,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'_AssoGene_filterMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
  
}

ActDerm=subset(ActDerm,ActDerm$pva.qval<=0.001)
SLE=subset(SLE,SLE$pva.qval<=0.001)
CD=subset(CD,CD$pva.qval<=0.001)
UC=subset(UC,UC$pva.qval<=0.001)
RA=subset(RA,RA$pva.qval<=0.001)

# number of genes whose expression was associated with multiple traits
autoimmune=matrix(0,5,1)
autoimmune[1,1]=length(unique(ActDerm$gene))
autoimmune[2,1]=length(unique(SLE$gene))
autoimmune[3,1]=length(unique(CD$gene))
autoimmune[4,1]=length(unique(UC$gene))
autoimmune[5,1]=length(unique(RA$gene))
# make node size of autoimmune disorder
rownames(autoimmune)<-c('ActDerm','SLE','CD','UC','RA')
write.table(autoimmune,'/Users/wenzhang/Desktop/Trait_AGene/autoimmuneNodeSizeMHC.txt',row.names = TRUE,col.names = FALSE,quote = FALSE)
autoimmune=read.table('/Users/wenzhang/Desktop/Trait_AGene/autoimmuneNodeSizeMHC.txt',header=FALSE,stringsAsFactors = FALSE)
rownames(autoimmune)<-autoimmune$V1

autoshared=matrix(0,dim(autoimmune)[1]*(dim(autoimmune)[1]-1)/2,1)
autoshared=data.table(autoshared)
write.table(autoshared,file='/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_sharedMHC.txt',quote=F,row.names=F,col.names = F,sep='\t')
autoshared=read.table('/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_sharedMHC.txt',header=FALSE,stringsAsFactors = FALSE)

## filter genes with high effects (up or downregulated or both)

for (i in 1:(dim(autoimmune)[1])){
    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'_AssoGene_filterMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    a=gene1$gene
    aa=paste('tmp1=',autoimmune[i,1],'[',autoimmune[i,1],'$gene %in% a,]',sep='')
    eval(parse(text = aa))
    assogene=tmp1[abs(tmp1$zscore)>=mean(abs(tmp1$zscore)),]
    write.table(assogene[,1:3],file=paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = TRUE)
}

autosharedindex=1
for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'_SharedGeneHigheffectMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    a=intersect(gene1$gene,gene2$gene)
    write.table(a,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,col.names = FALSE,row.names = FALSE)
    autoshared[autosharedindex,1]=length(a)
    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
    autosharedindex=autosharedindex+1
  }
}


#a=unique(datSig1$gene)
#for (i in 1:(dim(autoimmune)[1]-1)){
#  for (j in (i+1):dim(autoimmune)[1]){
#    if (file.info(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffect.txt',sep=''))$size !=0){
#    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffect.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
#    a=intersect(gene1[,1],gene2$gene)
#    write.table(a,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffect.txt',sep=''),quote=FALSE,col.names = FALSE,row.names = FALSE)
#    autoshared[autosharedindex,1]=length(a)
#    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
#    autosharedindex=autosharedindex+1}
    
#  }
#}


autoimmuneGeneTrait=matrix(0,sum(autoshared)*2,1)
autoimmuneGeneTrait=data.table(autoimmuneGeneTrait)
write.table(autoimmuneGeneTrait,file='/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_Gene_correlationMHC.txt',quote=F,row.names=F,col.names = F,sep='\t')
autoimmuneGeneTrait=read.table('/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_Gene_correlationMHC.txt',header=FALSE,stringsAsFactors = FALSE)

autogenetraitindex=1

for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    if (file.info(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''))$size !=0)
    {gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
    if (dim(gene1)[1]>0){
      autoimmuneGeneTrait[autogenetraitindex:(dim(gene1)[1]+autogenetraitindex-1),1]=autoimmune[i,1]
      autoimmuneGeneTrait[autogenetraitindex:(dim(gene1)[1]+autogenetraitindex-1),2]=gene1$V1
      autoimmuneGeneTrait[(dim(gene1)[1]+autogenetraitindex):(2*dim(gene1)[1]+autogenetraitindex-1),1]=autoimmune[j,1]
      autoimmuneGeneTrait[(dim(gene1)[1]+autogenetraitindex):(2*dim(gene1)[1]+autogenetraitindex-1),2]=gene1$V1
      autogenetraitindex=autogenetraitindex+2*dim(gene1)[1]
    }
}
  }
}

# set node size
autoimmuneGeneTrait[,3]=autoimmune[autoimmuneGeneTrait$V1,2]
# set trait color (corresponds to yellow)
autoimmuneGeneTrait[,4]=4
# set gene color (corresponds to grey)
autoimmuneGeneTrait[,5]=5
# set gene size (should be the same)
autoimmuneGeneTrait[,6]=2.5


allnumber=0
for(i in 1:dim(autoimmuneGeneTrait)[1]){
  aa=(datSigMHC[datSigMHC$trait==autoimmuneGeneTrait[i,1],])
  aa=aa[aa$gene==autoimmuneGeneTrait[i,2],]
  # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
  autoimmuneGeneTrait[i,7]=dim(aa)[1]
  # raw evaluation of up/down regulations
  upordow=(aa$zscore>=0)
  allnumber=0
  for (j in 1:length(upordow))
  {
    if (isTRUE(upordow[j])) allnumber=allnumber+1
  }
  # all positively correlated
  if (allnumber==length(upordow)) autoimmuneGeneTrait[i,8]=1
  # there are tissues where the gene is negatively correlated with the trait 
  if (allnumber<length(upordow)) {
    autoimmuneGeneTrait[i,8]=3
    # positively correlated cases are less than half --> negatively correlated
    if (allnumber<length(upordow)/2) autoimmuneGeneTrait[i,8]=2
    }
  
}

for(i in 1:dim(autoimmuneGeneTrait)[1]){
  aa=datSigMHC[datSigMHC$gene==autoimmuneGeneTrait[i,2],]
  autoimmuneGeneTrait[i,2]=aa[1,2]
}

# advanced evaluations of up/down regulations
autoimmuneGeneReg=autoimmuneGeneTrait[,1:2]
autoimmuneGeneReg[,3:17]=0
colnames(autoimmuneGeneReg)<-c('trait','gene',unique(datSig1$tissue))
for(i in 1:dim(autoimmuneGeneReg)[1]){
  aa=(datSigMHC[datSigMHC$trait==autoimmuneGeneReg[i,1],])
  aa=aa[aa$gene_name==autoimmuneGeneReg[i,2],]
  for (j in 3:16){
    tmp=aa[aa$tissue==colnames(autoimmuneGeneReg)[j],]
    if (dim(tmp)[1]>0){
    if (tmp$zscore<0){
      autoimmuneGeneReg[i,j]=1
    }
      if (tmp$zscore>0){
        autoimmuneGeneReg[i,j]=2
      }
    }
  }
  
  # the number counts of no, down and up regulations:
  a0=0
  a1=0
  a2=0
  for (j in 3:16){
    if (autoimmuneGeneReg[i,j]==0) a0=a0+1
    if (autoimmuneGeneReg[i,j]==1) a1=a1+1
    if (autoimmuneGeneReg[i,j]==2) a2=a2+1
  }
  if (a0==14) autoimmuneGeneReg[i,17]=0
  if (a0!=14){
    if (a1==a2) autoimmuneGeneReg[i,17]=0
    if (a1>a2) autoimmuneGeneReg[i,17]=1
    if (a1<a2)  autoimmuneGeneReg[i,17]=2
    }
  
  
 # # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
#  autoimmuneGeneTrait[i,7]=dim(aa)[1]
#  # raw evaluation of up/down regulations
 # upordow=(aa$zscore>=0)
#  allnumber=0
#  for (j in 1:length(upordow))
#  {
 #   if (isTRUE(upordow[j])) allnumber=allnumber+1
#  }
 # # all positively correlated
#  if (allnumber==length(upordow)) autoimmuneGeneTrait[i,8]=1
#  # there are tissues where the gene is negatively correlated with the trait 
#  if (allnumber<length(upordow)) {
#    autoimmuneGeneTrait[i,8]=3
#    # positively correlated cases are less than half --> negatively correlated
#    if (allnumber<length(upordow)/2) autoimmuneGeneTrait[i,8]=2
#  }
  
}

colnames(autoimmuneGeneReg)[17]<-'Regulation(1down/2up)'

autoimmuneGeneTrait[,8]<-autoimmuneGeneReg$`Regulation(1down/2up)`

autoimmuneGeneTrait=autoimmuneGeneTrait[,-5]

colnames(autoimmuneGeneTrait)<-c('trait','gene','node1Size','Node1color','Node2Color/Size','edgeWidth','Regulation(1down/2up)')

autoimmuneGeneTrait<-autoimmuneGeneTrait[autoimmuneGeneTrait$`Regulation(1down/2up)`!=0,]

write.table(autoimmuneGeneTrait,'/Users/wenzhang/Desktop/Trait_AGene/Autoimmune_Trait_geneMHC.txt',quote = FALSE,row.names = FALSE,col.names = TRUE)

####### #### Make neuropsychiatric disorder associated gene plots
## get subsets for each autoimmune disorder

datSigMHC1=datSigMHC[abs(datSigMHC$zscore)>=mean(abs(datSigMHC$zscore)),]

SCZ=subset(datSigMHC,datSigMHC$trait=='PGC2_SCZ')
Eduyear=subset(datSigMHC,datSigMHC$trait=='EduYear')
Neuroticism=subset(datSigMHC,datSigMHC$trait=='WellBeing_Neuroticism')
AD=subset(datSigMHC,datSigMHC$trait=='AD')
ADHC=subset(datSigMHC,datSigMHC$trait=='iPSYCH_ADHD_EUR')

SCZ=subset(SCZ,SCZ$pred_perf_qval<=0.001)
Eduyear=subset(Eduyear,Eduyear$pred_perf_qval<=0.001)
Neuroticism=subset(Neuroticism,Neuroticism$pred_perf_qval<=0.001)
AD=subset(AD,AD$pred_perf_qval<=0.001)
ADHC=subset(ADHC,ADHC$pred_perf_qval<=0.001)

SCZ=subset(SCZ,SCZ$pva.qval<=0.001)
Eduyear=subset(Eduyear,Eduyear$pva.qval<=0.001)
Neuroticism=subset(Neuroticism,Neuroticism$pva.qval<=0.001)
AD=subset(AD,AD$pva.qval<=0.001)
ADHC=subset(ADHC,ADHC$pva.qval<=0.001)

# number of genes whose expression was associated with multiple traits
neuropsy=matrix(0,5,1)
neuropsy[1,1]=length(unique(SCZ$gene))
neuropsy[2,1]=length(unique(Eduyear$gene))
neuropsy[3,1]=length(unique(Neuroticism$gene))
neuropsy[4,1]=length(unique(AD$gene))
neuropsy[5,1]=length(unique(ADHC$gene))
# make node size of autoimmune disorder
rownames(neuropsy)<-c('PGC2_SCZ','EduYear','WellBeing_Neuroticism','AD','iPSYCH_ADHD_EUR')
write.table(neuropsy,'/Users/wenzhang/Desktop/Trait_AGene/neuropsyNodeSizeMHC.txt',row.names = TRUE,col.names = FALSE,quote = FALSE)
neuropsy=read.table('/Users/wenzhang/Desktop/Trait_AGene/neuropsyNodeSizeMHC.txt',header=FALSE,stringsAsFactors = FALSE)
rownames(neuropsy)<-neuropsy$V1

autoshared=matrix(0,dim(neuropsy)[1]*(dim(neuropsy)[1]-1)/2,1)
autoshared=data.table(autoshared)
write.table(autoshared,file='/Users/wenzhang/Desktop/Trait_AGene/Neuropsy_sharedMHC.txt',quote=F,row.names=F,col.names = F,sep='\t')
autoshared=read.table('/Users/wenzhang/Desktop/Trait_AGene/Neuropsy_sharedMHC.txt',header=FALSE,stringsAsFactors = FALSE)

autosharedindex=1
for (i in 1:(dim(neuropsy)[1]-1)){
  for (j in (i+1):dim(neuropsy)[1]){
    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',neuropsy[i,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',neuropsy[j,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    a=intersect(gene1$gene,gene2$gene)
    autoshared[autosharedindex,1]=length(a)
    rownames(autoshared)[autosharedindex]=paste(neuropsy[i,1],neuropsy[j,1],sep = '+')
    autosharedindex=autosharedindex+1
  }
}

write.table(autoshared,file='/Users/wenzhang/Desktop/Trait_AGene/Neuropsy_shared_edgenumberMHC.txt',row.names = TRUE,col.names = FALSE,quote = FALSE)


for (i in 1:dim(neuropsy)[1]){
  dattmp=subset(datSigMHC,datSigMHC$trait==neuropsy[i,1])
  # for five autoimmune disorders, filtering the genes with low pred.perf.qvalue 
  dattmp=subset(dattmp,dattmp$pred_perf_qval<=0.001)
  dattmp=subset(dattmp,dattmp$pva.qval<=0.001)
  dattmp=dattmp[,1:2]
  write.table(dattmp,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',neuropsy[i,1],'_AssoGene_filterMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
  
}

SCZ=subset(SCZ,SCZ$pva.qval<=0.001)
AD=subset(AD,AD$pva.qval<=0.001)
ADHC=subset(ADHC,ADHC$pva.qval<=0.001)
Neuroticism=subset(Neuroticism,Neuroticism$pva.qval<=0.001)
Eduyear=subset(Eduyear,Eduyear$pva.qval<=0.001)

# number of genes whose expression was associated with multiple traits
neuropsy=matrix(0,5,1)
neuropsy[1,1]=length(unique(SCZ$gene))
neuropsy[2,1]=length(unique(Eduyear$gene))
neuropsy[3,1]=length(unique(Neuroticism$gene))
neuropsy[4,1]=length(unique(AD$gene))
neuropsy[5,1]=length(unique(ADHC$gene))
# make node size of autoimmune disorder
rownames(neuropsy)<-c('PGC2_SCZ','EduYear','WellBeing_Neuroticism','AD','iPSYCH_ADHD_EUR')
write.table(neuropsy,'/Users/wenzhang/Desktop/Trait_AGene/NeuropsyNodeSizeMHC.txt',row.names = TRUE,col.names = FALSE,quote = FALSE)
neuropsy=read.table('/Users/wenzhang/Desktop/Trait_AGene/NeuropsyNodeSizeMHC.txt',header=FALSE,stringsAsFactors = FALSE)
rownames(neuropsy)<-neuropsy$V1

autoshared=matrix(0,dim(neuropsy)[1]*(dim(neuropsy)[1]-1)/2,1)
autoshared=data.table(autoshared)
write.table(autoshared,file='/Users/wenzhang/Desktop/Trait_AGene/neuropsy_sharedMHC.txt',quote=F,row.names=F,col.names = F,sep='\t')
autoshared=read.table('/Users/wenzhang/Desktop/Trait_AGene/neuropsy_sharedMHC.txt',header=FALSE,stringsAsFactors = FALSE)

## filter genes with high effects (up or downregulated or both)
neuropsy[,1]=c('SCZ','Eduyear','Neuroticism','AD','ADHC')

for (i in 1:(dim(neuropsy)[1])){
  gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',rownames(neuropsy)[i],'_AssoGene_filterMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
  a=gene1$gene
  aa=paste('tmp1=',neuropsy[i,1],'[',neuropsy[i,1],'$gene %in% a,]',sep='')
  eval(parse(text = aa))
  assogene=tmp1[abs(tmp1$zscore)>=mean(abs(tmp1$zscore)),]
  write.table(assogene[,1:3],file=paste('/Users/wenzhang/Desktop/Trait_AGene/',rownames(neuropsy)[i],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = TRUE)
}


autosharedindex=1
for (i in 1:(dim(neuropsy)[1]-1)){
  for (j in (i+1):dim(neuropsy)[1]){
    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',rownames(neuropsy)[i],'_SharedGeneHigheffectMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',rownames(neuropsy)[j],'_SharedGeneHigheffectMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    a=intersect(gene1$gene,gene2$gene)
    write.table(a,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',rownames(neuropsy)[i],'+',rownames(neuropsy)[j],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,col.names = FALSE,row.names = FALSE)
    autoshared[autosharedindex,1]=length(a)
    rownames(autoshared)[autosharedindex]=paste(rownames(neuropsy)[i],rownames(neuropsy)[j],sep = '+')
    autosharedindex=autosharedindex+1
  }
}



neuroGeneTrait=matrix(0,sum(autoshared)*2,1)
neuroGeneTrait=data.table(neuroGeneTrait)
write.table(neuroGeneTrait,file='/Users/wenzhang/Desktop/Trait_AGene/Neuro_Gene_correlationMHC.txt',quote=F,row.names=F,col.names = F,sep='\t')
neuroGeneTrait=read.table('/Users/wenzhang/Desktop/Trait_AGene/Neuro_Gene_correlationMHC.txt',header=FALSE,stringsAsFactors = FALSE)

neurogenetraitindex=1

for (i in 1:(dim(neuropsy)[1]-1)){
  for (j in (i+1):dim(neuropsy)[1]){
    if (file.info(paste('/Users/wenzhang/Desktop/Trait_AGene/',rownames(neuropsy)[i],'+',rownames(neuropsy)[j],'_SharedGeneHigheffectMHC.txt',sep=''))$size !=0)
    {gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',rownames(neuropsy)[i],'+',rownames(neuropsy)[j],'_SharedGeneHigheffectMHC.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
    if (dim(gene1)[1]>0){
      neuroGeneTrait[neurogenetraitindex:(dim(gene1)[1]+neurogenetraitindex-1),1]=neuropsy[i,1]
      neuroGeneTrait[neurogenetraitindex:(dim(gene1)[1]+neurogenetraitindex-1),2]=gene1$V1
      neuroGeneTrait[(dim(gene1)[1]+neurogenetraitindex):(2*dim(gene1)[1]+neurogenetraitindex-1),1]=neuropsy[j,1]
      neuroGeneTrait[(dim(gene1)[1]+neurogenetraitindex):(2*dim(gene1)[1]+neurogenetraitindex-1),2]=gene1$V1
      neurogenetraitindex=neurogenetraitindex+2*dim(gene1)[1]
    }
    }
  }
}
rownames(neuropsy)<-neuropsy$V1
# set node size
neuroGeneTrait[,3]=neuropsy[neuroGeneTrait$V1,2]
# set trait color (corresponds to yellow)
neuroGeneTrait[,4]=4
# set gene color (corresponds to grey)
neuroGeneTrait[,5]=5
# set gene size (should be the same)
neuroGeneTrait[,6]=2.5

neuroGeneTrait[neuroGeneTrait$V1=='SCZ',1]='PGC2_SCZ'
neuroGeneTrait[neuroGeneTrait$V1=='Eduyear',1]='EduYear'
neuroGeneTrait[neuroGeneTrait$V1=='ADHC',1]='iPSYCH_ADHD_EUR'
#  'WellBeing_Neuroticism','AD','iPSYCH_ADHD_EUR

allnumber=0
for(i in 1:dim(neuroGeneTrait)[1]){
  aa=(datSigMHC[datSigMHC$trait==neuroGeneTrait[i,1],])
  aa=aa[aa$gene==neuroGeneTrait[i,2],]
  # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
  neuroGeneTrait[i,7]=dim(aa)[1]
  # raw evaluation of up/down regulations
  upordow=(aa$zscore>=0)
  allnumber=0
  for (j in 1:length(upordow))
  {
    if (isTRUE(upordow[j])) allnumber=allnumber+1
  }
  # all positively correlated
  if (allnumber==length(upordow)) neuroGeneTrait[i,8]=1
  # there are tissues where the gene is negatively correlated with the trait 
  if (allnumber<length(upordow)) {
    neuroGeneTrait[i,8]=3
    # positively correlated cases are less than half --> negatively correlated
    if (allnumber<length(upordow)/2) neuroGeneTrait[i,8]=2
  }
  
}

for(i in 1:dim(neuroGeneTrait)[1]){
  aa=datSigMHC[datSigMHC$gene==neuroGeneTrait[i,2],]
  neuroGeneTrait[i,2]=aa[1,2]
}

# advanced evaluations of up/down regulations
neuroGeneReg=neuroGeneTrait[,1:2]
neuroGeneReg[,3:17]=0
colnames(neuroGeneReg)<-c('trait','gene',unique(datSigMHC1$tissue))
for(i in 1:dim(neuroGeneReg)[1]){
  aa=(datSigMHC[datSigMHC$trait==neuroGeneReg[i,1],])
  aa=aa[aa$gene_name==neuroGeneReg[i,2],]
  for (j in 3:16){
    tmp=aa[aa$tissue==colnames(neuroGeneReg)[j],]
    if (dim(tmp)[1]>0){
      if (tmp$zscore<0){
        neuroGeneReg[i,j]=1
      }
      if (tmp$zscore>0){
        neuroGeneReg[i,j]=2
      }
    }
  }
  
  # the number counts of no, down and up regulations:
  a0=0
  a1=0
  a2=0
  for (j in 3:16){
    if (neuroGeneReg[i,j]==0) a0=a0+1
    if (neuroGeneReg[i,j]==1) a1=a1+1
    if (neuroGeneReg[i,j]==2) a2=a2+1
  }
  if (a0==14) neuroGeneReg[i,17]=0
  if (a0!=14){
    if (a1==a2) neuroGeneReg[i,17]=0
    if (a1>a2) neuroGeneReg[i,17]=1
    if (a1<a2)  neuroGeneReg[i,17]=2
  }
  
  
  
}

colnames(neuroGeneReg)[17]<-'Regulation(1down/2up)'

neuroGeneTrait[,8]<-neuroGeneReg$`Regulation(1down/2up)`

neuroGeneTrait=neuroGeneTrait[,-5]

colnames(neuroGeneTrait)<-c('trait','gene','node1Size','Node1color','Node2Color/Size','edgeWidth','Regulation(1down/2up)')

neuroGeneTrait<-neuroGeneTrait[neuroGeneTrait$`Regulation(1down/2up)`!=0,]

write.table(neuroGeneTrait,'/Users/wenzhang/Desktop/Trait_AGene/Neuro_Trait_geneMHC.txt',quote = FALSE,row.names = FALSE,col.names = TRUE)


#### for each trait, get number of genes that are up/down regulated

traitupdownGene=matrix(0,length(TRAIT),3)
rownames(traitupdownGene)=TRAIT
# second column: number of up regulated genes
# third column: number of down regulated genes
# fourth column: number of ambiguous regulated genes


for ( i in 1:dim(traitupdownGene)[1]){
  dattmp=datSigMHC[datSigMHC$trait==rownames(traitupdownGene)[i],]
  # genes that associated with the trait
  onlyassogene=unique(dattmp$gene)
  for(j in 1:length(onlyassogene)){
    # for each gene: judge if it is up/down/ambiguous expressional changed during the trait
    # get all tissues that the gene is associated with the trait
    setgene=dattmp[dattmp$gene==onlyassogene[j],]
    
    unordow=matrix(0,1,14)
    colnames(unordow)<-unique(datSigMHC$tissue)
    upchange=0
    downchange=0
    for (k in 1:14) {
      tmp=setgene[setgene$tissue==colnames(unordow)[k],]
      if (dim(tmp)[1]>0){
        if (tmp$zscore<0){ # down: set value 1
          unordow[k]=1
        }
        if (tmp$zscore>0){ # up: set value 2
          unordow[k]=2
        }
      }
    
    if (unordow[k]!=0) {
      if (unordow[k]==2) upchange=upchange+1
      if (unordow[k]==1) downchange=downchange+1    }
    
    }
    if (upchange>downchange) traitupdownGene[i,1]=traitupdownGene[i,1]+1
    if(upchange<downchange)  traitupdownGene[i,2]=traitupdownGene[i,2]+1
    if(upchange==downchange)  traitupdownGene[i,3]=traitupdownGene[i,3]+1
    
  }
}

indcategory=rep(0,58)
for (i in 1:length(Traitname[,3])){
  for (j in 1:length(rownames(traitupdownGene))){
    if (rownames(traitupdownGene)[j]==Traitname[i,1])
      indcategory[i]=j
  }
}

traitupdownGene=traitupdownGene[indcategory,]
write.table(traitupdownGene,'/Users/wenzhang/Desktop/Trait_AGene/TraitGeneNumber.txt',quote=FALSE,col.names = FALSE)
traitupdownGene=read.table('/Users/wenzhang/Desktop/Trait_AGene/TraitGeneNumber.txt',header=FALSE,stringsAsFactors = FALSE)
rownames(traitupdownGene)<-traitupdownGene$V1
traitupdownGene[,5]=Traitname[rownames(traitupdownGene),2]
traitupdownGene=traitupdownGene[,-1]

# Trait<-traitupdownGene[,4]
# Trait<-c(rep(Trait,each=3))
# Regulation<-c(rep(c("Up", "Down", "Ambiguous"), times = 58))
# Frequency<-rep(0,58*3)
# for(i in 1:58){
#   for (j in 1:3 ){
#     Frequency[(i-1)*3+j]=traitupdownGene[i,j]
#   }
# }
# Data <- data.frame(Trait, Regulation, Frequency)
# 
# mpdf('barchar')
#  ggplot(Data, aes(x = Trait, y = Frequency, fill = Regulation, label = Frequency)) +
#  geom_bar(stat = "identity") +
#   geom_text(size = 3, position = position_stack(vjust = 90))
#  dev.off()

 library(plotly)
 
 data <- data.frame(traitupdownGene[,4], traitupdownGene[,1], traitupdownGene[,2],traitupdownGene[,3])
 colnames(data)<-c('Trait','Up','Down','Ambiguous')
 #The default order will be alphabetized unless specified as below:
 data$Trait <- factor(data$Trait, levels = data[["Trait"]])
 mpdf("BarPlot", width=2+58*0.30, height=4+length(60)*0.30)
 mar.default = c(6,5,5,3) + 0.1
 par(mar = mar.default + c(10, 12, 0, 0)) #c(bottom, left, top, right)
plot_ly(data, x = ~Trait, y = ~Down, type = 'bar', name = 'Downregulation', marker = list(color = 'rgb(49,130,189)')) %>%
  add_trace(y = ~Ambiguous, name = 'Ambiguous', marker = list(color = 'rgb(4,154,4)')) %>%
  add_trace(y = ~Up,  name = 'Upregulation', marker = list(color = 'rgb(204,0,14)')) %>%
     layout(xaxis = list(title = "", tickangle = -65),
          yaxis = list(title = "Number of correlated genes"),font=t,
          margin = list(b = 225),
          barmode = 'group')
 dev.off()
 
 

#### plot shared gene correlations between 'Year of educations' and other traits 
for (i in 1:length(TRAIT)){
  dattmp=datSigMHC1[datSigMHC1$trait==TRAIT[i],]
  dattmp=dattmp[dattmp$pred_perf_qval<=0.01,]
  dattmp=dattmp[dattmp$pva.qval<=0.01,]
  dattmp=dattmp[1:2]
  write.table(dattmp,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',TRAIT[i],'_HighAssoGeneExpMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

#Edu1=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduYear_HighAssoGeneExp.txt',header=TRUE,stringsAsFactors = FALSE)
#Edu2=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduYears_HighAssoGeneExp.txt',header=TRUE,stringsAsFactors = FALSE)
#Edu=rbind(Edu1,Edu2)
Edu=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduYear_HighAssoGeneExpMHC.txt',header=TRUE,stringsAsFactors = FALSE)


suppressMessages(library(data.table)) 
DT <- data.table(Edu)
y<-unique(DT, by = "gene")
# keep the EduYear_AssoGene.txt as the final
write.table(y,file='/Users/wenzhang/Desktop/Trait_AGene/EduYear_HighAssoGeneExpMHC.txt',quote=F,row.names=F,sep='\t')

# all other traits
OtherTrait=TRAIT
OtherTrait=OtherTrait[-58]


EduYear=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduYear_HighAssoGeneExpMHC.txt',header=TRUE,stringsAsFactors = FALSE)
EduYearshared=OtherTrait[1:14]
shareindex=1
for (i in 1:length(OtherTrait)){
  other=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',OtherTrait[i],'_HighAssoGeneExpMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
  a<-intersect(EduYear$gene,other$gene)
  if (length(a)!=0){
    write.table(a,file=paste('/Users/wenzhang/Desktop/Trait_AGene/','EduYear_',OtherTrait[i],'_HighAssoGenesMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
    EduYearshared[shareindex]=OtherTrait[i]
    shareindex=shareindex+1
    }
}


# count numbers of genes other traits shared with EduYear
shareEduTraitcount=matrix(0,length(EduYearshared),1)
rownames(shareEduTraitcount)<-EduYearshared

for(i in 1:length(EduYearshared)){
  aa<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','EduYear_',EduYearshared[i],'_HighAssoGenesMHC.txt',sep=''),header=FALSE)
  shareEduTraitcount[i]=dim(aa)[1]
}

shareEduPlot=shareEduTraitcount
write.table(shareEduPlot,file='/Users/wenzhang/Desktop/Trait_AGene/EduyearplotTraitMHC.txt',quote=FALSE,row.names = FALSE,col.names = FALSE)
shareEduPlot=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduyearplotTraitMHC.txt',header=FALSE,stringsAsFactors = FALSE)
index1=1
for(i in 1:length(shareEduTraitcount)){
  #if(shareEduTraitcount[i]>=20) {
    shareEduPlot[index1,1]=shareEduTraitcount[i]
    rownames(shareEduPlot)[index1]<-rownames(shareEduTraitcount)[i]
    index1=index1+1
   # }
}

shareEduPlot[,2]=rownames(shareEduPlot)
shareEduPlot[,3]=0
shareEduPlot[,4]=0

EduYear=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduYear_HighAssoGeneExpMHC.txt',header=TRUE,stringsAsFactors = FALSE)
dattmp1=datSigMHC[datSigMHC$gene %in% EduYear$gene,]
dattmp1=dattmp1[dattmp1$pred_perf_qval<=0.01,]
dattmp1=dattmp1[dattmp1$pva.qval<=0.01,]
dattmp1=dattmp1[grep('EduYear', dattmp1$trait),]
dattmp2=dattmp1[dattmp1$trait=='EduYear',]


for(i in 1:dim(shareEduPlot)[1]){
  aa<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','EduYear_',shareEduPlot[i,2],'_HighAssoGenesMHC.txt',sep=''),header=FALSE)
  dattmp=datSigMHC[datSigMHC$gene %in% aa$V1,]
  dattmp=dattmp[dattmp$pred_perf_qval<=0.01,]
  dattmp=dattmp[dattmp$pva.qval<=0.01,]
  dattmp=dattmp[dattmp$trait==shareEduPlot[i,2],]
  intergene=intersect(dattmp$gene,dattmp1$gene)
  dattmpEdu=datSigMHC[datSigMHC$gene %in% aa$V1,]
  dattmpEdu=dattmpEdu[dattmpEdu$pred_perf_qval<=0.01,]
  dattmpEdu=dattmpEdu[dattmpEdu$pva.qval<=0.01,]
  dattmpEdu=dattmpEdu[dattmpEdu$trait=='EduYear',]
  ## loose evaluation of anticorrelation or correlation, get some values:
  ## number of 
  for(j in 1:length(intergene)){
    othergeneasso=dattmp[dattmp$gene==intergene[j],]
    Edugeneasso=dattmpEdu[dattmpEdu$gene==intergene[j],]
    unordowEdu=matrix(0,1,14)
    upordowother=matrix(0,1,14)
    colnames(unordowEdu)<-unique(datSigMHC1$tissue)
    colnames(upordowother)<-unique(datSigMHC1$tissue)
    anticorre=0
    corre=0
   for (k in 1:14) {
    tmp=othergeneasso[othergeneasso$tissue==colnames(upordowother)[k],]
    if (dim(tmp)[1]>0){
      if (tmp$zscore<0){
        upordowother[k]=1
      }
      if (tmp$zscore>0){
        upordowother[k]=2
      }
    }
    
    tmp=Edugeneasso[Edugeneasso$tissue==colnames(unordowEdu)[k],]
    if (dim(tmp)[1]>0){
      if (tmp$zscore<0){
        unordowEdu[k]=1
      }
      if (tmp$zscore>0){
        unordowEdu[k]=2
      }
    }
   
    # to determine if the gene is correlated or anti-correlated between Edu and other Trait
    # if in one tissue, one gene is only regulated during one trait, it can not be counted as anti-/correlated
    
    if (unordowEdu[k]!=0 & upordowother[k]!=0) {
      if (unordowEdu[k]!=upordowother[k]) anticorre=anticorre+1
      if (unordowEdu[k]==upordowother[k]) corre=corre+1    }
    
   }
   # third column: number of anticorrelated expressional genes
    # fourth column: number of correlated expressional genes
    if (anticorre>corre){ shareEduPlot[i,3]=shareEduPlot[i,3]+1  }
    if (anticorre<corre){ shareEduPlot[i,4]=shareEduPlot[i,4]+1  }
    
    
    }
    
    
    
  ## loose evaluation of anticorrelation or correlation, get some values:
  ## number of 
  #for(j in 1:length(intergene)){
  #  othergeneasso=dattmp[dattmp$gene==intergene[j],]
  #  upordowother=(othergeneasso$zscore>=0)
  #  allnumber=0
  #  for (j in 1:length(upordowother))
  #  {
  #    if (isTRUE(upordowother[j])) allnumber=allnumber+1
  #  }
  #  # all positively correlated
  #  if (allnumber==length(upordowother)) othercorr=1
  #  # there are tissues where the gene is negatively correlated with the trait 
  #  if (allnumber<=length(upordowother)) {
  #    othercorr=1
  #    # positively correlated cases are less than half --> negatively correlated
  #    if (allnumber<length(upordowother)/2) othercorr=2
  #  }
  #  
  #  eduYgeneasso=dattmp1[dattmp1$gene==intergene[j],]
  #  upordowedu=(eduYgeneasso$zscore>=0)
  #  allnumber=0
  #  for (j in 1:length(eduYgeneasso))
  #  {
  #    if (isTRUE(eduYgeneasso[j])) allnumber=allnumber+1
  #  }
  #  # all positively correlated
  #  if (allnumber==length(eduYgeneasso)) EduYcorr=1
  #  # there are tissues where the gene is negatively correlated with the trait 
  #  if (allnumber<=length(eduYgeneasso)) {
  #    EduYcorr=1
  #    # positively correlated cases are less than half --> negatively correlated
  #    if (allnumber<length(eduYgeneasso)/2) EduYcorr=2
  #  }
  #  if(EduYcorr==othercorr)  shareEduPlot[i,3]=shareEduPlot[i,3]+1
  #  if(EduYcorr!=othercorr)  shareEduPlot[i,4]=shareEduPlot[i,4]+1
    
  #}
  
}

shareEduPlot[,5]=Traitname[shareEduPlot[,2],2]

for(i in 1:dim(shareEduPlot)[1]){
 rownames(shareEduPlot)[i]<-paste(shareEduPlot[i,5],'(',shareEduPlot[i,4],'/',shareEduPlot[i,3],') ',sep='')
  
}


shareEduPlot[,2]='Years of education'
shareEduPlot[,5]=length(unique(EduYear$gene))
write.table(shareEduPlot,file ='/Users/wenzhang/Desktop/Trait_AGene/EduYear_PlotMHC.txt',quote=FALSE,col.names = FALSE,row.names = TRUE,sep='\t')

### To make plot of Antagonistic Pleiotropy between Schizophrenia and Other Traits

# To get upregulated and downregulated genes of Schizophrenia
dattmp=datSigMHC[datSigMHC$trait=='PGC2_SCZ',]
dattmp=dattmp[dattmp$pred_perf_qval<=0.001,]
dattmp=dattmp[dattmp$pva.qval<=0.001,]
SCZassogene=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]
SCZassogene=SCZassogene[,1:2]
write.table(SCZassogene,file = '/Users/wenzhang/Desktop/Trait_AGene/PGC2_SCZ_TopAssoGeneMHC.txt',quote=FALSE,row.names=FALSE,col.names=TRUE)

sczOthertrait=TRAIT
sczOthertrait=sczOthertrait[-35]
#sczOthertrait=sczOthertrait[-42]
#sczOthertrait=sczOthertrait[-56]
# # no need for Edu
# Edu1=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduYear_TopAssoGeneMHC.txt',header=TRUE,stringsAsFactors=FALSE)
# Edu2=read.table('/Users/wenzhang/Desktop/Trait_AGene/EduYears_TopAssoGene.txt',header=TRUE,stringsAsFactors=FALSE)
# Edu=rbind(Edu1,Edu2)
# write.table(Edu,file = '/Users/wenzhang/Desktop/Trait_AGene/EduYear_TopAssoGene.txt',quote=FALSE,row.names=FALSE,col.names=TRUE)

for (i in 1:length(TRAIT)){
  dattmp=datSigMHC[datSigMHC$trait==TRAIT[i],]
  dattmp=dattmp[dattmp$pred_perf_qval<=0.001,]
  dattmp=dattmp[dattmp$pva.qval<=0.001,]
  tmp1=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]
  tmp1=tmp1[1:2]
  write.table(tmp1,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',TRAIT[i],'_TopAssoGeneMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

for (i in 1:length(sczOthertrait)){
  other=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',sczOthertrait[i],'_TopAssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
  a<-intersect(SCZassogene$gene,other$gene)
  if (length(a)!=0){
    write.table(a,file=paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_',sczOthertrait[i],'_TopAssoGenesMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
    #EduYearshared[shareindex]=OtherTrait[i]
    #shareindex=shareindex+1
  }
}


# EduYear
dattmp=datSig[datSig$trait=='EduYear',]
dattmp=dattmp[dattmp$pred_perf_qval<=0.001,]
dattmp=dattmp[dattmp$pva.qval<=0.001,]
Eduassogene=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]

# IBD
dattmp=datSig[datSig$trait=='IBD',]
dattmp=dattmp[dattmp$pred_perf_qval<=0.001,]
dattmp=dattmp[dattmp$pva.qval<=0.001,]
IBDassogene=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]

# MI_ADD
dattmp=datSig[datSig$trait=='MI_ADD',]
dattmp=dattmp[dattmp$pred_perf_qval<=0.001,]
dattmp=dattmp[dattmp$pva.qval<=0.001,]
MI_ADDassogene=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]

# CAD_ADD
dattmp=datSig[datSig$trait=='CAD_ADD',]
dattmp=dattmp[dattmp$pred_perf_qval<=0.001,]
dattmp=dattmp[dattmp$pva.qval<=0.001,]
CAD_ADDassogene=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]

library(RCircos)
sczgeneall=NULL
#sczOthertrait=c('CD','EduYear','IBD','RA','UC','PGC2_BD')
sczOthertrait=c('CD','EduYear','IBD','PGC2_BD','UC','RA')

for(i in 1:length(sczOthertrait)){
  tmp1=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_',sczOthertrait[i],'_TopAssoGenes.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
  sczgeneall=union(sczgeneall,tmp1[,1])
  
}

CD_sczgene=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_CD_TopAssoGenes.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
EduYear_sczgene=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_EduYear_TopAssoGenes.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
IBD_sczgene=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_IBD_TopAssoGenes.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
RA_sczgene=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_RA_TopAssoGenes.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
UC_sczgene=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_UC_TopAssoGenes.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
PGC2_BD_sczgene=read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/','SCZ_PGC2_BD_TopAssoGenes.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)

dim(CD_sczgene)[1]+dim(EduYear_sczgene)[1]+dim(IBD_sczgene)[1]+dim(UC_sczgene)[1]+dim(PGC2_BD_sczgene)[1]+dim(RA_sczgene)[1]



## make table of Antagonistic pleiotropy between Schizophrenia and other traits
sczGeneTrait=matrix(0,56,1)
sczGeneTrait=data.table(sczGeneTrait)
write.table(sczGeneTrait,file='/Users/wenzhang/Desktop/Trait_AGene/SCZ_Gene_correlation.txt',quote=F,row.names=F,col.names = F,sep='\t')
SCZGeneTrait=read.table('/Users/wenzhang/Desktop/Trait_AGene/SCZ_Gene_correlation.txt',header=FALSE,stringsAsFactors = FALSE)

sczgenetraitindex=1

for (i in 1:(length(sczOthertrait)[1])){
    if (file.info(paste('/Users/wenzhang/Desktop/Trait_AGene/SCZ_',sczOthertrait[i],'_TopAssoGenes.txt',sep=''))$size!=0)
    {gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/SCZ_',sczOthertrait[i],'_TopAssoGenes.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
    if (dim(gene1)[1]>0){
      SCZGeneTrait[sczgenetraitindex:(dim(gene1)[1]+sczgenetraitindex-1),1]=gene1$V1
      SCZGeneTrait[sczgenetraitindex:(dim(gene1)[1]+sczgenetraitindex-1),2]='SCZ'
      SCZGeneTrait[sczgenetraitindex:(dim(gene1)[1]+sczgenetraitindex-1),3]=sczOthertrait[i]
      sczgenetraitindex=sczgenetraitindex+dim(gene1)[1]
    }
    }
  
}

allnumber=0
for(i in 1:dim(SCZGeneTrait)[1]){
  aa=(datSig[datSig$trait=='PGC2_SCZ',])
  aa=aa[aa$gene==SCZGeneTrait[i,1],]
  # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
  SCZGeneTrait[i,4]=dim(aa)[1]
  upordow=(aa$zscore>=0)
  allnumber=0
  for (j in 1:length(upordow))
  {
    if (isTRUE(upordow[j])) allnumber=allnumber+1
  }
  # all positively correlated
  if (allnumber==length(upordow)) SCZGeneTrait[i,5]=1 # if upregulated: V5=1
  # there are tissues where the gene is negatively correlated with the trait 
  if (allnumber<length(upordow)) {
    SCZGeneTrait[i,5]=1
    # positively correlated cases are less than half --> negatively correlated
    if (allnumber<length(upordow)/2) SCZGeneTrait[i,5]=2 # downregulated: V5=2
  }
  
  aa=(datSig[datSig$trait==SCZGeneTrait[i,3],])
  aa=aa[aa$gene==SCZGeneTrait[i,1],]
  # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
  SCZGeneTrait[i,6]=dim(aa)[1]
  upordow=(aa$zscore>=0)
  allnumber=0
  for (j in 1:length(upordow))
  {
    if (isTRUE(upordow[j])) allnumber=allnumber+1
  }
  # all positively correlated
  if (allnumber==length(upordow)) SCZGeneTrait[i,7]=1 # if upregulated: V7=1
  # there are tissues where the gene is negatively correlated with the trait 
  if (allnumber<length(upordow)) {
    SCZGeneTrait[i,7]=1
    # positively correlated cases are less than half --> negatively correlated
    if (allnumber<length(upordow)/2) SCZGeneTrait[i,7]=2 # downregulated: V7=2
  }
}

for(i in 1:dim(SCZGeneTrait)[1]){
  aa=datSig[datSig$gene==SCZGeneTrait[i,1],]
  SCZGeneTrait[i,8]=aa[1,2]
}

sczGeneTrait=SCZGeneTrait

sczGeneTrait[,2]=SCZGeneTrait[,8]
sczGeneTrait[,3]=SCZGeneTrait[,2]
sczGeneTrait[,4]=SCZGeneTrait[,5]
sczGeneTrait[,5]=SCZGeneTrait[,4]
sczGeneTrait[,6]=SCZGeneTrait[,3]
sczGeneTrait[,7]=SCZGeneTrait[,7]
sczGeneTrait[,8]=SCZGeneTrait[,6]

colnames(sczGeneTrait)=c('Gene','GeneName','Trait1','Direction1','NumberofTissue1','Trait2','Direction2','NumberofTissue2')

sczGeneTrait[,9]=0
for(i in 1:dim(sczGeneTrait)[1]){
  if(sczGeneTrait[i,4]!=sczGeneTrait[i,7])
    sczGeneTrait[i,9]=1
}

for(i in 3:9){
  if(sczGeneTrait[i,9]==1) sczGeneTrait[i,9]=0
  else sczGeneTrait[i,9]=1
}

colnames(sczGeneTrait)=c('Gene','GeneName','Trait1','Direction1','NumberofTissue1','Trait2','Direction2','NumberofTissue2','isAntagonisticPleiotropy')
sum(sczGeneTrait[,9])

SCZGeneTrait=sczGeneTrait[sczGeneTrait$isAntagonisticPleiotropy==1,]
df=data.frame(SCZGeneTrait)

df2=df[,1:3]
df2[,1]=df[,3]
df2[,2]=df[,6]
#df2[,3]=df[,3]
for (i in 1:dim(df)[1]){
  df2[i,3]=max(df[i,5],df[i,8])
}

write.table(df2,file = '/Users/wenzhang/Desktop/scz_circos.txt',quote=FALSE,row.names = FALSE,col.names = FALSE)
df2<-read.table('/Users/wenzhang/Desktop/scz_circos.txt',header=FALSE,stringsAsFactors = FALSE)
df<-data.frame(df2)
df[,3]=df[,3]/5
gap.degree=df[,3]
col = rand_color(nrow(df))
grid.col =c(SCZ = "green", CD = "red", EduYear = "black", IBD = "blue", PGC2_BD = "yellow", UC = "grey", RA = "orange")
mpdf("Circos")
#chordDiagramFromDataFrame(df)
chordDiagram(df, grid.col = grid.col, col = col,annotationTrack = c("grid","name"),preAllocateTracks = list(track.height = 0.1))

#chordDiagram(df,directional = 1)
circos.clear()

dev.off()




mpdf("Circos_test")

circos.par(points.overflow.warning = FALSE)
circos.initialize(factors = c('SCZ',sczOthertrait), xlim = c(0, 10))

circos.trackPlotRegion(ylim = c(0, 30))
circos.axis(sector.index = "SCZ",h=19,labels = df[,4],
            labels.facing = "reverse.clockwise")
circos.axis(sector.index = "CD",h=1,labels = df[1,4],
            labels.facing = "reverse.clockwise")
circos.axis(sector.index = "EduYear",h=3,labels = df[2:4,4],
            labels.facing = "reverse.clockwise")
circos.axis(sector.index = "IBD",h=1,labels = df[5,4],
            labels.facing = "reverse.clockwise")
circos.axis(sector.index = "PGC2_BD",h=3,labels = df[6:8,4],
            labels.facing = "reverse.clockwise")
circos.axis(sector.index = "UC",h=1,labels = df[9,4],
            labels.facing = "reverse.clockwise")
circos.axis(sector.index = "RA",h=10,labels = df[10:19,4],
            labels.facing = "reverse.clockwise")

chordDiagram(df, grid.col = grid.col, col = col,annotationTrack = c("grid","name"),preAllocateTracks = list(track.height = 0.1))


circos.clear()
dev.off()



library(circlize)
library(RCircos)
mat = matrix(1:18, 3, 6)
rownames(mat) = paste0("S", 1:3)
colnames(mat) = paste0("E", 1:6)

rn = rownames(mat)
cn = colnames(mat)
factors = c(rn, cn)
factors = factor(factors, levels = factors)
col_sum = apply(mat, 2, sum)
row_sum = apply(mat, 1, sum)
xlim = cbind(rep(0, length(factors)), c(row_sum, col_sum))

mpdf('circosexample')
par(mar = c(1, 1, 1, 1))
circos.par(cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = factors, xlim = xlim)
circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA,
                       bg.col = c("red", "green", "blue", rep("grey", 6)), track.height = 0.05,
                       panel.fun = function(x, y) {
                         sector.name = get.cell.meta.data("sector.index")
                         xlim = get.cell.meta.data("xlim")
                         circos.text(mean(xlim), 1.5, sector.name, adj = c(0.5, 0))
                       })

col = c("#FF000020", "#00FF0020", "#0000FF20")
for(i in seq_len(nrow(mat))) {
  for(j in seq_len(ncol(mat))) {
    circos.link(rn[i], c(sum(mat[i, seq_len(j-1)]), sum(mat[i, seq_len(j)])),
                cn[j], c(sum(mat[seq_len(i-1), j]), sum(mat[seq_len(i), j])),
                col = col[i], border = "white")
  }
}
circos.clear()
dev.off()

df = read.table(textConnection("
 brand_from model_from brand_to model_to
                               VOLVO        s80      BMW  5series
                               BMW    3series      BMW  3series
                               VOLVO        s60    VOLVO      s60
                               VOLVO        s60    VOLVO      s80
                               BMW    3series     AUDI       s4
                               AUDI         a4      BMW  3series
                               AUDI         a5     AUDI       a5
                               "), header = TRUE, stringsAsFactors = FALSE)

brand = c(structure(df$brand_from, names=df$model_from),
          structure(df$brand_to,names= df$model_to))
brand = brand[!duplicated(names(brand))]
brand = brand[order(brand, names(brand))]
brand_color = structure(2:4, names = unique(brand))
model_color = structure(2:8, names = names(brand))


mpdf('circosexample')
gap.degree = do.call("c", lapply(table(brand), function(i) c(rep(2, i-1), 8)))

circos.par(gap.degree = gap.degree)

chordDiagram(df[, c(2, 4)], order = names(brand), grid.col = model_color,
             directional = 1, annotationTrack = "grid", preAllocateTracks = list(
               list(track.height = 0.02))
)


dev.off()


library(circlize)
library(migest)
library(dplyr)

### Make data
m <- data.frame(order = 1:6,
                country = c("Ausralia", "India", "China", "Japan", "Thailand", "Malaysia"),
                V3 = c(1, 150000, 90000, 180000, 15000, 10000),
                V4 = c(35000, 1, 10000, 12000, 25000, 8000),
                V5 = c(10000, 7000, 1, 40000, 5000, 4000),
                V6 = c(7000, 8000, 175000, 1, 11000, 18000),
                V7 = c(70000, 30000, 22000, 120000, 1, 40000),
                V8 = c(60000, 90000, 110000, 14000, 30000, 1),
                r = c(255,255,255,153,51,51),
                g = c(51, 153, 255, 255, 255, 255),
                b = c(51, 51, 51, 51, 51, 153),
                stringsAsFactors = FALSE)
df1 <- m[, c(1,2, 9:11)]
m <- m[,-(1:2)]/1e04
m <- as.matrix(m[,c(1:6)])
dimnames(m) <- list(orig = df1$country, dest = df1$country)
#Sort order of data.frame and matrix for plotting in circos
df1 <- arrange(df1, order)
df1$country <- factor(df1$country, levels = df1$country)
m <- m[levels(df1$country),levels(df1$country)]


### Define ranges of circos sectors and their colors (both of the sectors and the links)
df1$xmin <- 0
df1$xmax <- rowSums(m) + colSums(m)
n <- nrow(df1)
df1$rcol<-rgb(df1$r, df1$g, df1$b, max = 255)
df1$lcol<-rgb(df1$r, df1$g, df1$b, alpha=200, max = 255)
mpdf('circosexample')
### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()

### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)

### Sector details
circos.initialize(factors = df1$country, xlim = cbind(df1$xmin, df1$xmax))

### Plot sectors
circos.trackPlotRegion(ylim = c(0, 1), factors = df1$country, track.height=0.1,
                       #panel.fun for each sector
                       panel.fun = function(x, y) {
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot country labels
                         circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                     col = df1$rcol[i], border=df1$rcol[i])
                         
                         #blank in part of main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3, 
                                     col = "white", border = "white")
                         
                         #white line all the way around
                         circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
                         
                         #plot axis
                         circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5), 
                                     minor.ticks=1, labels.away.percentage = 0.15)
                       })

### Plot links (inner part)
### Add sum values to df1, marking the x-position of the first links
### out (sum1) and in (sum2). Updated for further links in loop below.
df1$sum1 <- colSums(m)
df1$sum2 <- numeric(n)

### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
df2 <- cbind(as.data.frame(m),orig=rownames(m),  stringsAsFactors=FALSE)
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
               timevar="dest", time=rownames(m),  v.names = "m")
df2 <- arrange(df2,desc(m))

### Keep only the largest flows to avoid clutter
df2 <- subset(df2, m > quantile(m,0.6))

### Plot links
for(k in 1:nrow(df2)){
  #i,j reference of flow matrix
  i<-match(df2$orig[k],df1$country)
  j<-match(df2$dest[k],df1$country)
  
  #plot link
  circos.link(sector.index1=df1$country[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
              sector.index2=df1$country[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
              col = df1$lcol[i])
  
  #update sum1 and sum2 for use when plotting the next link
  df1$sum1[i] = df1$sum1[i] + abs(m[i, j])
  df1$sum2[j] = df1$sum2[j] + abs(m[i, j])
}


dev.off()


mpdf('NULL')
set.seed(123)
circos.initializeWithIdeogram(plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)
dev.off()


mpdf('NULL1')
set.seed(123)
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 10), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 10), ]

circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col="blue",  #col = rand_color(nrow(bed1),  transparency = 0.5),
                   border = 2)
dev.off()

circos.par(start.degree = -5, gap.after = c(rep(1, n_states-1), 10, rep(1, n_states-1), 10),
           cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)



mat1 =matrix(sample(20, 25, replace = TRUE), 5)
gap.degree =c(rep(1.2, 4), 1,rep(2, 4), 1)
mpdf("Circos1")
circos.clear()
circos.par(gap.degree = gap.degree, start.degree = -10/2)
chordDiagram(mat1, directional = 1, grid.col =rep(1:5, 2))
circos.clear()
dev.off()



df=df[,-1]
df[,2]=paste(df[,2],df[,1],sep = '_')

df[,5]=paste(df[,5],df[,1],sep = '_')

for (i in 1:dim(df)[1]){
  df[i,9]=max(df[i,4],df[i,7])
}
  
df2=df
df2=df2[,-1]
df2=df2[,-2]
df2=df2[,-2]
df2=df2[,-3]
df2=df2[,-3]
df2=df2[,-3]


Traitname=read.csv('/Users/wenzhang/Desktop/Trait_AGene/TraitNameCategory.csv',header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Traitname)<-Traitname$V1
Tissuename=read.csv('/Users/wenzhang/Desktop/Trait_AGene/TissueName.csv',header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Tissuename)<-Tissuename$V1



# Estimate significant and causal counts
SIGNIFICANT_COUNTS = gwasTissueSorter(table(dat[dat$isSignificant,c("tissue","trait")]))
CAUSAL_COUNTS = gwasTissueSorter(table(dat[dat$isCausal,c("eQTL","GWAS")]))
ALL_COUNTS = gwasTissueSorter(table(dat1[,c("eQTL","GWAS")]))

# Estimate normalized counts
NORMALIZED_SIGNIFICANT_COUNTS = (SIGNIFICANT_COUNTS/ALL_COUNTS[rownames(SIGNIFICANT_COUNTS),colnames(SIGNIFICANT_COUNTS)])*100
NORMALIZED_CAUSAL_COUNTS = (CAUSAL_COUNTS/ALL_COUNTS[rownames(CAUSAL_COUNTS),colnames(CAUSAL_COUNTS)])*100
if(length(unique(dat1$GWAS[dat1$isCausal]))>1){#only do this if we have 2+ traits with at least one causal var each  
  
  ## MAKE PLOT FOR SIGNIFICANT FINDINGS
  my_palette = colorRampPalette(c("blue","white","red"))(n = 1299)
  
  myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")),space="Lab")
  myPalette=colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),space="Lab")
  
  
  ## MAKE PLOT FOR CAUSAL FINDINGS
  my_palette = colorRampPalette(myColorTheme)(n = 299)
  
  #Plot HM plot
  mpdf("CAUSAL_COUNTS.NORMALIZED_RATIO", width=2+length(unique(dat1$GWAS))*0.30, height=4+length(unique(dat1$eQTL))*0.30)
  mar.default = c(6,5,5,3) + 0.1
  par(mar = mar.default + c(10, 12, 0, 0)) #c(bottom, left, top, right)
  labeledHeatmap(Matrix = apply(NORMALIZED_CAUSAL_COUNTS, 2, scale),
    xLabels = paste0(cleanGwasName(colnames(NORMALIZED_CAUSAL_COUNTS)), " (", gwasInfo[colnames(NORMALIZED_CAUSAL_COUNTS),"causalGenes"], ")"), 
    yLabels = paste0(cleanTissueName(rownames(NORMALIZED_CAUSAL_COUNTS)), " (", eqtlInfo[rownames(NORMALIZED_CAUSAL_COUNTS),"causalGenes"], ")"), 
    colorLabels = F,
    colors = my_palette,
    textMatrix = CAUSAL_COUNTS,
    setStdMargins = FALSE,
    cex.text = 0.65, 
    zlim = c(-3,3),
    naColor = "white",
    main = ""
  )
  dev.off()

  #Penalize small gwasses
  mpdf("CAUSAL_COUNTS.NORMALIZED_RATIO_PENALIZED", width=2+length(unique(dat1$GWAS))*0.30, height=4+length(unique(dat1$eQTL))*0.30)
  mar.default = c(6,5,5,3) + 0.1
  par(mar = mar.default + c(10, 12, 0, 0)) #c(bottom, left, top, right)
  labeledHeatmap(Matrix = t(t(apply(NORMALIZED_CAUSAL_COUNTS, 2, scale))*(gwasInfo[colnames(NORMALIZED_CAUSAL_COUNTS),"causalGenes"]/(10+gwasInfo[colnames(NORMALIZED_CAUSAL_COUNTS),"causalGenes"]))),
    xLabels = paste0(cleanGwasName(colnames(NORMALIZED_CAUSAL_COUNTS)), " (", gwasInfo[colnames(NORMALIZED_CAUSAL_COUNTS),"causalGenes"], ")"), 
    yLabels = paste0(cleanTissueName(rownames(NORMALIZED_CAUSAL_COUNTS)), " (", eqtlInfo[rownames(NORMALIZED_CAUSAL_COUNTS),"causalGenes"], ")"), 
    colorLabels = F, 
    colors = my_palette,
    textMatrix = CAUSAL_COUNTS,
    setStdMargins = FALSE,
    cex.text = 0.65, 
    zlim = c(-2.5,2.5),
    naColor = "white",
    main = "" 
  )
  dev.off()
 
}
###########################################################################


###########################################################################
## Make graphs showing tissue specificity of findings by trait
###########################################################################

dir.create(paste0(outDir,"/geneByTissueRaster"))


for(myGwas in unique(dat1$GWAS)){
  print(myGwas)
  if(is.na(numberOfTopGenesToPlot)){
    myGenes=sort(unique(dat1[dat1$isCausal & dat1$GWAS==myGwas,"updatedGene"]))
  }else{
    z=dat1[dat1$GWAS==myGwas & dat1$p_HET>HET_CUTOFF,]
    myGenes=z$updatedGene[order(z$p_SMR)]
    myGenes=myGenes[!duplicated(myGenes)][1:numberOfTopGenesToPlot]
    rm(z)
  }
  z=dat1[dat1$updatedGene %in% myGenes & dat1$GWAS==myGwas & !is.na(dat1$p_HET),]
  if(nrow(z)>0){
  
    #normalize (not currently used)
    overallMeanP_eQTL=mean(z$p_eQTL)
    for(myEqtlSet in unique(z$eQTL)){
      sel=z$eQTL==myEqtlSet
      z[sel,"logpfdrNorm"]=z$logpfdr[sel]/( mean(z$p_eQTL[sel])/overallMeanP_eQTL )
    }
  
    myCast=function(myCol){
      myMat=dcast(z[,c("updatedGene","eQTL",myCol)],updatedGene ~  eQTL, value.var=myCol)
      rownames(myMat)=myMat[,1];myMat[,1]=NULL
      myMat[myGenes,]
    }
    z_unmelt=myCast("logpfdr_signed")
    colOrder=TISSUE_ORDER[TISSUE_ORDER %in% colnames(z_unmelt)]
    z_unmelt=z_unmelt[,colOrder]
    z_fancyLabel=myCast("plotLabel")
    z_fancyLabel=z_fancyLabel[,colOrder]
  
    myPlotMatrix=z_unmelt
    mpdf(paste0("genes_",gwasInfo[myGwas,"CleanedNameFileSystemFriendly"]), width=6.5+0.18*ncol(z_unmelt), height=6+0.18*nrow(z_unmelt))
      mar.default = c(6,5,5,3) + 0.1
      par(mar = mar.default + c(15, 12, 0, 0))
      print(
      labeledHeatmap(Matrix = myPlotMatrix,
        xLabels = cleanTissueName(colnames(myPlotMatrix)),
        xLabelsAngle = 90,
        xLabelsAdj = 1,
        yLabels = rownames(myPlotMatrix),
        colorLabels = F,
        colors = colorRampPalette(myColorTheme)(n = 299),
        textMatrix = z_fancyLabel,
        setStdMargins = F,
        cex.text = 1,
        zlim = c(-5,5),
        naColor = "grey",
        main = cleanGwasName(myGwas)
      )
    );dev.off()


    #Raster version  
    png(paste0(outDir,"/geneByTissueRaster/rasterTraitByTissue_",gwasInfo[myGwas,"Abbreviation"],".png"), width=2*72*(6.5+0.18*ncol(z_unmelt)), height=2*72*(6.4+0.18*nrow(z_unmelt)),pointsize=2*11)
      mar.default = c(6,5,5,3) + 0.1
      par(mar = mar.default + c(15, 12, 0, 0))
      labeledHeatmap(Matrix = myPlotMatrix,
        xLabels = cleanTissueName(colnames(myPlotMatrix)),
        xLabelsAngle = 90,
        xLabelsAdj = 1,
        yLabels = rownames(myPlotMatrix),
        colorLabels = F,
        colors = colorRampPalette(myColorTheme)(n = 299),
        textMatrix = z_fancyLabel,
        setStdMargins = F,
        cex.text = 1,
        zlim = c(-5,5),
        naColor = "grey",
        main = cleanGwasName(myGwas),
        cex.main = 2 
      )
    dev.off()
  }

  rm(z,myPlotMatrix,z_unmelt,z_fancyLabel,colOrder,myCast,myGenes)
}
rm(myGwas)

###########################################################################



###########################################################################
## Uniqueness plots
###########################################################################


#############
#plot tissue specificity of genes per gwas

z=dat1[dat1$isCausal,]
z=aggregate(1:nrow(z), by=z[c("ProbeID","GWAS")], FUN=length)
z[z[,3]>10,3]=10
z=t(as.matrix(table(z[,2:3])))
z=as.data.frame.matrix(z)
z=data.frame(scale(z, center=F, scale=colSums(z)))
z$row = as.character(seq_len(nrow(z)))
rownames(z)=NULL
z = melt(z, id.vars = "row")
z$variable=as.character(z$variable)
z0=z
z$row=ordered(gsub("10","10+",as.character(z$row)),levels=c(as.character(1:9),"10+"))
z$variable=ordered(paste0(gwasInfo[z$variable,"CleanedName"]," (", gwasInfo[z$variable,"causalGenes"], ")"),levels=paste0(gwasInfo[,"CleanedName"]," (", gwasInfo[,"causalGenes"], ")"))


my_palette = colorRampPalette(c("#ccffcc", "#007700"))(n = 10)

mpdf("tissue_uniqueness");print( 
ggplot(z, aes(x=variable, y=value, fill=row)) + 
  geom_bar(stat="identity") +
  guides(fill = guide_legend(reverse=F)) +
  theme_bw() +
  theme_classic() +
  scale_fill_manual(values=my_palette,name="# of tissues\nwith eQTL") +
  ylab("Fraction of associated genes") +
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 1))
); dev.off()


z=z0
z$row=ordered(gsub("10","10+",as.character(z$row)),levels=c(as.character(1:9),"10+"))


tempNames=gwasInfo 
tempNames=tempNames[order(tempNames$CleanedName),]
tempNames=tempNames[order(tempNames$causalGenes),]
tempNames=paste0(tempNames[,"CleanedName"]," (", tempNames[,"causalGenes"], ")")

z$variable=ordered(paste0(gwasInfo[z$variable,"CleanedName"]," (", gwasInfo[z$variable,"causalGenes"], ")"),levels=tempNames)



mpdf("tissue_uniqueness_by_size");print(
ggplot(z, aes(x=variable, y=value, fill=row)) + 
  geom_bar(stat="identity") +
  guides(fill = guide_legend(reverse=F)) +
  theme_bw() +
  theme_classic() +
  scale_fill_manual(values=my_palette,name="# of tissues\nwith eQTL") +
  ylab("Fraction of associated genes") +
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 1))
);dev.off()

rm(z,z0)

#Uniqueness score (not that informative)
z=dat1[dat1$isCausal,]
z=ddply(z, c("ProbeID","GWAS"), summarise, x=length(ProbeID))
z$x=1/z$x
z=z[,2:3]

z=ddply(z,c("GWAS"), summarise, score = mean(x))

z3=z[match(gwasInfo[order(gwasInfo$causalGenes),"GWAS"],z$GWAS),]
z3=z3[!is.na(z3$GWAS),]

#see it ordered by GWAS causal count
z[match(gwasInfo[order(gwasInfo$causalGenes),"GWAS"],z$GWAS),]
mpdf("uniqueness-score");print(ggplot(z3,aes(x=ordered(GWAS,levels=z3$GWAS),y=score))+geom_bar(stat="identity")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)));dev.off()
###########################################################################


###########################################################################
## Compare tissues by causal genes they infer
###########################################################################

if(length(unique(dat1$GWAS[dat1$isCausal]))>2){ #only cluster if we have causal genes of multiple GWASes

  #################################
  #Tissue similarity causal
  tissueSimilarity=dat1[dat1$isCausal,"eQTL",drop=F]
  tissueSimilarity$probeAndGwas=apply(dat1[dat1$isCausal,c("ProbeID","GWAS")],1,paste0,collapse="|")
  tissueSimilarity$value=1
  
  unmelt_tissueSimilarity=dcast(tissueSimilarity,eQTL ~ probeAndGwas, value.var="value")
  unmelt_tissueSimilarity[is.na(unmelt_tissueSimilarity)]=0
  rownames(unmelt_tissueSimilarity)=unmelt_tissueSimilarity[,1]
  unmelt_tissueSimilarity=unmelt_tissueSimilarity[,-1]
  
  dist.mat=vegdist(unmelt_tissueSimilarity,method="jaccard",binary=T)
  
  #a:complete clustering
  myClustering=hclust(dist.mat) #, main="Tissue similarity causal")
  myClustering$labels=cleanTissueName(myClustering$labels)
  if(length(dist.mat)>2){ # don't cluster tissues using too little data
    mpdf("tissueClustering")
      plot(myClustering,main="",xlab="",sub="")
    dev.off()
  }
  
  #b:wardd clustering
  myClustering=hclust(dist.mat,method="ward.D") #, main="Tissue similarity causal")
  myClustering$labels=cleanTissueName(myClustering$labels)
  if(length(dist.mat)>2){ # don't cluster tissues using too little data
    mpdf("tissueClustering_wardd")
      plot(myClustering,main="",xlab="",sub="")
    dev.off()
  }
  rm(myClustering,tissueSimilarity,unmelt_tissueSimilarity,dist.mat)

  #################################
  #Tissue similarity nominal

  z=dat1$p_SMR<0.05 & dat1$p_HET > HET_CUTOFF & !is.na(dat1$p_HET)
  tissueSimilarity=dat1[z,"eQTL",drop=F]
  tissueSimilarity$probeAndGwas=apply(dat1[z,c("ProbeID","GWAS")],1,paste0,collapse="|")
  tissueSimilarity$value=1
  rm(z)
  
  unmelt_tissueSimilarity=dcast(tissueSimilarity,eQTL ~ probeAndGwas, value.var="value")
  unmelt_tissueSimilarity=as.data.frame(lapply(unmelt_tissueSimilarity, function(x){replace(x, is.na(x),0)}))
  rownames(unmelt_tissueSimilarity)=unmelt_tissueSimilarity[,1]
  unmelt_tissueSimilarity=unmelt_tissueSimilarity[,-1]
  
  dist.mat=vegdist(unmelt_tissueSimilarity,method="jaccard",binary=T)
  
  #a:complete
  myClustering=hclust(dist.mat)#, main="Tissue similarity nominal")
  myClustering$labels=cleanTissueName(myClustering$labels)
  mpdf("tissueClustering_nominal")
    plot(myClustering,main="",xlab="",sub="")
  dev.off()
  
  #b:wardd
  myClustering=hclust(dist.mat,method="ward.D")#, main="Tissue similarity nominal")
  myClustering$labels=cleanTissueName(myClustering$labels)
  mpdf("tissueClustering_nominal_wardd")
    plot(myClustering,main="",xlab="",sub="")
  dev.off()
  

  rm(myClustering,tissueSimilarity,unmelt_tissueSimilarity,dist.mat)
}
###########################################################################


###########################################################################
## Compare GWASes by causal genes
###########################################################################

if(length(unique(dat1$GWAS[dat1$isCausal]))>2){ #only cluster if we have causal genes of multiple GWASes
  
  #################################
  #gwas similarity causal

  gwasSimilarity=unique(dat1[dat1$isCausal,c("GWAS","ProbeID")])
  gwasSimilarity$value=1
  
  unmelt_gwasSimilarity=dcast(gwasSimilarity,GWAS ~ ProbeID, value.var="value")
  unmelt_gwasSimilarity[is.na(unmelt_gwasSimilarity)]=0
  rownames(unmelt_gwasSimilarity)=unmelt_gwasSimilarity[,1]
  unmelt_gwasSimilarity=unmelt_gwasSimilarity[,-1]
  
  dist.mat=vegdist(unmelt_gwasSimilarity,method="jaccard",binary=T)
  
  myClustering=hclust(dist.mat)
  myClustering$labels=cleanGwasName(myClustering$labels)
  mpdf("gwasClustering")
    plot(myClustering)
  dev.off()
  
  rm(myClustering,gwasSimilarity,unmelt_gwasSimilarity)


  #################################
  #gwas similarity nominal
  z=dat1$p_SMR<0.05 & dat1$p_HET > HET_CUTOFF & !is.na(dat1$p_HET)
  gwasSimilarity=unique(dat1[z,c("GWAS","ProbeID")])
  gwasSimilarity$value=1
  
  unmelt_gwasSimilarity=dcast(gwasSimilarity,GWAS ~ ProbeID, value.var="value")
  unmelt_gwasSimilarity=as.data.frame(lapply(unmelt_gwasSimilarity, function(x){replace(x, is.na(x),0)}))
  rownames(unmelt_gwasSimilarity)=unmelt_gwasSimilarity[,1]
  unmelt_gwasSimilarity=unmelt_gwasSimilarity[,-1]
  
  dist.mat=vegdist(unmelt_gwasSimilarity,method="jaccard",binary=T)
  
  #a:complete
  myClustering=hclust(dist.mat)
  myClustering$labels=cleanGwasName(myClustering$labels)
  if(length(dist.mat)>2){ # don't cluster uing one gwas
    mpdf("gwasClustering_nominal",height=10,width=15*nrow(gwasInfo)/57)
      plot(myClustering,main="",xlab="",sub="")
    dev.off()
  }
  
  #b:wardd
  myClustering=hclust(dist.mat,method="ward.D")
  myClustering$labels=cleanGwasName(myClustering$labels)
  if(length(dist.mat)>2){ # don't cluster uing one gwas
    mpdf("gwasClustering_nominal_wardd",height=10,width=15*nrow(gwasInfo)/57)
      plot(myClustering,main="",xlab="",sub="")
    dev.off()
  }
  
  rm(myClustering,gwasSimilarity,unmelt_gwasSimilarity,dist.mat)
}

###########################################################################
## Manhattan Plots
###########################################################################

z=dat1[dat1$p_HET>0.05 & !is.na(dat1$p_HET),c("GWAS","eQTL","ProbeID","SNP_Chr","SNP_bp","fdr_SMR")]
z=z[order(z$fdr_SMR),]
z$id=1:nrow(z)

myPlots=list()
myPlotNames=list()

#Max across all traits and all tissues
myName="ALL_ALL"
w=z
w=w$id[!duplicated(w$ProbeID)]
myPlots[[myName]]=w
myPlotNames[[myName]]=paste0("All traits - All tissues")

#All traits one tissue
for(i in unique(z$eQTL)){
  myName=paste0("ALL_", i)
  w=z[z$eQTL==i,]
  w=w$id[!duplicated(w$ProbeID)]
  myPlots[[myName]]=w
  myPlotNames[[myName]]=paste0("All traits - ", eqtlInfo[eqtlInfo$Sets==i,"CleanedName"])
}

#One trait all tissues
for(i in unique(z$GWAS)){
  myName=paste0(gwasInfo[i,"Abbreviation"],"_ALL")
  w=z[z$GWAS==i,]
  w=w$id[!duplicated(w$ProbeID)]
  myPlots[[myName]]=w
  myPlotNames[[myName]]=paste0(gwasInfo[i,"CleanedName"], " - All tissues")
}

#trait by tissue
for(i in unique(z$GWAS)){
  for(j in unique(z$eQTL)){
    myName=paste0(gwasInfo[i,"Abbreviation"], "_", j)
    w=z[z$GWAS==i&z$eQTL==j,]
    w=w$id[!duplicated(w$ProbeID)]
    myPlots[[myName]]=w
    myPlotNames[[myName]]=paste0(gwasInfo[i,"CleanedName"], " - ", eqtlInfo[eqtlInfo$Sets==j,"CleanedName"])
  }
}

dir.create(paste0(outDir,"/manhattanPlots"))

for(i in names(myPlots)){
  png(paste0(outDir,"/manhattanPlots/manhattan_",i,".png"),width = 2*800, height = 2*400, units = "px", pointsize = 2*12);print(
    manhattan(
      z[myPlots[[i]],], 
      chr = "SNP_Chr",
      bp = "SNP_bp",
      p = "fdr_SMR",
      snp="ProbeID",
      col = c("gray10", "gray60"),
      chrlabs = NULL,
      suggestiveline = F,
      genomewideline = -log10(0.05),
      highlight = NULL,
      logp = T,
      ylab=expression('-log'[10]*'(P-FDR)'),
      main=myPlotNames[[i]]
    )
 );dev.off()
}

rm(i,j,z,w,myPlots,myPlotNames)


#make MHC exclusion file:
MHC<-read.table('/Users/wenzhang/Desktop/Trait_AGene/makeMHC.txt',header = FALSE,stringsAsFactors = FALSE)
MHC=paste('cp getCKD-MHC.R get',MHC$V1,'-MHC.R',sep='')
write.table(MHC,file='/Users/wenzhang/Desktop/Trait_AGene/makeMHC.txt',quote=FALSE,row.names = FALSE,col.names = FALSE)



