###########################################################################
#   Make graphs of correlation and pleiotropy results
###########################################################################
# Folder to store intermediate files and annotations --
# First, copy the data files and annotations to this folder.
interFolder='/Users/panos/Documents/Publications/GENET/Analysis/'
# Comfortable output directory:
outDir='/Users/panos/Documents/Publications/GENET/Analysis/output/'

# Setup FDR cutoff for GENET and metaXcan
GENET_FDR_CUTOFF = 0.01
METAXCAN_FDR_CUTOFF = 0.01

options(stringsAsFactors = FALSE)
##  MetaXcan result data that is stored:
# The metaXcan result with removing MHC SNPs: metaXcanData. 
# The results are obtained through GENET and MetaXcan
load(paste(interFolder,'MetaXcan-MHC.RData',sep=''))
# chromosome 6 MHC region genes: chr6
load(paste(interFolder,'chr6MHC.RData',sep=''))

# 25 genes from MHC are in the MetaXcan results, so need to remove them

mhcgene<-intersect(chr6$gene_id,metaXcanData$gene)
inmhc<-(metaXcanData$gene %in% mhcgene)
datMHC=metaXcanData[!inmhc,]

## load trait annotation (trait names and the categories):
Traitname=read.csv(paste(interFolder,'TraitNameCategory.csv',sep=''),header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Traitname)<-Traitname$V1
# order the trait according to category
a=order(Traitname[unique(datMHC$trait),3])
# display the 58 traits that involved
Traitname[unique(datMHC$trait)[a],2]

# qvalue package to adjust association pvalues and filter with FDR 
suppressMessages(library(qvalue))

# set heatmap color panals
library(RColorBrewer)
myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")
my_palette=colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),space="Lab")

# to plot heatmaps
library(WGCNA)

# output pdf file function
mpdf=function(x,width=10,height=7)eval.parent(substitute({ pdf(paste0(outDir,"/plot_",gsub("(\\(|\\))","",gsub(" ","_",x)),".pdf"),width=width,height=height) }))

# filter out null pvalues: On rare occasions the LD structure as inferred from the reference panel is ill-posed. 
# On most cases there are no snps with matching alleles in the intersection of the prediction models and the GWAS.
match=is.na(datMHC$pvalue)
datMHC=datMHC[!match,]

# set fdr as significantly associated genes
# filter out genes with predictive performance qvalue > 1%
datMHC=datMHC[datMHC$pred_perf_qval<=GENET_FDR_CUTOFF,]

# Estimate FDR for metaXcan results
qobjwo <- qvalue(datMHC$pvalue, fdr.level = METAXCAN_FDR_CUTOFF)
datMHC$pva.qval=qobjwo$qvalues
datMHC$isSignificant=qobjwo$significant

# STARNET alleles are misset, zscore should be opposite
STARNETmis=grep('STARNET', datMHC$tissue)
datMHC[STARNETmis,3]=-datMHC[STARNETmis,3]
datMHC[STARNETmis,4]=-datMHC[STARNETmis,4]

# tissue annotations with full names
Tissuename=read.csv(paste(interFolder,'TissueName.csv',sep=''),header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Tissuename)<-Tissuename$V1
Tissuename$V2 = gsub('Fat', 'Adipose', Tissuename$V2)

# set of significant genes
datSigMHC=datMHC[datMHC$isSignificant,]

# number of significant genes -- to display in the heapmap text
datTabMHC=table(datMHC[datMHC$isSignificant,c("tissue","trait")])

# equivalently:
# number of significant genes
sigTabMHC=table(datSigMHC[,c("tissue","trait")])

# The number of all the genes that tested: use them to normalize the significant gene numbers and get enrichment scores
allTabMHC=table(datMHC[,c("tissue","trait")])
# normalize count numbers
NORMALIZED_SIGNIFICANT_COUNTS_MHC = (sigTabMHC/(allTabMHC[rownames(allTabMHC),colnames(allTabMHC)]))

# order the traits according to categories
indcategory=rep(0,58)
for (i in 1:length(Traitname[,3])){
  for (j in 1:length(colnames(datTabMHC))){
    if (colnames(datTabMHC)[j]==Traitname[i,1])
      indcategory[i]=j
  }
}

# scale normalized count by subtracting mean and dividing the standard deviaton (std)
matrix_normMHC=NORMALIZED_SIGNIFICANT_COUNTS_MHC
for(i in 1:14){
  for (j in 1:58){
    matrix_normMHC[i,j]=(NORMALIZED_SIGNIFICANT_COUNTS_MHC[i,j]-mean(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j]))/sd(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j])
  }
}

# unique gene names regarding to trait, tissue in the results
TRAIT=unique(datSigMHC$trait)
TISSUE=unique(datSigMHC$tissue)

# number of unique genes that contributed per trait/tissue
traitGenecountMHC=matrix(0,length(TRAIT),1)
tissueGenecountMHC=matrix(0,length(TISSUE),1)

rownames(traitGenecountMHC)<-TRAIT
rownames(tissueGenecountMHC)<-TISSUE

for (i in 1:length(TRAIT)){
  dattmp=subset(datSigMHC,datSigMHC$trait==TRAIT[i])
  traitGenecountMHC[i]=length(unique(dattmp$gene))
}

for (i in 1:length(TISSUE)){
  dattmp=subset(datSigMHC,datSigMHC$tissue==TISSUE[i])
  tissueGenecountMHC[i]=length(unique(dattmp$gene))
}

# reorder trait columns according to category
matrix_normMHC=matrix_normMHC[,indcategory]
datTabMHC=datTabMHC[,indcategory]

# plot heatmap of trait/tissue contributions of gene sets
# display tissue/trait name, followed by (number of genes in brackets)
colNamMHC<-colnames(matrix_normMHC)
colNam1MHC=Traitname[colNamMHC,2]
colNam1MHC<-paste(colNam1MHC," (",traitGenecountMHC[colNamMHC,],") ",sep='')
rowNamMHC<-rownames(matrix_normMHC)
rowNam1MHC=Tissuename[rowNamMHC,2]
rowNam1MHC<-paste(rowNam1MHC," (",tissueGenecountMHC[rowNamMHC,],") ",sep='')

# ******************************************************
#   Make the heatmap plot to show enrichment/depletion
# ******************************************************

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
               zlim = c(-2.5,2.5),
               naColor = "white",
               main = "")
dev.off()

# ******************************************************
# venn diagram to see gene contribution intersections
# ******************************************************
library(VennDiagram)
library(limma)

# get unique gene sets contributed by each tissue:
STARNET=datSigMHC[grep('STARNET', datSigMHC$tissue),]
CMC=datSigMHC[grep('CMC', datSigMHC$tissue),]
GTEx=datSigMHC[grep('GTEx', datSigMHC$tissue),]

# unique gene names for each of the tissue
STARNET=unique(STARNET$gene)
GTEx=unique(GTEx$gene)
CMC=unique(CMC$gene)

# make venn diagram of three gene sets
universe <- sort(unique(c(STARNET, GTEx, CMC)))
Counts <- matrix(0, nrow=length(universe), ncol=3)
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% STARNET
  Counts[i,2] <- universe[i] %in% GTEx
  Counts[i,3] <- universe[i] %in% CMC
}

# Name the columns with tissue names
colnames(Counts) <- c("STARNET","GTEx","CMC")

# Specify colors for the sets to make Figure 
cols<-c("Red", "Green", "Blue")
pdf(paste(outDir, 'venn.pdf', sep = ""))
vennDiagram(vennCounts(Counts), circle.col=cols)
dev.off()

######################################################################
## Cytoscape preparations for networks:
######################################################################
# get traits that will be plotted with sharing associations
# Only traits that have at least 50 associated genes will be plotted
plotTrait=subset(traitGenecountMHC,traitGenecountMHC[,1]>=50)
# store the intermediate files
write.table(plotTrait,file = paste(interFolder,'plotTraitNodeMHC.txt',sep=''),quote=F,row.names = TRUE,col.names = FALSE)
plotTrait=read.table(paste(interFolder,'plotTraitNodeMHC.txt',sep=''),header = FALSE,stringsAsFactors = FALSE)

# to get associated genes with respect to each trait
# and store them in the intermediate files 
for (i in 1:length(TRAIT)){
  dattmp=subset(datSigMHC,datSigMHC$trait==TRAIT[i])
  dattmp=dattmp[,1:2]
  write.table(dattmp,file=paste(interFolder,TRAIT[i],'_AssoGeneMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

# full names of the traits
plotTrait[,3]=Traitname[plotTrait[,1],3]

# order traits according to categories
plotTrait=plotTrait[order(plotTrait[,3]),]

# full names in the fourth column
plotTrait[,4]=Traitname[plotTrait[,1],2]

suppressMessages(library(data.table)) 

# store number/information of genes that shared between each trait-pair
sharednumber=matrix(0,dim(plotTrait)[1]*(dim(plotTrait)[1]-1)/2,1)
sharednumber=data.table(sharednumber)
write.table(sharednumber,file=paste(interFolder,'Trait_trait_sharedMHC.txt',sep=''),quote=F,row.names=F,col.names = F,sep='\t')
sharednumber=read.table(paste(interFolder,'Trait_trait_sharedMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)

# calculate possible trait-pairs and the number of genes that each trait-pair shared
sharedindex=1
for (i in 1:(dim(plotTrait)[1]-1)){
  for (j in (i+1):dim(plotTrait)[1]){
    gene1<-read.table(paste(interFolder,plotTrait[i,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste(interFolder,plotTrait[j,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    sharedgenes=intersect(gene1$gene,gene2$gene)
    # number of genes that two traits share
    sharednumber[sharedindex,1]=length(sharedgenes)
    # rownames indicate the trait pair: linked by '+' sign
    rownames(sharednumber)[sharedindex]=paste(plotTrait[i,1],plotTrait[j,1],sep = '+')
    sharedindex=sharedindex+1
  }
}

# *****************************************************************************
# Make the network file of trait category plot via correlated gene connections
# *****************************************************************************

suppressMessages(library(stringr))
sharedAsso = data.frame(str_split_fixed(as.character(rownames(sharednumber)), "\\+", 2), sharednumber$V1)
colnames(sharedAsso) = c("V1", "V2", "V3")

# set edge colors: for traits belonging to the same category, they are connected by green edge ( value 1)
# for traits in different categories, they are connected by yellow edge ( value 2 )
sharedAsso$edge=2
for (i in 1:dim(sharedAsso)[1]){
  if (Traitname[sharedAsso[i,1],3]==Traitname[sharedAsso[i,2],3]) { 
    sharedAsso[i,4]=1}
}

# set abbreviation of trait names as row names
rownames(plotTrait)<-plotTrait$V1 

# columns for node attributes
sharedAsso$node=1
sharedAsso$node2=1

for (i in 1:dim(sharedAsso)[1]){
  # node size
  sharedAsso[i,5]=plotTrait[sharedAsso[i,1],2]
  sharedAsso[i,6]=plotTrait[sharedAsso[i,2],2]
  # give node full names of the trait
  sharedAsso[i,1]=plotTrait[sharedAsso[i,1],4]
  sharedAsso[i,2]=plotTrait[sharedAsso[i,2],4]
}
colnames(sharedAsso)<-c('node1',  'node2',   'edgewidth',   'edgeColor',       'node1Size',       'node2Size')

# get Cytoscape importing file for trait association networks in terms of category 
# column1: node 1, column2: node 2, 
# edgewith indicate number of genes that both traits share
# node size is proportional to number of genes that associated with the trait
# the file is ready for importing to Cytoscape
write.table(sharedAsso,paste(interFolder,'CategoryPlotMHC.txt',sep=''),row.names = FALSE,quote = FALSE,col.names = TRUE,sep='\t')

# ****************************************************************************
# Following is to make autoimmune disorder gene-trait network
# ****************************************************************************

# get the subset of each autoimmune trait
ActDerm=subset(datSigMHC,datSigMHC$trait=='ActDerm')
SLE=subset(datSigMHC,datSigMHC$trait=='SLE')
CD=subset(datSigMHC,datSigMHC$trait=='CD')
UC=subset(datSigMHC,datSigMHC$trait=='UC')
RA=subset(datSigMHC,datSigMHC$trait=='RA')

# filter each subset with predictive performance qvalue <=0.5%
ActDerm=subset(ActDerm,ActDerm$pred_perf_qval<=0.005)
SLE=subset(SLE,SLE$pred_perf_qval<=0.005)
CD=subset(CD,CD$pred_perf_qval<=0.005)
UC=subset(UC,UC$pred_perf_qval<=0.005)
RA=subset(RA,RA$pred_perf_qval<=0.005)

# In addition, filter each subset with adjusted association pvalue (fdr) <=0.5%
ActDerm=subset(ActDerm,ActDerm$pva.qval<=0.005)
SLE=subset(SLE,SLE$pva.qval<=0.005)
CD=subset(CD,CD$pva.qval<=0.005)
UC=subset(UC,UC$pva.qval<=0.005)
RA=subset(RA,RA$pva.qval<=0.005)

# matrices to store numbers of genes whose expression was associated with multiple autoimmune traits
autoimmune=matrix(0,5,1)
autoimmune[1,1]=length(unique(ActDerm$gene))
autoimmune[2,1]=length(unique(SLE$gene))
autoimmune[3,1]=length(unique(CD$gene))
autoimmune[4,1]=length(unique(UC$gene))
autoimmune[5,1]=length(unique(RA$gene))

# node size of autoimmune disorder
rownames(autoimmune)<-c('ActDerm','SLE','CD','UC','RA')
write.table(autoimmune,paste(interFolder,'autoimmuneNodeSizeMHC.txt',sep=''),row.names = TRUE,col.names = FALSE,quote = FALSE)
autoimmune=read.table(paste(interFolder,'autoimmuneNodeSizeMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)
rownames(autoimmune)<-autoimmune$V1

# get shared genes of each autoimmune trait pairs
autoshared=matrix(0,dim(autoimmune)[1]*(dim(autoimmune)[1]-1)/2,1)
autoshared=data.table(autoshared)
write.table(autoshared,file=paste(interFolder,'Autoimmune_sharedMHC.txt',sep=''),quote=F,row.names=F,col.names = F,sep='\t')
autoshared=read.table(paste(interFolder,'Autoimmune_sharedMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)

# store the number of shared genes between each autoimmune trait pair
autosharedindex=1
for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    gene1<-read.table(paste(interFolder,autoimmune[i,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste(interFolder,autoimmune[j,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    autoimmuneshared=intersect(gene1$gene,gene2$gene)
    autoshared[autosharedindex,1]=length(autoimmuneshared)
    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
    autosharedindex=autosharedindex+1
  }
}

# edge width indicates number of genes that shared between each trait pair
write.table(autoshared,file=paste(interFolder,'Autoimmune_shared_edgenumberMHC.txt',sep=''),row.names = TRUE,col.names = FALSE,quote = FALSE)

# for each autoimmune trait, get subsets where genes have high predictive performance and significant associations
for (i in 1:dim(autoimmune)[1]){
  dattmp=subset(datSigMHC,datSigMHC$trait==autoimmune[i,1])
  # for five autoimmune disorders, filtering the genes with significant pred.perf.qvalue 0.5%
  dattmp=subset(dattmp,dattmp$pred_perf_qval<=0.005)
  dattmp=subset(dattmp,dattmp$pva.qval<=0.005)
  dattmp=dattmp[,1:2]
  write.table(dattmp,file=paste(interFolder,autoimmune[i,1],'_AssoGene_filterMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
}
rm(dattmp)

## filter genes with high effects (up or downregulated or ambiguous)
for (i in 1:(dim(autoimmune)[1])){
  gene1<-read.table(paste(interFolder,autoimmune[i,1],'_AssoGene_filterMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
  highassogene=gene1$gene
  # executive command line: 
  execucom=paste('tmp1=',autoimmune[i,1],'[',autoimmune[i,1],'$gene %in% highassogene,]',sep='')
  eval(parse(text = execucom))
  # filter out genes with low effects
  assogene=tmp1[abs(tmp1$zscore)>=mean(abs(tmp1$zscore)),]
  write.table(assogene[,1:3],file=paste(interFolder,autoimmune[i,1],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = TRUE)
}
rm(tmp1)

# shared genes of each autoimmune pair
autosharedindex=1
for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    gene1<-read.table(paste(interFolder,autoimmune[i,1],'_SharedGeneHigheffectMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste(interFolder,autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    sharedgenes=intersect(gene1$gene,gene2$gene)
    write.table(sharedgenes,file=paste(interFolder,autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,col.names = FALSE,row.names = FALSE)
    autoshared[autosharedindex,1]=length(sharedgenes)
    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
    autosharedindex=autosharedindex+1
  }
}

#************************
# Make the full file
#************************
# each row contains one autoimmune trait and one highly correlated gene
autoimmuneGeneTrait=matrix(0,sum(autoshared)*2,1)
autoimmuneGeneTrait=data.table(autoimmuneGeneTrait)
write.table(autoimmuneGeneTrait,file=paste(interFolder,'Autoimmune_Gene_correlationMHC.txt',sep=''),quote=F,row.names=F,col.names = F,sep='\t')
autoimmuneGeneTrait=read.table(paste(interFolder,'Autoimmune_Gene_correlationMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)

# index indicator
autogenetraitindex=1
for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    if (file.info(paste(interFolder,autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''))$size !=0)
    {sharedgenes<-read.table(paste(interFolder,autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
    if (dim(sharedgenes)[1]>0){ # column 1: trait, column 2: genes
      autoimmuneGeneTrait[autogenetraitindex:(dim(sharedgenes)[1]+autogenetraitindex-1),1]=autoimmune[i,1]
      autoimmuneGeneTrait[autogenetraitindex:(dim(sharedgenes)[1]+autogenetraitindex-1),2]=sharedgenes$V1
      # the gene is correlated with both of the traits
      autoimmuneGeneTrait[(dim(sharedgenes)[1]+autogenetraitindex):(2*dim(sharedgenes)[1]+autogenetraitindex-1),1]=autoimmune[j,1]
      autoimmuneGeneTrait[(dim(sharedgenes)[1]+autogenetraitindex):(2*dim(sharedgenes)[1]+autogenetraitindex-1),2]=sharedgenes$V1
      autogenetraitindex=autogenetraitindex+2*dim(sharedgenes)[1]
    }
    }
  }
}

# node size: number of genes that the trait associated
autoimmuneGeneTrait[,3]=autoimmune[autoimmuneGeneTrait$V1,2]
# trait color (corresponds to yellow)
autoimmuneGeneTrait[,4]=4
# gene color (corresponds to grey)
autoimmuneGeneTrait[,5]=5
# gene size (should be the same, anyvalue that are different from above values)
autoimmuneGeneTrait[,6]=2.5

# To get the edge width indicating the number of tissues that gene associated with traits
for(i in 1:dim(autoimmuneGeneTrait)[1]){
  # for each trait-gene pair, get the subset
  autosubset=(datSigMHC[datSigMHC$trait==autoimmuneGeneTrait[i,1],])
  autosubset=autosubset[autosubset$gene==autoimmuneGeneTrait[i,2],]
  # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
  autoimmuneGeneTrait[i,7]=dim(autosubset)[1]
}

# get the gene names and put it as the second column
for(i in 1:dim(autoimmuneGeneTrait)[1]){
  genesname=datSigMHC[datSigMHC$gene==autoimmuneGeneTrait[i,2],]
  autoimmuneGeneTrait[i,2]=genesname[1,2]
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# following advanced evaluations of up/down regulations are utilized
# store up/down regulation information
autoimmuneGeneReg=autoimmuneGeneTrait[,1:2]
# column 3 to column 16 indicate the up/down regulations in each of the 14 tissues
# column 17 indicate the judgements 
autoimmuneGeneReg[,3:17]=0
# last column name will be set later
colnames(autoimmuneGeneReg)<-c('trait','gene',unique(datSigMHC$tissue))

for(i in 1:dim(autoimmuneGeneReg)[1]){
  autosubset=(datSigMHC[datSigMHC$trait==autoimmuneGeneReg[i,1],])
  autosubset=autosubset[autosubset$gene_name==autoimmuneGeneReg[i,2],]
  for (j in 3:16){ # in each of the correlated tissue, judge if the gene is negative/positively related
    tmp=autosubset[autosubset$tissue==colnames(autoimmuneGeneReg)[j],]
    if (dim(tmp)[1]>0){ # in this tissue, the gene is correlated with the trait
      if (tmp$zscore<0){ # negatively: set 1
        autoimmuneGeneReg[i,j]=1
      }
      if (tmp$zscore>0){ # positively: set 2
        autoimmuneGeneReg[i,j]=2
      }
    }
  }
  # the number counts of ambiguous, down and up regulations:
  no0=0
  down1=0
  up2=0
  for (j in 3:16){
    if (autoimmuneGeneReg[i,j]==0) no0=no0+1
    if (autoimmuneGeneReg[i,j]==1) down1=down1+1
    if (autoimmuneGeneReg[i,j]==2) up2=up2+1
  }
  if (no0==14) autoimmuneGeneReg[i,17]=0 # neither up nor down regulated: ambiguous
  if (no0!=14){
    if (down1==up2) autoimmuneGeneReg[i,17]=0 # equal number of tissues that are up/down -- ambiguous 
    if (down1>up2) autoimmuneGeneReg[i,17]=1 # down regulation
    if (down1<up2)  autoimmuneGeneReg[i,17]=2 # up regulation
  } 
}

colnames(autoimmuneGeneReg)[17]<-'Regulation(1down/2up)'
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# set the last column to the importing file: indicator of regulations
autoimmuneGeneTrait[,8]<-autoimmuneGeneReg$`Regulation(1down/2up)`

# remove duplicated column
autoimmuneGeneTrait=autoimmuneGeneTrait[,-5]

# set the name to get attributions for each column (for cytoscape use)
colnames(autoimmuneGeneTrait)<-c('trait','gene','node1Size','Node1color','Node2Color/Size','edgeWidth','Regulation(1down/2up)')

# only include up or down regulated relations, discard ambiguous genes
autoimmuneGeneTrait<-autoimmuneGeneTrait[autoimmuneGeneTrait$`Regulation(1down/2up)`!=0,]

#############################################
# The file could be imported by Cytoscape
#############################################
write.table(autoimmuneGeneTrait,paste(interFolder,'Autoimmune_Trait_geneMHC.txt',sep=''),quote = FALSE,row.names = FALSE,col.names = TRUE)

###################################################################################
# Following is to make neuropsychiatric disorder associated gene plots analogously
###################################################################################

# each subset
PGC2_SCZ=subset(datSigMHC,datSigMHC$trait=='PGC2_SCZ')
EduYear=subset(datSigMHC,datSigMHC$trait=='EduYear')
WellBeing_Neuroticism=subset(datSigMHC,datSigMHC$trait=='WellBeing_Neuroticism')
AD=subset(datSigMHC,datSigMHC$trait=='AD')
iPSYCH_ADHD_EUR=subset(datSigMHC,datSigMHC$trait=='iPSYCH_ADHD_EUR')

# filter with 0.5% predictive performance qvalue
PGC2_SCZ=subset(PGC2_SCZ,PGC2_SCZ$pred_perf_qval<=0.005)
EduYear=subset(EduYear,EduYear$pred_perf_qval<=0.005)
WellBeing_Neuroticism=subset(WellBeing_Neuroticism,WellBeing_Neuroticism$pred_perf_qval<=0.005)
AD=subset(AD,AD$pred_perf_qval<=0.005)
iPSYCH_ADHD_EUR=subset(iPSYCH_ADHD_EUR,iPSYCH_ADHD_EUR$pred_perf_qval<=0.005)

# and 0.5% of fdr for association pvalue 
PGC2_SCZ=subset(PGC2_SCZ,PGC2_SCZ$pva.qval<=0.005)
EduYear=subset(EduYear,EduYear$pva.qval<=0.005)
WellBeing_Neuroticism=subset(WellBeing_Neuroticism,WellBeing_Neuroticism$pva.qval<=0.005)
AD=subset(AD,AD$pva.qval<=0.005)
iPSYCH_ADHD_EUR=subset(iPSYCH_ADHD_EUR,iPSYCH_ADHD_EUR$pva.qval<=0.005)

# number of genes whose expression was associated with multiple neuropsychiatric traits
neuropsy=matrix(0,5,1)
neuropsy[1,1]=length(unique(PGC2_SCZ$gene))
neuropsy[2,1]=length(unique(EduYear$gene))
neuropsy[3,1]=length(unique(WellBeing_Neuroticism$gene))
neuropsy[4,1]=length(unique(AD$gene))
neuropsy[5,1]=length(unique(iPSYCH_ADHD_EUR$gene))

# node size of neuropsychiatric disorder
rownames(neuropsy)<-c('PGC2_SCZ','EduYear','WellBeing_Neuroticism','AD','iPSYCH_ADHD_EUR')
write.table(neuropsy,paste(interFolder,'neuropsyNodeSizeMHC.txt',sep=''),row.names = TRUE,col.names = FALSE,quote = FALSE)
neuropsy=read.table(paste(interFolder,'neuropsyNodeSizeMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)
rownames(neuropsy)<-neuropsy$V1

neuroshared=matrix(0,dim(neuropsy)[1]*(dim(neuropsy)[1]-1)/2,1)
neuroshared=data.table(neuroshared)
write.table(neuroshared,file=paste(interFolder,'Neuropsy_sharedMHC.txt',sep=''),quote=F,row.names=F,col.names = F,sep='\t')
neuroshared=read.table(paste(interFolder,'Neuropsy_sharedMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)

# index indicator
neurosharedindex=1
for (i in 1:(dim(neuropsy)[1]-1)){
  for (j in (i+1):dim(neuropsy)[1]){
    gene1<-read.table(paste(interFolder,neuropsy[i,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    gene2<-read.table(paste(interFolder,neuropsy[j,1],'_AssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
    neurosharegene=intersect(gene1$gene,gene2$gene)
    write.table(neurosharegene,file=paste(interFolder,neuropsy[i,1],'+',neuropsy[j,1],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,col.names = FALSE,row.names = FALSE)
    neuroshared[neurosharedindex,1]=length(neurosharegene)
    rownames(neuroshared)[neurosharedindex]=paste(neuropsy[i,1],neuropsy[j,1],sep = '+')
    neurosharedindex=neurosharedindex+1
  }
}

# store the intermediate file for use
write.table(neuroshared,file=paste(interFolder,'Neuropsy_shared_edgenumberMHC.txt',sep=''),row.names = TRUE,col.names = FALSE,quote = FALSE)

# store neuropsychiatric trait correlated genes
for (i in 1:dim(neuropsy)[1]){
  dattmp=subset(datSigMHC,datSigMHC$trait==neuropsy[i,1])
  # for the neuropsychiatric disorders, filtering the genes with low pred.perf.qvalue 
  dattmp=subset(dattmp,dattmp$pred_perf_qval<=0.005)
  # and association qvalues
  dattmp=subset(dattmp,dattmp$pva.qval<=0.005)
  # get gene and gene names
  dattmp=dattmp[,1:2]
  write.table(dattmp,file=paste(interFolder,neuropsy[i,1],'_AssoGene_filterMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
}
rm(dattmp)

neuropsy[,1]=rownames(neuropsy)
# get number of genes that are highly correlated with the traits
for (i in 1:(dim(neuropsy)[1])){
  gene1<-read.table(paste(interFolder,rownames(neuropsy)[i],'_AssoGene_filterMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
  assogenes=gene1$gene
  # get command
  execucom=paste('tmp1=',neuropsy[i,1],'[',neuropsy[i,1],'$gene %in% assogenes,]',sep='')
  eval(parse(text = execucom))
  ## filter genes with high effects (up or downregulated or both)
  assogene=tmp1[abs(tmp1$zscore)>=mean(abs(tmp1$zscore)),]
  write.table(assogene[,1:3],file=paste(interFolder,rownames(neuropsy)[i],'_SharedGeneHigheffectMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = TRUE)
}
rm(gene1,tmp1)

# make the full file
neuroGeneTrait=matrix(0,sum(neuroshared)*2,1)
neuroGeneTrait=data.table(neuroGeneTrait)
write.table(neuroGeneTrait,file=paste(interFolder,'Neuro_Gene_correlationMHC.txt',sep=''),quote=F,row.names=F,col.names = F,sep='\t')
neuroGeneTrait=read.table(paste(interFolder,'Neuro_Gene_correlationMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)

# make first two columns: trait and genes that correlated
neurogenetraitindex=1

for (i in 1:(dim(neuropsy)[1]-1)){
  for (j in (i+1):dim(neuropsy)[1]){
    if (file.info(paste(interFolder,rownames(neuropsy)[i],'+',rownames(neuropsy)[j],'_SharedGeneHigheffectMHC.txt',sep=''))$size !=0)
    {gene1<-read.table(paste(interFolder,rownames(neuropsy)[i],'+',rownames(neuropsy)[j],'_SharedGeneHigheffectMHC.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
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
# set trait node size 
neuroGeneTrait[,3]=neuropsy[neuroGeneTrait$V1,2]
# set trait node color (corresponds to yellow)
neuroGeneTrait[,4]=4
# set gene node color (corresponds to grey)
neuroGeneTrait[,5]=5
# set gene node size (should be the same)
neuroGeneTrait[,6]=2.5

# set the trait name to abbreviations -- for convenience
neuroGeneTrait[neuroGeneTrait$V1=='SCZ',1]='PGC2_SCZ'
neuroGeneTrait[neuroGeneTrait$V1=='Eduyear',1]='EduYear'
neuroGeneTrait[neuroGeneTrait$V1=='ADHC',1]='iPSYCH_ADHD_EUR'
neuroGeneTrait[neuroGeneTrait$V1=='Neuroticism',1]='WellBeing_Neuroticism'

# get the gene names instead
for(i in 1:dim(neuroGeneTrait)[1]){
  geneset=datSigMHC[datSigMHC$gene==neuroGeneTrait[i,2],]
  neuroGeneTrait[i,2]=geneset[1,2]
}

for(i in 1:dim(neuroGeneTrait)[1]){
  sharedgenes=(datSigMHC[datSigMHC$trait==neuroGeneTrait[i,1],])
  sharedgenes=sharedgenes[sharedgenes$gene_name==neuroGeneTrait[i,2],]
  # edgewidth: number of traits (tissues via which) the gene is highly correlated with
  neuroGeneTrait[i,7]=dim(sharedgenes)[1]
  
}

# strict filtering of gene sets: absolute value of zscores >= average
datSigMHC1=datSigMHC[abs(datSigMHC$zscore)>=mean(abs(datSigMHC$zscore)),]

# advanced evaluations of up/down regulations
neuroGeneReg=neuroGeneTrait[,1:2]
# track up/down/ambiguous status in each of the 14 tissues
# column 3 to 17: indicate in each tissue, whether the gene is down/up/ambiguous regulated with the trait
neuroGeneReg[,3:17]=0
colnames(neuroGeneReg)<-c('trait','gene',unique(datSigMHC1$tissue))
for(i in 1:dim(neuroGeneReg)[1]){
  # get subsets that contain only the gene and the trait
  subsets=(datSigMHC[datSigMHC$trait==neuroGeneReg[i,1],])
  subsets=subsets[subsets$gene_name==neuroGeneReg[i,2],]
  for (j in 3:16){ # to get in each tissue, if there are associations between gene and trait
    tmp=subsets[subsets$tissue==colnames(neuroGeneReg)[j],]
    if (dim(tmp)[1]>0){ # there is association
      if (tmp$zscore<0){ # downregulation: set value 1
        neuroGeneReg[i,j]=1
      }
      if (tmp$zscore>0){ # upregulation: set value 2
        neuroGeneReg[i,j]=2
      }
    }
  }
  
  # the number counts of no, down and up regulations:
  no0=0
  down1=0
  up2=0
  for (j in 3:16){
    if (neuroGeneReg[i,j]==0) no0=no0+1 # if there is no association the tissue, will not count this tissue
    if (neuroGeneReg[i,j]==1) down1=down1+1 # downregulation
    if (neuroGeneReg[i,j]==2) up2=up2+1 # upregulation
  }
  if (no0==14) neuroGeneReg[i,17]=0 # no regulations: ambiguous
  if (no0!=14){
    if (down1==up2) neuroGeneReg[i,17]=0 # equal negative and positive effects : ambiguous
    if (down1>up2) neuroGeneReg[i,17]=1 # downregulation
    if (down1<up2)  neuroGeneReg[i,17]=2 # up regulation
  }
}

# set the column name indicating regulations
colnames(neuroGeneReg)[17]<-'Regulation(1down/2up)'

# add the regulation information in the import file
neuroGeneTrait[,8]<-neuroGeneReg$`Regulation(1down/2up)`
neuroGeneTrait=neuroGeneTrait[,-6]
# colnames for attributions
colnames(neuroGeneTrait)<-c('trait','gene','node1Size','Node1color','Node2Color/Size','edgeWidth','Regulation(1down/2up)')

# only include up/down regulated genes
neuroGeneTrait<-neuroGeneTrait[neuroGeneTrait$`Regulation(1down/2up)`!=0,]
# filter edges: number of associated tissues is >2
neuroGeneTrait<-neuroGeneTrait[neuroGeneTrait$edgeWidth>2,]

################################################
# Make the file ready for Cytoscape importing
################################################
write.table(neuroGeneTrait,paste(interFolder,'Neuro_Trait_geneMHC.txt',sep=''),quote = FALSE,row.names = FALSE,col.names = TRUE)

###############################################
# plot tissue specificity of genes per trait
###############################################
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(reshape2))
suppressMessages(library(plotly))

z=datSigMHC
z=aggregate(1:nrow(z), by=z[c("gene","trait")], FUN=length)
z[z[,3]>7,3]=7
z=t(as.matrix(table(z[,2:3])))
z=as.data.frame.matrix(z)
z=data.frame(scale(z, center=F, scale=colSums(z)))
z$row = as.character(seq_len(nrow(z)))
rownames(z)=NULL
z = melt(z, id.vars = "row")
z$variable=as.character(z$variable)
z0=z
z$row=ordered(gsub("7","7+",as.character(z$row)),levels=c(as.character(1:6),"7+"))
z$variable=paste0(Traitname[z$variable,2]," (", traitGenecountMHC[z$variable,1], ")") 

Traitname[,4]=traitGenecountMHC[Traitname$V1,1]
tempNames=Traitname
tempNames=tempNames[order(tempNames$V4),]
tempNames=paste0(tempNames[,2]," (", tempNames[,4], ")")

z$variable=ordered(z$variable,levels=tempNames)
my_palette = colorRampPalette(c("#e6e6fa", "#000090"))(n = 7)

mpdf("tissue_uniqueness");print( 
  ggplot(z, aes(x=variable, y=value, fill=row)) + 
    geom_bar(stat="identity") +
    guides(fill = guide_legend(reverse=F)) +
    theme_bw() +
    theme_classic() +
    scale_fill_manual(values=my_palette,name="Number of tissues") +
    ylab("Fraction of associated genes") +
    theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    scale_y_continuous(expand = c(0,0),limits = c(0, 1))
); dev.off()

###########################################################################################
# For genes associated with trait(s) coming from one tissue, explore the proportion
###########################################################################################

z=datSigMHC
z=aggregate(1:nrow(z), by=z[c("gene","trait")], FUN=length)
# genes associated with only one gene
zOneTissue=z[z$x==1,]

for (i in 1:dim(zOneTissue)[1]){
  dattmp=datSigMHC[datSigMHC$gene==zOneTissue[i,1],]
  dattmp=dattmp[dattmp$trait==zOneTissue[i,2],]
  zOneTissue[i,4]=dattmp$tissue
}

zOneTissue_list <- sapply(TRAIT,function(x) NULL)

t <- list(
  family = "Times New Roman",
  size = 30,
  color = 'black')

pdf(paste(outDir, "venn_unique_genes.pdf", sep = ""))
for (k in TRAIT){
  # genes that only correlate with schizophrenia and from only one tissue
  dat1=zOneTissue[zOneTissue$trait==k,]
  dat2=data.frame(Tissuename[,1:2], V3 = 0)
  # count the gene number in each tissue
  for(i in 1:dim(dat1)[1]){
    dat2[dat1[i,4],3]=dat2[dat1[i,4],3]+1
  }
  
  ### Tissue that majority of the 'one-tissue' genes come from 
  zOneTissue_list[[k]] = rbind(dat2[which.max(dat2$V3),], 
             data.frame(V1 = "Others", 
                        V2 = "Others", 
                        V3 = sum(dat2[-which.max(dat2$V3),"V3"])))
  
  ## output plots
  p <- plot_ly(zOneTissue_list[[k]],
               labels = ~paste(V2,' (',V3,') ',sep=''), values = ~V3, type = 'pie',
               textposition = 'inside',
               textinfo = 'label+percent',
               insidetextfont = list(color = '#FFFFFF',family='Times New Roman'),
               hoverinfo = 'text',
               text = ~paste(V2,' (',V3,') ',sep=''), 
               marker = list(colors = colors,
                             line = list(color = '#FFFFFF', width = 1)),
               #The 'pull' attribute can also be used to create space between the sectors
               showlegend = FALSE) %>%
    layout(title = k, 
           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),font=t)
  print(p)
}
dev.off()

##############################################################
# Trait specificity of genes per tissue
##############################################################
z=datSigMHC
z=aggregate(1:nrow(z), by=z[c("gene","tissue")], FUN=length)
z[z[,3]>6,3]=6
z=t(as.matrix(table(z[,2:3])))
z=as.data.frame.matrix(z)
z=data.frame(scale(z, center=F, scale=colSums(z)))
z$row = as.character(seq_len(nrow(z)))
rownames(z)=NULL
z = melt(z, id.vars = "row")
z$variable=as.character(z$variable)
z0=z
z$row=ordered(gsub("6","6+",as.character(z$row)),levels=c(as.character(1:5),"6+"))
z$variable=paste0(Tissuename[z$variable,2]," (", tissueGenecountMHC[z$variable,1], ")") 

Tissuename[,3]=tissueGenecountMHC[Tissuename$V1,1]
tempNames=Tissuename
tempNames=tempNames[order(tempNames$V3),]
tempNames=paste0(tempNames[,2]," (", tempNames[,3], ")")

z$variable=ordered(z$variable,levels=tempNames)

my_palette = colorRampPalette(c("#fffaf0", "#8b4513"))(n = 7)

mpdf("trait_uniqueness");print( 
  ggplot(z, aes(x=variable, y=value, fill=row)) + 
    geom_bar(stat="identity") +
    guides(fill = guide_legend(reverse=F)) +
    theme_bw() +
    theme_classic() +
    scale_fill_manual(values=my_palette,name="Number of traits") +
    ylab("Fraction of associated genes") +
    theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    scale_y_continuous(expand = c(0,0),limits = c(0, 1))
); dev.off()


###########################################################################
#          Compare tissues by associated genes they infer
###########################################################################

#################################
#Tissue similarity significant
#################################
tissueSimilarity=datMHC[datMHC$isSignificant,"tissue",drop=F]
tissueSimilarity$geneAndTrait=apply(datMHC[datMHC$isSignificant,c("gene","trait")],1,paste0,collapse="|")
tissueSimilarity$value=1

unmelt_tissueSimilarity=dcast(tissueSimilarity,tissue ~ geneAndTrait, value.var="value")
unmelt_tissueSimilarity[is.na(unmelt_tissueSimilarity)]=0
rownames(unmelt_tissueSimilarity)=unmelt_tissueSimilarity[,1]
unmelt_tissueSimilarity=unmelt_tissueSimilarity[,-1]

#source("http://bioconductor.org/biocLite.R")
#biocLite("vegan")
library(vegan)

dist.mat=vegdist(unmelt_tissueSimilarity,method="jaccard")#,binary=T)

#a:complete clustering
myClustering=hclust(dist.mat) #, main="Tissue similarity causal")
myClustering$labels=Tissuename[myClustering$labels,2]
if(length(dist.mat)>2){ # don't cluster tissues using too little data
  mpdf("tissueClustering")
  plot(myClustering,main="",xlab="",sub="")
  dev.off()
}

#b:wardd clustering
myClustering=hclust(dist.mat,method="ward.D") #, main="Tissue similarity causal")
myClustering$labels=Tissuename[myClustering$labels,2]
if(length(dist.mat)>2){ # don't cluster tissues using too little data
  mpdf("tissueClustering_wardd")
  plot(myClustering,main="",xlab="",sub="")
  dev.off()
}
rm(myClustering,tissueSimilarity,unmelt_tissueSimilarity,dist.mat)

#################################
#Tissue similarity nominal
#################################
z=datSigMHC$isSignificant
tissueSimilarity=datSigMHC[z,"tissue",drop=F]
tissueSimilarity$geneAndTrait=apply(datSigMHC[z,c("gene","trait")],1,paste0,collapse="|")
tissueSimilarity$value=1
rm(z)

unmelt_tissueSimilarity=dcast(tissueSimilarity,tissue ~ geneAndTrait, value.var="value")
unmelt_tissueSimilarity=as.data.frame(lapply(unmelt_tissueSimilarity, function(x){replace(x, is.na(x),0)}))
rownames(unmelt_tissueSimilarity)=unmelt_tissueSimilarity[,1]
unmelt_tissueSimilarity=unmelt_tissueSimilarity[,-1]

dist.mat=vegdist(unmelt_tissueSimilarity,method="jaccard")#,binary=T)

#a:complete
myClustering=hclust(dist.mat)#, main="Tissue similarity nominal")
myClustering$labels=Tissuename[myClustering$labels,2]
mpdf("tissueClustering_nominal")
plot(myClustering,main="",xlab="",sub="")
dev.off()

#b:wardd
myClustering=hclust(dist.mat,method="ward.D")#, main="Tissue similarity nominal")
myClustering$labels=Tissuename[myClustering$labels,2]
mpdf("tissueClustering_nominal_wardd")
plot(myClustering,main="",xlab="",sub="")
dev.off()

###########################################################################
#                Compare Traits by associated genes
###########################################################################

###trait similarity associations

traitSimilarity=unique(datSigMHC[datSigMHC$isSignificant,c("trait","gene")])
traitSimilarity$value=1

unmelt_traitSimilarity=dcast(traitSimilarity,trait ~ gene, value.var="value")
unmelt_traitSimilarity[is.na(unmelt_traitSimilarity)]=0
rownames(unmelt_traitSimilarity)=unmelt_traitSimilarity[,1]
unmelt_traitSimilarity=unmelt_traitSimilarity[,-1]

dist.mat=vegdist(unmelt_traitSimilarity,method="jaccard")#,binary=T)

myClustering=hclust(dist.mat)
myClustering$labels=Traitname[myClustering$labels,2]
mpdf("traitClustering")
plot(myClustering,main="",xlab="",sub="")
dev.off()

rm(myClustering,traitSimilarity,unmelt_traitSimilarity)

#################################
#gwas similarity nominal
z=datSigMHC$isSignificant
gwasSimilarity=unique(datSigMHC[z,c("trait","gene")])
gwasSimilarity$value=1

unmelt_gwasSimilarity=dcast(gwasSimilarity,trait ~ gene, value.var="value")
unmelt_gwasSimilarity=as.data.frame(lapply(unmelt_gwasSimilarity, function(x){replace(x, is.na(x),0)}))
rownames(unmelt_gwasSimilarity)=unmelt_gwasSimilarity[,1]
unmelt_gwasSimilarity=unmelt_gwasSimilarity[,-1]

dist.mat=vegdist(unmelt_gwasSimilarity,method="jaccard")#,binary=T)

#a:complete
myClustering=hclust(dist.mat)
myClustering$labels=Traitname[myClustering$labels,2]
if(length(dist.mat)>2){ # don't cluster uing one gwas
  mpdf("traitClustering_nominal",height=10,width=15*nrow(Traitname)/58)
  plot(myClustering,main="",xlab="",sub="")
  dev.off()
}

#b:wardd
myClustering=hclust(dist.mat,method="ward.D")
myClustering$labels=Traitname[myClustering$labels,2]
if(length(dist.mat)>2){ # don't cluster uing one gwas
  mpdf("traitClustering_nominal_wardd",height=10,width=15*nrow(Traitname)/58)
  plot(myClustering,main="",xlab="",sub="")
  dev.off()
}
rm(myClustering,gwasSimilarity,unmelt_gwasSimilarity,dist.mat)

################################################################################
#    Make the plot of up/down/ambiguous gene distributions for each trait
################################################################################

# for each trait, get number of genes that are up/down regulated or ambiguous 
traitupdownGene=matrix(0,length(TRAIT),3)
rownames(traitupdownGene)=TRAIT

# for this file (traitupdownGene):
# second column: number of up regulated genes
# third column: number of down regulated genes
# fourth column: number of ambiguous regulated genes
for ( i in 1:dim(traitupdownGene)[1]){
  # the subset of each trait
  dattmp=datSigMHC[datSigMHC$trait==rownames(traitupdownGene)[i],]
  # all the genes that associated with the trait
  onlyassogene=unique(dattmp$gene)
  for(j in 1:length(onlyassogene)){
    # for each gene: judge if it has up/down/ambiguous expressional change during the trait
    # get all tissues through which the gene is associated with the trait:
    setgene=dattmp[dattmp$gene==onlyassogene[j],]
    # for each gene, judge in every tissue to see if it is up/down regulated
    unordow=matrix(0,1,14)
    colnames(unordow)<-unique(datSigMHC$tissue)
    upchange=0 # number of tissues in which the gene is upregulated
    downchange=0 # number of tissues in which the gene is downregulated
    for (k in 1:14) {
      tmp=setgene[setgene$tissue==colnames(unordow)[k],]
      if (dim(tmp)[1]>0){ # in this tissue, the gene is associated with the trait
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
    # count number of regulated genes:
    if (upchange>downchange) traitupdownGene[i,1]=traitupdownGene[i,1]+1 # up
    if(upchange<downchange)  traitupdownGene[i,2]=traitupdownGene[i,2]+1 # down
    if(upchange==downchange)  traitupdownGene[i,3]=traitupdownGene[i,3]+1 # ambiguous
  }
}


# order traits according to category
indcategory=rep(0,58)
for (i in 1:length(Traitname[,3])){
  for (j in 1:length(rownames(traitupdownGene))){
    if (rownames(traitupdownGene)[j]==Traitname[i,1])
      indcategory[i]=j
  }
}
traitupdownGene=traitupdownGene[indcategory,]

write.table(traitupdownGene,paste(interFolder,'TraitGeneNumber.txt',sep=''),quote=FALSE,col.names = FALSE)
traitupdownGene=read.table(paste(interFolder,'TraitGeneNumber.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)
rownames(traitupdownGene)<-traitupdownGene$V1
# get the full name of trait 
traitupdownGene[,5]=Traitname[rownames(traitupdownGene),2]
# remove redundancy
traitupdownGene=traitupdownGene[,-1]

# make the bar plot for distributions of gene regulations regarding to each trait

# dateframe: trait name, number of genes that are up/down/ambiguous regulated
data <- data.frame(traitupdownGene[,4], traitupdownGene[,1], traitupdownGene[,2],traitupdownGene[,3])
colnames(data)<-c('Trait','Up','Down','Ambiguous')
#The default order will be alphabetized unless specified as below:
data$Trait <- factor(data$Trait, levels = data[["Trait"]])
## due to style, this plot can not be exported to pdf files: view in the 'Viewer' window of R
#mpdf("BarPlot", width=2+58*0.30, height=4+length(60)*0.30)
########## the following plot can only be viewed in the 'Viewer' subwindow
mar.default = c(6,5,5,3) + 0.1
t <- list(
  family = "sans serif",
  size = 14,
  color = 'black')

par(mar = mar.default + c(10, 12, 0, 0)) #c(bottom, left, top, right)
plot_ly(data, x = ~Trait, y = ~Down, type = 'bar', name = 'Downregulation', marker = list(color = 'rgb(49,130,189)')) %>%
  add_trace(y = ~Ambiguous, name = 'Ambiguous', marker = list(color = 'rgb(4,154,4)')) %>%
  add_trace(y = ~Up,  name = 'Upregulation', marker = list(color = 'rgb(204,0,14)')) %>%
  layout(xaxis = list(title = "", tickangle = -70,family='Times New Roman'),
         yaxis = list(title = "Number of correlated genes",family='Times New Roman'),font=t,
         margin = list(b = 225),
         barmode = 'group')
#dev.off()

##############################################################################################################
#  Following code is to make the correlated/anticorrelated plot between Years of Education and other traits
##############################################################################################################

## plot shared gene correlations between 'Year of educations' and other traits 
for (i in 1:length(TRAIT)){
  # Use high effect gene sets
  dattmp=datSigMHC[datSigMHC$trait==TRAIT[i],]
  dattmp=dattmp[1:2]
  write.table(dattmp,file=paste(interFolder,TRAIT[i],'_HighAssoGeneExpMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

# To do the correlation analysis, use only high effect genes
Edu=read.table(paste(interFolder,'EduYear_HighAssoGeneExpMHC.txt',sep=''),header=TRUE,stringsAsFactors = FALSE)

suppressMessages(library(data.table)) 
DT <- data.table(Edu)
# reduce duplicates
y<-unique(DT, by = "gene")
# keep the EduYear_AssoGene.txt as the final
write.table(y,file=paste(interFolder,'EduYear_HighAssoGeneExpMHC.txt',sep=''),quote=F,row.names=F,sep='\t')

# all other traits except EduYear
OtherTrait=TRAIT
OtherTrait=OtherTrait[-12]

EduYear=read.table(paste(interFolder,'EduYear_HighAssoGeneExpMHC.txt',sep=''),header=TRUE,stringsAsFactors = FALSE)
EduYearshared=OtherTrait
shareindex=1
for (i in 1:length(OtherTrait)){
  # genes highly correlated with other trait
  other=read.table(paste(interFolder,OtherTrait[i],'_HighAssoGeneExpMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
  # To check if there are shared genes between EduYear and other trait
  highshare<-intersect(EduYear$gene,other$gene)
  if (length(highshare)!=0){
    write.table(highshare,file=paste(interFolder,'EduYear_',OtherTrait[i],'_HighAssoGenesMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
    EduYearshared[shareindex]=OtherTrait[i]
    shareindex=shareindex+1
  }
}

## After checking results, there are 38 traits that correlated with EduYear, which are stored in the first 15 elements of EduYearshared
EduYearshared=EduYearshared[1:38]

# count numbers of genes that the 15 traits shared with EduYear
shareEduTraitcount=matrix(0,length(EduYearshared),1)
rownames(shareEduTraitcount)<-EduYearshared

for(i in 1:length(shareEduTraitcount)){
  sharedgenes<-read.table(paste(interFolder,'EduYear_',EduYearshared[i],'_HighAssoGenesMHC.txt',sep=''),header=FALSE)
  # get number of genes that shared: the node size of other traits in the plot
  shareEduTraitcount[i]=dim(sharedgenes)[1]
}

######### Make the plot file for Cytoscape import
shareEduPlot=shareEduTraitcount
write.table(shareEduPlot,file=paste(interFolder,'EduyearplotTraitMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
shareEduPlot=read.table(paste(interFolder,'EduyearplotTraitMHC.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)
index1=1
for(i in 1:length(shareEduTraitcount)){
  shareEduPlot[index1,1]=shareEduTraitcount[i]
  rownames(shareEduPlot)[index1]<-rownames(shareEduTraitcount)[i]
  index1=index1+1
}

shareEduPlot[,2]=rownames(shareEduPlot)
# number of correlated genes
shareEduPlot[,3]=0
# number of anti-correlated genes
shareEduPlot[,4]=0

# filter traits that shared at least 20 genes with EduYear
shareEduPlot=shareEduPlot[shareEduPlot$V1>=20,]

numbergene=0
uniquegene=NULL
for(i in 1:dim(shareEduPlot)[1]){
  sharedgenes<-read.table(paste(interFolder,'EduYear_',rownames(shareEduPlot)[i],'_HighAssoGenesMHC.txt',sep=''),header=FALSE)
  # get number of genes that shared: the node size of other traits in the plot
  uniquegene=union(uniquegene,sharedgenes$V1)
}

EduYear=read.table(paste(interFolder,'EduYear_HighAssoGeneExpMHC.txt',sep=''),header=TRUE,stringsAsFactors = FALSE)
dattmp1=datSigMHC[datSigMHC$gene %in% EduYear$gene,]
#dattmp1=dattmp1[dattmp1$pred_perf_qval<=0.05,]
#dattmpEduYear=dattmp1[dattmp1$pva.qval<=0.05,]
dattmpEduYear=dattmp1

for(i in 1:dim(shareEduPlot)[1]){
  sharedgenes<-read.table(paste(interFolder,'EduYear_',shareEduPlot[i,2],'_HighAssoGenesMHC.txt',sep=''),header=FALSE)
  #dattmp=datSigMHC[datSigMHC$gene %in% sharedgenes$V1,]
  #dattmp=dattmp[dattmp$pred_perf_qval<=0.05,]
  #dattmp=dattmp[dattmp$pva.qval<=0.05,]
  dattmp=datSigMHC[datSigMHC$trait==shareEduPlot[i,2],]
  #intergene=intersect(dattmp$gene,dattmpEduYear$gene)
  intergene=sharedgenes
  #dattmpEdu=datSigMHC[datSigMHC$gene %in% sharedgenes$V1,]
  #dattmpEdu=dattmpEdu[dattmpEdu$pred_perf_qval<=0.05,]
  #dattmpEdu=dattmpEdu[dattmpEdu$pva.qval<=0.05,]
  dattmpEdu=datSigMHC[datSigMHC$trait=='EduYear',]
  ## evaluation of anticorrelation or correlation, get numbers of genes that are anti-/correlated
  for(j in 1:dim(intergene)[1]){
    othergeneasso=dattmp[dattmp$gene==intergene[j,1],]
    Edugeneasso=dattmpEdu[dattmpEdu$gene==intergene[j,1],]
    # for each gene, judge if it is up/down regulated during EduYear and the other trait
    unordowEdu=matrix(0,1,14)
    upordowother=matrix(0,1,14)
    colnames(unordowEdu)<-unique(datSigMHC1$tissue)
    colnames(upordowother)<-unique(datSigMHC1$tissue)
    # compute anti-/correlated gene numbers
    anticorre=0
    corre=0
    for (k in 1:14) {
      tmp=othergeneasso[othergeneasso$tissue==colnames(upordowother)[k],]
      if (dim(tmp)[1]>0){
        if (tmp$zscore<0){ # down regulated: 1
          upordowother[k]=1
        }
        if (tmp$zscore>0){ # up regulated: 2
          upordowother[k]=2
        }
      }
      # for EduYear, do the same 
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
      # if in one tissue, one gene is only regulated during one trait (but not the other), it can not be counted as anti-/correlated
      
      if (unordowEdu[k]!=0 & upordowother[k]!=0) {
        if (unordowEdu[k]!=upordowother[k]) anticorre=anticorre+1
        if (unordowEdu[k]==upordowother[k]) corre=corre+1    }
      
    }
    # third column: number of anticorrelated expressional genes
    # fourth column: number of correlated expressional genes
    if (anticorre>corre){ shareEduPlot[i,3]=shareEduPlot[i,3]+1  }
    if (anticorre<corre){ shareEduPlot[i,4]=shareEduPlot[i,4]+1  }
  }
  
}

# full trait names:
shareEduPlot[,5]=Traitname[shareEduPlot[,2],2]

# set node names: trait + ( correlated #/ anticorrelated #)
for(i in 1:dim(shareEduPlot)[1]){
  rownames(shareEduPlot)[i]<-paste(shareEduPlot[i,5],'(',shareEduPlot[i,4],'/',shareEduPlot[i,3],') ',sep='')
}

shareEduPlot[,2]='Years of education'
shareEduPlot[,5]=length(uniquegene)

#  ++++++++++++++++++++++++++++++++++++++++++++++++++++
#     The file ready for importing cytoscape:
#  ++++++++++++++++++++++++++++++++++++++++++++++++++++
write.table(shareEduPlot,file =paste(interFolder,'EduYear_PlotMHC1.txt',sep=''),quote=FALSE,col.names = FALSE,row.names = TRUE,sep='\t')

###################################################################################################
#  To prepare files for making plots of Antagonistic Pleiotropy between One (SCZ) and Other Traits
###################################################################################################

# get upregulated and downregulated genes of Schizophrenia
datSCZ=datSigMHC[datSigMHC$trait=='PGC2_SCZ',]
# more strict cutoffs
datSCZ=datSCZ[datSCZ$pred_perf_qval<=0.005,]
datSCZ=datSCZ[datSCZ$pva.qval<=0.005,]
# highly correlated genes
SCZassogene=datSCZ[abs(datSCZ$zscore)>=mean(abs(datSCZ$zscore)),]
SCZassogene=SCZassogene[,1:2]
write.table(SCZassogene,file = paste(interFolder,'PGC2_SCZ_TopAssoGeneMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)

# all other traits
sczOthertrait=TRAIT
sczOthertrait=sczOthertrait[-36]

# get all highly correlated genes of other traits
for (i in 1:length(sczOthertrait)){
  dattmp=datSigMHC[datSigMHC$trait==sczOthertrait[i],]
  dattmp=dattmp[dattmp$pred_perf_qval<=0.005,]
  dattmp=dattmp[dattmp$pva.qval<=0.005,]
  tmp1=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]
  tmp1=tmp1[1:2]
  write.table(tmp1,file=paste(interFolder,sczOthertrait[i],'_TopAssoGeneMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

# get shared genes of SCZ and other trait
for (i in 1:length(sczOthertrait)){
  otherTrait=read.table(paste(interFolder,sczOthertrait[i],'_TopAssoGeneMHC.txt',sep=''),header=TRUE,stringsAsFactors=FALSE)
  sharedgenes<-intersect(SCZassogene$gene,otherTrait$gene)
  if (length(sharedgenes)!=0){
    write.table(sharedgenes,file=paste(interFolder,'SCZ_',sczOthertrait[i],'_TopAssoGenesMHC.txt',sep=''),quote=FALSE,row.names = FALSE,col.names = FALSE)
  }
}

## After checking outputs, there are several traits that could be plotted:
## several of them belong to neutral due to narrow pleiotropy
#sczothers=c('EduYear','GIANT_HIPadjBMI','GIANT_WCadjBMI','GIANT_WHRadjBMI','LIPIDS_HDL','LIPIDS_LDL','LIPIDS_TC','LIPIDS_TG','TAG_CPD','WellBeing_Neuroticism','iPSYCH_ADHD_EUR')
# only the following has non-neural Determine
# EduYear
dattmp=datSigMHC[datSigMHC$trait=='EduYear',]
dattmp=dattmp[dattmp$pred_perf_qval<=0.005,]
dattmp=dattmp[dattmp$pva.qval<=0.005,]
Eduassogene=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]

# iPSYCH_ADHD_EUR
dattmp=datSigMHC[datSigMHC$trait=='iPSYCH_ADHD_EUR',]
dattmp=dattmp[dattmp$pred_perf_qval<=0.005,]
dattmp=dattmp[dattmp$pva.qval<=0.005,]
iPSYCH_ADHD_EURassogene=dattmp[abs(dattmp$zscore)>=mean(abs(dattmp$zscore)),]

sczgeneall=NULL
sczOthertrait=c('EduYear','iPSYCH_ADHD_EUR')

for(i in 1:length(sczOthertrait)){
  tmp1=read.table(paste(interFolder,'SCZ_',sczOthertrait[i],'_TopAssoGenesMHC.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
  sczgeneall=union(sczgeneall,tmp1[,1])
  
}

EduYear_sczgene=read.table(paste(interFolder,'SCZ_EduYear_TopAssoGenesMHC.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)
iPSYCH_ADHD_EUR_sczgene=read.table(paste(interFolder,'SCZ_iPSYCH_ADHD_EUR_TopAssoGenesMHC.txt',sep=''),header = FALSE,stringsAsFactors  = FALSE)

# check how many genes (overlaps are counted twice!)
dim(EduYear_sczgene)[1]+dim(iPSYCH_ADHD_EUR_sczgene)[1]

## make table of Antagonistic pleiotropy information between schizophrenia and other traits
sczGeneTrait=matrix(0,13,1)
sczGeneTrait=data.table(sczGeneTrait)
write.table(sczGeneTrait,file=paste(interFolder,'SCZ_Gene_correlation.txt',sep=''),quote=F,row.names=F,col.names = F,sep='\t')
SCZGeneTrait=read.table(paste(interFolder,'SCZ_Gene_correlation.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)

sczgenetraitindex=1
for (i in 1:(length(sczOthertrait)[1])){
  if (file.info(paste(interFolder,'SCZ_',sczOthertrait[i],'_TopAssoGenesMHC.txt',sep=''))$size!=0)
  {gene1<-read.table(paste(interFolder,'SCZ_',sczOthertrait[i],'_TopAssoGenesMHC.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
  if (dim(gene1)[1]>0){
    SCZGeneTrait[sczgenetraitindex:(dim(gene1)[1]+sczgenetraitindex-1),1]=gene1$V1
    SCZGeneTrait[sczgenetraitindex:(dim(gene1)[1]+sczgenetraitindex-1),2]='PGC2_SCZ'
    SCZGeneTrait[sczgenetraitindex:(dim(gene1)[1]+sczgenetraitindex-1),3]=sczOthertrait[i]
    sczgenetraitindex=sczgenetraitindex+dim(gene1)[1]
  }
  }
}

for(i in 1:dim(SCZGeneTrait)[1]){
  # get subset first:
  sczset=(datSigMHC[datSigMHC$trait=='PGC2_SCZ',])
  sczset=sczset[sczset$gene==SCZGeneTrait[i,1],]
  # the fourth column: number of traits (tissues via which) the gene is highly correlated with SCZ
  SCZGeneTrait[i,4]=dim(sczset)[1]
  otherset=(datSigMHC[datSigMHC$trait==SCZGeneTrait[i,3],])
  otherset=otherset[otherset$gene==SCZGeneTrait[i,1],]
  # the fifth column: number of traits (tissues via which) the gene is highly correlated with the other trait
  SCZGeneTrait[i,5]=dim(otherset)[1]
}

# judge the regulation direction between each gene and SCZ/the other trait
sczGeneReg=SCZGeneTrait[,1:2]
otherGeneReg=SCZGeneTrait[,1:3]
# column 3 to column 16 indicate the up/down regulations in each of the 14 tissues
# column 17 indicate the final judgements 
sczGeneReg[,3:17]=0
# column name of the last will be set in the end
colnames(sczGeneReg)<-c('gene','trait',unique(datSigMHC$tissue))

for(i in 1:dim(sczGeneReg)[1]){
  sczsubset=(datSigMHC[datSigMHC$trait==sczGeneReg[i,2],])
  sczsubset=sczsubset[sczsubset$gene==sczGeneReg[i,1],]
  for (j in 3:16){ # in each of the correlated tissue, judge if the gene is negative/positively related
    tmp=sczsubset[sczsubset$tissue==colnames(sczGeneReg)[j],]
    if (dim(tmp)[1]>0){ # in this tissue, the gene is correlated with the trait
      if (tmp$zscore<0){ # negatively: set 1
        sczGeneReg[i,j]=1
      }
      if (tmp$zscore>0){ # positively: set 2
        sczGeneReg[i,j]=2
      }
    }
  }
  
  # the number counts of none, down and up regulations:
  no0=0
  down1=0
  up2=0
  for (j in 3:16){
    if (sczGeneReg[i,j]==0) no0=no0+1
    if (sczGeneReg[i,j]==1) down1=down1+1
    if (sczGeneReg[i,j]==2) up2=up2+1
  }
  if (no0==14) sczGeneReg[i,17]=0 # neither up nor down regulated
  if (no0!=14){
    if (down1==up2) sczGeneReg[i,17]=0 # ambiguous 
    if (down1>up2) sczGeneReg[i,17]=1 # down regulation
    if (down1<up2)  sczGeneReg[i,17]=2 # up regulation
  }
  
}
otherGeneReg[,2]=otherGeneReg[,3]
otherGeneReg[,3:17]=0
colnames(otherGeneReg)<-c('gene','trait',unique(datSigMHC$tissue))

for(i in 1:dim(otherGeneReg)[1]){
  othersubset=(datSigMHC[datSigMHC$trait==otherGeneReg[i,2],])
  othersubset=othersubset[othersubset$gene==otherGeneReg[i,1],]
  for (j in 3:16){ # in each of the correlated tissue, judge if the gene is negative/positively related
    tmp=othersubset[othersubset$tissue==colnames(otherGeneReg)[j],]
    if (dim(tmp)[1]>0){ # in this tissue, the gene is correlated with the trait
      if (tmp$zscore<0){ # negatively: set 1
        otherGeneReg[i,j]=1
      }
      if (tmp$zscore>0){ # positively: set 2
        otherGeneReg[i,j]=2
      }
    }
  }
  # the number counts of none, down and up regulations:
  no0=0
  down1=0
  up2=0
  for (j in 3:16){
    if (otherGeneReg[i,j]==0) no0=no0+1
    if (otherGeneReg[i,j]==1) down1=down1+1
    if (otherGeneReg[i,j]==2) up2=up2+1
  }
  if (no0==14) otherGeneReg[i,17]=0 # neither up nor down regulated
  if (no0!=14){
    if (down1==up2) otherGeneReg[i,17]=0 # ambiguous 
    if (down1>up2) otherGeneReg[i,17]=1 # down regulation
    if (down1<up2)  otherGeneReg[i,17]=2 # up regulation
  }
}

SCZGeneTrait[,6]=sczGeneReg[,17]
SCZGeneTrait[,7]=otherGeneReg[,17]
colnames(SCZGeneTrait)<-c('Gene','Trait/SCZ','Trait/other','NumberofAsso_SCZ','NumberofAsso_other','Regulation_SCZ(Down1/Up2)','Regulation_Other(Down1/Up2)')
sczGeneTrait=SCZGeneTrait[,1:4]
sczGeneTrait[,4]=0
colnames(sczGeneTrait)<-c('Gene','Trait/SCZ','Trait/other','isAntagonistic')
for(i in 1:dim(SCZGeneTrait)[1]){
  if(i<=12){
    if(SCZGeneTrait[i,6]!=0 & SCZGeneTrait[i,7]!=0){
      if(SCZGeneTrait[i,6]==SCZGeneTrait[i,7])
        sczGeneTrait[i,4]=1
    }
  }
  if (i==13){
    if(SCZGeneTrait[i,6]!=0 & SCZGeneTrait[i,7]!=0){
      if(SCZGeneTrait[i,6]!=SCZGeneTrait[i,7])
        sczGeneTrait[i,4]=1
    }
  }
}

for(i in 1:dim(sczGeneTrait)[1]){
  genesname=datSigMHC[datSigMHC$gene==sczGeneTrait[i,1],]
  sczGeneTrait[i,5]=genesname[1,2]
}

# Column 8: antagonistic or not
SCZGeneTrait[,8]=sczGeneTrait[,4]
colnames(SCZGeneTrait)[8]='isAntagonistic'

sczAntago=SCZGeneTrait[SCZGeneTrait$isAntagonistic==1,]
for(i in 1:dim(sczAntago)[1]){
  genesname=datSigMHC[datSigMHC$gene==sczAntago[i,1],]
  sczAntago[i,9]=genesname[1,2]
}

###########################################################################
## Make graphs showing tissue specificity of findings by trait: downstream
###########################################################################

dir.create(paste0(outDir,"/geneByTissueRaster"))
numberOfTopGenesToPlot=NA
needDetailTrait=c('EGG_BW','GIANT_HIP','GIANT_HIPadjBMI','HEIGHT','LIPIDS_LDL','LIPIDS_TC','MENARCHE','PGC2_SCZ','RA','UC','eGFRcrea')
noneedDetailTrait=setdiff(TRAIT,needDetailTrait)
#oldw <- getOption("warn")
#options(warn = -1)
for(myTrait in unique(datMHC$trait)){ #unique(datMHC$trait)){
  print(myTrait)
  if(is.na(numberOfTopGenesToPlot)){
    myGenes=sort(unique(datMHC[datMHC$isSignificant & datMHC$pred_perf_qval<=0.05 &datMHC$trait==myTrait,"gene_name"]))
  }else{
    z=datMHC[datMHC$trait==myTrait & datMHC$pva.qval<=0.05 &datMHC$pred_perf_qval<=0.05,]
    myGenes=z$gene_name[order(z$pva.qval)]
    myGenes=myGenes[!duplicated(myGenes)][1:numberOfTopGenesToPlot]
    rm(z)
  }
  z=datMHC[datMHC$gene_name %in% myGenes & datMHC$trait==myTrait & !is.na(datMHC$pva.qval),]
  z$PlotLabel[z$pvalue<0.05]="Â·"
  z$PlotLabel[z$pred_perf_pval>0.05]="-"
  z$PlotLabel[z$pva.qval<0.05]="#"
  
  if(nrow(z)>0){
    
    myCast=function(myCol){
      myMat=dcast(z[,c("gene_name","tissue",myCol)],gene_name ~  tissue, value.var=myCol)
      rownames(myMat)=myMat[,1];myMat[,1]=NULL
      myMat[myGenes,]
    }
    #z_unmelt=myCast("logpfdrNorm")
    z_unmelt=myCast("zscore")
    maxvalue=0
    minvalue=0
    for(i in 1:nrow(z_unmelt)){
      for (j in 1:ncol(z_unmelt)){
        if (dim(z[z$gene_name==rownames(z_unmelt)[i]&z$tissue==colnames(z_unmelt)[j],])[1]==0)
          z_unmelt[i,j]=NA
        if (dim(z[z$gene_name==rownames(z_unmelt)[i]&z$tissue==colnames(z_unmelt)[j],])[1]!=0)
        {tmp=z[z$gene_name==rownames(z_unmelt)[i]&z$tissue==colnames(z_unmelt)[j],]
        z_unmelt[i,j]=tmp[1,3]
        if (z_unmelt[i,j]>maxvalue) maxvalue=z_unmelt[i,j]
        if(z_unmelt[i,j]<minvalue) minvalue=z_unmelt[i,j]
        }
      }
    }
    # rescale to [-4.75,4.75]
    for(i in 1:nrow(z_unmelt)){
      for (j in 1:ncol(z_unmelt)){
        
        if (!is.na(z_unmelt[i,j]))
        {
          if(abs(maxvalue)>abs(minvalue)){
            if (z_unmelt[i,j]>0) z_unmelt[i,j]=z_unmelt[i,j]*4.75/maxvalue
            if(z_unmelt[i,j]<0) z_unmelt[i,j]=-z_unmelt[i,j]*4.75*(abs(minvalue/maxvalue))/minvalue
          }
          if(abs(maxvalue)<abs(minvalue)){
            if (z_unmelt[i,j]>0) z_unmelt[i,j]=z_unmelt[i,j]*4.75*(abs(maxvalue/minvalue))/maxvalue
            if(z_unmelt[i,j]<0) z_unmelt[i,j]=-z_unmelt[i,j]*4.75/minvalue
          }
        }
      }
    }
    
    colOrder=Tissuename[Tissuename$V1 %in% colnames(z_unmelt),1]
    z_unmelt=z_unmelt[,colOrder]
    z_fancyLabel=myCast("PlotLabel")
    z_fancyLabel=z_fancyLabel[,colOrder]
    for(i in 1:nrow(z_fancyLabel)){
      for (j in 1:ncol(z_fancyLabel)){
        if (dim(z[z$gene_name==rownames(z_fancyLabel)[i]&z$tissue==colnames(z_fancyLabel)[j],])[1]==0)
          z_fancyLabel[i,j]=''
        if (dim(z[z$gene_name==rownames(z_fancyLabel)[i]&z$tissue==colnames(z_fancyLabel)[j],])[1]!=0)
        {tmp=z[z$gene_name==rownames(z_fancyLabel)[i]&z$tissue==colnames(z_fancyLabel)[j],]
        z_fancyLabel[i,j]=tmp[1,17]
        }
      }
    }
    
    myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")
    my_palette=colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),space="Lab")
    
    myPlotMatrix=z_unmelt
    pdf(paste(outDir,'/geneByTissueRaster/',"genes_",Traitname[myTrait,1],".pdf",sep=''), width=6.5+0.18*ncol(z_unmelt), height=6+0.18*nrow(z_unmelt))
    mar.default = c(6,5,5,3) + 0.1
    par(mar = mar.default + c(15, 12, 0, 0))
    #print(
    labeledHeatmap(Matrix = myPlotMatrix,
                   xLabels = Tissuename[colnames(myPlotMatrix),2],
                   xLabelsAngle = 90,
                   xLabelsAdj = 1,
                   yLabels = rownames(myPlotMatrix),
                   colorLabels = F,
                   colors = myPalette2way(1000), #colorRampPalette(myColorTheme)(n = 299),
                   textMatrix = z_fancyLabel,
                   setStdMargins = F,
                   cex.text = 1,
                   zlim = c(-5,5),
                   naColor = "grey",
                   main = Traitname[myTrait,2]
    )
    #);
    dev.off()
    
    #png files  
    png(paste0(outDir,"/geneByTissueRaster/rasterTraitByTissue_",Traitname[myTrait,1],".png"), width=2*72*(6.5+0.18*ncol(z_unmelt)), height=2*72*(6.4+0.18*nrow(z_unmelt)),pointsize=2*11)
    mar.default = c(6,5,5,3) + 0.1
    par(mar = mar.default + c(15, 12, 0, 0))
    labeledHeatmap(Matrix = myPlotMatrix,
                   xLabels = Tissuename[colnames(myPlotMatrix),2],
                   xLabelsAngle = 90,
                   xLabelsAdj = 1,
                   yLabels = rownames(myPlotMatrix),
                   colorLabels = F,
                   colors = myPalette2way(1000), #colorRampPalette(myColorTheme)(n = 299),
                   textMatrix = z_fancyLabel,
                   setStdMargins = F,
                   cex.text = 1,
                   zlim = c(-5,5),
                   naColor = "grey",
                   main = Traitname[myTrait,2],
                   cex.main = 2 
    )
    dev.off()
  }
  
  rm(z,myPlotMatrix,z_unmelt,z_fancyLabel,colOrder,myCast,myGenes)
}
rm(myTrait)
#options(warn = oldw)
###########################################################################
# to plot genes that are interesting:
# association of the genes with traits per tissue:

dir.create(paste0(outDir,"/specialGenes"))
numberOfTopTraitsToPlot=NA
# provide a gene list:
genelist=c('ENSG00000178952','ENSG00000120088','ENSG00000161180','ENSG00000169592','ENSG00000197653','ENSG00000166949','ENSG00000140564')
for(myGene in genelist){ 
  print(myGene)
  if(is.na(numberOfTopTraitsToPlot)){
    myTraits=sort(unique(datMHC[datMHC$isSignificant & datMHC$pred_perf_qval<=0.05 &datMHC$gene==myGene,"trait"]))
  }else{
    z=datMHC[datMHC$gene==myGene & datMHC$pva.qval<=0.05 &datMHC$pred_perf_qval<=0.05,]
    myTraits=z$trait[order(z$pva.qval)]
    myTraits=myTraits[!duplicated(myTraits)][1:numberOfTopTraitsToPlot]
    rm(z)
  }
  z=datMHC[datMHC$trait %in% myTraits & datMHC$gene==myGene & !is.na(datMHC$pva.qval),]
  z$PlotLabel[z$pvalue<0.05]="Â·"
  z$PlotLabel[z$pred_perf_pval>0.05]="-"
  z$PlotLabel[z$pva.qval<0.05]="#"
  
  if(nrow(z)>0){
    
    myCast=function(myCol){
      myMat=dcast(z[,c("trait","tissue",myCol)],trait ~  tissue, value.var=myCol)
      rownames(myMat)=myMat[,1];myMat[,1]=NULL
      myMat[myTraits,]
    }
    z_unmelt=myCast("zscore")
    maxvalue=0
    minvalue=0
    for(i in 1:nrow(z_unmelt)){
      for (j in 1:ncol(z_unmelt)){
        if (dim(z[z$trait==rownames(z_unmelt)[i]&z$tissue==colnames(z_unmelt)[j],])[1]==0)
          z_unmelt[i,j]=NA
        if (dim(z[z$trait==rownames(z_unmelt)[i]&z$tissue==colnames(z_unmelt)[j],])[1]!=0)
        {tmp=z[z$trait==rownames(z_unmelt)[i]&z$tissue==colnames(z_unmelt)[j],]
        z_unmelt[i,j]=tmp[1,3]
        if (z_unmelt[i,j]>maxvalue) maxvalue=z_unmelt[i,j]
        if(z_unmelt[i,j]<minvalue) minvalue=z_unmelt[i,j]
        }
      }
    }
    # rescale to [-4.75,4.75]
    for(i in 1:nrow(z_unmelt)){
      for (j in 1:ncol(z_unmelt)){
        
        if (!is.na(z_unmelt[i,j]))
        {
          if(abs(maxvalue)>abs(minvalue)){
            if (z_unmelt[i,j]>0) z_unmelt[i,j]=z_unmelt[i,j]*4.75/maxvalue
            if(z_unmelt[i,j]<0) z_unmelt[i,j]=-z_unmelt[i,j]*4.75*abs(minvalue/maxvalue)/minvalue
          }
          if(abs(minvalue)>abs(maxvalue)){
            if (z_unmelt[i,j]>0) z_unmelt[i,j]=z_unmelt[i,j]*4.75*abs(maxvalue/minvalue)/maxvalue
            if(z_unmelt[i,j]<0) z_unmelt[i,j]=-z_unmelt[i,j]*4.75/minvalue
          }
        }
      }
    }
    
    colOrder=Tissuename[Tissuename$V1 %in% colnames(z_unmelt),1]
    z_unmelt=z_unmelt[,colOrder]
    z_fancyLabel=myCast("PlotLabel")
    z_fancyLabel=z_fancyLabel[,colOrder]
    # for(i in 1:nrow(z_fancyLabel)){
    #   for (j in 1:ncol(z_fancyLabel)){
    #     if (dim(z[z$trait==rownames(z_fancyLabel)[i]&z$tissue==colnames(z_fancyLabel)[j],])[1]==0)
    #       z_fancyLabel[i,j]=''
    #     if (dim(z[z$trait==rownames(z_fancyLabel)[i]&z$tissue==colnames(z_fancyLabel)[j],])[1]!=0)
    #     {tmp=z[z$trait==rownames(z_fancyLabel)[i]&z$tissue==colnames(z_fancyLabel)[j],]
    #     z_fancyLabel[i,j]=tmp[1,17]
    #     }
    #   }
    # }
    # 
    myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")
    my_palette=colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),space="Lab")
    
    myPlotMatrix=z_unmelt
    datmp=datMHC[datMHC$gene==myGene,]
    genenames=datmp[1,2]
    pdf(paste(outDir,'specialGenes/',genenames,".pdf",sep=''), width=6.5+0.18*ncol(z_unmelt), height=6+0.18*nrow(z_unmelt))
    mar.default = c(6,5,5,3) + 0.1
    par(mar = mar.default + c(15, 12, 0, 0))
    #print(
    
    labeledHeatmap(Matrix = myPlotMatrix,
                   xLabels = Tissuename[colnames(myPlotMatrix),2],
                   xLabelsAngle = 90,
                   xLabelsAdj = 1,
                   yLabels = Traitname[rownames(myPlotMatrix),2],
                   colorLabels = F,
                   colors = myPalette2way(1000), #colorRampPalette(myColorTheme)(n = 299),
                   textMatrix = z_fancyLabel,
                   setStdMargins = F,
                   cex.text = 1,
                   zlim = c(-5,5),
                   naColor = "grey",
                   main = paste0(genenames,' (',myGene,')')
    )
    #);
    dev.off()
    
    #png files  
    png(paste0(outDir,'specialGenes/',genenames,".png"), width=2*72*(6.5+0.18*ncol(z_unmelt)), height=2*72*(6.4+0.18*nrow(z_unmelt)),pointsize=2*11)
    mar.default = c(6,5,5,3) + 0.1
    par(mar = mar.default + c(15, 12, 0, 0))
    labeledHeatmap(Matrix = myPlotMatrix,
                   xLabels = Tissuename[colnames(myPlotMatrix),2],
                   xLabelsAngle = 90,
                   xLabelsAdj = 1,
                   yLabels = Traitname[rownames(myPlotMatrix),2],
                   colorLabels = F,
                   colors = myPalette2way(1000), #colorRampPalette(myColorTheme)(n = 299),
                   textMatrix = z_fancyLabel,
                   setStdMargins = F,
                   cex.text = 1,
                   zlim = c(-5,5),
                   naColor = "grey",
                   main = paste0(genenames,' (',myGene,')'),
                   cex.main = 1 
    )
    dev.off()
  }
  
  rm(z,myPlotMatrix,z_unmelt,z_fancyLabel,colOrder,myCast,myTraits)
}
rm(myGene)


### Correlation between PrediXcan and GENET data:


predixcanV=rep(0,14*58)
genetV=rep(0,14*58)

for (i in 1:14){
  for (j in 1:58){
    predixcanV[(i-1)*58+j]=matrixPredixcan[i,j]
    genetV[(i-1)*58+j]=matrixGenet[i,j]
  }
}

# res<-cor.test(predixcanV, genetV)
# info<-c(round(res$estimate^2,7),signif(res$p.value,8))
cor(predixcanV,genetV)


pdf('/Users/wenzhang/Desktop/correlationPrediXcan.pdf')
plot(predixcanV,genetV,pch = 16, cex = 0.8, col = rgb(0,10,200,70,maxColorValue=255), main = "Normalized values", xlab = "PrediXcan", ylab = "GENET")
abline(lm(predixcanV~genetV),col='red',lty=5)
dev.off()

#################################################################################
#     Check correlation between PrediXcan and GENET association results
#################################################################################

GENETSig=datSigMHC

####### load PrediXcan results:
interFolder1='/Users/wenzhang/Desktop/Trait_AGene1/'
##  MetaXcan result data (for PrediXcan) is stored:
# The metaXcan result with removing MHC SNPs: metaXcanData. 
load(paste(interFolder1,'Meta_PrediXcanResult.RData',sep=''))

# 24 genes from MHC are in the MetaXcan results, so remove them

mhcgene<-intersect(chr6$gene_id,MetaXcan_predixcanresult$gene)
inmhc<-(MetaXcan_predixcanresult$gene %in% mhcgene)
datMHCP=MetaXcan_predixcanresult[!inmhc,]

# filter out null pavlues
match=is.na(datMHCP$pvalue)
datMHCP=datMHCP[!match,]

# set fdr<=0.01 as significantly associated genes
qobjwo <- qvalue(datMHCP$pvalue, fdr.level = 0.01)
datMHCP$pva.qval=qobjwo$qvalues
datMHCP$isSignificant=qobjwo$significant
# filter out genes with predictive performance qvalue > 1%
datMHCP=datMHCP[datMHCP$pred_perf_qval<=0.01,]

# STARNET alleles are misset, zscore should be opposite
STARNETmis=grep('STARNET', datMHCP$tissue)
datMHCP[STARNETmis,3]=-datMHCP[STARNETmis,3]

# set of significant genes
PrediXcanSig=datMHCP[datMHCP$isSignificant,]

PrediXcanSig$resultSet='PrediXcan'
GENETSig$resultSet='GENET'

datSig=rbind(PrediXcanSig,GENETSig)

# folder to store correlation plots
dir.create(paste0(outDir,"/PGcorrelation_abs"))
dir.create(paste0(outDir,"/GENETOnly"))
dir.create(paste0(outDir,"/PrediXcanOnly"))
# source("http://bioconductor.org/biocLite.R")
biocLite("maptools")
# biocLite("plotrix")


library(dplyr)
library(maptools) 
library(plotrix)

for(trait in Traitname$V1){
  for (tissue in Tissuename$V1){
    dattmp=datSig[datSig$trait==trait &datSig$tissue== tissue,]
    
    if(dim(dattmp)[1]!=0){
      dattmp$color='grey'
      dattmp$zscore=abs(dattmp$zscore)
      for (checkgene in unique(dattmp$gene)){
        if(dim(dattmp[dattmp$gene==checkgene,])[1]!=2){
          if(dattmp[dattmp$gene==checkgene,17]=='GENET') dattmp[dattmp$gene==checkgene,18]='red'
          if(dattmp[dattmp$gene==checkgene,17]=='PrediXcan') dattmp[dattmp$gene==checkgene,18]='blue'
          
        }
      }
      
    }
    datplot=dattmp[dattmp$color=='grey',]
    
    
    if(dim(datplot)[1]!=0){
      plotgene=unique(datplot$gene)
      predixcanV=matrix(0,length(plotgene),2)
      genetV=matrix(0,length(plotgene),2)
      
      for(i in 1:length(plotgene)){
        a=datplot[datplot$gene==plotgene[i],]
        predixcanV[i,1]=a[a$resultSet=='PrediXcan',3]
        predixcanV[i,2]=a[a$resultSet=='PrediXcan',2]
        genetV[i,1]=a[a$resultSet=='GENET',3]
        genetV[i,2]=a[a$resultSet=='GENET',2]
      }
    }
    
    write.table(predixcanV,file = paste(interFolder,'predixtemp.txt',sep=''),quote=F,row.names = FALSE,col.names = FALSE)
    write.table(genetV,file = paste(interFolder,'genettemp.txt',sep=''),quote=F,row.names = FALSE,col.names = FALSE)
    predixcanV=read.table(paste(interFolder,'predixtemp.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)
    genetV=read.table(paste(interFolder,'genettemp.txt',sep=''),header=FALSE,stringsAsFactors = FALSE)
    
    datplotG=dattmp[dattmp$color=='red',3]
    plotGtext=dattmp[dattmp$color=='red',2]
    datplotGx=rep(0,length(datplotG))
    datplotP=dattmp[dattmp$color=='blue',3]
    plotPtext=dattmp[dattmp$color=='blue',2]
    datplotPy=rep(0,length(datplotP))
    maxx=max(predixcanV[,1])
    maxx=max(maxx,datplotGx,datplotP)
    maxy=max(genetV[,1])
    maxy=max(maxy,datplotG,datplotPy)
    pdf(paste(outDir,'PGcorrelation_abs/',Tissuename[tissue,1],Traitname[trait,1],".pdf",sep=''))
    
    l=NA
    g=NA
    p=NA
    if(dim(datplot)[1]!=0){
      
      plot(predixcanV[,1],genetV[,1],xlim=c(-1,maxx+0.1), ylim=c(-1,maxy+0.1),pch = 16, cex = 0.8, col = rgb(125,127,130,210,maxColorValue=255), main = paste("|z-scores| (",Traitname[trait,1], " from ", Tissuename[tissue,1],")\n correlation R= ", round(cor(predixcanV[,1],genetV[,1]),6)), xlab = "PrediXcan", ylab = "GENET")
      l=NULL
    }
    if(length(datplotG)!=0) { if (is.null(l)) {
      points(datplotGx,datplotG,pch = 16, cex = 0.8,col = rgb(250,27,30,210,maxColorValue=255)) 
      #require(plotrix)
      #thigmophobe.labels(datplotGx, datplotG, labels = plotGtext, cex=0.7, offset=0.5)
      if (length(datplotG)<=5) pointLabel(datplotGx, datplotG, labels = paste("  ", plotGtext, "  ", sep=""), cex=0.7)
      if(length(datplotG)>5)
      {
        plotGtext=plotGtext[order(datplotG,decreasing = TRUE)]
        datplotG=datplotG[order(datplotG,decreasing = TRUE)]
        pointLabel(datplotGx[1:5], datplotG[1:5], labels = paste("  ", plotGtext[1:5], "  ", sep=""), cex=0.7)
      }
      
    }
      if (!is.null(l)) {plot(datplotGx,datplotG,pch = 16, cex = 0.8,col = rgb(250,27,30,210,maxColorValue=255),main = paste0('|z-scores| (',Traitname[trait,1], ' from ', Tissuename[tissue,1],')'), xlab = "PrediXcan", ylab = "GENET")
        if (length(datplotG)<=5) pointLabel(datplotGx, datplotG, labels = paste("  ", plotGtext, "  ", sep=""), cex=0.7)
        if(length(datplotG)>5)
        {
          plotGtext=plotGtext[order(datplotG,decreasing = TRUE)]
          datplotG=datplotG[order(datplotG,decreasing = TRUE)]
          pointLabel(datplotGx[1:5], datplotG[1:5], labels = paste("  ", plotGtext[1:5], "  ", sep=""), cex=0.7)
        }
        g=NULL}}
    if(length(datplotP)!=0) {
      if(is.null(g) | is.null(l)){
        points(datplotP,datplotPy,pch = 16, cex = 0.8,col = rgb(20,27,130,210,maxColorValue=255))
        if (length(datplotP)<=5) pointLabel(datplotP, datplotPy, labels = paste("  ", plotPtext, "    ", sep="   "), cex=0.7,las=3,offset=0.75)
        if(length(datplotP)>5)
        {
          plotPtext=plotPtext[order(datplotPy,decreasing = TRUE)]
          datplotPy=datplotPy[order(datplotPy,decreasing = TRUE)]
          pointLabel(datplotP[1:5], datplotPy[1:5], labels = paste("  ", plotPtext[1:5], "    ", sep="   "), cex=0.7,las=3,offset=0.75)
        }}
      if (!is.null(g) & !is.null(l))  {plot(datplotP,datplotPy,pch = 16, cex = 0.8,col = rgb(20,27,130,210,maxColorValue=255),main = paste0('|z-scores| (',Traitname[trait,1], ' from ', Tissuename[tissue,1],')'), xlab = "PrediXcan", ylab = "GENET")
        if (length(datplotP)<=5) pointLabel(datplotP, datplotPy, labels = paste("  ", plotPtext, "    ", sep="   "), cex=0.7)
        if(length(datplotP)>5)
        {
          plotPtext=plotPtext[order(datplotPy,decreasing = TRUE)]
          datplotPy=datplotPy[order(datplotPy,decreasing = TRUE)]
          pointLabel(datplotP[1:5], datplotPy[1:5], labels = paste("  ", plotPtext[1:5], "    ", sep="   "), cex=0.7,las=3,offset=0.75)
        }
        p=NULL
      }}
    xm=lm(predixcanV[,1]~genetV[,1])
    if(!is.na(xm$coefficients[1]) & !is.na(xm$coefficients[2])) {
      if (is.null(l)) abline(lm(predixcanV[,1]~genetV[,1]),col='red',lty=5)
    }
    
    dev.off()
    if(dim(datplot)[1]==0 & length(datplotG)==0 & length(datplotP)==0){system(paste('rm ',paste(outDir,'PGcorrelation_abs/',Tissuename[tissue,1],Traitname[trait,1],".pdf",sep=''),sep=" "))}
    
  }
}



##############################################################################
# Make bar plots for gene regulations identified solely by GENET/PrediXcan
##############################################################################

# bar plots for GENET uniquely identified genes

library(ggplot2)
library(plyr)
library(reshape2)

for(i in 1:dim(Traitname)[1]){
  for (j in 1:dim(Tissuename)[1]){
    
    dattmp=datSig[datSig$trait==Traitname[i,1] &datSig$tissue== Tissuename[j,1],]
    genetgene=dattmp[dattmp$resultSet=='GENET',1:2]
    predixcangene=dattmp[dattmp$resultSet=='PrediXcan',1:2]
    GenetGeneOnly=setdiff(unique(genetgene$gene),unique(predixcangene$gene))
    if(length(GenetGeneOnly)>0){
      
      dataGenet=dattmp[dattmp$resultSet=='GENET' & dattmp$gene %in% GenetGeneOnly,]
      dataGenet=dataGenet[,1:3]
      if (dim(dataGenet)[1]>0){
        
        for(ii in 1:dim(dataGenet)[1]){
          if (dataGenet[ii,3]>0) dataGenet[ii,1]='Up'
          if (dataGenet[ii,3]<0) dataGenet[ii,1]='Down'}
        colnames(dataGenet)=c('Regulation','Gene','zScore')
        dataGenet=dataGenet[t(order(dataGenet$zScore)),]
        temGene=dataGenet$Gene
        dataGenet$Gene=ordered(dataGenet$Gene,level=temGene)
        dataGenet<-data.frame(dataGenet)
        pdf(paste(outDir,'GENETOnly/',Tissuename[j,1],Traitname[i,1],".pdf",sep=''))
        p<-ggplot(dataGenet, aes(x=Gene, y=zScore, fill=Regulation)) +
          geom_bar(stat="identity",width=0.65)+theme_minimal()+
          scale_fill_manual(values=c(rgb(13, 33, 136,250,maxColorValue=255),rgb(168, 21, 36,250,maxColorValue=255)))+
          ylab("zScore") +
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+ 
          ggtitle(paste("GENET uniquely identified genes that are significantly associated with\n",Traitname[i,2],' (from ',Tissuename[j,2],' )',sep=''))
        print(p)
        dev.off() }}
  }
}


## PrediXcan uniquely identified genes

for(i in 1:dim(Traitname)[1]){
  for (j in 1:dim(Tissuename)[1]){
    
    dattmp=datSig[datSig$trait==Traitname[i,1] &datSig$tissue== Tissuename[j,1],]
    genetgene=dattmp[dattmp$resultSet=='GENET',1:2]
    predixcangene=dattmp[dattmp$resultSet=='PrediXcan',1:2]
    PrediXcanGeneOnly=setdiff(unique(predixcangene$gene),unique(genetgene$gene))
    if(length(PrediXcanGeneOnly)>0){
      
      dataPredixcan=dattmp[dattmp$resultSet=='PrediXcan' & dattmp$gene %in% PrediXcanGeneOnly,]
      dataPredixcan=dataPredixcan[,1:3]
      if (dim(dataPredixcan)[1]>0){
        
        for(ii in 1:dim(dataPredixcan)[1]){
          if (dataPredixcan[ii,3]>0) dataPredixcan[ii,1]='Up'
          if (dataPredixcan[ii,3]<0) dataPredixcan[ii,1]='Down'}
        colnames(dataPredixcan)=c('Regulation','Gene','zScore')
        dataPredixcan=dataPredixcan[t(order(dataPredixcan$zScore)),]
        temGene=dataPredixcan$Gene
        dataPredixcan$Gene=ordered(dataPredixcan$Gene,level=temGene)
        dataPredixcan<-data.frame(dataPredixcan)
        pdf(paste(outDir,'PrediXcanOnly/',Tissuename[j,1],Traitname[i,1],".pdf",sep=''))
        p<-ggplot(dataPredixcan, aes(x=Gene, y=zScore, fill=Regulation)) +
          geom_bar(stat="identity",width=0.65)+theme_minimal()+
          scale_fill_manual(values=c(rgb(13, 33, 136,250,maxColorValue=255),rgb(168, 21, 36,250,maxColorValue=255)))+
          ylab("zScore") +
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+ 
          ggtitle(paste("PrediXcan uniquely identified genes that are significantly associated with\n",Traitname[i,2],' (from ',Tissuename[j,2],' )',sep=''))
        print(p)
        dev.off() }}
  }
}
