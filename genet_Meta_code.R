################
# GENET / META # 
################
# Requirements: around 12GB of RAM or a large enough swap to compensate.

# How to run script under linux
## Rscript [location of script]/genet_Meta_code.R args[1] args[2] args [3] args [4]
### e.g. Rscript genet_Meta_code.R /Users/wenzhang/Desktop/Trait_AGene/ /Users/wenzhang/Desktop/ 0.01 0.01

# If you plan to reuse the code please review parts of the code that have #!CAUTION!#

# Arguments
args <- commandArgs(trailingOnly = TRUE)
## args[1]: set interFolder: directory containing annotations and data files e.g. '/Users/wenzhang/Desktop/Trait_AGene/'
if (is.na(args[1])) args[1] = '/Users/wenzhang/Desktop/Trait_AGene/'; interFolder<-args[1] # folder with input data and annotations must exist
setwd (interFolder)
tempDir <- paste0(interFolder,"temp/") ; if (file.exists(tempDir) == FALSE) dir.create(tempDir)
## args[2]: set outDir: output directory e.g. '/Users/wenzhang/Desktop/'
if (is.na(args[2])) args[2] = '/Users/wenzhang/Desktop/'; outDir<-args[2] # output directory
if (file.exists(outDir) == FALSE) dir.create(outDir)
## args [3]: false discovery rate (FDR) for GENET; default :0.01
if (is.na(args[3])) args[3] = 0.01; GENET_FDR_CUTOFF <- as.numeric(args[3])
## args [4]: FDR for metaXcan; default: 0.01
if (is.na(args[4])) args[4] = 0.01; METAXCAN_FDR_CUTOFF <- as.numeric(args[4])

# Search for "Wen needs help here"

###############################################################################
#  Dependencies, sripts, global default options, parameters, folder checking  #
###############################################################################
# Setting the rest of the environment directories and folders
scriptsDir <- paste0(interFolder, "R/")
dataDir <- paste0(interFolder, "data/")
figsDir <- paste0(interFolder, "figs/")
resourcesDir <- paste0(interFolder, "resources/")
source(paste0(scriptsDir,"0.01_create_directories.R"))
gv_create_directories(interFolder)
rm(gv_create_directories)

load(paste(interFolder,'chr6MHC.RData',sep=''))
# INPUT FILES
# Check if you have all the necessary files in the data directory
required.files <- c("chr6MHC.RData", # list of genes in chr6 MHC
                    "genet.RData", 
                    "Meta_PrediXcanResult.RData", # MetaXcan results from ENet algorithm
                    "MetaXcan-MHC.RData", # MetaXcan results from WENet algorithm
                    "softpanel.csv", # contains gene lists for traits from SoftPanel
                    # http://www.isb.pku.edu.cn/SoftPanel/index.html. 
                    # Gene lists were manually acquired as described 
                    # in the header of the file
                    "TissueName.csv",
                    "TraitNameCategory.csv")
source(paste0(scriptsDir,"0.02_directory_checker.R"))
missing.files <- gv_directory_checker(required.files, dataDir)
if (length(missing.files) >= 1) {stop(paste0(missing.files, "are missing"))}
rm (gv_directory_checker, required.files, missing.files)

# Installing and loading packages
cran.packages <- c( "RColorBrewer", "ggrepel", "VennDiagram", "data.table",
                    "stringr", "ggplot2", "gridExtra", "plyr", "reshape2",
                    "vegan", "dplyr","maptools", "plotrix", "stringr",
                    "GenABEL", # for lambda estimation
                    "rgeos" # required for maptools
) 
bioc.packages<- c("qvalue", "WGCNA", "limma", "vegan")
source(paste0(scriptsDir,"0.03_install_packages.R"))
gv_install_packages(cran.packages, bioc.packages)
rm(cran.packages, bioc.packages, gv_install_packages)

## Global default options
options(stringsAsFactors = FALSE) # Every data frame you create after executing that line will not auto-convert to factors unless explicitely told to do so

## Color Palettes
# grey color-blind compatible palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
# black color-blind comaptible palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

## OtherScripts
housecleaning <- function() rm(list =ls()[!ls() %in% critical.object.list])
housecleaning_except <- function(vector.with.extra.items) 
  rm(list =ls()[!ls() %in% c(critical.object.list, vector.with.extra.items)])
mpdf=function(x,folder=figsDir,width=10,height=7)
  eval.parent(substitute({ pdf(paste0(folder,"/plot_",gsub("(\\(|\\))","",gsub(" ","_",x)),".pdf"),
                               width=width,height=height) })) # output pdf file function
mpdf_nodingbats=function(x,folder=figsDir, subfolder=character(), width=10, height=7) {
  eval.parent(substitute({
    if (file.exists(paste0(folder, subfolder)) == FALSE) dir.create(paste0(folder, subfolder))
    pdf(paste0(folder, subfolder, "/plot_",gsub("(\\(|\\))","",gsub(" ","_",x)),".pdf"), width=width, height=height, useDingbats=FALSE) 
  })) }
make_safe <- function(x) gsub("(\\(|\\))","",gsub(" ","_", gsub("/","_",x))) # function to make names compatible with paths
percent <- function(x, digits = 2, format = "f", ...) { paste0(formatC(100 * x, format = format, digits = digits, ...), "%") } # percent function
# Add label function:
#' @param xfrac The fraction over from the left side.
#' @param yfrac The fraction down from the top.
#' @param label The text to label with.
#' @param pos Position to pass to text()
#' @param ... Anything extra to pass to text(), e.g. cex, col.
add_label <- function(xfrac, yfrac, label, pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

###########################################################################
#   Load and prepare data, annotation files
###########################################################################
# require("qvalue")
critical.object.list <- c(ls(), # variables and functions as above
                          "datMHC", "datMHCP", # WENet and ENet data respectively
                          "PrediXcanSig", "datSigMHC", # significant data only
                          "critical.object.list", # this list
                          "Tissuename", "Traitname") # Annotation tables
# this list contains data objects that are needed throughout the script.

message("Loading and processing data")
# Load list of chr6 MHC genes (695 total, ~20 of them are in the MetaXcan results)
# load(paste(interFolder,'chr6MHC.RData',sep='')) # type: dataframe / name: chr6

###############################################################
# Load MetaXcan results from ENet (PrediXcan) and WENet (GENET)
datMHCP <- NA; datMHC <- NA; PrediXcanSig <- NA; datSigMHC <- NA
preparedata.metaxcan <- function(RData.file, raw, filtered.data, significantdata) {
  # filter MHC genes out of metaxcan results and subset significant genes
  # Args:
  #   RData.file: metaXcan output
  #   raw: data element contained in RData file
  #   mhcfiltereddata: object name that will hold the filtered data with MHC genes
  #   significantdata: object name that will hold the significant data
  load(paste0(interFolder, RData.file)) # Load MetaXcan results
  # Remove entries for genes that belong to the MHC region
  all.data <- raw
  mhcgene<-intersect(chr6$gene_id,all.data$gene)
  inmhc<-(all.data$gene %in% mhcgene)
  filtered.data <- all.data[!inmhc,]
  # filter out null pavlues
  match=is.na(filtered.data$pvalue)
  filtered.data=filtered.data[!match,]
  # set FDR as significantly associated genes
  # filter out genes with predictive performance qvalue > 1%
  qobjwo <- qvalue(filtered.data$pvalue, fdr.level = METAXCAN_FDR_CUTOFF)
  filtered.data$pva.qval=qobjwo$qvalues
  filtered.data$isSignificant=qobjwo$significant
  # filter out genes with predictive performance qvalue > 1%
  # use same cutoff as GENET
  filtered.data=filtered.data[filtered.data$pred_perf_qval <= GENET_FDR_CUTOFF,]
 
  # set of significant genes
  significantdata=filtered.data[filtered.data$isSignificant,]
  #assign(paste0(significantdata), significant.data, envir = .GlobalEnv)
  #housecleaning_except("chr6")
  invisible(gc())
}


preparedata.metaxcan("Meta_PrediXcanResult.RData", MetaXcan_predixcanresult, datMHCP, PrediXcanSig)
preparedata.metaxcan("MetaXcan-MHC.RData", metaXcanData, datMHC, datSigMHC)

# failed to use above function, so just load it and do the preparations

load(paste0(interFolder, "Meta_PrediXcanResult.RData")) # Load MetaXcan results
# Remove entries for genes that belong to the MHC region
mhcgene<-intersect(chr6$gene_id,MetaXcan_predixcanresult$gene)
inmhc<-(MetaXcan_predixcanresult$gene %in% mhcgene)
datMHCP <- MetaXcan_predixcanresult[!inmhc,]
# filter out null pavlues
match=is.na(datMHCP$pvalue)
datMHCP=datMHCP[!match,]
# set FDR as significantly associated genes
# filter out genes with predictive performance qvalue > 1%
qobjwo <- qvalue(datMHCP$pvalue, fdr.level = METAXCAN_FDR_CUTOFF)
datMHCP$pva.qval=qobjwo$qvalues
datMHCP$isSignificant=qobjwo$significant
# filter out genes with predictive performance qvalue > 1%
# use same cutoff as GENET
datMHCP=datMHCP[datMHCP$pred_perf_qval <= GENET_FDR_CUTOFF,]
PrediXcanSig=datMHCP[datMHCP$isSignificant,]


load(paste0(interFolder, "MetaXcan-MHC.RData")) # Load MetaXcan results
# Remove entries for genes that belong to the MHC region
mhcgene<-intersect(chr6$gene_id,metaXcanData$gene)
inmhc<-(metaXcanData$gene %in% mhcgene)
datMHC <- metaXcanData[!inmhc,]
# filter out null pavlues
match=is.na(datMHC$pvalue)
datMHC=datMHC[!match,]
# set FDR as significantly associated genes
# filter out genes with predictive performance qvalue > 1%
qobjwo <- qvalue(datMHC$pvalue, fdr.level = METAXCAN_FDR_CUTOFF)
datMHC$pva.qval=qobjwo$qvalues
datMHC$isSignificant=qobjwo$significant
# filter out genes with predictive performance qvalue > 1%
# use same cutoff as GENET
datMHC=datMHC[datMHC$pred_perf_qval <= GENET_FDR_CUTOFF,]
datSigMHC=datMHC[datMHC$isSignificant,]


rm(chr6)

# load trait annotation (trait names and the categories):
Traitname=read.csv(paste(interFolder,'TraitNameCategory.csv',sep=''),header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Traitname)<-Traitname$V1 # to be able to call by rowname for the axis titles
a=order(Traitname[unique(datMHC$trait),3]) # order the trait according to trait category
# Traitname[unique(datMHC$trait)[a],2] # display the 58 traits that involved

## Load tissue annotation file with full names
Tissuename=read.csv(paste(interFolder,'TissueName.csv',sep=''),header=FALSE,stringsAsFactors = FALSE,sep='|')
rownames(Tissuename)<-Tissuename$V1 # to be able to call by rowname for the axis titles
message("Done")

###########################################################################
#   Make graphs of correlation and pleiotropy results
###########################################################################
# library(RColorBrewer) # to set heatmap color panels
# library(WGCNA) # to plot heatmaps
# Figure 5B

message("Generating heatmap")
## set heatmap color panels
myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")
my_palette=colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),space="Lab")

## Prepare tables
datSigMHC=datMHC[datMHC$isSignificant,] # set of significant genes
datTabMHC=table(datMHC[datMHC$isSignificant,c("tissue","trait")]) # number of significant genes per tissue/trait to display in the heapmap text
sigTabMHC=table(datSigMHC[,c("tissue","trait")]) # equivalently: number of significant genes
allTabMHC=table(datMHC[,c("tissue","trait")]) # The number of all the genes that tested: use them to normalize the significant gene numbers and get enrichment scores
NORMALIZED_SIGNIFICANT_COUNTS_MHC = (sigTabMHC/(allTabMHC[rownames(allTabMHC),colnames(allTabMHC)])) # normalize count numbers

## Find order of appearance for the traits in the data file when compared to trait table
indcategory=rep(0,58) # 58 distinct traits - unique(datMHC$trait) OR unique(Traitname$V1)
for (i in 1:length(Traitname[,3])){
  for (j in 1:length(colnames(datTabMHC))){
    if (colnames(datTabMHC)[j]==Traitname[i,1])
      indcategory[i]=j
  }
} # asigns each column from table to a row from the trait table
## Find order of appearance for the tissues in the data file when compared to tissue table
indQTL=rep(0,14) # 14 distinct tissues - unique(datMHC$tissue) OR unique(Tissuename$V1)
for (i in 1:length(Tissuename[,1])){
  for (j in 1:length(rownames(datTabMHC))){
    if (rownames(datTabMHC)[j]==Tissuename[i,1])
      indQTL[i]=j
  }
} # assigns each row from the table to a row from the tissue table

# scale normalized count by subtracting mean and dividing the standard deviaton (std). This will be used for the block color of the heatmap
matrix_normMHC=NORMALIZED_SIGNIFICANT_COUNTS_MHC
for(i in 1:14){
  for (j in 1:58){
    matrix_normMHC[i,j]=(NORMALIZED_SIGNIFICANT_COUNTS_MHC[i,j]-mean(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j]))/sd(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j])
  }
}

# This could be done faster with table function - IMPROVEMENT
# Generate two matrices (one for traits and one for tissues) of significant gene counts. This will be used for the count annotation in the heatmap blocks
TRAIT=unique(datSigMHC$trait) # list of traits with the order they appear in parent data file
TISSUE=unique(datSigMHC$tissue) # list of tissues with the order they appear in parent data file
traitGenecountMHC=matrix(0,length(TRAIT),1) # make a matrix with zeros with and nrow of number of traits and ncol=1
tissueGenecountMHC=matrix(0,length(TISSUE),1) # make a matrix with zeros with nrow of number of tissues and ncol=1
rownames(traitGenecountMHC)<-TRAIT # set row names
rownames(tissueGenecountMHC)<-TISSUE # set row names
for (i in 1:length(TRAIT)){
  dattmp=subset(datSigMHC,datSigMHC$trait==TRAIT[i])
  traitGenecountMHC[i]=length(unique(dattmp$gene))
} # Calculate number of counts per trait from significant data subset and populate trait matrix with it
for (i in 1:length(TISSUE)){
  dattmp=subset(datSigMHC,datSigMHC$tissue==TISSUE[i])
  tissueGenecountMHC[i]=length(unique(dattmp$gene))
} # Calculate number of counts per tissue from significant data subset and populate tissue matrix with it
# reorder traits and tissues in the matrix and significant counts df to reflect order in trait and tissue annotation files
matrix_normMHC=matrix_normMHC[,indcategory]
datTabMHC=datTabMHC[,indcategory]
matrix_normMHC=matrix_normMHC[indQTL,]
datTabMHC=datTabMHC[indQTL,]

# plot heatmap of trait/tissue contributions of gene sets
# display tissue/trait name, followed by (number of genes in brackets)
colNamMHC<-colnames(matrix_normMHC)
colNam1MHC=Traitname[colNamMHC,2]
colNam1MHC<-paste(colNam1MHC," (",traitGenecountMHC[colNamMHC,],") ",sep='')
rowNamMHC<-rownames(matrix_normMHC)
rowNam1MHC=Tissuename[rowNamMHC,2]
rowNam1MHC<-paste(rowNam1MHC," (",tissueGenecountMHC[rowNamMHC,],") ",sep='')

#   Make the heatmap plot to show enrichment/depletion
mpdf("SIGNIFICANT_COUNTS.NORMALIZED_RATIO_countNumber1-MHC1", width=3+length(unique(datMHC$trait))*0.30, height=4.5+length(unique(datMHC$tissue))*0.30)
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
message("Done")


###########################################################################
# venn diagram to see gene contribution intersections
###########################################################################

# Figure 5A
message("Generating Venn Diagram")
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
forvenn = list(STARNET = STARNET, GTEx = GTEx, CMC=CMC) # Venn diagram reads list files
# generate venn diagram
x<-"tissue types"
plot.new()
p1.venn.diagram<-venn.diagram(
  x = forvenn,
  filename = NULL,
  col = "transparent",
  fill = c("#E69F00", "#56B4E9", "#009E73"),
  alpha = 0.5,
  label.col = c("gray33", "white", "gray33", "white", "white", "white", "gray33"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("#000000", "#000000", "#000000"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 0
)
library(grDevices)
mpdf('Venn_diagram', width = 8, height = 8)
grid.draw(p1.venn.diagram) # To save
dev.off()
message("Done")


###########################################################################
# Generate tables for Cytoscape networks and Up/down/ambiquous graph
###########################################################################
## Uses gene names instead of ENSEMBL

# Figure 7 (autoimmune) and Figure S10 (Neuropsychiatric)
## Based on 1% FDR: edge Width, and Regulation (1 = down = blue / 2 = up = Red)
## Based on 0.5% FDR: Node1 size
## Based on 0.5% FDR and higher than mean z score:  selects for Node2 (genes) to be calculated

# For figure S9 generation (uses the 1% threshold) - gene count for up/down/ambiguous zscore for each trait
# Direction is judged based on zscore; if equal number of tissues show negative and positive zscores then the direction is considered ambiguous for that gene

cbPalette.upnodown <- cbPalette[c(6,4,7)] # color pallete for barplot down, ambi, up
library(data.table)

cyto.trait.temp <- sort(unique(Traitname$V3))
up.down.matrix <- Traitname[with(Traitname, order (Traitname$V3, Traitname$V2)), ]  # this is for figure S9
up.down.matrix$up <- 0 ; up.down.matrix$down <- 0; up.down.matrix$ambi <- 0
for (j in 1:length(cyto.trait.temp)) {
  trait.name <- cyto.trait.temp[j] # get trait category name  
  cyto.temp <- sort(Traitname[Traitname$V3 == trait.name , 1]) # get trait category trait names
  trait.matrix <- matrix(0,length(cyto.temp),1) # matrices to store numbers of genes whose expression was associated with multiple traits
  rownames(trait.matrix) <- cyto.temp # row names for matrix
  for (i in 1: length(cyto.temp)) {
    cyto.temp.1 <- subset(datSigMHC, datSigMHC$trait == cyto.temp[i]) # FDR 1% will be used for edge Width
    # filter each subset with predictive performance qvalue <= 0.5% and adjusted association pvalue (fdr) <= 0.5%
    cyto.temp.2 <- cyto.temp.1[cyto.temp.1$pred_perf_qval <= 0.01 & cyto.temp.1$pva.qval <= 0.01, ] # will be used for node size
    cyto.temp.3 <- cyto.temp.2[abs(cyto.temp.2$zscore) >= mean(abs(cyto.temp.2$zscore)), ] # will be used for gene nodes
    trait.matrix[i,1] = length(unique(cyto.temp.2$gene))
    # Calculate number of up/down/ambiguous regulated genes per trait 
    if (nrow(cyto.temp.1) == 0) {} else { # if there are no genes move on...
      for (k in 1:length(unique(cyto.temp.1$gene))) { # for each unique gene...
        no0=0; down1=0; up2=0 # to keep track of regulation direction per tissue
        all.tissues <- unique(cyto.temp.1$tissue) # how many tissues
        for (m in 1:length(all.tissues)) { # for every tissue
          reg.temp <- cyto.temp.1[cyto.temp.1$gene == unique(cyto.temp.1$gene)[k] & cyto.temp.1$tissue == unique(cyto.temp.1$tissue)[m] , "zscore"] # get zscore for specific tissue and gene
          ifelse (reg.temp > 0, up2 <- up2+1, ifelse(reg.temp < 0, down1 <- down1 + 1, no0 <- no0 + 1))
        } # finds in how many tissues it is up or downregulated
        ifelse (up2>down1, up.down.matrix[cyto.temp[i], "up"] <- up.down.matrix[cyto.temp[i], "up"] + 1, 
                ifelse(down1>up2, up.down.matrix[cyto.temp[i], "down"] <- up.down.matrix[cyto.temp[i], "down"] + 1,
                       up.down.matrix[cyto.temp[i], "ambi"] <- up.down.matrix[cyto.temp[i], "ambi"] + 1)) # keeps score for each gene
      }
    } # calculation of up/down/ambiguous regulated genes per trait done.
    assign(cyto.temp[i], cyto.temp.1, pos = .GlobalEnv) # data with 1% threshold
    assign(paste0(cyto.temp[i],".filtered"), cyto.temp.2, pos = .GlobalEnv) # data with 0.5% threshold
    assign(paste0(cyto.temp[i],".higheffect"), cyto.temp.3, pos = .GlobalEnv) # data with 0.5% threshold and high effect
    assign(paste0(trait.name,".matrix"), trait.matrix, pos = .GlobalEnv) # get matrix out (has the 0.5% threshold)
  } # this loop creates a subset table for all traits within trait category and names it as trait + creates a matrix with common genes
  write.table(get(paste0(trait.name,".matrix")), paste0(tempDir, "cyto_", make_safe(trait.name),'NodeSizeMHC.txt'),row.names = TRUE,col.names = FALSE,quote = FALSE) # node size is based on 0.5% threshold
  # Pairwise comparisons section starts here
  #!CAUTION!# START
  if (trait.name == "Autoimmune") { cyto.temp <- cyto.temp[!grepl("IBD", cyto.temp)] } # remove IBD if autoimmune, makes network graph too heavy
  if (trait.name == "Neuropsychiatric") { cyto.temp <- c("PGC2_SCZ", "EduYear", "WellBeing_Neuroticism", "AD", "iPSYCH_ADHD_EUR") } # custom selection for publication
  #!CAUTION!# END
  combo.temp <- combn(cyto.temp, 2) # generate list of pairwise comparisons needed within trait category.
  autoshared <- data.table(matrix(0, nrow = dim(combo.temp)[2],1)) # setting as data.table is necessary at this step do not remove it or change to data.frame
  autoshared$V2 <- NA
  autoshared <- cbind (name = NA, autoshared) # to add a column with name, data.tables do not support row names very well.
  autoshared.filtered <- autoshared # pairwise comparison w 0.5% threshold
  autoshared.higheffect <- autoshared # pairwise comparison w zscore threshold will be used for
  trait.gene.map <- matrix(0, 0, 7) # final cytoscape output
  colnames(trait.gene.map) <- c("trait", "gene", "Node1Size", "Node1Color", "Node2Color/Size", "edgeWidth","Regulation(1down/2up)")
  # Count the number of shared genes between all pairwise comparisons with different thresholds as defined above
  for (i in 1:dim(combo.temp)[2]){ 
    trait.1 <- get(paste0(combo.temp[1, i])) # loads trait 1 gene table of pairwise comparison 1% threshold  
    trait.2 <- get(paste0(combo.temp[2, i])) # loads trait 1 gene table of pairwise comparison 1% threshold
    trait.1.filtered <- get(paste0(combo.temp[1, i],".filtered")) # loads trait 1 gene table of pairwise comparison 0.5% threshold
    trait.2.filtered <- get(paste0(combo.temp[2, i],".filtered")) # loads trait 1 gene table of pairwise comparison 0.5% threshold
    trait.1.higheffect <- get(paste0(combo.temp[1, i],".higheffect")) # loads trait 1 gene table of pairwise comparison 0.5% & zscore threshold
    trait.2.higheffect <- get(paste0(combo.temp[2, i],".higheffect")) # loads trait 1 gene table of pairwise comparison 0.5% & zscore threshold
    autoshared$name[i] <- paste0(as.character(combo.temp[1, i]),"+",as.character(combo.temp[2, i]))
    autoshared$V1[i] <- length(intersect(trait.1$gene, trait.2$gene))
    autoshared$V2[i] <- paste(unlist(intersect(trait.1$gene, trait.2$gene)), collapse=",")
    assign(paste0(trait.name,".pairwise"), autoshared, pos = .GlobalEnv) # pairwise comparison w 1% threshold, used for edge number?
    autoshared.filtered$name[i] <- paste0(as.character(combo.temp[1, i]),"+",as.character(combo.temp[2, i]))
    autoshared.filtered$V1[i] <- length(intersect(trait.1.filtered$gene, trait.2.filtered$gene))
    autoshared.filtered$V2[i] <- paste(unlist(intersect(trait.1.filtered$gene, trait.2.filtered$gene)), collapse=",")
    assign(paste0(trait.name,".pairwise.filtered"), autoshared.filtered, pos = .GlobalEnv) # pairwise comparison w 0.5% threshold
    autoshared.higheffect$name[i] <- paste0(as.character(combo.temp[1, i]),"+",as.character(combo.temp[2, i]))
    autoshared.higheffect$V1[i] <- length(intersect(trait.1.higheffect$gene, trait.2.higheffect$gene))
    autoshared.higheffect$V2[i] <- paste(unlist(intersect(trait.1.higheffect$gene, trait.2.higheffect$gene)), collapse=",")
    assign(paste0(trait.name,".pairwise.higheffect"), autoshared.higheffect, pos = .GlobalEnv)
    # Final cytoscape file
    intersect.higheffect <- intersect(trait.1.higheffect$gene_name, trait.2.higheffect$gene_name) # !!! be careful here what if it 0
    trait.gene.temp <- matrix(0, 2*length(intersect.higheffect), 7) # high z score genes
    if (nrow(trait.gene.temp) == 0) { 
    } else {
      for (k in 1:(2*length(intersect.higheffect))) {
        if (k <= length(intersect.higheffect)) {
          gene.list <- get(paste0(as.character(combo.temp[1, i]))) # gets trait list genes (1% FDR)
          trait.gene.temp [k, 1] <- paste0(as.character(combo.temp[1, i]))
          trait.gene.temp [k, 2] <- intersect.higheffect[k]
          trait.gene.temp [k, 3] <- trait.matrix[row.names(trait.matrix) == paste0(as.character(combo.temp[1, i])), 1]
          trait.gene.temp [k, 6] <- nrow(gene.list[gene.list$gene_name == intersect.higheffect[k], ])
          gene.list <- gene.list[gene.list$gene_name == intersect.higheffect[k], ] # subset for trait
        } else {
          gene.list <- get(paste0(as.character(combo.temp[2, i]))) # gets trait list genes (1% FDR)
          trait.gene.temp [k, 1] <- paste0(as.character(combo.temp[2, i]))
          l = k - length(intersect.higheffect) # for the other trait
          trait.gene.temp [k, 2] <- intersect.higheffect[l]
          trait.gene.temp [k, 3] <- trait.matrix[row.names(trait.matrix) == paste0(as.character(combo.temp[2, i])), 1]
          trait.gene.temp [k, 6] <- nrow(gene.list[gene.list$gene_name == intersect.higheffect[l], ])
          gene.list <- gene.list[gene.list$gene_name == intersect.higheffect[l], ] # subset for trait
        }
        trait.gene.temp [k, 4] <- 4 # corresponds to yellow
        trait.gene.temp [k, 5] <- 2.5 # gray is 5 and 2.5 is size
        no0=0; down1=0; up2=0 # to keep track of regulation direction per tissue
        all.tissues <- unique(datSigMHC$tissue)
        # We will check down or upregulation of the specific gene in every tissue
        for (m in 1:length(all.tissues)) {
          avg.zscore <- mean(gene.list[gene.list$tissue == all.tissues[m],][,3])
          ifelse (avg.zscore > 0, up2 <- up2+1, ifelse(avg.zscore < 0, down1 <- down1 + 1, no0 <- no0 + 1))
        }
        trait.gene.temp [k, 7] <- ifelse (no0 == 14 | down1 == up2, 0, ifelse(down1>up2, 1, 2))
      } # populates trait.gene.temp
      # assign(paste0(trait.name,".nodes.MHC"), trait.gene.map, pos = .GlobalEnv)
    }
    trait.gene.map <- rbind(trait.gene.map, trait.gene.temp) # includes resultws in trait.gene.map
  }
  trait.gene.map <- trait.gene.map[trait.gene.map[, 7]!=0, ] # only include up or down regulated relations
  write.table(trait.gene.map, paste0(outDir, make_safe(trait.name), '_Trait_geneMHC.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE)
  autoshared <- get(paste0(trait.name,".pairwise"))
  autoshared.filtered <- get(paste0(trait.name,".pairwise.filtered"))
  autoshared.higheffect <- get(paste0(trait.name,".pairwise.higheffect")) # will be used in cytoscape for figure generation
  write.table(autoshared, file=paste(tempDir, "cyto_", make_safe(unique(Traitname$V3)[j]),'_shared_edgenumberMHC.txt',sep=''),row.names = FALSE,col.names = FALSE,quote = FALSE) # 1% threshold
  write.table(autoshared.filtered, file=paste(tempDir,"cyto_", make_safe(unique(Traitname$V3)[j]),'_shared_edgenumberMHC.filtered.txt',sep=''),row.names = FALSE,col.names = FALSE,quote = FALSE) #0.5% threshold
  write.table(autoshared.higheffect, file=paste(tempDir,"cyto_", make_safe(unique(Traitname$V3)[j]),'_Sshared_edgenumberMHC.higheffect.txt.txt',sep=''),row.names = FALSE,col.names = FALSE,quote = FALSE) #0.5% threshold
}
# Generate Figure S9
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(reshape2))
suppressMessages(library(plotly))

colnames(up.down.matrix) <- c("Trait.Code", "Trait", "Trait.Category", "Upregulated", "Downregulated", "Ambiguous")

up.down.plot <- as.data.table(up.down.matrix)
up.down.plot$Trait <- factor(up.down.plot$Trait, levels = up.down.plot[["Trait"]]) # so that traits of same categories cluster together
up.down.plot <- melt(up.down.plot, 
                     id.vars = c("Trait", "Trait.Category"),
                     measure.vars = c("Downregulated", "Ambiguous", "Upregulated"),
                     variable.name = "Direction", 
                     value.name = "Number.of.Genes")
# up.down.plot <- up.down.plot[!up.down.plot$Number.of.Genes == 0, ] # remove 0 values
mpdf("Barplot_up-down-ambi_per_trait", width=8,height=10)
ggplot(data=up.down.plot, aes(x=Trait, y=Number.of.Genes)) +
  geom_bar(aes(fill=Direction), stat="identity", position = "dodge") + # without the stat the graph collapses
  theme_minimal() +
  theme(legend.position = "top", 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_y_continuous(expand = c(0, 0)) + # expand = c(0, 0) is added to save space
  # scale_y_continuous(trans="log10", expand = c(0, 0)) + # expand = c(0, 0) is added to save space
  coord_flip() +
  labs(x = "Trait", y = "Number of correlated genes") +
  scale_fill_manual(values = cbPalette.upnodown)
dev.off()

# Generate Figure 6 
## Only traits that have at least 50 associated genes will be plotted in the CategoryPlotMHC
## Output is node names and sizes (1% FDR), edge width (1% FDR), edge color (are they in the same trait category?)
plot.trait <- data.table(subset(traitGenecountMHC, traitGenecountMHC[, 1] >= 50), keep.rownames = TRUE) # traitGenecountMHC was generated in the heatmap section
plot.trait <- sort(plot.trait$rn)
combo.temp.2 <- combn(plot.trait, 2) # generate list of pairwise comparisons needed within trait category.
sharedAsso <- matrix(0, dim(combo.temp.2)[2], 6) # final cytoscape output for Fig. 6
colnames(sharedAsso)<-c('node1short',  'node2short',   'edgewidth',   'edgeColor',       'node1Size',       'node2Size')
for (i in 1:dim(combo.temp.2)[2]) {
  sharedAsso[i,1] <- combo.temp.2[1,i]
  sharedAsso[i,2] <- combo.temp.2[2,i]
  trait.1 <- get(paste0(combo.temp.2[1, i])) # loads trait 1 gene table of pairwise comparison 1% threshold  
  trait.2 <- get(paste0(combo.temp.2[2, i])) # loads trait 1 gene table of pairwise comparison 1% threshold
  sharedAsso[i,3] <- length(intersect(trait.1$gene, trait.2$gene)) # number of shared genes
  sharedAsso[i,4] <- ifelse (Traitname[Traitname$V1 == combo.temp.2[1,i], 3] == Traitname[Traitname$V1 == combo.temp.2[2,i], 3],
                             1, # green (1) if they belong to same category
                             2) # yellow (2) if they don't
  sharedAsso[i,5] <- length(unique(trait.1$gene)) # node 1 size: # of unique genes
  sharedAsso[i,6] <- length(unique(trait.2$gene)) # node 2 size: # of unique genes
  #assign(sharedAsso, sharedAsso, pos = .GlobalEnv)
}
annot.trait <- Traitname[,1:2] # we will replace short names with full names in the final file
colnames(annot.trait) <- c("node2short","node2") ; sharedAsso <- merge(annot.trait,sharedAsso, by="node2short")
colnames(annot.trait) <- c("node1short","node1") ; sharedAsso <- merge(annot.trait,sharedAsso, by="node1short")
sharedAsso<-sharedAsso[, -grep("short", colnames(sharedAsso))]
write.table(sharedAsso,paste(outDir,'CategoryPlotMHC.txt',sep=""),row.names = FALSE,quote = FALSE,col.names = TRUE,sep='\t')

# Generate Figure S11 / EduYear is used for the graph/ For EduYear: i=15 (and all other traits).
# Only genes with more than 15 correlated genes will be mapped in the end
# For the correlation analysis dataset with 1% filtering is used
# Each count of correlated or anticorrelated corresponds to a specific gene and represents the majority of direction across tissues.
all.traits <- Traitname$V1
for (i in 1:length(all.traits)) { # for each main trait...
  # message(paste0("Generating (anti-)correlation cytoscape table for ", all.traits[i])) # for diagnostics
  cor.anticor <- as.data.frame(matrix(0, length(all.traits)-1, 8))
  cor.anticor [, 1] <- all.traits[i]
  cor.anticor [, 8] <- 165 # node size for cytoscape
  #rownames(cor.anticor) <- cor.anticor[, 1]
  trait.1.genes <- get (paste0(all.traits[i])) # get gene table for trait 1
  other.traits <- all.traits[-i]
  cor.anticor [,2] <- other.traits # populate all pairwise comparisons
  colnames(cor.anticor) <- c("Trait.1", "Trait.2", "Common.genes", "Correlated", "Anticorrelated", "Ambiguous", "Genes.present.in.different.tissues", "Node.Size")
  for (k in 1:length(other.traits)) { # for each other trait...
    trait.2.genes <- get(paste0(other.traits[k]))
    common.genes <- intersect(trait.1.genes$gene, trait.2.genes$gene)
    if (length(common.genes) == 0) { # if there are no common genes conserve CPU cycles
    } else { # calculate correlated and anticorrelated genes
      cor.anticor [k, 3] <- length(common.genes) # number of common genes shared by trait
      for (m in 1:length(common.genes)) { # for each gene...
        common.tissues <- intersect(trait.1.genes[trait.1.genes$gene == common.genes[m], "tissue"], trait.2.genes[trait.2.genes$gene == common.genes[m], "tissue"])
        if (length(common.tissues) == 0) {
          cor.anticor[k, 7] <- cor.anticor[k, 7] + 1
        } else {
          correlated = 0; anticorrelated = 0 # to keep score per tissue
          for (n in 1:length(common.tissues)) { # for each common tissue...
            ifelse (trait.1.genes[trait.1.genes$gene == common.genes[m] & trait.1.genes$tissue == common.tissues[n], "zscore"] * trait.2.genes[trait.2.genes$gene == common.genes[m] & trait.2.genes$tissue == common.tissues[n], "zscore"] > 0,
                    correlated <- correlated + 1,
                    anticorrelated <- anticorrelated + 1)
          } # loop for each tissue closes
          ifelse (correlated > anticorrelated,
                  cor.anticor[k, 4] <- cor.anticor[k, 4] + 1,
                  ifelse (anticorrelated > correlated,
                          cor.anticor[k, 5] <- cor.anticor[k, 5] + 1,
                          cor.anticor[k, 6] <- cor.anticor[k, 6] + 1))
        } # if statement for each tissue closes
      } # loop for each gene closes
    } # if statement for each gene closes 
  } # loop for other traits closes
  #!CAUTION!# Start - Remove all correlated traits with fewer of 15 common genes
  cor.anticor <- cor.anticor[cor.anticor$Common.genes >= 15, ]
  #!CAUTION!# End
  if (nrow(cor.anticor) == 0) { # if needed to compensate for aggressive filtering
    message("Trait: ", paste0(all.traits[i], " did not have any (anti-)correlated genes with other traits. No table file was generated"))
  } else { # saving data loop start
    cor.anticor[, 1] <- Traitname[cor.anticor[, 1], 2]  
    cor.anticor[, 2] <- paste0(Traitname[cor.anticor[, 2], 2], " (", cor.anticor[, 4], "/", cor.anticor[, 5], ")") 
    write.table(cor.anticor,paste(outDir,"plot_Cor-Anti-Cor_cytoscape_", make_safe(all.traits[i]), '_PlotMHC1.txt',sep=""),row.names = FALSE,quote = FALSE,col.names = TRUE,sep='\t')
  } # saving data loop end
} # loop for main trait closes


###############################################
# plot tissue specificity of genes per trait
###############################################
## Figure S7

z=datSigMHC
z=aggregate(1:nrow(z), by=z[c("gene","trait")], FUN=length) # aggregates gene counts per trait. Counts represent number of tissues. 
z[z[,3]>7,3]=7 # sets max = 7 for plot generation
z=t(as.matrix(table(z[,2:3]))) # creates table of gene counts per trait
z=as.data.frame.matrix(z)
z=data.frame(scale(z, center=F, scale=colSums(z))) # convert number of genes to percent within trait
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
# Figure S8

suppressMessages(library(stringr))
suppressMessages(library(plyr))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
z <- datSigMHC
z <- aggregate(1:nrow(z), by=z[c("gene", "trait", "tissue")], FUN=length)
z.one.tissue <- z[!(duplicated(z[c("gene","trait")]) | duplicated(z[c("gene","trait")], fromLast = TRUE)), ] # keep only unique gene/trait combinations
trait.list.pie <- unique(z$trait)
z.one.tissue <- aggregate(1:nrow(z.one.tissue), by=z.one.tissue[c("trait","tissue")], FUN=length)
z.one.tissue <- z.one.tissue[order(z.one.tissue$trait, -z.one.tissue$x), ]
for (i in 1:length(trait.list.pie)) {
  e <- z.one.tissue[z.one.tissue$trait == trait.list.pie[i], ]
  e$count.all <- sum(e$x)
  e$x <- e$x/sum(e$x)
  e$tissue <- paste(unlist(e[e$x == max(e$x), 2]), collapse=' & ')
  e <- e[1, ]
  z.one.tissue <- z.one.tissue[!z.one.tissue$trait == trait.list.pie[i], ]
  z.one.tissue <- rbind.fill(z.one.tissue, e) # requires plyr package
} # calculate percents of top tissues contributing most unique genes per trait
z.one.tissue$trait.category <- Traitname[z.one.tissue$trait, 3]
z.one.tissue$trait <- Traitname[z.one.tissue$trait, 2]
z.one.tissue <- z.one.tissue[order(z.one.tissue$trait.category,z.one.tissue$trait, decreasing = FALSE), ]
z.one.tissue$graph.label <- paste0(z.one.tissue$tissue, " (", percent(z.one.tissue$x), ")")
z.one.tissue$trait <- paste0(z.one.tissue$trait, " (", z.one.tissue$count.all, ")")
z.one.tissue$trait <- factor(z.one.tissue$trait, levels = z.one.tissue$trait) # make trait an ordered factor so that ggplot does not reorder it
caption <- paste(Tissuename$V1, Tissuename$V2, sep = ": ", collapse = "; ")
caption <- str_wrap(caption, width = 140) # requires stringr
mpdf("Major_tissue_contributors_of_unique_genes_per_trait", width=10,height=8)
ggplot(data=z.one.tissue, aes(x=trait, y=x, fill=trait.category)) +
  geom_bar(stat="identity") + # without that the graph collapses
  theme_minimal() +
  theme(legend.position = "top", 
        legend.text=element_text(size=8), 
        legend.key.size = unit(0.3, "cm"),
        plot.caption = element_text(size=7, hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_y_continuous(limits = c(0, 1.24), labels = scales::percent, breaks=c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) + # expand = c(0, 0) is added to save space
  geom_text(aes(label=graph.label), hjust=0, vjust=0.5, color="black", size=2.4) +
  coord_flip() +
  labs(x = "Trait", y = "Percent of unique genes coming from highest contributing tissue(s)", fill='Trait Category', caption = caption) +
  scale_fill_manual(values = cbPalette)
dev.off()

###################################################################################################
#  To prepare files for making plots of Antagonistic Pleiotropy between One (SCZ) and Other Traits
###################################################################################################

# get upregulated and downregulated genes of Schizophrenia
datSCZ=datSigMHC[datSigMHC$trait=='PGC2_SCZ',]
# more strict cutoffs
datSCZ=datSCZ[datSCZ$pred_perf_qval<=0.01,]
datSCZ=datSCZ[datSCZ$pva.qval<=0.01,]
# highly correlated genes
SCZassogene=datSCZ[abs(datSCZ$zscore)>=mean(abs(datSCZ$zscore)),]
SCZassogene=SCZassogene[,1:2]
write.table(SCZassogene,file = paste(interFolder,'PGC2_SCZ_TopAssoGeneMHC.txt',sep=''),quote=FALSE,row.names=FALSE,col.names=TRUE)

# all other traits
sczOthertrait=TRAIT
sczOthertrait=sczOthertrait[-38]

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

###########################################################################
## Make graphs showing tissue specificity of findings by trait: downstream
###########################################################################

dir.create(paste0(outDir,"/geneByTissueRaster"))
numberOfTopGenesToPlot=NA
#needDetailTrait=c('EGG_BW','GIANT_HIP','GIANT_HIPadjBMI','HEIGHT','LIPIDS_LDL','LIPIDS_TC','MENARCHE','PGC2_SCZ','RA','UC','eGFRcrea')
#noneedDetailTrait=setdiff(TRAIT,needDetailTrait)
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
  z$PlotLabel[z$pvalue<0.05]="·"
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
  z$PlotLabel[z$pvalue<0.05]="·"
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


#################################################################################
#     Check correlation between PrediXcan and GENET association results         #
#################################################################################
# only in significant results

GENETSig=datSigMHC
PrediXcanSig$resultSet='PrediXcan'
GENETSig$resultSet='GENET'
datSig=rbind(PrediXcanSig,GENETSig)

# folder to store correlation plots
dir.create(paste0(figsDir,"/PGcorrelation_abs"))
dir.create(paste0(figsDir,"/GENETOnly"))
dir.create(paste0(figsDir,"/PrediXcanOnly"))
# source("https://urldefense.proofpoint.com/v2/url?u=http-3A__bioconductor.org_biocLite.R&d=DwIFbw&c=shNJtf5dKgNcPZ6Yh64b-A&r=Er1GsV62JoD1xPAgXExBQrCtgoKB3MOzI2fvEhlKEtM&m=-UKsSEUnDAnJ5w60oqdQEneElZodgdDPGc5-JBGRTWE&s=o-yhiLOOFoD69kP95KVBUQtJ7fKJOyXN72mAntG3HbE&e=")
require(maptools)

for(trait in Traitname$V1){
  for (tissue in Tissuename$V1){
    dattmp=datSig[datSig$trait==trait & datSig$tissue== tissue,]
    if(dim(dattmp)[1]!=0){
      dattmp$color='grey' # grey for existing in both
      dattmp$zscore=as.numeric(abs(dattmp$zscore))
      for (checkgene in unique(dattmp$gene)){
        if(dim(dattmp[dattmp$gene==checkgene,])[1]!=2){
          if(dattmp[dattmp$gene==checkgene,17]=='GENET') dattmp[dattmp$gene==checkgene,18]='red' # red for GENET only
          if(dattmp[dattmp$gene==checkgene,17]=='PrediXcan') dattmp[dattmp$gene==checkgene,18]='blue' # blue for PrediXcan only
        }
      }
    } # set colors for datasets
    datplot=dattmp[dattmp$color=='grey',]
    
    if(dim(datplot)[1]!=0){
      plotgene=unique(datplot$gene) # unique genes
      predixcanV=data.frame(matrix(0,length(plotgene),2)) # separate matrix for ENet results
      genetV=data.frame(matrix(0,length(plotgene),2)) # separate matrix for WENet results
      
      
      for(i in 1:length(plotgene)){ # for each gene populate the matrices
        a=datplot[datplot$gene==plotgene[i],]
        predixcanV[i,1]=as.numeric(a[a$resultSet=='PrediXcan',3]) # zscore
        predixcanV[i,2]=a[a$resultSet=='PrediXcan',2] # gene name
        genetV[i,1]=as.numeric(a[a$resultSet=='GENET',3]) # zscore
        genetV[i,2]=a[a$resultSet=='GENET',2] # gene name
      }
    }
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
    #mpdf_nodingbats(paste0(Tissuename[tissue,1],"-",Traitname[trait,1]),subfolder = "PGcorrelation_abs")
    #pdf(paste(outDir,'PGcorrelation_abs/',Tissuename[tissue,1],Traitname[trait,1],".pdf",sep=''))
    # png files are easier for website
    png(paste(figsDir,'PGcorrelation_abs/',Tissuename[tissue,1],Traitname[trait,1],".png",sep=''),width = 2*600, height = 2*600, pointsize = 2*12)
    l=NA
    g=NA
    p=NA
    if(dim(datplot)[1]!=0){
      plot(predixcanV[,1],genetV[,1],xlim=c(-1,maxx+0.1), ylim=c(-1,maxy+0.1),
           pch = 16, cex = 0.8, col = rgb(125,127,130,210,maxColorValue=255), 
           main = paste("|z-scores| (", Traitname[trait,1], " from ", Tissuename[tissue,1],
                        ")\n correlation R= ", round(cor(predixcanV[,1],genetV[,1]),6)), 
           xlab = "PrediXcan", ylab = "GENET")
      l=NULL
    }
    if(length(datplotG)!=0) { if (is.null(l)) {
      points(datplotGx,datplotG,pch = 16, cex = 0.8,col = rgb(250,27,30,210,maxColorValue=255))
      #require(plotrix)
      #thigmophobe.labels(datplotGx, datplotG, labels = plotGtext, cex=0.7, offset=0.5)
      if (length(datplotG)<=3) 
        pointLabel(datplotGx, datplotG, labels = paste("  ", plotGtext, "  ", sep=""), cex=0.7)
      if(length(datplotG)>3) {
        plotGtext=plotGtext[order(datplotG,decreasing = TRUE)]
        datplotG=datplotG[order(datplotG,decreasing = TRUE)]
        pointLabel(datplotGx[1:3], datplotG[1:3], labels = paste("  ", plotGtext[1:3], "  ", sep=""), cex=0.7)
      }
      
    }
      if (!is.null(l)) {
        plot(datplotGx,datplotG,pch = 16, cex = 0.8,col = rgb(250,27,30,210,maxColorValue=255),
             main = paste0('|z-scores| (',Traitname[trait,1], ' from ', Tissuename[tissue,1],')'),
             xlab = "PrediXcan", ylab = "GENET")
        if (length(datplotG)<=3) 
          pointLabel(datplotGx, datplotG, labels = paste("  ", plotGtext, "  ", sep=""), cex=0.7)
        if(length(datplotG)>3)
        {
          plotGtext=plotGtext[order(datplotG,decreasing = TRUE)]
          datplotG=datplotG[order(datplotG,decreasing = TRUE)]
          pointLabel(datplotGx[1:3], datplotG[1:3], labels = paste("  ", plotGtext[1:3], "  ", sep=""), cex=0.7)
        }
        g=NULL
      }
    }
    if(length(datplotP)!=0) {
      if(is.null(g) | is.null(l)) {
        points(datplotP,datplotPy,pch = 16, cex = 0.8,col = rgb(20,27,130,210,maxColorValue=255))
        if (length(datplotP)<=3) 
          pointLabel(datplotP, datplotPy, labels = paste("  ", plotPtext, "    ", sep="   "),
                     cex=0.7,las=3,offset=0.75)
        if(length(datplotP)>3) {
          plotPtext=plotPtext[order(datplotPy,decreasing = TRUE)]
          datplotPy=datplotPy[order(datplotPy,decreasing = TRUE)]
          pointLabel(datplotP[1:3], datplotPy[1:3], labels = paste("  ", plotPtext[1:3], "    ", sep="   "), cex=0.7,las=3,offset=0.75)
        }}
      if (!is.null(g) & !is.null(l))  {
        plot(datplotP,datplotPy,pch = 16, cex = 0.8,col = rgb(20,27,130,210,maxColorValue=255),
             main = paste0('|z-scores| (',Traitname[trait,1], ' from ', Tissuename[tissue,1],')'), 
             xlab = "PrediXcan", ylab = "GENET")
        if (length(datplotP)<=3) 
          pointLabel(datplotP, datplotPy, labels = paste("  ", plotPtext, "    ", sep="   "), cex=0.7)
        if(length(datplotP)>3)
        {
          plotPtext=plotPtext[order(datplotPy,decreasing = TRUE)]
          datplotPy=datplotPy[order(datplotPy,decreasing = TRUE)]
          pointLabel(datplotP[1:3], datplotPy[1:3], 
                     labels = paste("  ", plotPtext[1:3], "    ", sep="   "), 
                     cex=0.7,las=3,offset=0.75)
        }
        p=NULL
      }
    }
    xm=lm(predixcanV[,1]~genetV[,1])
    if(!is.na(xm$coefficients[1]) & !is.na(xm$coefficients[2])) {
      if (is.null(l)) abline(lm(predixcanV[,1]~genetV[,1]),col='red',lty=5)
    }
    dev.off()
    # make warning pictures for those that miss the correlations this works only when script is run via linux
     if(dim(datplot)[1]==0 & length(datplotG)==0 & length(datplotP)==0){
     system(paste('cp ',paste(figsDir,'Warning/Warning.png',sep=''),' ',paste(figsDir,'PGcorrelation_abs/',Tissuename[tissue,1],Traitname[trait,1],".png",sep=''),sep=" "))}
    
  }
}

##############################################################################
# Make bar plots for gene regulations identified solely by GENET/PrediXcan   #
##############################################################################

# bar plots for GENET uniquely identified genes

require(reshape2)
require(ggplot2)

for(i in 1:dim(Traitname)[1]){
  for (j in 1:dim(Tissuename)[1]){
    dattmp=datSig[datSig$trait==Traitname[i,1] &datSig$tissue== Tissuename[j,1],]
    genetgene=dattmp[dattmp$resultSet=='GENET',1:2]
    predixcangene=dattmp[dattmp$resultSet=='PrediXcan',1:2]
    GenetGeneOnly=setdiff(unique(genetgene$gene),unique(predixcangene$gene))
    if(length(GenetGeneOnly)>0){
      dataGenet=dattmp[dattmp$resultSet=='GENET' & dattmp$gene %in% GenetGeneOnly,]
      dataGenet=dataGenet[,1:3]
      n_occur <- data.frame(table(dataGenet[,2]))
      if(nrow(n_occur[n_occur$Freq > 1, ])>=1) { # some genes have the same gene name
        require(dplyr)
        dataGenet <- dataGenet %>%
          group_by(gene_name) %>%
          dplyr::summarise(zscore= mean(zscore, na.rm=TRUE))
        dataGenet <- cbind(Regulation=NA, dataGenet) # adding a column in the beginning
      }
      if (dim(dataGenet)[1]>0){
        for(ii in 1:dim(dataGenet)[1]){
          if (dataGenet[ii,3]>0) dataGenet[ii,1]='Up'
          if (dataGenet[ii,3]<0) dataGenet[ii,1]='Down'}
        colnames(dataGenet)=c('Regulation','Gene','zScore')
        dataGenet=dataGenet[t(order(dataGenet$zScore)),]
        temGene=dataGenet$Gene
        dataGenet$Gene=ordered(dataGenet$Gene,level=temGene)
        dataGenet<-data.frame(dataGenet)
        png(paste(figsDir,'GENETOnly/',Tissuename[j,1],Traitname[i,1],".png",sep=''),width = 2*400, height = 2*400,  res=95, pointsize = 2*12)
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
          ggtitle(paste("GENET uniquely identified genes that are significantly associated with\n",Traitname[i,2],' (from ',Tissuename[j,2],')\nGENET performance q-value<=0.01, adjusted association p-value<=0.01',sep=''))
        print(p)
        dev.off()
      }
    }
    if(length(GenetGeneOnly)==0){
      system(paste('cp ',paste(figsDir,'Warning/Warning1.png',sep=''),' ',paste(figsDir,'GENETOnly/',Tissuename[j,1],Traitname[i,1],".png",sep=''),sep=" "))
    }
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
        png(paste(figsDir,'PrediXcanOnly/',Tissuename[j,1],Traitname[i,1],".png",sep=''),width = 2*400, height = 2*400,  res=95, pointsize = 2*12)
        #mpdf_nodingbats(paste0(Tissuename[j,1],"-",Traitname[i,1]), subfolder = "PrediXcanOnly")
        #pdf(paste(outDir,'PrediXcanOnly/',Tissuename[j,1],Traitname[i,1],".pdf",sep=''))
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
          ggtitle(paste("PrediXcan uniquely identified genes that are significantly associated with\n",Traitname[i,2],' (from ',Tissuename[j,2],')\nPrediXcan performance q-value<=0.01, adjusted association p-value<=0.01',sep=''))
        print(p)
        dev.off() }}
    if(length(PrediXcanGeneOnly)==0){
      system(paste('cp ',paste(figsDir,'Warning/Warning1.png',sep=''),' ',paste(figsDir,'PrediXcanOnly/',Tissuename[j,1],Traitname[i,1],".png",sep=''),sep=" "))
    }
    
  }
}


## make database file:
dataMHC=datMHC
dataMHC$trait_id=dataMHC$trait
dataMHC$tissue_id=dataMHC$tissue

for (trait in Traitname$V1)
{ dataMHC[dataMHC$trait==trait,14]=Traitname[trait,2]}

for (tissue in Tissuename$V1){ 
  dataMHC[dataMHC$tissue==tissue,13]=Tissuename[tissue,2]}

write.table(dataMHC,file='/Users/wenzhang/Desktop/GENET_meta.csv',quote=F,row.names=F,col.names = TRUE,sep='|')
save(dataMHC,file='/Users/wenzhang/Desktop/GENET_meta.RData')

genes=unique(dataMHC$gene)
gene=dataMHC[1:length(gene),1:2]

for(i in 1:length(genes)){
  dattmp=dataMHC[dataMHC$gene==genes[i],]
  dattmp=dattmp[1,1:2]
  gene[i,]=dattmp
}


##############################################################################
# CV qqplots comparing WENET and ENET
##############################################################################

# Plotting functions
## nn is the sample size, number of individuals used to compute correlation.
## needs correlation vector as input.
## nullcorvec generates a random sample from correlation distributions, under the null hypothesis of 0 correlation using Fisher's approximation.
# blue: rgb(20,27,230,200)
#abline(0,1,col='grey',lty=5)

qqR2<-function(corvec1,nn1,corvec2,nn2) # used for PNG files and pulled graph
{ set.seed(12345)
  mm1 <- length(corvec1)
  nullcorvec1 = tanh(rnorm(mm1)/sqrt(nn1-3)) ## null correlation vector
  mm2 <- length(corvec2)
  nullcorvec2 = tanh(rnorm(mm2)/sqrt(nn2-3)) ## null correlation vector
  corvec2[is.na(corvec2)]=0
  qqplot(nullcorvec2^2,corvec2^2, xlab=expression("Expected R"^"2"), ylab=expression(" Observed CV R"^"2 "),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.8,col = rgb(30,37,30,220,maxColorValue=255)); abline(0,1,col='grey',lty=1)
  points(sort(nullcorvec1^2),sort(corvec1^2),xlab=expression("Expected R"^"2"), ylab=expression(" Observed CV R"^"2 "),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.8,col = rgb(250,27,30,210,maxColorValue=255))
  legend(0,0.8, bty="n",legend = c("WENet", "ENet"), pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
}

qqR2_w_title<-function(corvec1,nn1,corvec2,nn2) {
  set.seed(12345)
  mm1 <- length(corvec1)
  nullcorvec1 = tanh(rnorm(mm1)/sqrt(nn1-3)) ## null correlation vector
  mm2 <- length(corvec2)
  nullcorvec2 = tanh(rnorm(mm2)/sqrt(nn2-3)) ## null correlation vector
  corvec2[is.na(corvec2)]=0
  qqplot(nullcorvec2^2, corvec2^2, 
         xlab=expression("Expected R"^"2"), 
         ylab=expression(" Observed CV R"^"2 "),
         ylim=c(0,0.925),
         xlim=c(0,0.0016),
         main=paste0(cv.table[i, "name"]), # Added title
         pch = 16,
         bty="n", # no box around the figure
         cex = 0.3,
         col = rgb(30,37,30,220,maxColorValue=255)); abline(0,1,col='grey',lty=1)
  points(sort(nullcorvec1^2), sort(corvec1^2), 
         xlab=expression("Expected R"^"2"), 
         ylab=expression(" Observed CV R"^"2 "),
         ylim=c(0,0.925),
         xlim=c(0,0.0016),
         pch = 16,
         bty="n",
         cex = 0.3,
         col = rgb(250,27,30,210,maxColorValue=255))
  # legend(0,0.8, bty="n",legend = c("WENet", "ENet"), pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
}

# Folder and file lists
# Create a data.frame with file names and tissue names
cv.folder.list <- c("CMC", "GTEX", "STARNET")
cv.table <- data.frame(folder=character(), name=character(), prior=character(), woprior=character()) # table will hold CV comparison table for GENET and PrediXcan
for (i in 1:length(cv.folder.list)) { # get needed file names for CV calculation
  cv.file.list <- list.files(path = paste0(interFolder, cv.folder.list[i]))
  #!CAUTION!# START selecting specific files: this is due to file names in the folder
  if (cv.folder.list[i] == "CMC") {
    cv.file.list <- cv.file.list[grep("CMC466_good", cv.file.list)]
    cv.file.list <- cv.file.list[-grep("maf5", cv.file.list)]
  }
  if (cv.folder.list[i] == "STARNET") {
    cv.file.list <- cv.file.list[grep("STARNET_", cv.file.list)]
    cv.file.list <- cv.file.list[-grep("Bezier", cv.file.list)]
    cv.file.list <- cv.file.list[-grep("(?i)Aor_", cv.file.list)] # use the Aorta ?i ignores case
    cv.file.list <- cv.file.list[-grep("MAM_noprior", cv.file.list)] # there is also a woprior
  }
  if (cv.folder.list[i] == "GTEX") { 
    cv.file.list <- cv.file.list[-grep("TW_", cv.file.list)]
    cv.file.list <- cv.file.list[-grep("[123].txt", cv.file.list)]
    cv.file.list <- cv.file.list[-grep("1_", cv.file.list)]
    cv.file.list <- cv.file.list[-grep("(?i)Aor", cv.file.list)] # too complicated to grep
    cv.file.list <- c("GTEX_Aor_noprior3.txt", "GTEX_Aor_prior3.txt", cv.file.list) # add manually
  }
  #!CAUTION!# END
  cv.file.list.prior <- cv.file.list[grep("_prior3{0,1}.txt", cv.file.list)]
  cv.file.list.woprior <- cv.file.list[grep("_[wn]oprior3{0,1}.txt", cv.file.list)]
  for (j in 1:length(cv.file.list.prior)) {
    if (sub("_[wn]{0,1}o{0,1}prior3{0,1}.txt", "", cv.file.list.prior[j]) !=sub("_[wn]{0,1}o{0,1}prior3{0,1}.txt", "", cv.file.list.woprior[j]) ) {
      message("There is something wrong with txt file table pairing priors and no priors. Troubleshooting is required.")
    }
  }
  # #!CAUTION!# START put a grep regular expression here to get the name and capitalize everything
  cv.table.temp <- cbind(folder=cv.folder.list[i], name="", prior = cv.file.list.prior, woprior = cv.file.list.woprior)
  if (cv.folder.list[i] == "CMC") { # generate name that will be used for output files
    cv.table.temp[, "name"] <- "CMC"
  } else {
    cv.table.temp[, "name"] <- toupper(sub("_[wn]{0,1}o{0,1}prior3{0,1}.txt", "", cv.table.temp[, "prior"]))
    cv.table.temp[, "name"] <- sub("BLOOD", "BLD", cv.table.temp[, "name"])
    cv.table.temp[, "name"] <- sub("LIVER", "LIV", cv.table.temp[, "name"])
    cv.table.temp[, "name"] <- sub("AORTA", "AOR", cv.table.temp[, "name"])
  } #!CAUTION!# END
  cv.table <- rbind(cv.table, cv.table.temp)
}

# Generate graph composite, per tissue, and pngs for website
mpdf("Figure S1 - cross validation WENET vs ENET", width = 7, height = 10)
par(mfrow = c(5, 3))  # 5 rows and 3 columns
par(mar=c(4.2,4.2,2,1)) # set proper margins bottom, left, top and right margins respectively of the plot region in number of lines of text.
for (i in 1:nrow(cv.table)) {
  rm(cv.wo, cv.prior)
  cv.wo <- read.table (paste0(interFolder, cv.table[i, "folder"], "/", cv.table[i, "woprior"]), header = TRUE)
  cv.prior <- read.table (paste0(interFolder, cv.table[i, "folder"], "/", cv.table[i, "prior"]), header = TRUE)
  qqR2_w_title(sqrt(cv.prior[, "R2"]),length(cv.prior[, "R2"]),sqrt(cv.wo[, "R2"]),length(cv.wo[, "R2"]))
  if (i == nrow(cv.table)) { # Also plot the legend 
    plot.new() # alternatively # plot(1, type="n", axes=FALSE, xlab="", ylab="")
    legend(x="center", bty="n", title = "Cross Validation", legend = c("WENet","ENet"),pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))  
  }
}
dev.off()

# Generate crossvalidation plot containing all measurements and individual png files
cv.wo.pooled <- setNames(data.frame(matrix(ncol=7, nrow=0)), c("gene", "genename", "R2", "adjustedR2", "n.snps", "pval", "rmse"))
cv.prior.pooled <- cv.wo.pooled 
R2Dir <- paste0(outDir,"CV") ; if (file.exists(R2Dir) == FALSE) dir.create(R2Dir)
for (i in 1:nrow(cv.table)) {
  rm(cv.wo, cv.prior)
  # message(paste0("Now doing: ", cv.table[i, 2]))
  cv.wo <- read.table (paste0(interFolder, cv.table[i, "folder"], "/", cv.table[i, "woprior"]), header = TRUE)
  cv.prior <- read.table (paste0(interFolder, cv.table[i, "folder"], "/", cv.table[i, "prior"]), header = TRUE)
  # generate the png files
  png(paste0(R2Dir,"/",cv.table[i,2],".png"), width = 2*600, height = 2*600, units = "px", pointsize = 2*12)
  qqR2(sqrt(cv.prior[,"R2"]),length(cv.prior[,"R2"]),sqrt(cv.wo[,"R2"]),length(cv.wo[,"R2"]))
  dev.off()
  # generate the pdf files
  pdf(paste0(R2Dir,"/",cv.table[i,2],".pdf"))
  qqR2(sqrt(cv.prior[,"R2"]),length(cv.prior[,"R2"]),sqrt(cv.wo[,"R2"]),length(cv.wo[,"R2"]))
  dev.off()
  cv.wo.pooled <- rbind.fill(cv.wo.pooled, cv.wo) # some tables have more columns
  cv.prior.pooled <- rbind.fill(cv.prior.pooled, cv.prior) # some tables have more columns
}


##############################################################################
# predictive correlation R2 comparing WENET and ENET
##############################################################################
# Brain Area abbreviations are based on "The Atlas of the Human Brain" by Jürgen K. Mai, George Paxinos, Thomas Voss
# list of abbreviations can be found http://www.thehumanbrain.info/database/nomenclature.php except for ACC

#!CAUTION!# START
cv.folder.list <- c("R2comp")
cv.table.R2 <- data.frame(folder= "R2comp",
                          name= c(
                            "CMC_HBCC",
                            "Starnet_GtexSKLM", "Starnet_GtexVAF", "Starnet_GtexSF", "Starnet_GtexAOR", "Starnet_GtexLIV","Starnet_GtexMAM", "Starnet_GtexBLD",
                            "GTEx_STARNETLIV", "GTEx_STARNETMAM", "GTEx_STARNETSF", "GTEx_STARNETVAF", "GTEx_STARNETBLD", "GTEx_STARNETAOR", "GTEx_STARNETSKLM",
                            "CMC_GTExAnterior_cingulate_cortex", "CMC_GTExCaudate_basal_ganglia", "CMC_GTExCerebellar_Hemisphere",
                            "CMC_GTExPutamen_basal_ganglia", "CMC_GTExNucleus_accumbens_basal_ganglia", "CMC_GTExHippocampus",
                            "CMC_GTExFC", "CMC_GTExCortex", "CMC_GTExCerebellum", "CMC_GTExHypothalamus"
                          ),
                          short.name = c(
                            "CMC -> HBCC",
                            "SKLM: Starnet -> Gtex", "VAF: Starnet -> Gtex", "SF: Starnet -> Gtex", "AOR: Starnet -> Gtex", "LIV: Starnet -> Gtex","MAM: Starnet -> Gtex", "BLD: Starnet -> Gtex",
                            "LIV: GTEx -> STARNET", "MAM: GTEx -> STARNET", "SF: GTEx -> STARNET", "VAF: GTEx -> STARNET", "BLD: GTEx -> STARNET", "AOR: GTEx -> STARNET", "SKLM: GTEx -> STARNET",
                            "CMC -> GTEx_ACC", "CMC -> GTEx_Caudate (BG)", "CMC -> GTEx_CH",
                            "CMC -> GTEx_Putamen (BG)", "CMC -> GTEx_Ac (BG)", "CMC -> GTEx_Hippocampus",
                            "CMC -> GTEx_FC", "CMC -> GTEx_Cortex", "CMC -> GTEx_Cerebellum", "CMC -> GTEx_Hypothalamus"
                          ),
                          prior= c("cmcprior_Hbcc_R2.txt",
                                   
                                   "Starprior_GtexSKLM_R2.txt", "Starprior_GtexVAF_R2.txt", "Starprior_GtexSF_R2.txt", 
                                   "Starprior_GtexAOR_R2.txt", "Starprior_GtexLiver_R2.txt", "Starprior_GtexMAM_R2.txt", 
                                   "Starprior_GtexBlood_R2.txt",
                                   
                                   "Gtexprior_StarLiver_R2.txt", "Gtexprior_StarMam_R2.txt", "Gtexprior_StarSF_R2.txt", 
                                   "Gtexprior_StarVAF_R2.txt", "Gtexprior_StarBlood_R2.txt", "Gtexprior_StarAor_R2.txt",
                                   "Gtexprior_StarSKLM_R2.txt",
                                   
                                   "cmcprior_GtexBrainAnterior_cingulate_cortex_R2.txt", "cmcprior_GtexBrainCaudate_basal_ganglia_R2.txt", 
                                   "cmcprior_GtexBrainCerebellar_Hemisphere_R2.txt", "cmcprior_GtexBrainPutamen_basal_ganglia_R2.txt", 
                                   "cmcprior_GtexBrainNucleus_accumbens_basal_ganglia_R2.txt", "cmcprior_GtexBrainHippocampus_R2.txt",
                                   "cmcprior_GtexBrainFC_R2.txt", "cmcprior_GtexBrainCortex_R2.txt", "cmcprior_GtexBrainCerebellum_R2.txt",
                                   "cmcprior_GtexBrainHypothalamus_R2.txt"
                          ),
                          woprior=c("cmcwo_Hbcc_R2.txt",
                                    
                                    "Starnoprior_GtexSKLM_R2.txt", "Starnoprior_GtexVAF_R2.txt", "Starnoprior_GtexSF_R2.txt", 
                                    "Starnoprior_GtexAOR_R2.txt", "Starnoprior_GtexLiver_R2.txt", "Starnoprior_GtexMAM_R2.txt", 
                                    "Starnoprior_GtexBlood_R2.txt",
                                    
                                    "Gtexnoprior_StarLiver_R2.txt", "Gtexnoprior_StarMAM_R2.txt", "Gtexnoprior_StarSF_R2.txt", 
                                    "Gtexnoprior_StarVAF_R2.txt", "Gtexnoprior_StarBlood_R2.txt", "Gtexnoprior_StarAor_R2.txt",
                                    "Gtexnoprior_StarSKLM_R2.txt",
                                    
                                    "cmcnoprior_GtexBrainAnterior_cingulate_cortex_R2.txt", "cmcnoprior_GtexBrainCaudate_basal_ganglia_R2.txt", 
                                    "cmcnoprior_GtexBrainCerebellar_Hemisphere_R2.txt", "cmcnoprior_GtexBrainPutamen_basal_ganglia_R2.txt", 
                                    "cmcnoprior_GtexBrainNucleus_accumbens_basal_ganglia_R2.txt", "cmcnoprior_GtexBrainHippocampus_R2.txt",
                                    "cmcnoprior_GtexBrainFC_R2.txt", "cmcnoprior_GtexBrainCortex_R2.txt", "cmcnoprior_GtexBrainCerebellum_R2.txt",
                                    "cmcnoprior_GtexBrainHypothalamus_R2.txt"
                          ) # table will hold CV comparison table for GENET and PrediXcan
)
R2.file.list <- list.files(path = paste0(interFolder, cv.folder.list))
needed.files <- c(as.character(cv.table.R2$prior), as.character(cv.table.R2$woprior))
# Check whether files are in the folder
for (i in 1:length(needed.files)) {
  if (needed.files[i] %in% R2.file.list) { } else {
    message("Caution! File: ", needed.files[i], " in folder ", paste0(interFolder, cv.folder.list), " does not exist!")
  }
}
#!CAUTION!# END

# Generate graph composite, per tissue, and pngs for website
qqR2cor_w_title<-function(corvec1,nn1,corvec2,nn2)
{
  set.seed(12345)
  corvec2=corvec2[!is.na(corvec2)]
  corvec1=corvec1[!is.na(corvec1)]
  mm1 <- length(corvec1)
  nullcorvec1 = tanh(rnorm(mm1)/sqrt(nn1-3)) ## null correlation vector
  mm2 <- length(corvec2)
  nullcorvec2 = tanh(rnorm(mm2)/sqrt(nn2-3)) ## null correlation vector
  qqplot(nullcorvec2^2,corvec2^2, main=figure.title, xlab=expression("Expected R"^"2"), ylab=expression("Observed correlation R"^"2"),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.3,col = rgb(30,37,30,220,maxColorValue=255)); abline(0,1,col='grey',lty=1)
  points(sort(nullcorvec1^2),sort(corvec1^2),xlab=expression("Expected R"^"2"), ylab=expression("Observed correlation R"^"2"),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.3,col = rgb(230,27,30,200,maxColorValue=255))
  #legend(0,0.8, bty="n",legend = c("WENet", "ENet"), pch = 16, col = c(rgb(230,27,30,200,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
}

R2.nonbrain <- cv.table.R2[2:15,]
rownames(R2.nonbrain[order(R2.nonbrain$short.name), ])
drawing.batch.1 <- rownames(R2.nonbrain[order(R2.nonbrain$short.name), ]) # non brain list / Figure S3
drawing.batch.2 <- c(1, 16:nrow(cv.table.R2))
#drawing.batch.2 <- c(1,18,20,23,21,19,22,16,17,25,24) # brain list / Figure S2

# Figure S3
mpdf("Figure S3 - predictive correlation WENET vs ENET non brain", width = 7, height = 10)
par(mfrow = c(5, 3))  # 5 rows and 3 columns
par(mar=c(4.2,4.2,2,1)) # set proper margins bottom, left, top and right margins respectively of the plot region in number of lines of text.
for (i in 1:length(drawing.batch.1)) {
  element.list <- ls()
  if ("cv.wo" %in% element.list) {rm(cv.wo, cv.prior)}
  figure.title <- cv.table.R2[drawing.batch.1[i], "short.name"]
  cv.wo <- read.table (paste0(interFolder, "R2comp", "/", cv.table.R2[drawing.batch.1[i], "woprior"]), header = TRUE)
  cv.prior <- read.table (paste0(interFolder, "R2comp", "/", cv.table.R2[drawing.batch.1[i],  "prior"]), header = TRUE)
  qqR2cor_w_title(sqrt(cv.prior[, "R2"]),length(cv.prior[, "R2"]),sqrt(cv.wo[, "R2"]),length(cv.wo[, "R2"]))
  if (i == length(drawing.batch.1)) { # Also plot the legend 
    plot.new()
    legend(x="center", bty="n", title = "Predictive Correlation", legend = c("WENet","ENet"),pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))  
  }
}
dev.off()

# Figure S2
mpdf("Figure S2 - predictive correlation WENET vs ENET brain", width = 7, height = 10)
par(mfrow = c(5, 3))  # 5 rows and 3 columns
par(mar=c(4.2,4.2,2,1)) # set proper margins bottom, left, top and right margins respectively of the plot region in number of lines of text.
for (i in 1:length(drawing.batch.2)) {
  element.list <- ls()
  if ("cv.wo" %in% element.list) {rm(cv.wo, cv.prior)}
  figure.title <- cv.table.R2[drawing.batch.2[i], "short.name"]
  cv.wo <- read.table (paste0(interFolder, "R2comp", "/", cv.table.R2[drawing.batch.2[i], "woprior"]), header = TRUE)
  cv.prior <- read.table (paste0(interFolder, "R2comp", "/", cv.table.R2[drawing.batch.2[i],  "prior"]), header = TRUE)
  qqR2cor_w_title(sqrt(cv.prior[, "R2"]),length(cv.prior[, "R2"]),sqrt(cv.wo[, "R2"]),length(cv.wo[, "R2"]))
  if (i == length(drawing.batch.2)) { # Also plot the legend 
    plot.new()
    legend(x="center", bty="n", title = "Predictive Correlation", legend = c("WENet","ENet"),pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))  
  }
}
dev.off()

# Generate crossvalidation plot containing all measurements and individual png files
# Figure 2B, Figure S2 panels, Figure S3 panels

# For pdfs and pngs
qqR2cor<-function(corvec1,nn1,corvec2,nn2)
{ set.seed(12345)
  corvec2=corvec2[!is.na(corvec2)]
  corvec1=corvec1[!is.na(corvec1)]
  mm1 <- length(corvec1)
  nullcorvec1 = tanh(rnorm(mm1)/sqrt(nn1-3)) ## null correlation vector
  mm2 <- length(corvec2)
  nullcorvec2 = tanh(rnorm(mm2)/sqrt(nn2-3)) ## null correlation vector
  qqplot(nullcorvec2^2,corvec2^2, xlab=expression("Expected R"^"2"), ylab=expression("Observed correlation R"^"2"),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.8,col = rgb(30,37,30,220,maxColorValue=255)); abline(0,1,col='grey',lty=1)
  points(sort(nullcorvec1^2),sort(corvec1^2),xlab=expression("Expected R"^"2"), ylab=expression("Observed correlation R"^"2"),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.8,col = rgb(230,27,30,200,maxColorValue=255))
  legend(0,0.8, bty="n",legend = c("WENet", "ENet"), pch = 16, col = c(rgb(230,27,30,200,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
}

# Figure 2B

R2.wo.pooled <- setNames(data.frame(matrix(ncol=7, nrow=0)), c("gene", "genename", "R2", "adjustedR2", "n.snps", "pval", "rmse"))
R2.prior.pooled <- R2.wo.pooled 
R2Dir <- paste0(outDir,"R2") ; if (file.exists(R2Dir) == FALSE) dir.create(R2Dir)
for (i in 1:nrow(cv.table.R2)) {
  element.list <- ls()
  if ("cv.wo" %in% element.list) {rm(cv.wo, cv.prior)}
  # message(paste0("Now doing: ", cv.table.R2[i, 2]))
  cv.wo <- read.table (paste0(interFolder, "R2comp", "/", cv.table.R2[i, "woprior"]), header = TRUE)
  cv.prior <- read.table (paste0(interFolder, "R2comp", "/", cv.table.R2[i, "prior"]), header = TRUE)
  # generate the png files
  png(paste0(R2Dir,"/",cv.table.R2[i,2],".png"), width = 2*600, height = 2*600, units = "px", pointsize = 2*12)
  qqR2cor(sqrt(cv.prior[,"R2"]),length(cv.prior[,"R2"]),sqrt(cv.wo[,"R2"]),length(cv.wo[,"R2"]))
  dev.off()
  # generate the pdf files
  pdf(paste0(R2Dir,"/",cv.table.R2[i,2],".pdf"), useDingbats=FALSE)
  qqR2cor(sqrt(cv.prior[,"R2"]),length(cv.prior[,"R2"]),sqrt(cv.wo[,"R2"]),length(cv.wo[,"R2"]))
  dev.off()
  R2.wo.pooled <- rbind.fill(R2.wo.pooled, cv.wo) # some tables have more columns
  R2.prior.pooled <- rbind.fill(R2.prior.pooled, cv.prior) # some tables have more columns
}


# FIGURE 2
# cv.prior.pooled, cv.wo.pooled, R2.prior.pooled, R2.wo.pooled files generated from Figures S1,S2,S3 are needed to generate figure 2

# Figure 2a plot function
qqR2_pooled<-function(corvec1,nn1,corvec2,nn2)
{
  set.seed(12345)
  mm1 <- length(corvec1)
  nullcorvec1 = tanh(rnorm(mm1)/sqrt(nn1-3)) ## null correlation vector
  mm2 <- length(corvec2)
  nullcorvec2 = tanh(rnorm(mm2)/sqrt(nn2-3)) ## null correlation vector
  corvec2[is.na(corvec2)]=0
  qqplot(nullcorvec2^2,corvec2^2, xlab=expression("Expected R"^"2"), ylab=expression(" Observed CV R"^"2 "),ylim=c(0,0.9),xlim=c(0,0.00018),pch = 16,bty="n",cex = 0.5,col = rgb(30,37,30,220,maxColorValue=255)); abline(0,1,col='grey',lty=1)
  points(sort(nullcorvec1^2),sort(corvec1^2),xlab=expression("Expected R"^"2"), ylab=expression(" Observed CV R"^"2 "),ylim=c(0,0.9),xlim=c(0,0.00018),pch = 16,bty="n",cex = 0.5,col = rgb(250,27,30,210,maxColorValue=255))
  #legend(0,0.8, bty="n",legend = c("WENet", "ENet"), pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
}

# Figure 2b plot function
qqR2_pooled_R2<-function(corvec1,nn1,corvec2,nn2)
{ set.seed(12345)
  corvec2=corvec2[!is.na(corvec2)]
  corvec1=corvec1[!is.na(corvec1)]
  mm1 <- length(corvec1)
  nullcorvec1 = tanh(rnorm(mm1)/sqrt(nn1-3)) ## null correlation vector
  mm2 <- length(corvec2)
  nullcorvec2 = tanh(rnorm(mm2)/sqrt(nn2-3)) ## null correlation vector
  qqplot(nullcorvec2^2,corvec2^2, xlab=expression("Expected R"^"2"), ylab=expression("Observed correlation R"^"2 "),ylim=c(0,0.9),xlim=c(0,0.00018),pch = 16,bty="n",cex = 0.5,col = rgb(30,37,30,220,maxColorValue=255)); abline(0,1,col='grey',lty=1)
  points(sort(nullcorvec1^2),sort(corvec1^2),xlab=expression("Expected R"^"2"), ylab=expression("Observed correlation R"^"2 "),ylim=c(0,0.9),xlim=c(0,0.00018),pch = 16,bty="n",cex = 0.5,col = rgb(250,27,30,210,maxColorValue=255))
  # legend(x = "bottomright" , bty="n",legend = c("WENet", "ENet"), pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
}

mpdf ("Figure 2 - cross validation WENET vs ENET pooled", height = 4 ) 
par(mfrow = c(1, 2))  # 1 rows and 2 columns
par(mar=c(4.2,4.2,2,1)) # set proper margins bottom, left, top and right margins respectively of the plot region in number of lines of text.
qqR2_pooled(sqrt(cv.prior.pooled[,"R2"]),length(cv.prior.pooled[,"R2"]),sqrt(cv.wo.pooled[,"R2"]),length(cv.wo.pooled[,"R2"]))
add_label(0.02, 0.07, "Crossvalidation")
qqR2_pooled_R2(sqrt(R2.prior.pooled[,"R2"]),length(R2.prior.pooled[,"R2"]),sqrt(R2.wo.pooled[,"R2"]),length(R2.wo.pooled[,"R2"]))
add_label (0.02, 0.07, "Predictive Correlation")
legend(x = "bottomright" , bty="n",legend = c("WENet", "ENet"), pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
dev.off()

##############################################################################
# Clinvar, OMIM and softpanel enrichment enrichment.
##############################################################################
# Plots only for more than 5 clinvar traits.
# requires: dplyr and ggplot2

############ !!!!!!!! UPDATE THIS PART !!!!!!!!! ##########

# source("./clinvar.R") # run the script to generate the clinvar annotation file (gets loaded as clinvar.genes)
# More information for the clinvar file can be found in the script above.
load("./Resources/clinvar.omim.clinical.RData") # check this one where is the df?
# Data generated by script above are:
# clinvar.genes.list.all                list of all clinvar genes
# clinvar.genes.list.pathogenic         list of clinvar genes that most likely are pathogenic
# clinvar.sig                           above list with phenotype data
# clinvar.genes.list.omim               list of clinvar genes with OMIM IDs
# clinvar.genes.omim.clinicalSynopsis   above list with clinical synopsis from OMIM

# calculate lambda
gv_lambda_pvalue <- function(pvalues.atomic.vector) {
  set.seed(12345)
  chisq <- qchisq(1-pvalues.atomic.vector, 1)
  # chisq <- qchisq(pvalues.atomic.vector, 1, lower.tail=FALSE) # same as
  lambda <- median(chisq) / qchisq(0.5, 1)
  return (lambda)
}
### CHECK ###
gv_lambda_pvalue_v2 <-function(pvalues.atomic.vector) {
  set.seed(12345)
  chisq <- qchisq(pvalues.atomic.vector, 1, lower.tail=FALSE)
  lambda <- median(chisq) / qchisq(0.5, 1)
  return (lambda)
}

###
# QQplots for unique traits #
###
# Only if >=3 genes are identified
# Only if >=3 points to be plotted
# Update to include values that are not statistically significant



#NOTES
## Traits except for trait groups as below:
## Dyslipidemia Group: "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia"
## BMI group contains "BMI","BMI, waist adjusted","BMI, hip adjusted","BMI, waist/hip ratio adjusted","Waist circumference", "Hip circumference", "Waist/hip ratio" : "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther"
## Childhood BMI and childhood obesity: "thisisaCScategory_growthWeight; thisisaCScategory_growthOther
## Bone mineral density group: bone mineral density; osteopetrosis; osteoporosis

# Prepare annotation table
Traitname$clinVar <- c( # for clinvar.sig pattern matching
  "Schizophrenia",                           "Alzheimer",                               "Autism",               
  "Attention Deficit Hyperactivity disorder","Anxiety",                                 "Anxiety",                 
  "Smoking",                                 "Smoking",                                 "Smoking",                
  "Smoking",                                 "Extraversion",                            "Neuroticism",                            
  "Depressive symptoms",                     "Subjective well being",                   "Years of education",                     
  "Crohn",                                   "Ulcerative colitis",                      "Inflammatory bowel disease",             
  "Atopic dermatitis",                       "Rheumatoid arthritis",                    "Systemic lupus erythematosus",           
  "Coronary artery disease; coronary heart disease","Myocardial infarction",            "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia",               
  "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia", "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia", "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia",                          
  "Diabetes mellitus type 2; Insulin-resistant diabetes mellitus; type II diabetes","fasting glucose; hyperglycemia; hypoglycemia","Glycated hemoglobin",                    
  "insulinem",                             "Insulin resistance",                        "Beta-cell function, HOMA",               
  "Glucose, fasting 2",                      "insulinem",                               "overweight; obese; obesity; cachexia; build",                                    
  "overweight; obese; obesity; cachexia; build", "overweight; obese; obesity; cachexia; build", "overweight; obese; obesity; cachexia; build",          
  "overweight; obese; obesity; cachexia; build", "overweight; obese; obesity; cachexia; build", "overweight; obese; obesity; cachexia; build",                        
  "Glucose, 2hr tolerance test",             "dwarfism; gigantism",                     "Childhood obesity",                      
  "Childhood BMI",                           "Birth weight",                            "microcephaly; macrocephaly",                    
  "Birth length",                            "Menarche, age at",                        "Menopause, age at",                      
  "osteoporosis; osteopetrosis",           "osteoporosis; osteopetrosis",           "osteoporosis; osteopetrosis",            
  "Kidney, eGFR",                            "Kidney, urinary albumin/creatinine",      "microalbuminuria",               
  "chronic kidney disease") 
Traitname$OMIM <- c(
  "Schizophrenia",                           "Alzheimer",                               "Autism",               
  "Attention Deficit Hyperactivity disorder; ADHD","Anxiety",                           "Anxiety",                 
  "Smoking",                                 "Smoking",                                 "Smoking",                
  "Smoking",                                 "Extraversion",                            "HP:0000716",                            
  "HP:0000716",                              "Subjective well being",                   "Years of education",                     
  "Crohn",                                   "Ulcerative colitis",                      "Inflammatory bowel disease",             
  "Atopic dermatitis",                       "Rheumatoid arthritis",                    "Systemic lupus erythematosus",           
  "Coronary artery disease",                 "Myocardial infarction",                   "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia",               
  "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia", "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia", "lipidemia; lipoprotein; cholesterolemia; hypertriglyceridemia",                          
  "Noninsulin-dependent diabetes mellitus; resistant diabetes mellitus; Diabetes mellitus type 2","fasting glucose; hyperglycemia; hypoglycemia","Glycated hemoglobin",                    
  "Insulinoma; insulinemia",                 "Insulin resistance",                      "Beta cell",               
  "fasting glucose; glycem",                 "insulinem",                               "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther",                                    
  "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther", "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther", "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther",          
  "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther", "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther", "overweight; obese; obesity; cachexia; build; thisisaCScategory_growthWeight; thisisaCScategory_growthOther",                        
  "Glucose, 2hr tolerance test",             "thisisaCScategory_growthHeight",         "thisisaCScategory_growthWeight; thisisaCScategory_growthOther",                      
  "thisisaCScategory_growthWeight; thisisaCScategory_growthOther",            "birth weight; birthweight",               "microcephaly; macrocephaly; small head; large head; Head circumference ",                    
  "Birth length",                            "Menarche, age at",                        "Menopause, age at",                      
  "Bone mineral density; osteoporosis; osteopetrosis",           "Bone mineral density; osteoporosis; osteopetrosis",           "Bone mineral density; osteoporosis; osteopetrosis",            
  "glomerular filtration rate",                            "Kidney, urinary albumin/creatinine",      "microalbuminuria",               
  "chronic kidney disease") 
Traitname$Softpanel <- Traitname$V1
softpanel <- read.csv(paste0(resourcesDir, "/softpanel.csv"), fill = TRUE, comment.char = "#", skip=2)

require(ggplot2)
Traitname$lambda <- NA
Traitname$lambda.clinvar <- NA
Traitname$lambda.omimCS <- NA
Traitname$lambda.SoftPanel <- NA
all.trait.clinvar.genes <- character() # will hold clinvar genes present in our traits
all.trait.omimCS.genes <- character() # will hold all genes that are associated with the OMIM clinical synopsis of our traits
all.trait.SoftPanel.genes <- character() # will hold all genes that are associated with keywords for our traits

# clinvar.sig                           above list with phenotype data
# clinvar.genes.list.omim               list of clinvar genes with OMIM IDs
# clinvar.genes.omim.clinicalSynopsis   above list with clinical synopsis from OMIM

#max.y.value <- max(-log10(datMHC$pvalue)) + max(-log10(datMHC$pvalue))*0.05
for (i in 1:length(Traitname$V1)) {
  df <- subset(datMHC, datMHC$trait == Traitname$V1[i])
  max.y.value <- max(-log10(df$pvalue))
  df <- df[order(df$pvalue), ]
  N  <- length(df$pvalue)
  df$observed <- -log10(df$pvalue)
  df$expected <- -log10(ppoints(N))
  number.entries <- nrow(df); number.unique.genes <- length(unique(df$gene_name))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  lambda.1 <- gv_lambda_pvalue(df$pvalue)
  label.1 <- paste0("All (", format(lambda.1, digits=3), "; ", format(number.entries, big.mark = ","), "; ", format(number.unique.genes, big.mark = ","), ")")
  Traitname$lambda[i] <- lambda.1
  # Calculate clinVAR
  clinvar.trait <- Traitname[i, "clinVar"]
  if (!is.na(clinvar.trait)) {
    clinvar.trait <- unlist(strsplit(clinvar.trait, split="; "))
    for (j in 1:length(clinvar.trait))  { # get clinvar and OMIM genes for trait
      trait.clinvar.genes.temp <- unique(clinvar.sig[grep(clinvar.trait[j], clinvar.sig$PhenotypeList, ignore.case=TRUE),][, "GeneSymbol"])
      if (j==1) {
        trait.clinvar.genes <- trait.clinvar.genes.temp
      } else {
        trait.clinvar.genes <- c(trait.clinvar.genes, trait.clinvar.genes.temp)
      }
      trait.clinvar.genes <- unique(trait.clinvar.genes)
      if (length(all.trait.clinvar.genes)!=0) {
        all.trait.clinvar.genes <- c(all.trait.clinvar.genes, trait.clinvar.genes)
      } else {all.trait.clinvar.genes <- trait.clinvar.genes}
    }
  }
  if (length(trait.clinvar.genes)>=2) {
    df.clinvar <- df[df$gene_name %in% trait.clinvar.genes, ]
    number.entries.clinvar <- nrow(df.clinvar); number.unique.genes.clinvar <- length(unique(df.clinvar$gene_name))
    if (number.entries.clinvar>=2){
      df.clinvar$expected <- -log10(ppoints(length(df.clinvar$pvalue)))
      lambda.2 <- gv_lambda_pvalue(df.clinvar$pvalue)
      Traitname$lambda.clinvar[i] <- lambda.2
      label.2 <- paste0("clinVar (", format(lambda.2, digits=3), "; ", format(number.entries.clinvar, big.mark = ","), "; ", format(number.unique.genes.clinvar, big.mark = ","), ")") 
    } else {df.clinvar<-NA;number.entries.clinvar<-NA;lambda.2<-NA;label.2<-NA}
  } else {df.clinvar<-NA;number.entries.clinvar<-NA;lambda.2<-NA;label.2<-NA}
  # OMIM using clinicalSynopsis
  omim.trait <- Traitname[i, "OMIM"]
  collist <- colnames(clinvar.genes.omim.clinicalSynopsis)
  collist <- collist[5:length(collist)] # we will need this for the grep
  if (!is.na(omim.trait)) {
    omim.trait <- unlist(strsplit(omim.trait, split="; "))
    for (j in 1:length(omim.trait))  { # get clinvar and OMIM genes for trait
      if (any(grep("thisisaCScategory_", omim.trait[j]))) {
        omim.trait[j] <- sub("thisisaCScategory_","",omim.trait[j])
        trait.omim.genes.temp <- clinvar.genes.omim.clinicalSynopsis[!is.na(clinvar.genes.omim.clinicalSynopsis[, omim.trait[j]]), "GeneSymbol"]
      } else { # fast way to check the whole dataframe
        selection <- apply(clinvar.genes.omim.clinicalSynopsis[,collist], 1, function(row) length(grep(omim.trait[j], row))>0)
        trait.omim.genes.temp <- unique(clinvar.genes.omim.clinicalSynopsis[selection,][, "GeneSymbol"])
      }
      if (j==1) {
        trait.omim.genes <- trait.omim.genes.temp
      } else {
        trait.omim.genes <- c(trait.omim.genes, trait.omim.genes.temp)
      }
      trait.omim.genes <- unique(trait.omim.genes)
      if (length(all.trait.omimCS.genes)!=0) {
        all.trait.omimCS.genes <- c(all.trait.omimCS.genes, trait.omim.genes)
      } else {all.trait.omimCS.genes <- trait.omim.genes}
    }
  }
  if (length(trait.omim.genes)>=2) {
    df.omim <- df[df$gene_name %in% trait.omim.genes, ]
    number.entries.omim <- nrow(df.omim); number.unique.genes.omim <- length(unique(df.omim$gene_name))
    if (number.entries.omim>=2){
      df.omim$expected <- -log10(ppoints(length(df.omim$pvalue)))
      lambda.3 <- gv_lambda_pvalue(df.omim$pvalue)
      Traitname$lambda.omimCS[i] <- lambda.3
      label.3 <- paste0("OMIM CS (", format(lambda.3, digits=3), "; ", number.entries.omim, "; ", number.unique.genes.omim, ")")
    } else {df.omim<-NA;number.entries.omim<-NA;lambda.3<-NA;label.3<-NA}
  } else {df.omim<-NA;number.entries.omim<-NA;lambda.3<-NA;label.3<-NA}
  # SoftPanel gene list
  softpanel.trait <- Traitname[i, "Softpanel"]
  if (!is.na(softpanel.trait)) {
    trait.softpanel.genes <- as.character(softpanel[, softpanel.trait])
    trait.softpanel.genes <- trait.softpanel.genes[trait.softpanel.genes != ""]
    if (length(all.trait.SoftPanel.genes)!=0) {
      all.trait.SoftPanel.genes <- c(all.trait.SoftPanel.genes, trait.softpanel.genes)
    } else {all.trait.SoftPanel.genes <- trait.softpanel.genes}
  }
  if (length(trait.softpanel.genes)>=2) {
    df.softpanel <- df[df$gene_name %in% trait.softpanel.genes, ]
    number.entries.softpanel <- nrow(df.softpanel); number.unique.genes.softpanel <- length(unique(df.softpanel$gene_name))
    if (number.entries.softpanel>=2){
      df.softpanel$expected <- -log10(ppoints(length(df.softpanel$pvalue)))
      lambda.4 <- gv_lambda_pvalue(df.softpanel$pvalue)
      Traitname$lambda.SoftPanel[i] <- lambda.4
      label.4 <- paste0("SoftPanel (", format(lambda.4, digits=3), "; ", number.entries.softpanel, "; ", number.unique.genes.softpanel, ")")
    } else {df.softpanel<-NA;number.entries.softpanel<-NA;lambda.4<-NA;label.4<-NA}
  } else {df.softpanel<-NA;number.entries.softpanel<-NA;lambda.4<-NA;label.4<-NA}
  # construct label
  labels <- label.1
  color.values <- c(cbPalette[1])
  shape.values <- c(16)
  if (!is.na(number.entries.clinvar)) {
    labels <- c(labels, label.2) 
    color.values <- c(color.values, cbPalette[4])
    shape.values <- c(shape.values, 17)}
  if (!is.na(number.entries.omim)) {
    labels <- c(labels, label.3) 
    color.values <- c(color.values, cbPalette[6])
    shape.values <- c(shape.values, 18)}
  if (!is.na(number.entries.softpanel)) {
    labels <- c(labels, label.4) 
    color.values <- c(color.values, cbPalette[7])
    shape.values <- c(shape.values, 15)}
  plot.this <- ggplot() +   geom_point(data = df, aes(df$expected, df$observed, colour = "custom.1", shape="custom.1"), size = 3)
  if (!is.na(number.entries.clinvar)) {
    plot.this <- plot.this + # add the clinvar data
      geom_point(data = df.clinvar, aes(df.clinvar$expected, df.clinvar$observed, colour = "custom.2", shape="custom.2"), size = 3)
  }
  if (!is.na(number.entries.omim)) {
    plot.this <- plot.this + # add the omim data
      geom_point(data = df.omim, aes(df.omim$expected, df.omim$observed, colour = "custom.3", shape="custom.3"), size = 3)
  } 
  if (!is.na(number.entries.softpanel)) {plot.this <- plot.this + # add the omim data
    geom_point(data = df.softpanel, aes(df.softpanel$expected, df.softpanel$observed, colour = "custom.4", shape="custom.4"), size = 3)} 
  plot.this <- plot.this + 
    scale_colour_manual(name = expression(paste("Gene Set (",lambda,"; entries; genes)")), values = color.values, labels = labels) +
    #scale_colour_manual(name = expression(paste("Gene Set (",lambda,"; entries; genes)")), values =cbPalette[c(1, 4, 6, 7)], labels = labels) +
    scale_shape_manual(name = expression(paste("Gene Set (",lambda,"; entries; genes)")), values = shape.values, labels = labels) + 
    #scale_shape_manual(name = expression(paste("Gene Set (",lambda,"; entries; genes)")), values = c(16, 17, 18, 15), labels = labels) + 
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_minimal() +
    theme(legend.position="top")+
    coord_cartesian(ylim = c(0, max.y.value+(max.y.value*0.05))) +
    annotate("text", x = 0.05, y = max.y.value+(max.y.value*0.05), label = paste0(Traitname$V2[i]), vjust = "inward", hjust = "inward")
  mpdf_nodingbats(make_safe(paste0("qqplot.geneset.",Traitname$V2[i])))
  print(plot.this)
  dev.off()
}

all.trait.clinvar.genes <- unique(all.trait.clinvar.genes)
all.trait.omimCS.genes <- unique(all.trait.omimCS.genes)
all.trait.SoftPanel.genes <- unique(all.trait.SoftPanel.genes)

# ALL TRAITS but enrichment only for genes that have been identified for our traits as above
gv_ggqqplot_enrich_relevant <- function(df, ci = 0.95) {
  require(ggplot2)
  df <- df[order(df$pvalue), ]
  N  <- length(df$pvalue)
  df$observed <- -log10(df$pvalue)
  df$expected <- -log10(ppoints(N))
  number.entries <- nrow(df); number.unique.genes <- length(unique(df$gene_name))
  df.clinvar <- df[df$gene_name %in% all.trait.clinvar.genes, ]
  df.clinvar$expected <- -log10(ppoints(length(df.clinvar$pvalue)))
  number.entries.clinvar <- nrow(df.clinvar); number.unique.genes.clinvar <- length(unique(df.clinvar$gene_name))
  df.omim <- df[df$gene_name %in% all.trait.omimCS.genes, ]
  df.omim$expected <- -log10(ppoints(length(df.omim$pvalue)))
  number.entries.omim <- nrow(df.omim); number.unique.genes.omim <- length(unique(df.omim$gene_name))
  df.softpanel <- df[df$gene_name %in% all.trait.SoftPanel.genes, ]
  df.softpanel$expected <- -log10(ppoints(length(df.softpanel$pvalue)))
  number.entries.softpanel <- nrow(df.softpanel); number.unique.genes.softpanel <- length(unique(df.softpanel$gene_name))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  lambda.1 <- gv_lambda_pvalue(df$pvalue) 
  lambda.2 <- gv_lambda_pvalue(df.clinvar$pvalue) 
  lambda.3 <- gv_lambda_pvalue(df.omim$pvalue)
  lambda.4 <- gv_lambda_pvalue(df.softpanel$pvalue)
  labels <- c(paste0("All (", format(lambda.1, digits=3), "; ", format(number.entries, big.mark = ","), "; ", format(number.unique.genes, big.mark = ","), ")"),
              paste0("clinVar (", format(lambda.2, digits=3), "; ", format(number.entries.clinvar, big.mark = ","), "; ", format(number.unique.genes.clinvar, big.mark = ","), ")"),
              paste0("OMIM CS (", format(lambda.3, digits=3), "; ", format(number.entries.omim, big.mark = ","), "; ", format(number.unique.genes.omim, big.mark = ","), ")"),
              paste0("SoftPanel (", format(lambda.4, digits=3), "; ", format(number.entries.softpanel, big.mark = ","), "; ", format(number.unique.genes.softpanel, big.mark = ","), ")"))
  p.enrich <- ggplot() +
    geom_point(data = df, aes(df$expected, df$observed, colour = "custom.1", shape="custom.1"), size = 3) +
    geom_point(data = df.clinvar, aes(df.clinvar$expected, df.clinvar$observed, colour = "custom.2", shape="custom.2"), size = 3) +
    geom_point(data = df.omim, aes(df.omim$expected, df.omim$observed, colour = "custom.3", shape="custom.3"), size = 3) +
    geom_point(data = df.softpanel, aes(df.softpanel$expected, df.softpanel$observed, colour = "custom.4", shape="custom.4"), size = 3) +
    # Within aes() colour and shape are assigned to a string to enable custom label
    scale_colour_manual(name = expression(paste("Gene Set (",lambda,"; entries; genes)")), values =cbPalette[c(1, 4, 6, 7)], labels = labels) +
    scale_shape_manual(name = expression(paste("Gene Set (",lambda,"; entries; genes)")), values = c(16, 17, 18, 15), labels = labels) +
    # names and labels have to be the same and repeated, otherwise two legends are formed
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_minimal() +
    theme(legend.position="top") +
    guides(color=guide_legend(nrow=2, byrow=TRUE),
           shape=guide_legend(nrow=2, byrow=TRUE)) # using guides to wrap legend
  return(p.enrich)
}

require(ggplot2)
ggsave(paste0(outDir,"clinvar.OMIM.softpanel.alltraits.QQplot.png"), gv_ggqqplot_enrich_relevant(datMHC), height = 7, width = 10, dpi = 600)
ggsave(paste0(outDir,"clinvar.OMIM.softpanel.alltraits.QQplot.pdf"), gv_ggqqplot_enrich_relevant(datMHC), height = 7, width = 10)

# Create dataframes for ENet vs. WENet p values in different gene sets and calculating density for plotting
WENet <- datMHC; WENet <- WENet[, c(2,13,14,5)]
WENet$pvalue <- -log10(WENet$pvalue)
names(WENet) <- c("gene.name", "tissue", "trait", "WENet.pvalue")
ENet <- datMHCP; ENet <- ENet[, c(2,13,14,5)] 
ENet$pvalue <- -log10(ENet$pvalue)
names(ENet) <- c("gene.name", "tissue", "trait", "ENet.pvalue")
require(dplyr)
WENvsEN <- inner_join(WENet, ENet)
rm("clinicalSynopsis.trait.list", "clinvar.genes.list.all", "clinvar.genes.list.omim",
   "clinvar.genes.list.pathogenic", "clinvar.genes.omim.clinicalSynopsis", "clinvar.sig",                       
   "clinvar.trait", "collist", "color.values", "critical.object.list", "df", "df.clinvar",                        
   "df.omim",  "df.softpanel", "i", "j", "label.1", "label.2", "label.3", 
   "label.4", "labels", "lambda.1" , "lambda.2", "lambda.3", "lambda.4", "log10Pe", 
   "log10Po", "max.y.value", "N", "number.entries", "number.entries.clinvar", 
   "number.entries.omim", "number.entries.softpanel", "number.unique.genes", 
   "number.unique.genes.clinvar", "number.unique.genes.omim", "number.unique.genes.softpanel",
   "omim.trait", "plot.this", "selection", "shape.values", "softpanel", "softpanel.trait",
   "trait.clinvar.genes", "trait.clinvar.genes.temp", "trait.omim.genes", "trait.omim.genes.temp",            
   "trait.softpanel.genes", "WENet", "ENet"); gc()
WENet.max <- max(WENvsEN$WENet.pvalue)*1.15
ENet.max <- max(WENvsEN$ENet.pvalue)*1.15
ifelse(WENet.max > ENet.max, ENet.max <- WENet.max, WENet.max <- ENet.max)
all.noted.genes.enrichment <- unique(c(all.trait.clinvar.genes, all.trait.omimCS.genes, all.trait.SoftPanel.genes))
WENvsEN.notlisted <- WENvsEN[!WENvsEN$gene.name %in% all.noted.genes.enrichment, ]
WENvsEN.clinvar <- WENvsEN[WENvsEN$gene.name %in% all.trait.clinvar.genes, ]
WENvsEN.omimcs <- WENvsEN[WENvsEN$gene.name %in% all.trait.omimCS.genes, ]
WENvsEN.softpanel <- WENvsEN[WENvsEN$gene.name %in% all.trait.SoftPanel.genes, ]
rm("WENvsEN")
# calculate distance for labelling important genes
WENvsEN.notlisted$distance <- abs(WENvsEN.notlisted$WENet.pvalue-WENvsEN.notlisted$ENet.pvalue) 
WENvsEN.clinvar$distance <- abs(WENvsEN.clinvar$WENet.pvalue-WENvsEN.clinvar$ENet.pvalue)
WENvsEN.omimcs$distance <- abs(WENvsEN.omimcs$WENet.pvalue-WENvsEN.omimcs$ENet.pvalue)
WENvsEN.softpanel$distance <- abs(WENvsEN.softpanel$WENet.pvalue-WENvsEN.softpanel$ENet.pvalue)

log10WENet <- expression(paste("WENet -log"[10], plain(P)))
log10ENet <- expression(paste("ENet -log"[10], plain(P)))
require("ggrepel") # required for geom_text_repel; caution, data positions may be slightly off

p.clinvar <- ggplot(data=WENvsEN.clinvar,aes(ENet.pvalue, WENet.pvalue))+
  geom_point(colour=cbbPalette[4], shape=2, size=3, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  geom_text_repel(data=WENvsEN.clinvar, aes(label=ifelse(distance>=min(dplyr::top_n(WENvsEN.clinvar, 3, distance)[,"distance"]), 
                                                         paste0(as.character(gene.name)," (",as.character(trait)," - ",as.character(tissue),")"), 
                                                         '')), show.legend = FALSE) +
  xlab(log10ENet) + ylab(log10WENet) + theme_minimal() + xlim(0, ENet.max) + ylim(0, WENet.max) +
  annotate("text", x = 0.05, y = WENet.max, label = "Relevant clinVar genes", vjust = "inward", hjust = "inward")

p.omim <- ggplot(data=WENvsEN.omimcs,aes(ENet.pvalue, WENet.pvalue))+
  geom_point(colour=cbbPalette[6], shape=5, size=3, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  geom_text_repel(data=WENvsEN.omimcs, aes(label=ifelse(distance>=min(dplyr::top_n(WENvsEN.omimcs, 3, distance)[,"distance"]), 
                                                        paste0(as.character(gene.name)," (",as.character(trait)," - ",as.character(tissue),")"), 
                                                        '')), show.legend = FALSE) + 
  xlab(log10ENet) + ylab(log10WENet) + theme_minimal() +  xlim(0, ENet.max) + ylim(0, WENet.max) +
  annotate("text", x = 0.05, y = WENet.max,  label = "Relevant OMIM CS genes", vjust = "inward", hjust = "inward")

p.softpanel <- ggplot(data=WENvsEN.softpanel,aes(ENet.pvalue, WENet.pvalue))+
  geom_point(colour=cbbPalette[7], shape=0, size=3, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  geom_text_repel(data=WENvsEN.softpanel, aes(label=ifelse(distance>=min(dplyr::top_n(WENvsEN.softpanel, 3, distance)[,"distance"]), 
                                                           paste0(as.character(gene.name)," (",as.character(trait)," - ",as.character(tissue),")"), 
                                                           '')), show.legend = FALSE) +
  xlab(log10ENet) + ylab(log10WENet) + theme_minimal() +  xlim(0, ENet.max) + ylim(0, WENet.max) +
  annotate("text", x = 0.05, y = WENet.max,  label = "Relevant SoftPanel genes", vjust = "inward", hjust = "inward")

p.rest <- ggplot(data=WENvsEN.notlisted, aes(ENet.pvalue, WENet.pvalue))+
  geom_point(colour=cbPalette[1], shape=1, size=3, show.legend = FALSE) + # 
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + #
  geom_text_repel(data=WENvsEN.notlisted, aes(label=ifelse(distance>=min(dplyr::top_n(WENvsEN.notlisted, 3, distance)[,"distance"]),  
                                                           paste0(as.character(gene.name)," (",as.character(trait)," - ",as.character(tissue),")"), 
                                                           '')), show.legend = FALSE) +
  xlab(log10ENet) + ylab(log10WENet) + theme_minimal() +  xlim(0, ENet.max) + ylim(0, WENet.max) +
  annotate("text", x = 0.05, y = WENet.max,  label = "Rest of the genes", vjust = "inward", hjust = "inward")

require("gridExtra")
g<-arrangeGrob(p.clinvar, p.omim, p.softpanel, p.rest, ncol = 2)
ggsave(paste0(outDir,"plot_pvalue_enrichments_WENetvs.ENet.repel.png"), g, height = 7, width = 10, dpi = 600)
ggsave(paste0(outDir,"plot_pvalue_enrichments_WENetvs.ENet.repel.pdf"), g, height = 7, width = 10)
suppressMessages(rm("ENet.max", "g", "log10ENet", "log10WENet",
                    "p.clinvar", "p.omim", "p.rest", "p.softpanel", "p.softpanel.pvalues", "WENet.max",
                    "WENvsEN", "WENvsEN.clinvar", "WENvsEN.notlisted", "WENvsEN.omimcs", "WENvsEN.softpanel"))

qq_test<-function(corvec1,nn1,corvec2,nn2,pvalue1,pvalue2)
{
  genet.lambda=gv_lambda_pvalue(pvalue1)
  genet.lambda=round(genet.lambda,4)
  predixcan.lambda=gv_lambda_pvalue(pvalue2)
  predixcan.lambda=round(predixcan.lambda)
  set.seed(12345)
  ## nn is the sample size, number of individuals used to compute correlation.
  ## needs correlation vector as input.
  ## nullcorvec generates a random sample from correlation distributions, under the null hypothesis of 0 correlation using Fisher's approximation.
  mm1 <- length(corvec1)
  nullcorvec1 = tanh(rnorm(mm1)/sqrt(nn1-3)) ## null correlation vector
  
  mm2 <- length(corvec2)
  nullcorvec2 = tanh(rnorm(mm2)/sqrt(nn2-3)) ## null correlation vector
  corvec2[is.na(corvec2)]=0
  corvec1[is.na(corvec1)]=0
  qqplot(nullcorvec2^2,corvec2^2, xlab=expression("Expected R"^"2"), ylab=expression(" Observed CV R"^"2 "),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.8,col = rgb(30,37,30,220,maxColorValue=255)); abline(0,1,col='grey',lty=1)
  
  points(sort(nullcorvec1^2),sort(corvec1^2),xlab=expression("Expected R"^"2"), ylab=expression(" Observed CV R"^"2 "),ylim=c(0,0.925),xlim=c(0,0.0016),pch = 16,bty="n",cex = 0.8,col = rgb(250,27,30,210,maxColorValue=255))
  #plot(nullcorvec^2,corvec^2, xlab=expression("Expected R"^"2"), ylab=expression("Observed predictive R"^"2"),ylim=c(0,0.825),pch = 16,cex = 0.8,col = rgb(250,27,30,210,maxColorValue=255))
  #abline(0,1,col='grey',lty=5)
  # blue: rgb(20,27,230,200)
  legend(0,0.8, bty="n",legend = c(expression(paste(paste("GENET (",lambda,"="),genet.lambda,")")), expression(paste(paste("PrediXcan (",lambda,"="),predixcan.lambda,")"))), pch = 16, col = c(rgb(250,27,30,210,maxColorValue=255), rgb(30,37,30,220,maxColorValue=255)))
  
}

plot.new()
genet.lambda <- 1.1111111111
predixcan.lambda <- 1.33333333333
legend(0,0.8, bty="n",legend = substitute(paste("GENET (",lambda,"=",genet.lambda.label,")"),
                                          list (genet.lambda.label = format(genet.lambda, digits=3))), 
       pch = 16, col = rgb(250,27,30,210,maxColorValue=255))
legend(0,0.75, bty="n",legend = substitute(paste("PrediXcan (",lambda,"=", predixcan.lambda.label,")"),
                                           list (predixcan.lambda.label = format(predixcan.lambda, digits=3))), 
       pch = 16, col = rgb(30,37,30,220,maxColorValue=255))


require(GenABEL.data)
data(srdta)
pex <- summary(gtdata(srdta))[,"Pexact"]

gv_lambda_pvalue_v2(pex)

gv_lambda_pvalue(pex)
estlambda(pex, method="median", filter = FALSE, proportion = 1, plot=FALSE)
