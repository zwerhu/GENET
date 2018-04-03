# Scripts to prepare files and run qtlBHM pipeline to get priors

###########
# eQTL dataset and condition
# condition means specific tissue for the dataset expressions
# dataset="CMC" 
# condition="DLPFC"
# condition="ACC"

# dataset="HBCC"
# condition="DLPFC"

# dataset="liver_Schadt"
# condition="Liver"
# dataset="GTEX"
# condition="LIV"
dataset="STARNET"
# condition="MAM"
# condition = "LIV"
# condition = "Blood" 
condition = "SKLM"  # for VAF and SF, annotation (subset_name) is "FAT_ADIP_NUC"
#condition = "SF"

###########
# Annotations that employed
##annot="ATACSeq.REMCannot"
#annot="REMC.auxHMM"
#annot = "REMC.coreHMM"
annot = "REMC.coreHMM.imputed"
#annot = "ATACSeq_QLFFIT"
#annot="TFADD10.ATACSeq"
##annot="TFADD10_CORRECTED.ATACSeq"
#annot = "TF.ATACSeq"
#annot="PEAKS.TF.ATACSeq"

# THIS IS ONLY FOR REMC
# REMC DLPFC (use the auxHMM here)
subset = TRUE # subset annotations in REMC
#subset_name= "BRN_DL_PRFRNTL_CRTX"
#subset_name="LIV_ADLT" # annotation name 
#subset_name="BLD"
#subset_name="FAT_ADIP_NUC" # for VAF and SF
subset_name="MUS_PSOAS" # for SKLM
#subset_name= "LIV_ADLT"
if (!grepl("REMC", annot)) subset=FALSE
###################
# Combined Models
# considering only the annotations with significant weights in the single model (i.e. with Weight +/- (1.96*Stderr) not crossing zero), both enriched and repressed, in the same combined model
#exclude_Quies = TRUE #always
# Model
singleAnnot=TRUE # First, submit each annotation separately one at a time # if FALSE it is a combined model
## Do this after running the single annotations and have weights for each 
# singleAnnot=FALSE: This creates a combined bed file for the annotations that resulted significant in the single annotation models (so run this after singleAnnot=TRUE generating all the outputs)

significant = TRUE # this should always be true
significant_enriched = TRUE  # Combine in a single model only annotations that have a positive weight and a 95% CI that does not cross zero

set = FALSE # if have a specific set of annotations you want to look at
set_names = c("ADD10_NEU_CORRECTED_TF_motif_121", "ADD10_NEU_CORRECTED_TF_motif_315", "ADD10_NEU_CORRECTED_TF_motif_727", "ADD10_NEU_CORRECTED_TF_motif_2675")
if (significant_enriched) significant = TRUE
###################
# input
base.in="/sc/orga/projects/epigenAD/coloc"
eqtlFolderIn=paste(base.in, '/data/', dataset, '/eQTLs/matrixEQTL/', condition, sep = '')
#if (dataset=="HBCC") eqtlFolderIn="/sc/orga/projects/epigenAD/coloc/data/HBCC/eQTLs/solly/formatted/"
pathANNOT ="/sc/orga/projects/roussp01a/Claudia_TMP/data/" # Path of annotations 
RDataFile=paste(pathANNOT, annot, ".Rdata", sep="")

# Setting output folders
base.dir = "/sc/orga/projects/roussp01a/Wen/QTLbhm/qtlBHM/"
# Split by eQTL datasets (because comparisons are going to be based within each eQTL data)
eqtlFolder = paste(base.dir, "eQTLs/", dataset, "/", sep="")
   if (!file.exists (eqtlFolder)) dir.create(eqtlFolder, recursive = T)
#oFolder=paste(base.out, dataset, "_", condition, "_", annot, "/", sep="")
annotMainFolder=paste(base.dir, "annotations/", annot, "/", sep="")
if (subset) annotFolder=paste(annotMainFolder, subset_name, "/", sep="") # for REMC, i.e. LIV_ADLT
if (singleAnnot) {  annotFolder = paste(annotFolder, "/single/", sep="") }
if (!singleAnnot) { annotFolder = paste(annotFolder, "/combined/", sep="") }
    if (!file.exists (annotFolder)) dir.create(annotFolder, recursive = T)

oMainFolder = paste(base.dir, dataset, "_", condition, "_", annot, "/", sep="")
if (singleAnnot) {  oFolder = paste(oMainFolder, "/single/", sep="") }

if (!singleAnnot) { oFolder = paste(oMainFolder, "/combined/", sep="") }
if (!singleAnnot & !significant & !set) { oFolder = paste(oFolder, "/exclude_Quies/", sep="") }
if (!singleAnnot & significant & !significant_enriched & !set) { oFolder = paste(oFolder, "/significant/", sep="") }
if (!singleAnnot & significant_enriched & !set) { oFolder = paste(oFolder, "/significant_enriched/", sep="") }
if (!singleAnnot & set & !significant_enriched & !significant) { oFolder = paste(oFolder, "/set/", sep="") }

   if (!file.exists (oFolder)) dir.create(oFolder, recursive = T)

   if (!file.exists (paste(oFolder, "log/", sep=""))) dir.create(paste(oFolder, "log/", sep=""), recursive = T)
   if (!file.exists (paste(oFolder, "scripts/", sep=""))) dir.create(paste(oFolder, "scripts/", sep=""), recursive = T)
   if (!file.exists (paste(oFolder, "results/", sep=""))) dir.create(paste(oFolder, "results/", sep=""), recursive = T)
   #if (!file.exists (paste(eqtlFolderOut, "data/", sep=""))) dir.create(paste(eqtlFolderOut, "data/", sep=""), recursive = T)


# eQTL annotation is the same for all subdirectories
OutFileStats=paste(eqtlFolder, condition, "_statistics.txt", sep="")
# OutFileAnnot=paste(annotFolder, "annotation.bed", sep="")
# Must also include name of the eqtl dataset and condition because different eqtl/annotations will have different significant combined models
if (!singleAnnot & !significant & !significant_enriched ) (OutFileAnnot=paste(annotFolder, dataset, "_", condition, "_exclude_Quies.bed", sep=""))
if (!singleAnnot & significant & !significant_enriched ) (OutFileAnnot=paste(annotFolder, dataset, "_", condition, "_significant.bed", sep=""))
if (!singleAnnot & significant_enriched ) (OutFileAnnot=paste(annotFolder, dataset, "_", condition, "_significant_enriched.bed", sep=""))
if (!singleAnnot & set) (OutFileAnnot=paste(annotFolder, "set.bed", sep=""))

# if single have a list of files
if (singleAnnot)  OutFileAnnotList =list.files(annotFolder, full.names=T)

source("/sc/orga/projects/roussp01a/Wen/qtlBHM/qtlbhm/functions_pipeline.R")
##################
# the gzipped file containing the test statistics (effect sizes and standard errors)
# Create one cis eQTL file
prepare_eqtl = TRUE

# prepare eQTL file
if (prepare_eqtl) {
OutFileTEMP=paste(eqtlFolder, "temp", sep="")
allFiles=paste(eqtlFolderIn, "/chr", 1:22, "_all_pval.tab", sep="")
#eqtlFolderIn="/sc/orga/projects/epigenAD/coloc/data/HBCC/eQTLs/"
if (all(file.exists(allFiles))==TRUE) {
    #final = as.character(read.table(allFiles[1], nrows = 1, header = FALSE, sep =' ', stringsAsFactors = FALSE))
    # header should not be included 
    #system(paste("head -1 ", allFiles[1], " > ", OutFile, sep=""))
##    system(paste("tail -n +2 -q ", paste(allFiles, collapse=" "), " >> ", OutFileTEMP, sep=""))

    combineCHR(base.folder=eqtlFolderIn, prefix = "chr", suffix="_all_pval.tab", OutFile=OutFileTEMP)
##    message("Final input file in ", OutFileTEMP)
}

if (file.exists(OutFileStats)) file.remove(OutFileStats)
# find columns in summary file
# find EnsemblID column
locus.col = grep("ProbeID", names(read.table(OutFileTEMP, header = T, nrow = 1, stringsAsFactors=F))) 
# find chr, pos, effect, SE columns
colsAll = sniff(file = OutFileTEMP)
colsAll = colsAll[[1]][names(colsAll[[1]]) %in% c("CHR", "POS", "BETA", "SE")]
# Locus              Variant         EffectSize  StandardError
#if (!all(c("snp", "chr", "pos", "se", "pval", "effect")   %in% names(colsAll[[1]]))) message("Some columns are not in the data")
# awk -F $'\t' -v OFS='\t' '{gsub(/"/, "", $2);gsub(/"/, "", $7); print $2, "chr"$7"."$8, $3, $12}' /sc/orga/projects/epigenAD/qtlBHMCMC/DLPFC/test_stat_file
#cmd = "awk -F $'\t' -v OFS='\t' '{gsub(/\"/, \"\", $2);gsub(/\"/, \"\", $7); print $2, \"chr\"$7\".\"$8, $3, $12}'"
#cmd = paste("awk -v probe=", locus.col, " -v chr=", colsAll[["chr"]], " -v pos=", colsAll[["pos"]], " -v effect=", colsAll[["effect"]], " -v se=", colsAll[["se"]],  " 'NR>1{OFS=\"\t\"; print $probe, \"chr\"$chr\".\"$pos, $effect, $se}'", sep="")

cmd = paste("awk -v probe=", locus.col, " -v chr=", colsAll[["CHR"]], " -v pos=", colsAll[["POS"]], " -v effect=", colsAll[["BETA"]], " -v se=", colsAll[["SE"]],  " 'NR>1{OFS=\"\t\"; gsub(/\"/, \"\", $probe); gsub(/\"/, \"\", $chr); print $probe, \"chr\"$chr\".\"$pos, $effect, $se}'", sep="")

cmd = paste(cmd, " ", OutFileTEMP, " > ",  OutFileStats, sep="")

# Delete all lines for which the kth field is 'NA'
# cmd = paste("awk -v probe=", locus.col, " -v chr=", colsAll[["chr"]], " -v pos=", colsAll[["pos"]], " -v effect=", colsAll[["effect"]], " -v se=", colsAll[["se"]],  " 'NR>1 && $", colsAll[["se"]], " != \"NA\"{OFS=\"\t\"; gsub(/\"/, \"\", $probe); gsub(/\"/, \"\", $chr); print $probe, \"chr\"$chr\".\"$pos, $effect, $se}'", sep="")

print(cmd)
system(cmd)

system(paste("rm", OutFileTEMP, sep=" "))
#system(paste("mv", OutFileTEMP, OutFileStats, sep=" "))
#system(paste("gzip", OutFileStats, sep=" "))
#OutFileStats=paste(OutFileStats, ".gz", sep="")
#gz1 = gzfile(OutFileStats,"w")
#dput(OutFileStats, gz1)
#close(gz1)
#paste("gzip -r ", oFolder, "/data/", sep="")

message("eQTL file saved in ", OutFileStats)
}


##################
# the gzipped file(s) containing genomic and/or variant annotations
prepare_annot=FALSE
# if combined model, usually want to re-do the bed files
if (!singleAnnot) prepare_annot = TRUE

if (prepare_annot) {
  source("/sc/orga/projects/roussp01a/Wen/qtlBHM/qtlbhm/functions_pipeline.R")
  options(scipen=999) # Turn off scientific notation
  listANNOT=get(load(RDataFile))
  print(listANNOT)

#if (annot=="REMC.auxHMM") {
if (grepl("REMC", annot)) {
if (subset) {
   outlist = list()
   #nameList = names(listANNOT)[grep("BRN", names(listANNOT))]   
   nameList = names(listANNOT)[grep(subset_name, names(listANNOT))]   
   for (name in nameList) {
        listANNOTsubset = listANNOT[[name]]
        outlist[[length(outlist)+1]] <- listANNOTsubset
   }
  names(outlist) = nameList
  listANNOT=outlist
}

# bedList = GRListtoBED(grList=listANNOT)
# GRListtoBEDforIntervals function: make the GenomicRanges and BED files exactly same:
bedList = GRListtoBEDforIntervals(grList=listANNOT)
# need to also keep id for REMC
# keep only columns seqnames    starts      ends id
bedList = lapply(bedList, function(x) x[(names(x) %in% c("seqnames", "starts", "ends", "id"))])
# add AnnotationLabel
bedList <- mapply(cbind, bedList, "AnnotationLabel"=names(bedList), SIMPLIFY=F)
# If there is more than one annotation then stop
if (length(names(bedList))>1) stop("There is more than one annotation selected with this name: choose only one")

bed= as.data.frame(bedList[[1]])
# for blood
#bed= as.data.frame(bedList[["BLD_CD14_PC"]])
#bed= as.data.frame(bedList[["BRN_DL_PRFRNTL_CRTX"]])
# bed <- do.call("rbind", bedList)
# for CMC, select the DLPFC
#bed = bed[bed$AnnotationLabel=="BRN_DL_PRFRNTL_CRTX",]
bed = bed[order(bed$id, bed$seqnames, bed$starts),]
bed$AnnotationLabel = paste(bed$AnnotationLabel, bed$id, sep="_")
bed=bed[,-4]
}

if (class(listANNOT)=="GRangesList") { # For ATACSeq.REMCannot
bedList = GRListtoBED(grList=listANNOT)
# keep only columns seqnames    starts      ends
bedList = lapply(bedList, function(x) x[(names(x) %in% c("seqnames", "starts", "ends", "id"))])
# add AnnotationLabel
bedList <- mapply(cbind, bedList, "AnnotationLabel"=names(bedList), SIMPLIFY=F)

# rbind it into one file
bed <- do.call("rbind", bedList)
# Add another division: the non-significant peaks, i.e. [ total - (peaks in neuorns+glia) ] 
bed$id= as.character(bed$id)
bed_NONE = bed[grep("TOTAL", bed$AnnotationLabel),] 
bed_NONE = bed_NONE[!(bed_NONE$id %in% bed$id[grep("NON-NEURONAL|NEURONAL", bed$AnnotationLabel)]),]
bed_NONE$AnnotationLabel = gsub("TOTAL", "NONE", bed_NONE$AnnotationLabel)
bed = rbind.data.frame(bed, bed_NONE)
bed=bed[,-4]
}

#if (annot!="REMC.auxHMM" & class(listANNOT)=="list") { # For TFs for ex.
if (!grepl("REMC", annot) & class(listANNOT)=="list") { # For TFs for ex.
   nameList = names(listANNOT)
   bed = data.frame()
   for (name in nameList) {
       bedList = GRListtoBED(grList=listANNOT[[name]])
       bedList = lapply(bedList, function(x) x[(names(x) %in% c("seqnames", "starts", "ends"))])
       bedList <- mapply(cbind, bedList, "AnnotationLabel"=paste(name, "_", names(bedList), sep=""), SIMPLIFY=F)
       bed_subset <- do.call("rbind", bedList)
       bed = rbind.data.frame(bed, bed_subset)
   }
}

if (!singleAnnot) {
  # some annotation files are: 
  # filenames = list.files(paste("/sc/orga/projects/epigenAD/qtlBHM/DLPFC_REMC.auxHMM/LIV/", "/results/", sep=""), full.names=T, pattern="_annotations.txt$")
  #filenames = list.files(paste(oMainFolder, "prior/", sep=""), full.names=T, pattern="_annotations.txt$") # if don't use $ it includes both: *_variant_annotations.txt.gz and *_annotations.txt
  filenames = list.files(oMainFolder, full.names=T, pattern="_annotations.txt$")
  if (grepl("REMC", annot)) filenames = filenames[grep(subset_name, filenames)]
  if (length(filenames)==0) stop("No results to use")
  if (length(filenames)!=length(list.files(paste(annotMainFolder, subset_name, "/single/", sep="")))) warning("The output folder does not contain results for all the single models")

  annotToKeep = do.call("rbind", lapply(filenames, read.table, header = TRUE, stringsAsFactors=FALSE))
  annotToKeep$UpperCI = annotToKeep$Weight - (1.96*annotToKeep$Stderr)
  annotToKeep$LowerCI = annotToKeep$Weight + (1.96*annotToKeep$Stderr)
  message("These annotations will not be included in the model:")

  #if (exclude_Quies) {
  matches = !grepl("Quies", annotToKeep$Annotation)
  print(annotToKeep[!matches,])
  annotToKeep = annotToKeep[matches,]
    
  if (significant) {
      matches = !apply(annotToKeep[,c("UpperCI", "LowerCI")],1,function(x){prod(sign(x))<=0})
      print(annotToKeep[!matches,])
      annotToKeep = annotToKeep[matches,]

    if (significant_enriched) {
        matches = annotToKeep$Weight>0
        print(annotToKeep[!matches,])
        annotToKeep = annotToKeep[matches,]
    }
  }

  if (set) {
      matches = annotToKeep$Annotation %in% set_names
      #if (length(annotToKeep[matches,])==0) stop("The files in set names do not match files in directory")
      if (length(annotToKeep$Annotation[matches])!=length(set_names)) stop("Some annotations in set names do not match")
      annotToKeep = annotToKeep[matches,]
  }

  bed=bed[bed$AnnotationLabel %in% annotToKeep$Annotation,]  
  bed$AnnotationLabel=as.character(bed$AnnotationLabel)
  print(names(table(bed$AnnotationLabel)))
}

if (singleAnnot) {
  # Split by AnnotationLabel and save for single annotations
  dataByAnnot = split(bed, f=as.factor(bed$AnnotationLabel))
  lapply(names(dataByAnnot), function(x){write.table(dataByAnnot[[x]], row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t", file = paste(annotFolder, "annotation.bed_", x, sep = ""))})

#gz1 = gzfile("df1.gz","w")
#dput(df1, gz1)
#close(gz1)
#paste("gzip -r ", oFolder, "/data/", sep="")
  OutFileAnnotList =list.files(annotFolder, full.names=T) 

} else {
write.table(bed, file=OutFileAnnot, row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")

#system(paste("gzip", OutFileAnnot, sep=" "))
#OutFileAnnot=paste(OutFileAnnot, ".gz", sep="")

message("Annotation file saved in ", OutFileAnnot)
}

}


##################
# general command:
# module load python py_packages
# python /hpc/users/zhangw17/infer_causal_variants.py /sc/orga/projects/epigenAD/qtlBHMCMC/DLPFC/statistics.txt.gz /sc/orga/projects/epigenAD/qtlBHMCMC/DLPFC/annotation.bed.gz

# submit to minerva hpc:

submitRun=TRUE
interactive = FALSE
account_name = "acc_roussp01a"
#account_name = "acc_epigenAD"
#account_name = "acc_psychgen" #slow
#account_name = "acc_psychimage"
# This failed for liver shadt data but seemed to work for other data:
#           #BSUB -R 'rusage[mem=8000]'
#           #BSUB -R span[ptile=4]
#           #BSUB -m bode
#           #BSUB -q premium
#          #BSUB -R 'span[hosts=1]'

if (singleAnnot) {
    for (f in 1:length(OutFileAnnotList)) {
    scriptname=paste(oFolder, "/scripts/Run_", basename(OutFileAnnotList[f]), ".sh", sep="")
    #if (dataset=="liver_Schadt"|dataset=="HBCC") {
    if (dataset=="liver_Schadt") {
    write(file=scriptname, paste("#!/bin/bash
          #BSUB -J qtlBHMrun.", basename(OutFileAnnotList[f]), "
          #BSUB -q premium
          #BSUB -P ", account_name, "
          #BSUB -R 'rusage[mem=12000]'
          #BSUB -n 4
          #BSUB -R 'span[ptile=4]'
          #BSUB -W 50:00
          #BSUB -L /bin/bash
          #BSUB -m manda
          #BSUB -oo ", oFolder, "/log/qtlbhm.", basename(OutFileAnnotList[f]), ".out
          #BSUB -eo ", oFolder, "/log/qtlbhm.", basename(OutFileAnnotList[f]), ".err", sep=""), append=F)
    } else {
    write(file=scriptname, paste("#!/bin/bash
          #BSUB -J qtlBHMrun.", basename(OutFileAnnotList[f]), "
          #BSUB -q premium
          #BSUB -P ", account_name, "
          #BSUB -R 'rusage[mem=12000]'
          #BSUB -n 4
          #BSUB -R 'span[ptile=4]'
          #BSUB -W 50:00
          #BSUB -L /bin/bash
          #BSUB -m manda
          #BSUB -oo ", oFolder, "/log/qtlbhm.", basename(OutFileAnnotList[f]), ".out
          #BSUB -eo ", oFolder, "/log/qtlbhm.", basename(OutFileAnnotList[f]), ".err", sep=""), append=F)
       } 
      if (!interactive) {
           write(file=scriptname, "module load python py_packages", append=T)
          } else {
          write(file=scriptname, "module load python py_packages", append=F)
          }

           write(file=scriptname, paste("cd ", oFolder, "results/", sep=""), append=T)
           write(file=scriptname, paste("python /hpc/users/zhangw17/infer_causal_variants.py ",  OutFileStats, " ", OutFileAnnotList[f], " --output_prefix /sc/orga/projects/roussp01a/Wen/prior/STARNET_SKLM/", basename(OutFileAnnotList[f]), sep=""), append=T)


       if (submitRun) {
          message("Submit script: ", scriptname)
          if (!interactive) {
           system(paste("bsub < ", scriptname, sep=""))
          } else {
           system(paste("nohup bash ", scriptname, " > ", oFolder, "/log/qtlbhm.", basename(OutFileAnnotList[f]), ".out 2> ", oFolder, "/log/qtlbhm.", basename(OutFileAnnotList[f]), ".err &", sep=""))
         }
       }
}
} 

if (!singleAnnot) {
    outfile = paste("combined_", annot, "_", gsub(".bed", "", basename(OutFileAnnot)), sep="")
    if (subset)  outfile = paste("combined_", annot, "_", subset_name, "_", gsub(".bed", "", basename(OutFileAnnot)), sep="")
    scriptname = paste(oFolder, "/scripts/", outfile, "_Run.sh", sep="")
    if (dataset=="liver_Schadt") {
    write(file=scriptname, paste("#!/bin/bash
          #BSUB -J qtlBHMrun", outfile, "
          #BSUB -q premium
          #BSUB -P ", account_name, "
          #BSUB -R 'rusage[mem=20000]'
          #BSUB -n 10
          #BSUB -R 'span[ptile=4]'
          #BSUB -W 78:00 
          #BSUB -L /bin/bash
          #BSUB -m manda
          #BSUB -oo ", oFolder, "/log/", outfile, "_Run.out
          #BSUB -eo ", oFolder, "/log/", outfile, "_Run.err", sep=""), append=F)
    } else {
    write(file=scriptname, paste("#!/bin/bash
          #BSUB -J qtlBHMrun", outfile, "
          #BSUB -q premium
          #BSUB -P ", account_name, "
          #BSUB -n 4
          #BSUB -R span[ptile=4]
          #BSUB -R 'rusage[mem=20000]'
          #BSUB -W 78:00 
          #BSUB -L /bin/bash
          #BSUB -m manda
          #BSUB -oo ", oFolder, "/log/", outfile, "_Run.out
          #BSUB -eo ", oFolder, "/log/", outfile, "_Run.err", sep=""), append=F)
    }
    write(file=scriptname, "module load python py_packages", append=T)
    write(file=scriptname, paste("cd ", oFolder, "results/", sep=""), append=T)
    #write(file=scriptname, paste("cd ", oFolder, "/", sep=""), append=T)
    write(file=scriptname, paste("python /hpc/users/zhangw17/infer_causal_variants.py ",  OutFileStats, " ", OutFileAnnot, " --output_prefix /sc/orga/projects/roussp01a/Wen/prior/STARNET_SKLM/", outfile, sep=""), append=T)
  if (submitRun) {
    system(paste("bsub < ", scriptname, sep=""))
  }
}
