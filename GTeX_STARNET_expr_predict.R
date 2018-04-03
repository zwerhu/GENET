# compare observed and predicted expressions of Gtex Liver/Aorta using STARNET PredictDB with/without priors
# Predicted gene expression without prior (STARNET_Liver data)
# the file names maybe changed due to several runnings
ex1<-read.table('/sc/orga/projects/roussp01a/Wen/convert/predicted_expression_Liver_no.txt',header=TRUE,stringsAsFactors = FALSE,row.names='IID')
ex1=ex1[,-1]
obsex<-read.table('/sc/orga/projects/roussp01a/Wen/convert/Liver_Gtex_geneexpression.txt',header=TRUE,stringsAsFactors=FALSE,row.names='gene_id')

obsex=obsex[,-1]
obsex=obsex[,-1]
obsex=obsex[,-1]
obsex<-t(obsex)
ex1<-t(ex1)
write.table(ex1,file="predicted_expression_Liver_no1.txt",quote=F)
ex2<-read.table('predicted_expression_Liver_no1.txt',header=TRUE,stringsAsFactors = FALSE)
ex2<-t(ex2)
commonGenes1<-intersect(colnames(obsex),colnames(ex2))
resultsmat1 <- matrix(0,ncol=2,nrow=length(commonGenes1))
colnames(resultsmat1) <- c('R2','p')
commGene1=commonGenes1
gene_annot <- readRDS('/sc/orga/projects/roussp01a/Wen/PredictDB/PredictDBPipeline/data/intermediate/annotations/gene_annotation/gencode.v25.genes.parsed.RDS')
for(j in 1:length(commonGenes1)){ commGene1[j]= gene_annot[commonGenes1[j],3]}
rownames(resultsmat1) <- commGene1

ind<-intersect(rownames(obsex),rownames(ex2))
obx1=obsex[ind,]
ex3=ex2[ind,]

for(j in 1:length(commonGenes1)){
gene<-commonGenes1[j]
res<-cor.test(obx1[,gene], ex3[,gene])
info<-c(round(res$estimate^2,7),signif(res$p.value,8))
resultsmat1[j,]<-info
}
finalres<-resultsmat1[order(resultsmat1[,1],decreasing=T),]
finalres[1:20,]

# compare observed and predicted expressions of Gtex Aorta using STARNET PredictDB and vice versa
# Predicted gene expression without prior (STARNET_Liver data)
ex1<-read.table('/sc/orga/projects/roussp01a/Wen/predicted_expression_Aor_noprior.txt',header=TRUE,stringsAsFactors = FALSE,row.names='IID')
ex1=ex1[,-1]
obsex<-read.table('/sc/orga/projects/roussp01a/Wen/Aor_Gtex_expression.txt',header=TRUE,stringsAsFactors=FALSE,row.names='gene_id')

obsex=obsex[,-1]
obsex=obsex[,-1]
obsex=obsex[,-1]
obsex<-t(obsex)
ex1<-t(ex1)
write.table(ex1,file="convert/predicted_expression_Aor_no1.txt",quote=F)
ex2<-read.table('convert/predicted_expression_Aor_no1.txt',header=TRUE,stringsAsFactors = FALSE)
ex2<-t(ex2)
commonGenes1<-intersect(colnames(obsex),colnames(ex2))
resultsmat1 <- matrix(0,ncol=2,nrow=length(commonGenes1))
colnames(resultsmat1) <- c('R2','p')
commGene1=commonGenes1
gene_annot <- readRDS('/sc/orga/projects/roussp01a/Wen/PredictDB/PredictDBPipeline/data/intermediate/annotations/gene_annotation/gencode.v25.genes.parsed.RDS')
for(j in 1:length(commonGenes1)){ commGene1[j]= gene_annot[commonGenes1[j],3]}
rownames(resultsmat1) <- commGene1

ind<-intersect(rownames(obsex),rownames(ex2))
obx1=obsex[ind,]
ex3=ex2[ind,]

for(j in 1:length(commonGenes1)){
gene<-commonGenes1[j]
res<-cor.test(obx1[,gene], ex3[,gene])
info<-c(round(res$estimate^2,7),signif(res$p.value,8))
resultsmat1[j,]<-info
}
finalno<-resultsmat1[order(resultsmat1[,1],decreasing=T),]
finalno[1:20,]

write.table(finalno,file="/sc/orga/projects/roussp01a/Wen/nopriorLiv_R2.txt",quote=F)

# with prior:

ex1<-read.table('/sc/orga/projects/roussp01a/Wen/predicted_expression_Aor_noprior.txt',header=TRUE,stringsAsFactors = FALSE,row.names='IID')
ex1=ex1[,-1]
obsex<-read.table('/sc/orga/projects/roussp01a/Wen/Aor_Gtex_expression.txt',header=TRUE,stringsAsFactors=FALSE,row.names='gene_id')

obsex=obsex[,-1]
obsex=obsex[,-1]
obsex=obsex[,-1]
obsex<-t(obsex)
ex1<-t(ex1)
write.table(ex1,file="convert/predicted_expression_Aor_prior1.txt",quote=F)
ex2<-read.table('convert/predicted_expression_Aor_prior1.txt',header=TRUE,stringsAsFactors = FALSE)
ex2<-t(ex2)
commonGenes1<-intersect(colnames(obsex),colnames(ex2))
resultsmat1 <- matrix(0,ncol=2,nrow=length(commonGenes1))
colnames(resultsmat1) <- c('R2','p')
commGene1=commonGenes1
gene_annot <- readRDS('/sc/orga/projects/roussp01a/Wen/PredictDB/PredictDBPipeline/data/intermediate/annotations/gene_annotation/gencode.v25.genes.parsed.RDS')
for(j in 1:length(commonGenes1)){ commGene1[j]= gene_annot[commonGenes1[j],3]}
rownames(resultsmat1) <- commGene1

ind<-intersect(rownames(obsex),rownames(ex2))
obx1=obsex[ind,]
ex3=ex2[ind,]

for(j in 1:length(commonGenes1)){
gene<-commonGenes1[j]
res<-cor.test(obx1[,gene], ex3[,gene])
info<-c(round(res$estimate^2,7),signif(res$p.value,8))
resultsmat1[j,]<-info
}
finalpri<-resultsmat1[order(resultsmat1[,1],decreasing=T),]
finalpri[1:20,]
