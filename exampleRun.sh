#!/bin/bash

#BSUB -J acc_roussp01a
#BSUB -P acc_roussp01a
#BSUB -q premium
#BSUB -n 4
#BSUB -R span[ptile=4]
#BSUB -R rusage[mem=12000]
#BSUB -W 27:00
#BSUB -m manda
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

Rscript ./create_model.R CMC467 ../../data/intermediate/expression_phenotypes/CMC.expr.RDS ../../data/intermediate/genotypes/CMC.new.snps.chr22.txt ../../data/intermediate/annotations/gene_annotation/gencode.v27.unified.chr22.RDS ../../data/intermediate/annotations/snp_annotation/CMC.snps.anno.chr22.RDS /sc/orga/projects/roussp01a/Wen/DLPF_annoprior/DLPFCBRN_Anno_prior_rsid.chr22.RDS 10 0.5 ../../data/intermediate/model_by_chr/ 22 Ref 1e6 > chr22.out


