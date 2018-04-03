# function that runs the model
# To be consistent with PrediXcan pipeline, parameter names keep the same
# snp_prior_RDS file is needed to store the priors, which is rescaled to penalty factors

argv <- commandArgs(trailingOnly = TRUE)
source("GENET_Tissue_Wide_CV_elasticNet_penalty.R")

study <- argv[1]
expression_RDS <- argv[2]
geno_file <- argv[3]
gene_annot_RDS <- argv[4]
snp_annot_RDS <- argv[5]
snp_prior_RDS <- argv[6]
n_k_folds <- as.numeric(argv[7])
alpha <- as.numeric(argv[8])
out_dir <- argv[9]
chrom <- argv[10]
snpset <- argv[11]
window <- as.numeric(argv[12])

TW_CV_model(expression_RDS, geno_file, gene_annot_RDS, snp_annot_RDS,snp_prior_RDS,5,n_k_folds, alpha, out_dir, study, chrom, snpset, window)
