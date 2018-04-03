# GENET pipeline
 

This repository provides GENET (From GENotype and Epigenome to Transcriptome) scripts and files.
 
 For the scripts of training models, it is better to run them on the high performance cluster. It will
 take about several hours to complete the process. 

The scripts are for the GENET pipeline impelmentation and analysis results. 
Running examples are enclosed. 

All files are generated from the GENET approach which incorporate epigenomic information.
 
++++++++++++++++++++

Quick start:

   Please check the file RUNexample.txt (https://github.com/zwerhu/GENET/blob/master/RUNexample.txt)


++++++++++++++++++++

References: 
 
1. For GENET online database: http://icahn.mssm.edu/genet 

2. For PrediXcan pipeline: https://github.com/hakyim/PrediXcan

3. For qtlBHM pipeline: https://github.com/rajanil/qtlBHM

4. For qtlBHM file preparation and SNP prior generation: https://bitbucket.org/clagiamba/qtlbhm

5. For expression and genetic data preparation to run eQTL association: https://bitbucket.org/clagiamba/eqtl/overview

The cross validation (CV) training process using prior information is included in the script:

GENET_CV_elasticNet_penalty.R

or check the files after unzip the example.zip in example folder (refer to RUNexample.txt)

which incorporates qtlBHM priors as penalty factors for elastic net machine learning approach. The mapping function to project
priors to penalty factors is embedded into the script.

One example of running the GENET pipeline: exampleRun.sh, which shows the training process for chromosome 22. Cases for other chromosomes are the same, which could be execuated in parallel.

For the full version of example running process, please see RUNexample.txt.

The scripts provide the framework of integrating epigenomic information into the PredictDB training and predictions. A data-driven equation for mapping that converts priors to penalty factors is incorporated and the predictive performance (CV R2) and prediction correlations are improved.

The rescalings are approximated through optimal Bezier approach of second order with three control points.

++++++++++++++++++++

STEPS:

1. Prepare genotype, gene expressins, gene/snp annotations with the same format of PrediXcan pipeline.
   
2. Calculate qtlBHM priors using prepare.R and corresponding eQTL as well as roadmap epigenome mapping data.

3. Run the GENET pipeline.

4. Compile intermediate output files and get the predictors. 
