
## Load R functions and libraries
source('./scripts/rix_qtl_mapping_functions.r')

## Read phenotype data
cleaned_data_dir = '/Users/mooneymi/Documents/SIG/Mapping/WNV/16-Jul-2016_WNV_Gale_weight'
pheno_dir = '/Users/mooneymi/Documents/SIG/Mapping/phenotypes'
mapping_dir = '/Users/mooneymi/Documents/SIG/Mapping'

## The phenotype will be D12 weight change percentage
pheno = read.delim(file.path(cleaned_data_dir, 'd12_weight_pheno_table_GaleOnly_90lines_Oas1b.txt'), header=T, as.is=T)
dim(pheno)

## View the phenotype distribution accross the mapping population
library(lattice)
dotplot(reorder(pheno[,'UW_Line'], pheno[,'D12_Percentage'], mean, na.rm=T) ~ 
        pheno[,'D12_Percentage'] | pheno[,'Virus'], 
        panel = function(x,y,...) {panel.dotplot(x,y,...); panel.abline(v=0, col.line="red")}, 
        pch=19, ylab='UW Line', xlab="D12 Weight Change Percentage")

## Sort pheno dataframe and set rownames
pheno = pheno[with(pheno, order(Mating, RIX_ID)),]
rownames(pheno) = pheno$ID

## Add sex column
pheno$sex = 'M'

## Create covariate dataframe (must include sex)
covar = data.frame(sex = as.numeric(pheno$sex == 'M'))
rownames(covar) = pheno$ID

## Get IDs and Matings for each sample
samples = pheno$ID
matings = unlist(lapply(strsplit(samples, '_'), function(x) {x[1]}))
matings = unique(matings)

## Read strain ID mapping file (with Oas1b status for WNV mapping analyses)
## Note: use the Mx1 file (linked at top) for Flu analyses
strain_map = read.delim(file.path('./data', 'Oas1b_status_recoded.txt'), header=T, as.is=T, sep='\t')
head(strain_map)

## Check that matings match the strain mapping file
setdiff(matings, strain_map$Mating)

## Update covar with Oas1b status (only for WNV analyses)
covar$Mating = pheno$Mating
covar$Oas1b = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b[strain_map$Mating == x] else NA})
covar$Oas1b_High = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_High[strain_map$Mating == x] else NA})
covar$Oas1b_Mod = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_Mod[strain_map$Mating == x] else NA})
covar$Oas1b_Low = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_Low[strain_map$Mating == x] else NA})
pheno$Oas1b = covar$Oas1b
pheno$Oas1b_High = covar$Oas1b_High
pheno$Oas1b_Mod = covar$Oas1b_Mod
pheno$Oas1b_Low = covar$Oas1b_Low

## Load universal model probabilities (loads a model.probs object containing all RIX lines)
data_dir = './xdata'
load(file.path(data_dir, 'rix_universal_model_prob_males_27-Jun-2016.rda'))

## Create model.probs array specific to the mapping population
model.probs = model.probs[pheno$Mating, , ]
dimnames(model.probs)[[1]] = pheno$ID

## Check model.probs object
dim(model.probs)
names(dimnames(model.probs))
dim(model.probs)[1] == dim(pheno)[1]
model.probs[1,,1:5]

## Fix very small genotype probabilities
model.probs[model.probs < 0.005] = 1e-20

## Check model.probs object
model.probs[1,,1:5]

## Create kinship probability matrix
K = kinship.probs(model.probs)

## Check kinship matrix
K[1:5, 1:5]

## First get marker positions
marker_pos = read.csv(file.path(data_dir, 'CC001-Uncb38V01.csv'), as.is=T)
marker_pos = marker_pos[,1:3]
marker_pos$position_cM = NA
head(marker_pos)

## Run QTL scan
qtl = scanone(pheno=pheno, pheno.col='D12_Percentage', probs=model.probs, K=K, 
              addcovar=covar[, c('sex'), drop=F], snps=marker_pos)

## Run permutations to calculate significance threshold
perms = scanone.perm(pheno=pheno, pheno.col='D12_Percentage', probs=model.probs,
                     addcovar=covar[, c('sex'), drop=F], snps=marker_pos, nperm = 1000)

## Plot QTL results
plot(qtl, sig.thr = 6.75, main = 'Day 12 Weight Change Percentage')
abline(h=6.31, col='red', lty=2)

## Identify the Bayes Credible Interval for the peak on chromosome 5
interval = bayesint_v2(qtl, chr = 5)
interval

## Create founder effects (coefficient) plot for peak interval
coefplot_v2(qtl, chr=5, start=interval[1,3], end=interval[3,3])
