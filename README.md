# FINDOR (beta release)
Functionally Informed Novel Discovery of Risk Loci
## Description
This tool is designed to improve GWAS power for polygenic traits. For details of methodology please see [Kichaev et al. (bioRxiv)](https://www.biorxiv.org/content/early/2017/11/20/222265).
## Required data
The core of the methodology relies on stratifying SNPs int bins of predicted chi square statistics. FINDOR uses the BaselineLD model from [Gazal et al (2016 Nat Genet)](https://www.nature.com/articles/ng.3954) for prediction.
Download LDscores for the 75 annotation model of BaselineLD model [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v1.1_ldscores.tgz)
## Running FINDOR
1. Run LD-score regression with the BaselineLD model on your GWAS data to get annotation effect size estimates. [See LD-score regression github](https://github.com/bulik/ldsc)
2. Run FINDOR on the entire GWAS to get re-weighted pvalues. Three inputs required:

	A. Full GWAS data set. Requires `N, SNP, Z` columns.
 
	B. BaselineLD model LD scores. Can be gzipped. 

	C. `.results` file from an application of LD score regression with the BaselineLD model on GWAS data. 

To access details on usage flags:
`python FINDOR.py --help`

An example execution comand would look like:

`python FINDOR.py --ref-ld-chr "$PATH_TO_LDSCORES"/baselineLD. \

		--gwas-data "$PATH_TO_GWAS_DATA"/gwas.data \

	        --regression-results "$PATH_TO_GWAS_DATA"/gwas.data.results\  

		--out "$PATH_TO_GWAS_DATA"/gwas.data..reweighted`
