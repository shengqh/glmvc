GLMVC: Somatic Mutation Calling Using Both DNA and RNAseq Data

[Build from source](https://github.com/shengqh/glmvc/wiki/Build-from-source)

[Download binary file](https://github.com/shengqh/glmvc/releases)

[Requirement](https://github.com/shengqh/glmvc/wiki/Requirement)

[Configuration file](https://github.com/shengqh/glmvc/wiki/Configuration-file)

[Usage](https://github.com/shengqh/glmvc/wiki/Usage)

<a name="Changes"/>
#Changes
- 2015/12/07 Version 1.3.11
 1. Bugfix: add use_zero_minor_allele_strategy will cause filter processor exclude all candidates.
- 2015/12/07 Version 1.3.10
 1. Enhanced: add use_zero_minor_allele_strategy to pileup and filter. The candidate without minor allele detected at normal sample will be filtered by loose criteria.
- 2015/11/30 Version 1.3.9
 1. Enhanced: candidate position will not be used in mpileup for validation. Using "-l candidates.bed" option will cause different BAQ calibration result comparing to using all reads.
- 2015/11/24 Version 1.3.8
 1. Enhanced: read depth will be used as filter in call mpileup to accelerate the analysis but not in validate mpileup to report accurate result.
- 2015/11/18 Version 1.3.7
 1. Bug fix: validate: error when reading validation result without minor allele detected.
- 2015/11/17 Version 1.3.6
 1. Bug fix: remove a column distance_in_gene_range from result if gtf is provided for distance annotation
- 2015/11/16 Version 1.3.5
 1. Bug fix: table function will miss the first file if the file list doesn't have header
 2. Enhanced: table function will include refGeneAAChange in result.
- 2015/11/10 Version 1.3.4
 1. Enhanced: add option glm_use_raw_pvalue for filter
 2. New feature: summarize and build somatic mutation table
- 2015/11/05 Version 1.3.3
 1. Enhanced: add option exclude_bed to exclude specific ranges of genome in bed format
- 2015/11/01 Version 1.3.2
 1. Bug fixed: Parsing error of msm file name
- 2015/11/01 Version 1.3.1
 1. Bug fixed: Compatible with different chromosome names such like 4_JH584292_random in mm10 database
