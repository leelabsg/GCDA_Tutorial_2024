# Practice Session #3: Two-sample Mendelian Randomization (Nov 11, 2024)
In this session, we will learn how to perform two Sample MR with publicly available GWAS datasets. \
The overall process in this session is mainly based on the [TwoSampleMR in R vignette](https://mrcieu.github.io/TwoSampleMR/). \
This document was created on October 28, 2024 and the following contents were tested on LeeLab servers.

#### Two Sample MR

### 0. Installing twoSampleMR in R using Conda 
To avoid conflicts in dependecies, we will create new environment for two sample MR (twoMR) in this session.
First, lets install R and devtools in new conda environment (Already installed as twoMR under leelabguest)
```
conda create -n twoMR r-base r-essentials r-devtools r-remotes r-gmp -c conda-forge -c R 
conda activate twoMR ; R --no-save
```
TwoSampleMR pacakge is installed from github 
```
remotes::install_github("MRCIEU/TwoSampleMR")
```
Next, Load the necessary packages
```
library(TwoSampleMR)
library(ggplot2)
```

### 1. Two Sample Mandelian Randomization Overview
The workflow for performing MR is as follows:

#### a. Select instruments for the exposure (perform LD clumping if necessary)
#### b. Extract the instruments from the IEU GWAS database for the outcomes of interest
- Or instead, you can either
    1) Create your own GWAS summary statistics
    2) Find the correct GWAS summary statistics for your MR design 
#### c. Harmonise the effect sizes for the instruments on the exposures and the outcomes to be each for the same reference allele
#### d. Perform MR analysis, sensitivity analyses, create plots, compile reports

![Twosamplemr_Diagram](https://user-images.githubusercontent.com/22064612/235571598-61784c94-bb0e-44ce-887d-c73a0593dbe6.png)
### 2. IEU GWAS database
The [IEU](https://gwas.mrcieu.ac.uk/) GWAS database (IGD) contains complete GWAS summary statistics from a large number of studies. You can browse them through the website. \
<img width="746" alt="image" src="https://user-images.githubusercontent.com/22064612/235642892-b283511b-6e2e-410a-9f7b-3c0611b325f7.png">

From the package, available_outcomes() will list up outcomes obtainable from IEU. 
```
ao <- available_outcomes(options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/'))
head(ao)
```
### 3. Exposure Data Example
In this part of session, we will follow the example from the package using formatted summary statistics. \
We want to discover causal relationship between BMI and CHD. \
BMI comes from [GIANT](https://gwas.mrcieu.ac.uk/datasets/ieu-a-2/) Consortium \
CHD is from [CARDIoGRAMplusC4D](https://gwas.mrcieu.ac.uk/datasets/ieu-a-7/). 
```
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')
head(bmi_exp_dat)
```
The output from this function is a new data frame with standardized column names: 
- SNP
- exposure
- beta.exposure
- se.exposure
- effect_allele.exposure
- other_allele.exposure
- eaf.exposure
- mr_keep.exposure
- pval.exposure
- pval_origin.exposure
- id.exposure
- data_source.exposure
- units.exposure
- gene.exposure
- samplesize.exposure
### 4. Outcome Data Example
Following similar procedure as Exposure data acquirment, let's bring CHD summary statistics.
By taking snps argument, we can scope the common variant list with BMI  
```
chd_out_dat <- extract_outcome_data(
    snps = bmi_exp_dat$SNP,
    outcomes = 'ieu-a-7'
)
head(chd_out_dat)
```
### 5. Harmonizing Data
Harmonizing data step is required because there are mismatches between alleles in exposure and outcome. Recent GWASs typically present the effects of a SNP in reference to the allele on the forward strand but not necessarily.
![image](https://user-images.githubusercontent.com/22064612/235645818-0c25fa91-7ba8-4dcf-9066-406e2db0eb9e.png)

```
dat <- harmonise_data(
    exposure_dat = bmi_exp_dat, 
    outcome_dat = chd_out_dat
)
```
There are three options to harmonising the data.

1. Assume all alleles are presented on the forward strand 
2. Try to infer the forward strand alleles using allele frequency information 
3. Correct the strand for non-palindromic SNPs, but drop all palindromic SNPs \
By default, the harmonise_data function uses option 2, but this can be modified using the action argument, e.g. harmonise_data(exposure_dat, outcome_dat, action = 3).
### 6. Perform MR
Before we perform MR, lets find out which method is available in the package
```
mr_method_list()
```
Let's try MR with Egger, two sample maximum likelihood,simple median, and IVW. You can add list of methods or methods set by default.

```
res <- mr(dat, method_list = c("mr_egger_regression", "mr_two_sample_ml","mr_simple_median","mr_ivw"))
res
```
Some Methods can perform sensitivity test which are implemented in mr_heterogeneity()
```
mr_heterogeneity(dat)
```
Intercept term in MR egger regression indicate the presence of Horizontal Pleiotropy.
```
mr_pleiotropy_test(dat)
```
Single SNP MR can be also performed by following.
```
res_single <- mr_singlesnp(dat)
head(res_single)
```
To see if a single snp is driving the association, we can perform leave-one-out MR.
In this case, we can see that there is no single variant that driving the association between exposure and outcome.
```
res_loo <- mr_leaveoneout(dat)
head(res_loo)
```

Visualizing MR analysis with plot can be performed as below:
Let's save to pdf using ggsave() from ggplot2
```
res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
print(length(p1))
ggplot2::ggsave(p1[[1]], file = "/data/home/leelabguest/MR_tutorial/IEU_example_Scatter.png", width = 7, height = 7)
```
Download the plot to local computer.
```
scp leelabguest@147.47.200.131:/data/home/leelabguest/MR_tutorial/IEU_example_Scatter.png ./
scp -P 22555 leelabguest@147.47.200.131:/data/home/leelabguest/MR_tutorial/IEU_example_Scatter.png ./
```
![filename](https://user-images.githubusercontent.com/22064612/235652335-b57f0381-b820-47e0-a07f-da09564f9c0b.png)

### 7. Real Data Practice
Let's try two-sample MR between HDL exposure and Diabetes outcome 
#### [HDL](https://gwas.mrcieu.ac.uk/datasets/ieu-a-299/) Exposure data from IEU GWAS
```
HDL_exp_dat <- extract_instruments(outcomes = 'ieu-a-299')
head(HDL_exp_dat)
```

#### Diabetes outcome data downloaded from [Pheweb](https://koges.leelabsg.org/pheno/KoGES_DM)
The downloaded summary statistics on Diabetes Malitious from KoGES pheweb is located in ~/GCDA/MR_tutorial/data \
It is saved as tab-delimited format(.tsv), we can take a look. \
manual outcome data frame should include snps, beta, effect, se, alleles, eaf column. \
Since we have case af and control af, we can generate eaf.

```
library(data.table)
library(dplyr)
outcome_DM = fread('/data/home/leelabguest/MR_tutorial/data/phenocode-KoGES_DM.tsv')
head(outcome_DM)
outcome_DM = outcome_DM %>% mutate(eaf = (case_af*5083+control_af*67127)/72210,n = 72210,pheno = "DM") 
```
![image](https://github.com/chkang96/GCDA/blob/0a3eec28848d0fa1b7dd2c0ae8753f8ee9a28be0/Screenshot%202024-11-03%20at%208.27.28%20PM.png)

SNP columns
```
outcome_dat <- format_data(
    type="outcome",
    phenotype_col="pheno",
    dat = data.frame(outcome_DM),
    snps = HDL_exp_dat$SNP,
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    eaf_col = "eaf",
    samplesize_col = "n"
)
```
#### Harmonizing data step & perform MR
```
dat <- harmonise_data(
    exposure_dat = HDL_exp_dat, 
    outcome_dat = outcome_dat
)
res = mr(dat)
res
p1 <- mr_scatter_plot(res, dat)
ggplot2::ggsave(p1[[1]], file = "/data/home/leelabguest/MR_tutorial/results/HDL_T2D_Scatter.png", width = 7, height = 7)
```

Heterogeneity test indicates there is a presence of heterogeneity between exposure and outcome, \
which is reasonable since the population group is different (Global vs KOR)
```
mr_heterogeneity(dat)
```
There is a slight evidence of horizontal pleiotropy 
```
mr_pleiotropy_test(dat)
```
Single SNP MR can be also performed by following. We can visualize the results via forest and funnel plots
```
res_single <- mr_singlesnp(dat)
head(res_single)
p2 <- mr_forest_plot(res_single)
ggplot2::ggsave(p2[[1]], file = "/data/home/leelabguest/MR_tutorial/results/HDL_T2D_Forest.png", width = 7, height = 7)
p3 <- mr_funnel_plot(res_single)
ggplot2::ggsave(p3[[1]], file = "/data/home/leelabguest/MR_tutorial/results/HDL_T2D_Funnel.png", width = 7, height = 7)
```
To see if a single snp is driving the association, we can perform leave-one-out MR.
In this case, we can see that there is no single variant that driving the association between exposure and outcome.
```
res_loo <- mr_leaveoneout(dat)
head(res_loo) ; summary(res_loo$p)
p4 <- mr_leaveoneout_plot(res_loo)
ggplot2::ggsave(p4[[1]], file = "/data/home/leelabguest/MR_tutorial/results/HDL_T2D_LOO.png", width = 7, height = 7)
```

### 8. Real Data Practice #2
Let's try another example, this time when both the exposure and outcome summary statistics are from published sources. 
It is important to note that the source of your data fits the design of your MR experiment. The 2nd practice is from an MR design regarding use of Statins on the PLA2G7 protein expression level. 

Here, the exposure SNPs come from a MR study of statins and cancer. Since statins work directly on the HMGCR gene region in order to lower LDL levels, we use genome-wide significant SNPs (p < 5E-8) in the HMGCR gene region from the source summary statistic as "proxies" for statin use. 

Key point: When involving drug targets, it was important to identify SNPs from the target gene region associated with the exposure phenotype (LDL levels). 

Exposure source: https://elifesciences.org/articles/57191 \
Outcome source: UKB-PPP Olink Proteomics data

Same as before, we read in exposure and outcome files and then harmonize and then view outputs
```
file_exp<-'/data/home/leelabguest/MR_tutorial/data/Statin_GWAS_0806.txt' 
ldl_exp_dat_org <- read_exposure_data( 
    filename = file_exp,
    sep = "\t",
    snp_col = "pos",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "Freq.A1.1000G.EUR",
    pval_col = "P-value"
)
```
```
file_protein<-'/data/home/leelabguest/MR_tutorial/data/PLA2G7_outcome_SNP.txt'
outcome_dat <- read_outcome_data(
    filename = file_protein,
    sep = ",",
    snp_col = "POS_19",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P",
    phenotype_col = "FILE_Index"
)
```
```
dat <- harmonise_data(
    exposure_dat = ldl_exp_dat_org, 
    outcome_dat = outcome_dat
)

res <- mr(dat, method_list = c("mr_weighted_median", "mr_egger_regression", "mr_ivw"))
```
Optional P-value adjustments:
```
res$p.BH<-p.adjust(res$pval, method = "BH")
```

### 9. Tips

1) Use genome-wide significant SNPs (p < 5E-8) in your exposure/outcome files.
2) Consider LD clumping when you have too many SNPs. Most MR studies typically use from 10 to a few hundred SNPs. LD clumping and all the rest of the work we did today is available here: https://mrcieu.github.io/TwoSampleMR/index.html
3) Conduct many sensitivity analyses to give you a better idea about the MR results.


