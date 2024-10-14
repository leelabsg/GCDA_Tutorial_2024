# Practice Session #3: Polygenic Risk Score (November 12, 2024)

In this session, we are going to construct polygenic risk score using PRS-CS. \
References : [PRS-CS github](https://github.com/getian107/PRScs), [PRS-CS paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6467998/). \
The data we are going to use are already preprocessed or downloaded.

### 0. Log in to leelabguest
``` 
ssh leelabguest@147.47.200.131 -p 22555
```

### 1. Connect CPU
``` 
ssh leelabsg11
``` 

### 2. Activate conda environment
``` 
conda activate python_3
``` 

### 4. Make directory for practice session in your directory
``` 
mkdir /data/GCDA/usr/YOUR_DIRECTORY/practice_4 
``` 

### 5. Run PRScs 
``` 
python /home/n1/leelabguest/GCDA/3_PRS/PRScs/PRScs.py \
--ref_dir=/home/n1/leelabguest/GCDA/3_PRS/data/reference/ldblk_1kg_eas \
--bim_prefix=/home/n1/leelabguest/GCDA/3_PRS/data/plink/sample \
--sst_file=/home/n1/leelabguest/GCDA/3_PRS/data/summary_stat/sumstats.txt \
--n_gwas=177618 \
--out_dir=/data/GCDA/usr/YOUR_DIRECTORY/practice_4/prscs
``` 

### 6. Merge chr1 - chr22 beta files into one file 
``` 
for i in {1..22}; do cat "/data/GCDA/usr/YOUR_DIRECTORY/practice_4/prscs_pst_eff_a1_b0.5_phiauto_chr$i.txt" >> /data/GCDA/usr/YOUR_DIRECTORY/practice_4/prscs_chr1-22.txt; done
``` 

### 7. Calculate PRS using plink 
``` 
/home/n1/leelabguest/plink \
--bfile /home/n1/leelabguest/GCDA/3_PRS/data/plink/sample \
--score /data/GCDA/usr/YOUR_DIRECTORY/practice_4/prscs_chr1-22.txt 2 4 6 \
--out /data/GCDA/usr/YOUR_DIRECTORY/practice_4/score
``` 