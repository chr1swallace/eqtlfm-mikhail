#!/usr/bin/env bash

mkdir ~/rds/hpc-work/eqtlfm-mikhail
cd ~/rds/hpc-work/eqtlfm-mikhail

# login-mrc-bsu: eqtlfm-mikhail $ ftp ftp2.babraham.ac.uk
# Connected to ftp2.babraham.ac.uk (149.155.133.3).
# 220 (vsFTPd 2.0.5)
# Name (ftp2.babraham.ac.uk:cew54): ftpusr40
# 331 Please specify the password.
# Password:
# 230 Login successful.
# ftp> ls
# 227 Entering Passive Mode (149,155,133,3,118,40)
# 150 Here comes the directory listing.
# -rw-r--r--    1 542      543      1073566849 Feb 26 09:43 Merged_chromosomes_all_crms_359_LCLs.vcf
# -rw-r--r--    1 542      543      192248663 Feb 25 22:53 quantile_normalised_peer_residuals_LCLs.txt
# 226 Directory send OK.
# ftp> mget *

for i in `seq 1 22`; do
    wget https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GEUVADIS.chr${i}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz
done

# Subset the chr vcf files to individuals with genotype data, select variants with MAF > 0.01, exclude multiallelic snps
head -n1 quantile_normalised_peer_residuals_LCLs.txt | tr '\t' '\n' | egrep -v 'discard|gene_id' | sort > Samples1

cat Merged_chromosomes_all_crms_359_LCLs.vcf | head -n 300 | grep HG00096 | tr '\t' '\n' | egrep 'NA|HG' | sort> Samples2
comm -12 Samples1 Samples2 > Samples
wc Samples*

# VCF="Merged_chromosomes_all_crms_359_LCLs.vcf"
rm f*.sh
for i in `seq 22 -1 1`; do
    # VCF=GEUVADIS.chr${i}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz
    FROM="../eqtlfm/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    TO="jGeno-${i}.vcf.gz"
    # [[ -f $VCF ]] ||  echo "wget https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GEUVADIS.chr${i}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz" >> f${i}.sh
    echo "/scratch/wallace/local/bin/bcftools view -S Samples $FROM -Ou | /scratch/wallace/local/bin/bcftools view --min-af 0.01:minor --max-alleles 2 --min-alleles 2 -Oz -o $TO" >> f${i}.sh
    echo "/scratch/wallace/local/bin/bcftools index $TO" >> f${i}.sh
    if [ ! -f $TO ]; then
	echo "bash f${i}.sh > f${i}.log 2>&1" >> f.sh
    fi
done
qlines.rb -r f.sh


## convert to plink for twas
for i in `seq 1 22`; do
    # VCF=GEUVADIS.chr${i}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz
    FROM="jGeno-${i}.vcf.gz"
    TO="/rds/user/cew54/rds-cew54-wallace-share/Projects/twas/geuvadis/eur01-${i}"
    if [ -f $FROM -a ! -f ${TO}.bed ]; then
	echo "chromosome $i"
       plink --vcf $FROM \
	     --keep-allele-order \
	     --vcf-idspace-to _ \
	     --const-fid \
	     --allow-extra-chr 0 \
	     --split-x b37 no-fail \
	     --make-bed \
	     --out $TO
    fi
done

##  example code to convert to plink - don't run
bcftools view --regions 10:101316095-101316180 Geno.vcf.gz -Ob |\
plink --bcf /dev/stdin \
     --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr 0 \
    --split-x b37 no-fail \
    --make-bed \
    --out geno



