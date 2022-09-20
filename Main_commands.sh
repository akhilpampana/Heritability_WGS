##############################################################################################################################################################
# Script used to run GCTA/plink/bcftools for the analyses in the GREML-LDMS WGS 2022 paper						    		     #	
# It contains the most important commands for performing GREML-LDMS from WGS										     #	
# The Online Methods describe with greater details all the QC steps performed										     #
# Some steps are added from the original github page													     #
# Additional R scripts are also available in this git													     #
##############################################################################################################################################################

#### Variables
ncpu=12 
sample_list="phenotypes/NTproBNP_14k.tsv" #List of samples with phenotypic information available

###############################################################################################################################################################
# 			Create BED files and perform first filtering on SVM variants and other QC metrics - submitted as slurm job			      #
###############################################################################################################################################################
###### Create BED files and perform first filtering on SVM variants and other QC metrics - submitted as slurm job

for i in  {1..22} ; do
    	bcftools view -m2 -M2  -Ou  -i'FILTER="PASS"' --threads "$ncpu" -S "$sample_list" --force-samples freeze.10a.chr"$i".pass_only.gtonly.minDP10.vcf.gz | 
	bcftools  annotate --threads "$ncpu" -Ob -I +'%CHROM:%POS:%REF:%ALT' > heritability/freeze10.14k.chr"$i".pass.bcf ;
	plink --bcf gwas/heritability/freeze10.14k.chr"$i".pass.bcf --maf 0.0001 --allow-extra-chr --keep-allele-order --geno 0.05 --hwe 0.000001 --threads ${ncpu} -make-bed --out freeze10.14k.chr"$i".0.0001 ;
done

###############################################################################################################################################################
# 								Allele Frequency Calculations								      #
###############################################################################################################################################################

for i in {1..22}; do
plink2 --bfile freeze10.14k.chr${i}.0.0001 --freq --out freeze10.14k.chr${i}
done

###############################################################################################################################################################
# 							Subset variants as per freeze10 pca's generation encore	    					      #
###############################################################################################################################################################
module load PLINK
for i in {1..22}; do
plink2 --bfile originial/freeze10.14k.chr${i}.0.0001 --extract variants to filter to --make-bed  --out plink_format/prunned_list_included_in_encore_pcs_generation/freeze10.14k.chr${i}.pruned
done

module load PLINK
for i in {1..22}; do
plink2 --bfile freeze10.14k.chr${i}.pruned --set-all-var-ids @:#:'$r':'$a'  --new-id-max-allele-len 1000 --max-alleles 2 --make-bed --out freeze10.14k.chr${i}.pruned
done

ls -l | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile freeze10.14k.chr1.pruned --merge-list merge --make-bed --out freeze10.14k.hqp_encore

module load KING/2.1.2-foss-2016a
#king -b freeze10.14k.hqp_encore.bed --unrelated --degree 2
king -b freeze10.14k.hqp_encore.bed --kinship    


###############################################################################################################################################################
# 					Extraction and seperation on unrelated vs related using kinship score < 0.025 (king.kin)                 	      #
###############################################################################################################################################################
 # IN house R script used

###############################################################################################################################################################
# 								Subset to unrelated individuals								      #
###############################################################################################################################################################

module load PLINK/1.90-foss-2016a
for i in {1..22}; do
plink --bfile freeze10.14k.chr${i}.0.0001_var --keep /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/plink_format/prunned_list_included_in_encore_pcs_generation/unrelated_kinship_0.025.txt --make-bed --out freeze10.14k.chr${i}.0.0001_var_unrel
plink --bfile freeze10.14k.chr${i}.0.0001_var_unrel --freq --out freeze10.14k.chr${i}.0.0001_var_unrel
done


###############################################################################################################################################################
# 		Subset Variants to 4 MAF bins [0.05],[0.01,0.05),[0.001,0.01),[0.0001,0.001) based on allele frequences of unrelated individuals	      #																			      #
###############################################################################################################################################################
#module load R => R
require(data.table)
category1 = data.frame() ## [0.0001,0.001)
category2 = data.frame() ## [0.001,0.01)
category3 = data.frame() ## [0.01,0.05)
category4 = data.frame() ## [0.05]
for(i in 1:22) {
	tmp = fread(paste0("freeze10.14k.chr",i,".0.0001_var_unrel.frq"))
	tmp = tmp[which(tmp$MAF >= 0.0001),]
	tmp1 = tmp[which(tmp$MAF >= 0.0001 & tmp$MAF < 0.001),]
	category1 = rbind(category1,tmp1)

	tmp1 = tmp[which(tmp$MAF >= 0.001 & tmp$MAF < 0.01),]
	category2 = rbind(category2,tmp1)
	
	tmp1 = tmp[which(tmp$MAF >= 0.01 & tmp$MAF < 0.05),]
	category3 = rbind(category3,tmp1)
	
	tmp1 = tmp[which(tmp$MAF >= 0.05),]
	category4 = rbind(category4,tmp1)
	
	}

save(category1,category2,category3,category4,file="allele_frequence_cutoff_4cat_unrel_09162022.rda")

write.table(category1,file="category1_0.0001_0.001.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category2,file="category2_0.001_0.01.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category3,file="category3_0.01_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category4,file="category4_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)

###############################################################################################################################################################
# 							   category1 - [0.0001,0.001)						       			      #
###############################################################################################################################################################

for i in  {1..22} ; do 
plink2 --bfile plink_format/original/freeze10.14k.chr${i}.0.0001_var_unrel --extract plink_format/original/category1_0.0001_0.001.csv --make-bed --out subset_for_h2_calc/category1_0.0001_0.001/unrelated/category1_chr${i}  ; 
done


cd subset_for_h2_calc/category1_0.0001_0.001/unrelated/
ls -l | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile category1_chr1 --merge-list merge --make-bed --out cat1_chr_all

### Remove variants based on geno 0.05 , mind 0.05, phwe < 1e-6 
for i in {1..22}; do
plink --bfile category1_chr${i} --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out category1_chr${i}_qc
done

### Prunning to generate high quality variants
for i in {1..22}; do
plink --bfile category1_chr${i}_qc --indep-pairwise 50 5 0.1 --out cat1_${i}
done

## subset to prunned variants
for i in {1..22}; do
plink --bfile category1_chr${i}_qc --extract cat1_${i}.prune.in --make-bed --out cat1_chr${i}_hqp
done

## LD score calculation to do 
for i in {1..22}; do
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat1_chr${i}_hqp --ld-score-region 200 --out  cat1_chr${i}_hqp --thread-num 100
done

### Categories based on quartiles as suggested in the paper  ## remove monomorphic variants before cuting to snps based on quartiles - code from gcta tutorial - to do 

require(data.table)
category1 = data.frame() 
category2 = data.frame() 
category3 = data.frame() 
category4 = data.frame() 
for(i in 1:22) {
lds_seg = read.table(paste0("cat2_chr",i,"_hqp.score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = as.data.frame(lds_seg$SNP[lb1])
lb2_snp = as.data.frame(lds_seg$SNP[lb2])
lb3_snp =  as.data.frame(lds_seg$SNP[lb3])
lb4_snp =  as.data.frame(lds_seg$SNP[lb4])

category1 = rbind(category1,lb1_snp)
category2 = rbind(category2,lb2_snp)
category3 = rbind(category3,lb3_snp)
category4 = rbind(category4,lb4_snp)
}

write.table(category1, "snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(category2, "snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(category3, "snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(category4, "snp_group4.txt", row.names=F, quote=F, col.names=F)


###############################################################################################################################################################
# 							   category2 - [0.001,0.01)						       			      #
###############################################################################################################################################################

for i in  {1..22} ; do 
plink2 --bfile plink_format/original/freeze10.14k.chr${i}.0.0001_var_unrel --extract plink_format/original/category2_0.001_0.01.csv --make-bed --out subset_for_h2_calc/category2_0.001_0.01/unrelated/category2_chr${i}  ; 
done


cd subset_for_h2_calc/category2_0.001_0.01/unrelated/
ls -l | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile category2_chr1 --merge-list merge --make-bed --out cat2_chr_all

### Remove variants based on geno 0.05 , mind 0.05, phwe < 1e-6 
plink --bfile cat2_chr_all --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out cat2_chr

### Prunning to generate high quality variants
plink --bfile cat2_chr --indep-pairwise 50 5 0.1 --out cat2

## subset to prunned variants
plink --bfile cat2_chr --extract cat2.prune.in --make-bed --out cat2_chr_hqp

## LD score calculation to do 
for i in {1..22}; do
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_chr${i}_hqp --ld-score-region 200 --out  cat2_chr${i}_hqp --thread-num 100
done


### Categories based on quartiles as suggested in the paper  ## remove monomorphic variants before cuting to snps based on quartiles - code from gcta tutorial

require(data.table)
category1 = data.frame() 
category2 = data.frame() 
category3 = data.frame() 
category4 = data.frame() 
for(i in 1:22) {
lds_seg = read.table(paste0("cat2_chr",i,"_hqp.score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

category1 = rbind(category1,lb1_snp)
category2 = rbind(category2,lb2_snp)
category3 = rbind(category3,lb3_snp)
category4 = rbind(category4,lb4_snp)
}

write.table(category1, "snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(category2, "snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(category3, "snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(category4, "snp_group4.txt", row.names=F, quote=F, col.names=F)

### SUBSET TO 4 QUARTILES TO GENERATE GRM 
module load PLINK/1.90-foss-2016a
for i in {1..22}; do
plink --bfile cat2_chr${i}_hqp --extract snp_group1.txt --make-bed --out quartiles/cat2_chr${i}_hqp_q1
plink --bfile cat2_chr${i}_hqp --extract snp_group2.txt --make-bed --out quartiles/cat2_chr${i}_hqp_q2
plink --bfile cat2_chr${i}_hqp --extract snp_group3.txt --make-bed --out quartiles/cat2_chr${i}_hqp_q3
plink --bfile cat2_chr${i}_hqp --extract snp_group4.txt --make-bed --out quartiles/cat2_chr${i}_hqp_q4
done

### MERGE TO ONE MEGA FILE
ls -l | grep q1 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat2_chr1_hqp_q1 --merge-list merge --make-bed --out cat2_hqp_q1

ls -l | grep q2 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat2_chr1_hqp_q2 --merge-list merge --make-bed --out cat2_hqp_q2

ls -l | grep q3 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat2_chr1_hqp_q3 --merge-list merge --make-bed --out cat2_hqp_q3

ls -l | grep q4 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat2_chr1_hqp_q4 --merge-list merge --make-bed --out cat2_hqp_q4


### GRM CREATION
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q1  --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q1
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q2  --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q2
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q3  --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q3
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q4  --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q4

### GRM cutoff 
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q1 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q1_0.05
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q2 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q2_0.05
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q3 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q3_0.05
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q4 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q4_0.05

###############################################################################################################################################################
# 							   category3 - [0.01,0.05)						       			      #
###############################################################################################################################################################
for i in  {1..22} ; do 
plink2 --bfile plink_format/original/freeze10.14k.chr${i}.0.0001_var_unrel --extract plink_format/original/category3_0.01_0.05.csv  --make-bed --out subset_for_h2_calc/category3_0.01_0.05/unrelated/category3_chr${i}  ; 
done

cd subset_for_h2_calc/category3_0.01_0.05/unrelated
ls -l | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile category3_chr1 --merge-list merge --make-bed --out cat3_chr_all

### Remove variants based on geno 0.05 , mind 0.05, phwe < 1e-6 
plink --bfile cat3_chr_all --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out cat3_chr

### Prunning to generate high quality variants
plink --bfile cat3_chr --indep-pairwise 50 5 0.1 --out cat3

## subset to prunned variants
plink --bfile cat3_chr --extract cat3.prune.in --make-bed --out cat3_chr_all_hqp

## LD score calculation
../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat3_chr_all_hqp --ld-score-region 200 --out cat3_chr_all_hqp_unrel --thread-num 100

### Categories based on quartiles as suggested in the paper  ## remove monomorphic variants before cuting to snps based on quartiles - code from gcta tutorial

lds_seg = read.table("cat3_chr_all_hqp_unrel.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "snp_group4.txt", row.names=F, quote=F, col.names=F)


### GRM CREATION
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat3_chr_all_hqp  --extract snp_group1.txt --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q1
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat3_chr_all_hqp  --extract snp_group2.txt --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q2
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat3_chr_all_hqp  --extract snp_group3.txt --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q3
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat3_chr_all_hqp  --extract snp_group4.txt --make-grm-alg 1 --thread-num 4 --out cat3_hqp_q4

### GRM cutoff 
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q1 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q1_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q2 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q2_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q3 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q3_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q4 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q4_0.05


###############################################################################################################################################################
# 							 	  category4 - [0.05]						       			      #
###############################################################################################################################################################
for i in  {1..22} ; do 
plink2 --bfile plink_format/original/freeze10.14k.chr${i}.0.0001_var_unrel --extract plink_format/original/category4_0.05.csv --make-bed --out subset_for_h2_calc/category4_0.05/unrelated/category4_chr${i}  ; 
done

cd subset_for_h2_calc/category4_0.05/unrelated/
ls -l | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile category40_chr1 --merge-list merge --make-bed --out cat4_chr_all

### Remove variants based on geno 0.05 , mind 0.05, phwe < 1e-6 
plink --bfile cat4_chr_all --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out cat4_chr

### Prunning to generate high quality variants
plink --bfile cat4_chr --indep-pairwise 50 5 0.1 --out cat4

### subset to prunned variants
plink --bfile cat4_chr --extract cat4.prune.in --make-bed --out cat4_chr_hqp

### LD score calculation
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile cat4_chr_hqp --ld-score-region 200 --out chr_all_unrel --thread-num 100

### Categories based on quartiles as suggested in the paper  ## remove monomorphic variants before cuting to snps based on quartiles - code from gcta tutorial

lds_seg = read.table("cat4_chr_all_hqp_unrel.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, "snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "snp_group4.txt", row.names=F, quote=F, col.names=F)

### GRM CREATION
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat4_chr_hqp --extract snp_group1.txt --make-grm-alg 1 --thread-num 4 --out cat4_hqp_q1
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat4_chr_hqp --extract snp_group2.txt --make-grm-alg 1 --thread-num 4 --out cat4_hqp_q2
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat4_chr_hqp --extract snp_group3.txt --make-grm-alg 1 --thread-num 4 --out cat4_hqp_q3
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat4_chr_hqp --extract snp_group4.txt --make-grm-alg 1 --thread-num 4 --out cat4_hqp_q4

### GRM cutoff 
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q1 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q1_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q2 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q2_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q3 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q3_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q4 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q4_0.05

###############################################################################################################################################################
# 							 		GRM CREATION						       			      #
###############################################################################################################################################################


#Construct a GRM, extracting a list of variants (see 01_LD_bins.R to create LD bins from a bed of a MAF range)
#GRM constructed by part (one part computed per script on an array job)


### task array category 1 - script in server
softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  \
--bfile heritability/subset_for_h2_calc/category2_0.001_0.01/cat2_chr_all_hqp \
--extract subset_for_h2_calc/category2_0.001_0.01/parts/cat2_chr${SLURM_ARRAY_TASK_ID}_hqp.bim \
--make-grm\
--out subset_for_h2_calc/category2_0.001_0.01/cat2_chr${SLURM_ARRAY_TASK_ID}_hqp
--make-grm-alg 1 \


### task array category 2 - script in server
#category2_0.001_0.01/parts/
## subset in r
require(data.table)
data = fread("cat2_chr_all_hqp.bim")
for(i in 1:22){
	data1 = data[which(data$V1 == i),]
	write.table(data1[,"V2"],file=paste0("cat2_chr",i,"_hqp.bim"),sep="\t",col.names=F,dec=".",quote=F,row.names=F)
	}

softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
--bfile heritability/subset_for_h2_calc/category2_0.001_0.01/cat2_chr_all_hqp \
--extract subset_for_h2_calc/category2_0.001_0.01/parts/cat2_chr${SLURM_ARRAY_TASK_ID}_hqp.bim \
--make-grm \
--out subset_for_h2_calc/category2_0.001_0.01/cat2_chr${SLURM_ARRAY_TASK_ID}_hqp \
--make-grm-alg 1

#Merge all GRM parts together

cat cat2_chr*_hqp.grm.id > cat2_hqp.grm.id
cat cat2_chr*_hqp.grm.bin > cat2_hqp.grm.bin
cat cat2_chr*_hqp.grm.N.bin >  cat2_hqp.grm.N.bin

## Cat3
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat3_chr_all_hqp --make-grm-alg 1 --thread-num 4 --out test
## Cat4
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat4_chr_all_hqp --make-grm-alg 1 --thread-num 4 --out test



#Merge all GRM parts together

cat cat2_chr*_hqp.grm.id   > test.grm.id
cat cat2_chr*_hqp.grm.bin > test.grm.bin
cat cat2_chr*_hqp.grm.N.bin > test.grm.N.bin


cat ${GRM_out}.part_99_*.grm.id > ${GRM_out}.grm.id
cat ${GRM_out}.part_99_*.grm.bin > ${GRM_out}.grm.bin
cat ${GRM_out}.part_99_*.grm.N.bin > ${GRM_out}.grm.N.bin

###### Extract unrelated samples ######
#Here we keep a relatedness threshold of 0.05

GCTA \
	--grm ${GRM_out} \
	--grm-cutoff 0.05 \
	--make-grm \
	--out ${GRM_out}_unrelated

for i in {1..22}; do
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm  cat2_chr${i}_hqp --grm-cutoff 0.05 --make-grm --out category4_0.05/cat2_chr${i}_hqp_unrelated
done

###### Create a file containing multiple GRMs in a directory (need full path) ######

for i in *unrelated.grm.bin ; do readlink -f "$i"  | cut -d'.' -f1-2 >>  path; done



###############################################################################################################################################################
# 							 		Compute PCs						       			      #
###############################################################################################################################################################

##cat1
../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm category1_0.0001_0.001/test --pca 20 --threads 10 --out category1_0.0001_0.001/test
##cat2
../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm category2_0.001_0.01/cat2_hqp --pca 20 --threads 10 --out category2_0.001_0.01/test
##cat3
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm  cat3_chr_hqp --pca 20 --threads 10 --out  cat3_chr_hqp
##cat4
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat4_chr_hqp --pca 20 --threads 10 --out cat4_chr_hqp



###############################################################################################################################################################
# 							 		REML							       			      #
###############################################################################################################################################################
#Run unconstrained REML from multiple GRMs (for GREML-LDMS), fitting PCs as quantitative covariates
#REML-no-lrt doesn't calculate the reduced model

## phenotype subset 6 different phenotypes
cat phenotypes/Combined_4cohorts_NTproBNP_08222022.tsv | awk ' { print 0,"\t",$1,"\t" $10 } ' > phenotypes/Combined_4cohorts_NTproBNP_08222022.phen

GCTA \
	--reml \
	--mgrm ${mgrm_file_path} \
	--reml-no-constrain \
	--pheno ${phenotype_file} \
	--out ${REML_output_file} \
	--thread-num ${ncpu} \
	--qcovar ${PCA_out} \
	--reml-no-lrt

### cat 1

### cat 2

### cat 3
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm test --reml-no-constrain --pheno phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar category3_0.01_0.05/test.eigenvec --out test --thread-num 10 --reml-no-lrt
### cat 4
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm test --reml-no-constrain --pheno phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar category4_0.05/test.eigenvec --out test --thread-num 10 --reml-no-lrt



###### Computing IBD ###### Extra
#Computing IBD segments. For computational reasons, it is recomended to prune/select common SNPs for faster runtime.
#module load KING/2.1.2-foss-2016a

#king -b ../category3_0.01_0.05/cat3_chr_all_hqp.bed --ibdseg --prefix ../category3_0.01_0.05/cat3_chr_all_hqp --cpus 4 --seglength 3
#king -b cat4_chr_all_hqp.bed --ibdseg --prefix cat4_chr_all_hqp --cpus 4 --seglength 3


#king -b ${BED_file_merged_QC}.bed  --ibdseg --prefix ${BED_file_merged_QC_out_name} --cpus ${ncpu} --seglength 3







###############################################################################
# 	Clumping and thresholding to get variants for conditional analysis    #
###############################################################################
module load PLINK ## update rsids to chr:pos:ref:alt
for i in {1..22}; do 
plink --bfile  ../../heritability/plink_format/originial/freeze10.14k.chr${i}.0.0001 --set-all-var-ids @:#:'$r':'$a'  --new-id-max-allele-len 1000  --make-bed --out ../../heritability/plink_format/originial/freeze10.14k.chr${i}.0.0001_var ; 
done


module load PLINK/1.90-foss-2016a
for i in {1..22}; do
plink --bfile ../../heritability/plink_format/originial/freeze10.14k.chr${i}.0.0001_var --clump overall_gwas_09162022.csv  --clump-p1 5e-9	 --clump-p2 0.01   --clump-r2 0.60 --clump-kb 250  --out conditional_analysis/snps_for_conditioning/freeze10.14k.chr${i}.0.0001
done




