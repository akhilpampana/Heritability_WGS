
##############################################################################################################################################################
# Script used to run GCTA/plink/bcftools for the analyses in the GREML-LDMS WGS 2022 paper		        			    		     #	
# It contains the most important commands for performing GREML-LDMS from imputed dataset                                  				     #	
# The Online Methods describe with greater details all the QC steps performed										     #
# Some steps are added from the original github page													     #
# Additional R scripts are also available in this git													     #
##############################################################################################################################################################

#### Variables
ncpu=12 
sample_list="/mnt/project/NTproBNP_OLINK/hqp_ntprobnp.fam" #List of samples with phenotypic information available

###############################################################################################################################################################
# 			Create BED files and perform first filtering on SVM variants and other QC metrics - submitted as slurm job			      #
###############################################################################################################################################################
###### Create BED files and perform first filtering on SVM variants and other QC metrics - submitted as slurm job
## Use more memory for chromosome 1-5

for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bgen /mnt/project/${input_path}/chr${i}.bgen ref-first --sample chr2_topmed_imputed.txt --maf 0.0001 --allow-extra-chr --keep-allele-order  --keep /mnt/project/SBP/Null_Models/FEMALE/females_final.fam --geno 0.05 --hwe 0.000001 --threads 20 --make-bed --out chr${i}_hqp" --destination "${path}/Subset_from_BGEN/" --tag "subset_from_bgen" --name "Subset from bgen: chr${i}"  --instance-type "mem1_ssd1_v2_x8"; 
done  

for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bgen /mnt/project/Bulk/Imputation/Imputation\ from\ genotype\ \(TOPmed\)/ukb21007_c${i}_b0_v1.bgen ref-first --sample /mnt/project/NTproBNP_OLINK/chr2_topmed_imputed.txt --maf 0.0001 --allow-extra-chr --keep-allele-order  --keep /mnt/project/SBP/Null_Models/FEMALE/females_final.fam --geno 0.05 --hwe 0.000001 --threads 20 --make-bed --out female_chr${i}_hqp" --destination "Test 3.0:SBP/Heritability/" --tag "subset_from_bgen" --name "Subset from bgen: chr${i}"  --instance-type "mem1_ssd1_v2_x72"; 
done  


for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bgen /mnt/project/Bulk/Imputation/Imputation\ from\ genotype\ \(TOPmed\)/ukb21007_c${i}_b0_v1.bgen ref-first --sample /mnt/project/NTproBNP_OLINK/chr2_topmed_imputed.txt --maf 0.0001 --allow-extra-chr --keep-allele-order  --keep /mnt/project/SBP/Null_Models/MALE/males_final.fam --geno 0.05 --hwe 0.000001 --threads 20 --make-bed --out male_chr${i}_hqp" --destination "Test 3.0:SBP/Heritability/" --tag "subset_from_bgen" --name "Male Subset from bgen: chr${i}"  --instance-type "mem1_ssd1_v2_x"; 
done  


for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bgen /mnt/project/${input_path}/chr${i}.bgen ref-first --sample chr2_topmed_imputed.txt --maf 0.0001 --allow-extra-chr --keep-allele-order  --keep /mnt/project/NTproBNP_OLINK/hqp_ntprobnp.fam --geno 0.05 --hwe 0.000001 --threads 20 --make-bed --out chr${i}_hqp" --destination "${path}/Subset_from_BGEN/" --tag "subset_from_bgen" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  


for i in  1 2 3 4 5 6 7 8 9 10 11 12 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bgen /mnt/project/Bulk/Imputation/Imputation\ from\ genotype\ \(TOPmed\)/ukb21007_c${i}_b0_v1.bgen ref-first --sample /mnt/project/NTproBNP_OLINK/chr2_topmed_imputed.txt --maf 0.0001 --allow-extra-chr --keep-allele-order  --keep /mnt/project/SBP/Null_Models/FEMALE/females_final.fam --geno 0.05 --hwe 0.000001 --threads 20 --make-bed --out female_chr${i}_hqp" --destination "Test 3.0:SBP/Heritability/" --tag "subset_from_bgen" --name "Subset from bgen: chr${i}"  --instance-type "mem1_ssd1_v2_x72"; 
done  


for i in  1 2 3 4 5 6 7 8 9 10 11 12 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bgen /mnt/project/Bulk/Imputation/Imputation\ from\ genotype\ \(TOPmed\)/ukb21007_c${i}_b0_v1.bgen ref-first --sample /mnt/project/NTproBNP_OLINK/chr2_topmed_imputed.txt --maf 0.0001 --allow-extra-chr --keep-allele-order  --keep /mnt/project/SBP/Null_Models/MALE/males_final.fam --geno 0.05 --hwe 0.000001 --threads 20 --make-bed --out male_chr${i}_hqp" --destination "Test 3.0:SBP/Heritability/" --tag "subset_from_bgen" --name "Male Subset from bgen: chr${i}"  --instance-type "mem1_ssd1_v2_x72"; 
done  




###############################################################################################################################################################
# 								Allele Frequency Calculations								      #
###############################################################################################################################################################
## Subset to unrelated and qc filtering 
for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/SBP/Heritability/male_chr${i}_hqp  --keep /mnt/project/SBP/Heritability/unrelated_ukbiobank_reported.tsv --set-all-var-ids @:#:'$r':'$a' --geno 0.05 --mind 0.05 --hwe 1e-6  --make-bed --out male_chr${i}_unrel" --destination "/mnt/project/SBP/Heritability/unrelated_qced/" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  

for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/SBP/Heritability/female_chr${i}_hqp  --keep /mnt/project/SBP/Heritability/unrelated_ukbiobank_reported.tsv --set-all-var-ids @:#:'$r':'$a' --geno 0.05 --mind 0.05 --hwe 1e-6  --make-bed --out female_chr${i}_unrel" --destination "/mnt/project/SBP/Heritability/unrelated_qced" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  



###############################################################################################################################################################
# 								Filtering to unrelated and genotype missingness, phwe filter				      #
###############################################################################################################################################################
## Subset to unrelated and qc filtering 
for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/SBP/Heritability/male_chr${i}_hqp  --keep /mnt/project/SBP/Heritability/unrelated_ukbiobank_reported.tsv  --set-all-var-ids "@:#:'$r':'$a'" --new-id-max-allele-len 1600 --geno 0.05 --mind 0.05 --hwe 1e-6  --make-bed --out male_chr${i}_unrel" --destination "/mnt/project/SBP/Heritability/unrelated_qced/" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  

for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/SBP/Heritability/female_chr${i}_hqp  --keep /mnt/project/SBP/Heritability/unrelated_ukbiobank_reported.tsv --set-all-var-ids "@:#:'$r':'$a'" --geno 0.05 --mind 0.05 --hwe 1e-6  --make-bed --out female_chr${i}_unrel" --destination "/mnt/project/SBP/Heritability/unrelated_qced" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  




###############################################################################################################################################################
# 							 LD prunning 											      #
###############################################################################################################################################################
## Subset to unrelated and qc filtering 10
#1 2 3 4
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink --bfile /mnt/project/mnt/project/SBP/Heritability/unrelated_qced/male_chr${i}_unrel  --indep-pairwise 50 5 0.1  --out male_chr${i}_prunned" --destination "Test 3.0:/SBP/Heritability/unrelated_qced/" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  
  
#1 2 
for i in  3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink --bfile /mnt/project/mnt/project/SBP/Heritability/unrelated_qced/female_chr${i}_unrel  --indep-pairwise 50 5 0.1 --out female_chr${i}_prunned" --destination "Test 3.0:/SBP/Heritability/unrelated_qced/" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  

for i in  1 2; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink --bfile /mnt/project/mnt/project/SBP/Heritability/unrelated_qced/female_chr${i}_unrel  --indep-pairwise 50 5 0.1 --out female_chr${i}_prunned" --destination "Test 3.0:/SBP/Heritability/unrelated_qced/" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  

###############################################################################################################################################################
# 							  			Subset to pruned set							      #
###############################################################################################################################################################
## Subset to unrelated and qc filtering 10
#1 2 3 4
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/mnt/project/SBP/Heritability/unrelated_qced/male_chr${i}_unrel  --extract  /mnt/project/SBP/Heritability/unrelated_qced/male_chr${i}_prunned.prune.in --make-bed --out male_chr${i}_prunned" --destination "Test 3.0:/SBP/Heritability/unrelated_qced/" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  
  
for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/mnt/project/SBP/Heritability/unrelated_qced/female_chr${i}_unrel   --extract  /mnt/project/SBP/Heritability/unrelated_qced/female_chr${i}_prunned.prune.in --make-bed --out female_chr${i}_prunned" --destination "Test 3.0:/SBP/Heritability/unrelated_qced/" --tag "qc_filter" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  


##################################################################################################################################################################
# Merging and AF groups and LD score groups
####################################################################################################################################################################

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/SBP/Heritability/unrelated_qced/male_chr${i}_prunned --extract /mnt/project/SBP/Heritability/unrelated_qced/PRS_Variants/male_posterior_prob_03292023.tsv --make-bed --out male_chr${i}_prunned_prs" --destination "Test 3.0:/SBP/Heritability/unrelated_qced/PRS_Variants/" --tag "subset_to_prs" --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  
  
for i in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; 
do  
echo "y" |  dx run app-swiss-army-knife -icmd="plink2 --bfile /mnt/project/SBP/Heritability/unrelated_qced/female_chr${i}_prunned  --extract /mnt/project/SBP/Heritability/unrelated_qced/PRS_Variants/female_posterior_prob_03292023.tsv --make-bed --out female_chr${i}_prunned_prs" --destination "Test 3.0:/SBP/Heritability/unrelated_qced/PRS_Variants/" --tag "subset_to_prs"  --name "Subset from bgen: chr${i}"  --instance-type "mem1_hdd1_v2_x96"; 
done  
## Old Code
###############################################################################################################################################################
# 							Subset variants as per freeze10 pca's generation encore	    					      #
###############################################################################################################################################################
## CHange from here
module load PLINK
for i in {1..22}; do
plink2 --bfile originial/freeze10.14k.chr${i}.0.0001 --extract variants to filter to --make-bed  --out plink_format/prunned_list_included_in_encore_pcs_generation/freeze10.14k.chr${i}.pruned
done

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
./plink2 --bfile /mnt/project/SBP/Heritability/female_chr${i}_hqp --keep unrelated_ukbiobank_reported.tsv --set-all-var-ids @:#:'$r':'$a' --geno 0.05 --mind 0.05 --hwe 1e-6  --make-bed --out female_unrel_chr${i}
done

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
./plink2 --bfile /mnt/project/SBP/Heritability/male_chr${i}_hqp --keep unrelated_ukbiobank_reported.tsv --set-all-var-ids @:#:'$r':'$a' --geno 0.05 --mind 0.05 --hwe 1e-6  --make-bed --out male_unrel_chr${i}
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


#### 
tmp1 = tmp[which(tmp$AF >= 0.0001 & tmp$MAF < 0.001),]
category1 = rbind(category1,tmp1)

tmp1 = tmp[which(tmp$AF >= 0.001 & tmp$AF < 0.01),]
category2 = rbind(category2,tmp1)
	
tmp1 = tmp[which(tmp$AF >= 0.01 & tmp$AF < 0.05),]
category3 = rbind(category3,tmp1)
	
tmp1 = tmp[which(tmp$AF >= 0.05),]
category4 = rbind(category4,tmp1)

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
lds_seg = read.table(paste0("cat1_chr",i,"_hqp.score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
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


### SUBSET TO 4 QUARTILES TO GENERATE GRM 
module load PLINK/1.90-foss-2016a
for i in {1..22}; do
plink --bfile cat1_chr${i}_hqp --extract snp_group1.txt --make-bed --out quartiles/cat1_chr${i}_hqp_q1
plink --bfile cat1_chr${i}_hqp --extract snp_group2.txt --make-bed --out quartiles/cat1_chr${i}_hqp_q2
plink --bfile cat1_chr${i}_hqp --extract snp_group3.txt --make-bed --out quartiles/cat1_chr${i}_hqp_q3
plink --bfile cat1_chr${i}_hqp --extract snp_group4.txt --make-bed --out quartiles/cat1_chr${i}_hqp_q4
done


### MERGE TO ONE MEGA FILE
ls -l | grep q1 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat1_chr1_hqp_q1 --merge-list merge --make-bed --out cat1_hqp_q1

ls -l | grep q2 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat1_chr1_hqp_q2 --merge-list merge --make-bed --out cat1_hqp_q2

ls -l | grep q3 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat1_chr1_hqp_q3 --merge-list merge --make-bed --out cat1_hqp_q3

ls -l | grep q4 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
plink --bfile cat1_chr1_hqp_q4 --merge-list merge --make-bed --out cat1_hqp_q4


### GRM CREATION
### Job submissions for grms creation in 5 parts


### merge files together

cat cat1_hqp_q2.part_5_*.grm.id > cat1_hqp_q1.grm.id
cat cat1_hqp_q2.part_5_*.grm.bin > cat1_hqp_q1.grm.bin
cat cat1_hqp_q2.part_5_*.grm.N.bin > cat1_hqp_q1.grm.N.bin

cat cat1_hqp_q2.part_5_*.grm.id > cat1_hqp_q2.grm.id
cat cat1_hqp_q2.part_5_*.grm.bin > cat1_hqp_q2.grm.bin
cat cat1_hqp_q2.part_5_*.grm.N.bin > cat1_hqp_q2.grm.N.bin

cat cat1_hqp_q3.part_5_*.grm.id > cat1_hqp_q3.grm.id
cat cat1_hqp_q3.part_5_*.grm.bin > cat1_hqp_q3.grm.bin
cat cat1_hqp_q3.part_5_*.grm.N.bin > cat1_hqp_q3.grm.N.bin

cat cat1_hqp_q4.part_5_*.grm.id > cat1_hqp_q4.grm.id
cat cat1_hqp_q4.part_5_*.grm.bin > cat1_hqp_q4.grm.bin
cat cat1_hqp_q4.part_5_*.grm.N.bin > cat1_hqp_q4.grm.N.bin

				OR
				
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --mgrm Q3/merge --make-grm --out Q3/cat1_hqp_q3

### GRM cutoff
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1 --grm-cutoff 0.05 --make-grm --threads 10 --out cat1_hqp_q1_0.05
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2 --grm-cutoff 0.05 --make-grm --threads 10 --out cat1_hqp_q2_0.05
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3 --grm-cutoff 0.05 --make-grm --threads 10 --out cat1_hqp_q3_0.05
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4 --grm-cutoff 0.05 --make-grm --threads 10 --out cat1_hqp_q4_0.05

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


### GRM CREATION - job submission
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q1  --make-grm-alg 1 --thread-num 4 --out cat2_hqp_q1
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q2  --make-grm-alg 1 --thread-num 4 --out cat2_hqp_q2
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q3  --make-grm-alg 1 --thread-num 4 --out cat2_hqp_q3
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat2_hqp_q4  --make-grm-alg 1 --thread-num 4 --out cat2_hqp_q4

### GRM cutoff 
#../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q1 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q1_0.05
#../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q2 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q2_0.05
#../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q3 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q3_0.05
#../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q4 --grm-cutoff 0.05 --make-grm --out cat2_hqp_q4_0.05

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
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q1 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q1_0.05
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q2 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q2_0.05
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q3 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q3_0.05
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q4 --grm-cutoff 0.05 --make-grm --out cat3_hqp_q4_0.05


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
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q1 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q1_0.05
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q2 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q2_0.05
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q3 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q3_0.05
#../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q4 --grm-cutoff 0.05 --make-grm --out cat4_hqp_q4_0.05



###############################################################################################################################################################
# 							 		Compute PCs						       			      #
###############################################################################################################################################################

##cat1
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --pca 20 --threads 10 --out cat1_hqp_q1_0.05
#../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --pca 20 --threads 10 --out cat1_hqp_q2_0.05
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --pca 20 --threads 10 --out cat1_hqp_q3_0.05
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --pca 20 --threads 10 --out cat1_hqp_q4_0.05


##cat2
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q1_0.05 --pca 20 --threads 10 --out cat2_hqp_q1_0.05
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q2_0.05 --pca 20 --threads 10 --out cat2_hqp_q2_0.05
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q3_0.05 --pca 20 --threads 10 --out cat2_hqp_q3_0.05
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q4_0.05 --pca 20 --threads 10 --out cat2_hqp_q4_0.05
##cat3
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q1_0.05 --pca 20 --threads 10 --out cat3_hqp_q1_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q2_0.05 --pca 20 --threads 10 --out cat3_hqp_q2_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q3_0.05 --pca 20 --threads 10 --out cat3_hqp_q3_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q4_0.05 --pca 20 --threads 10 --out cat3_hqp_q4_0.05

##cat4
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q1_0.05 --pca 20 --threads 10 --out cat4_hqp_q1_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q2_0.05 --pca 20 --threads 10 --out cat4_hqp_q2_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q3_0.05 --pca 20 --threads 10 --out cat4_hqp_q3_0.05
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q4_0.05 --pca 20 --threads 10 --out cat4_hqp_q4_0.05


###############################################################################################################################################################
# 							 		REML							       			      #
###############################################################################################################################################################
#Run unconstrained REML from multiple GRMs (for GREML-LDMS), fitting PCs as quantitative covariates
#REML-no-lrt doesn't calculate the reduced model - repeat the multi ethnic step with all other categories

## phenotype subset 6 different phenotypes
cat phenotypes/Combined_4cohorts_NTproBNP_08222022.tsv | awk ' { print 0,"\t",$1,"\t" $10 } ' > phenotypes/Combined_4cohorts_NTproBNP_08222022.phen

### cat 1
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat1_hqp_q1_0.05.eigenvec --out cat1_hqp_q1_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat1_hqp_q2_0.05.eigenvec --out cat1_hqp_q2_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat1_hqp_q3_0.05.eigenvec --out cat1_hqp_q3_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat1_hqp_q4_0.05.eigenvec --out cat1_hqp_q4_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000

../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/black_NTproBNP_08222022.phen --qcovar cat1_hqp_q1_0.05.eigenvec --out cat1_hqp_q1_0.05_afr --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/black_NTproBNP_08222022.phen --qcovar cat1_hqp_q2_0.05.eigenvec --out cat1_hqp_q2_0.05_afr --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/black_NTproBNP_08222022.phen --qcovar cat1_hqp_q3_0.05.eigenvec --out cat1_hqp_q3_0.05_afr --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/black_NTproBNP_08222022.phen --qcovar cat1_hqp_q4_0.05.eigenvec --out cat1_hqp_q4_0.05_afr --thread-num 10 --reml-no-lrt --reml-maxit 10000


../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/white_NTproBNP_08222022.phen --qcovar cat1_hqp_q1_0.05.eigenvec --out cat1_hqp_q1_0.05_eur --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/white_NTproBNP_08222022.phen --qcovar cat1_hqp_q2_0.05.eigenvec --out cat1_hqp_q2_0.05_eur --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/white_NTproBNP_08222022.phen --qcovar cat1_hqp_q3_0.05.eigenvec --out cat1_hqp_q3_0.05_eur --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/white_NTproBNP_08222022.phen --qcovar cat1_hqp_q4_0.05.eigenvec --out cat1_hqp_q4_0.05_eur --thread-num 10 --reml-no-lrt --reml-maxit 10000

../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/others_NTproBNP_08222022.phen --qcovar cat1_hqp_q1_0.05.eigenvec --out cat1_hqp_q1_0.05_others --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/others_NTproBNP_08222022.phen --qcovar cat1_hqp_q2_0.05.eigenvec --out cat1_hqp_q2_0.05_others --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/others_NTproBNP_08222022.phen --qcovar cat1_hqp_q3_0.05.eigenvec --out cat1_hqp_q3_0.05_others --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/others_NTproBNP_08222022.phen --qcovar cat1_hqp_q4_0.05.eigenvec --out cat1_hqp_q4_0.05_others --thread-num 10 --reml-no-lrt --reml-maxit 10000

../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/male_NTproBNP_08222022.phen --qcovar cat1_hqp_q1_0.05.eigenvec --out cat1_hqp_q1_0.05_male --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/male_NTproBNP_08222022.phen --qcovar cat1_hqp_q2_0.05.eigenvec --out cat1_hqp_q2_0.05_male --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/male_NTproBNP_08222022.phen --qcovar cat1_hqp_q3_0.05.eigenvec --out cat1_hqp_q3_0.05_male --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/male_NTproBNP_08222022.phen --qcovar cat1_hqp_q4_0.05.eigenvec --out cat1_hqp_q4_0.05_male --thread-num 10 --reml-no-lrt --reml-maxit 10000

../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/female_NTproBNP_08222022.phen --qcovar cat1_hqp_q1_0.05.eigenvec --out cat1_hqp_q1_0.05_female --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/female_NTproBNP_08222022.phen --qcovar cat1_hqp_q2_0.05.eigenvec --out cat1_hqp_q2_0.05_female --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/female_NTproBNP_08222022.phen --qcovar cat1_hqp_q3_0.05.eigenvec --out cat1_hqp_q3_0.05_female --thread-num 10 --reml-no-lrt --reml-maxit 10000
../gcta-1.94.1-linux-kernel-2-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../../../phenotypes/female_NTproBNP_08222022.phen --qcovar cat1_hqp_q4_0.05.eigenvec --out cat1_hqp_q4_0.05_female --thread-num 10 --reml-no-lrt --reml-maxit 10000


### cat 2
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat2_hqp_q1_0.05.eigenvec --out cat2_hqp_q1_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat2_hqp_q2_0.05.eigenvec --out cat2_hqp_q2_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat2_hqp_q3_0.05.eigenvec --out cat2_hqp_q3_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat2_hqp_q4_0.05.eigenvec --out cat2_hqp_q4_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000


### cat 3
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat3_hqp_q1_0.05.eigenvec --out cat3_hqp_q1_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat3_hqp_q2_0.05.eigenvec --out cat3_hqp_q2_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat3_hqp_q3_0.05.eigenvec --out cat3_hqp_q3_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat3_hqp_q4_0.05.eigenvec --out cat3_hqp_q4_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000

### cat 4
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q1_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat4_hqp_q1_0.05.eigenvec --out cat4_hqp_q1_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q2_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat4_hqp_q2_0.05.eigenvec --out cat4_hqp_q2_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q3_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat4_hqp_q3_0.05.eigenvec --out cat4_hqp_q3_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000
../../../../../softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat4_hqp_q4_0.05 --reml-no-constrain --pheno ../../../../../phenotypes/Combined_4cohorts_NTproBNP_08222022.phen --qcovar cat4_hqp_q4_0.05.eigenvec --out cat4_hqp_q4_0.05 --thread-num 10 --reml-no-lrt --reml-maxit 10000



######################################################################################################################################################
# 							Clumping and thresholding to get variants for conditional analysis     			     #
######################################################################################################################################################
module load PLINK ## update rsids to chr:pos:ref:alt
for i in {1..22}; do 
plink --bfile  ../../heritability/plink_format/originial/freeze10.14k.chr${i}.0.0001 --set-all-var-ids @:#:'$r':'$a'  --new-id-max-allele-len 1000  --make-bed --out ../../heritability/plink_format/originial/freeze10.14k.chr${i}.0.0001_var ; 
done

### 45 variants among 5 cytobands where identified
module load PLINK/1.90-foss-2016a
for i in {1..22}; do
plink --bfile ../../heritability/plink_format/originial/freeze10.14k.chr${i}.0.0001_var --clump overall_gwas_09162022.csv  --clump-p1 5e-9	 --clump-p2 0.01   --clump-r2 0.60 --clump-kb 500  --out conditional_analysis/snps_for_conditioning/freeze10.14k.chr${i}.0.0001
done

### Subset to each loci and do suggestive significance fitering +/- 500kb per loci (ran in R)
module load PLINK/1.90-foss-2016a 

load("overall_loci_09202022.rda")
for(i in 1:5) {
  system(paste0("plink --bfile ../heritability/plink_format/original/freeze10.14k.chr",df$chr[i],".0.0001_var --chr ",df$chr[i]," --to-bp ",df$pos_500kb[i]," --from-bp ",df$neg_500kb[i]," --make-bed --out /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/coloc/chr",df$chr[i],"loci",i,"")) 
}


plink --bfile chr1loci1 --clump ../conditional_analysis/snps_for_conditioning/overall_gwas_09162022.csv --clump-p1 5e-9 --clump-p2 5e-7 --clump-r2 0.60 --clump-kb 500 --out chr1loci1
plink --bfile chr4loci2 --clump ../conditional_analysis/snps_for_conditioning/overall_gwas_09162022.csv --clump-p1 5e-9 --clump-p2 5e-7 --clump-r2 0.60 --clump-kb 500 --out chr4loci2
plink --bfile chr8loci3 --clump ../conditional_analysis/snps_for_conditioning/overall_gwas_09162022.csv --clump-p1 5e-9 --clump-p2 5e-7 --clump-r2 0.60 --clump-kb 500 --out chr8loci3
plink --bfile chr8loci4 --clump ../conditional_analysis/snps_for_conditioning/overall_gwas_09162022.csv --clump-p1 5e-9 --clump-p2 5e-7 --clump-r2 0.60 --clump-kb 500 --out chr8loci4
plink --bfile chr12loci5 --clump ../conditional_analysis/snps_for_conditioning/overall_gwas_09162022.csv --clump-p1 5e-9 --clump-p2 5e-7 --clump-r2 0.60 --clump-kb 500 --out chr12loci5

#### Clumping and thresholding for all dataset p < 5e-7 
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
plink	--bfile ../../../heritability/plink_format/original/rsid/freeze10.14k.chr${i}.0.0001  --write-snplist --out snplist #_chr${i}
cat snplist.snplist | sort | uniq -d > duplicated_snps.snplist
plink --bfile ../../../heritability/plink_format/original/rsid/freeze10.14k.chr${i}.0.0001 --exclude duplicated_snps.snplist --make-bed --out freeze10.14k.chr${i}.0.0001
plink --bfile freeze10.14k.chr${i}.0.0001 --clump gwas_variants_12162022.csv --clump-kb 500 --clump-p1 5e-7 --clump-p2 0.05 --clump-r2 0.60 --out group3_chr${i} 
done

### New Version
module load PLINK/1.90-foss-2016a
for i in 14 15 16 17 18 19 20 21 22 
do
plink --bfile chr${i}.pass_only_all --clump TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/gwas_orginal/tophits_0.05_12272022.csv --clump-kb 500 --clump-p1 5e-7 --clump-p2 0.05 --clump-r2 0.60 --out chr${i}_5e07_0.05_500 
done


for i in 16
do
plink --bfile chr${i}.pass_only_all --clump TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/gwas_orginal/tophits_0.05_12272022.csv --clump-kb 500 --clump-p1 5e-7 --clump-p2 0.05 --clump-r2 0.60 --out chr${i}_5e07_0.05_500 
done
