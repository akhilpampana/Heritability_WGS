for i in  {1..22} ; do 
plink2 --bfile freeze10.14k.chr${i}.0.0001_var_unrel --keep ../../../../../phenotypes/black_NTproBNP_08222022.phen --make-bed --out black/chr${i}  ; 
plink2 --bfile freeze10.14k.chr${i}.0.0001_var_unrel --keep ../../../../../phenotypes/white_NTproBNP_08222022.phen --make-bed --out white/chr${i}  ; 
plink2 --bfile freeze10.14k.chr${i}.0.0001_var_unrel --keep ../../../../../phenotypes/male_NTproBNP_08222022.phen --make-bed --out male/chr${i}   ; 
plink2 --bfile freeze10.14k.chr${i}.0.0001_var_unrel --keep ../../../../../phenotypes/female_NTproBNP_08222022.phen --make-bed --out male/chr${i}  ; 
done


for i in  {1..22} ; do 
plink2 --bfile black/chr${i}  --freq --out black/chr${i}  ; 
plink2 --bfile  white/chr${i} --freq --out white/chr${i}  ; 
plink2 --bfile male/chr${i} --freq  --out male/chr${i}   ; 
plink2 --bfile male/chr${i} --freq  --out female/chr${i}  ; 
done




###############################################################################################################################################################
# 							   category1 - [0.0001,0.001)						       			      #
###############################################################################################################################################################

for i in  {1..22} ; do 
plink2 --bfile freeze10.14k.chr${i}.0.0001_var_unrel --extract plink_format/original/category1_0.0001_0.001.csv --make-bed --out subset_for_h2_calc/category1_0.0001_0.001/unrelated/category1_chr${i}  ; 
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

