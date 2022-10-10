### ancestry subsets
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


### Variants subset to maf bins

#black
require(data.table)
category1 = data.frame() ## [0.0001,0.001)
category2 = data.frame() ## [0.001,0.01)
category3 = data.frame() ## [0.01,0.05)
category4 = data.frame() ## [0.05]
for(i in 1:22) {
	tmp = fread(paste0("black/chr",i,".afreq"))
	tmp = tmp[which(!(tmp$ALT_FREQS == 0)),]
	tmp = tmp[which(tmp$ALT_FREQS >= 0.0001),]
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.0001 & tmp$ALT_FREQS < 0.001),]
	category1 = rbind(category1,tmp1)

	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.001 & tmp$ALT_FREQS < 0.01),]
	category2 = rbind(category2,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.01 & tmp$ALT_FREQS < 0.05),]
	category3 = rbind(category3,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.05),]
	category4 = rbind(category4,tmp1)
	
	}

save(category1,category2,category3,category4,file="black/allele_frequence_cutoff_4cat_unrel_10062022.rda")

write.table(category1,file="black/category1_0.0001_0.001.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category2,file="black/category2_0.001_0.01.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category3,file="black/category3_0.01_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category4,file="black/category4_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)



#white
require(data.table)
category1 = data.frame() ## [0.0001,0.001)
category2 = data.frame() ## [0.001,0.01)
category3 = data.frame() ## [0.01,0.05)
category4 = data.frame() ## [0.05]
for(i in 1:22) {
	tmp = fread(paste0("white/chr",i,".afreq"))
	tmp = tmp[which(!(tmp$ALT_FREQS == 0)),]
	tmp = tmp[which(tmp$ALT_FREQS >= 0.0001),]
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.0001 & tmp$ALT_FREQS < 0.001),]
	category1 = rbind(category1,tmp1)

	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.001 & tmp$ALT_FREQS < 0.01),]
	category2 = rbind(category2,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.01 & tmp$ALT_FREQS < 0.05),]
	category3 = rbind(category3,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.05),]
	category4 = rbind(category4,tmp1)
	
	}

save(category1,category2,category3,category4,file="white/allele_frequence_cutoff_4cat_unrel_10062022.rda")

write.table(category1,file="white/category1_0.0001_0.001.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category2,file="white/category2_0.001_0.01.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category3,file="white/category3_0.01_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category4,file="white//category4_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)



#male
require(data.table)
category1 = data.frame() ## [0.0001,0.001)
category2 = data.frame() ## [0.001,0.01)
category3 = data.frame() ## [0.01,0.05)
category4 = data.frame() ## [0.05]
for(i in 1:22) {
	tmp = fread(paste0("male/chr",i,".afreq"))
	tmp = tmp[which(!(tmp$ALT_FREQS == 0)),]
	tmp = tmp[which(tmp$ALT_FREQS >= 0.0001),]
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.0001 & tmp$ALT_FREQS < 0.001),]
	category1 = rbind(category1,tmp1)

	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.001 & tmp$ALT_FREQS < 0.01),]
	category2 = rbind(category2,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.01 & tmp$ALT_FREQS < 0.05),]
	category3 = rbind(category3,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.05),]
	category4 = rbind(category4,tmp1)
	
	}

save(category1,category2,category3,category4,file="male/allele_frequence_cutoff_4cat_unrel_10062022.rda")

write.table(category1,file="male/category1_0.0001_0.001.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category2,file="male/category2_0.001_0.01.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category3,file="male/category3_0.01_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category4,file="male/category4_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)



#female
require(data.table)
category1 = data.frame() ## [0.0001,0.001)
category2 = data.frame() ## [0.001,0.01)
category3 = data.frame() ## [0.01,0.05)
category4 = data.frame() ## [0.05]
for(i in 1:22) {
	tmp = fread(paste0("female/chr",i,".afreq"))
	tmp = tmp[which(!(tmp$ALT_FREQS == 0)),]
	tmp = tmp[which(tmp$ALT_FREQS >= 0.0001),]
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.0001 & tmp$ALT_FREQS < 0.001),]
	category1 = rbind(category1,tmp1)

	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.001 & tmp$ALT_FREQS < 0.01),]
	category2 = rbind(category2,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.01 & tmp$ALT_FREQS < 0.05),]
	category3 = rbind(category3,tmp1)
	
	tmp1 = tmp[which(tmp$ALT_FREQS >= 0.05),]
	category4 = rbind(category4,tmp1)
	
	}

save(category1,category2,category3,category4,file="female/allele_frequence_cutoff_4cat_unrel_10062022.rda")

write.table(category1,file="female/category1_0.0001_0.001.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category2,file="female/category2_0.001_0.01.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category3,file="female/category3_0.01_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category4,file="female/category4_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)


#### QC
### Remove variants based on geno 0.05 , mind 0.05, phwe < 1e-6 
for i in {1..22}; do
plink2 --bfile black/chr${i} --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out black/chr${i}_qc
plink2 --bfile white/chr${i} --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out white/chr${i}_qc
plink2 --bfile male/chr${i} --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out male/chr${i}_qc
plink2 --bfile female/chr${i} --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out female/chr${i}_qc
done

### Prunning to generate high quality variants
for i in {1..22}; do
plink2 --bfile black/chr${i}_qc --extract black/category1_0.0001_0.001.csv --indep-pairwise 50 5 0.1 --out black/prune/cat1_${i}
plink2 --bfile black/chr${i}_qc --extract black/category2_0.001_0.01.csv --indep-pairwise 50 5 0.1 --out black/prune/cat2_${i}
plink2 --bfile black/chr${i}_qc --extract black/category3_0.01_0.05.csv  --indep-pairwise 50 5 0.1 --out black/prune/cat3_${i}
plink2 --bfile black/chr${i}_qc --extract black/category4_0.05.csv --indep-pairwise 50 5 0.1 --out black/prune/cat4_${i}
done

for i in {1..22}; do
plink2 --bfile white/chr${i}_qc --extract white/category1_0.0001_0.001.csv --indep-pairwise 50 5 0.1 --out white/prune/cat1_${i}
plink2 --bfile white/chr${i}_qc --extract white/category2_0.001_0.01.csv --indep-pairwise 50 5 0.1 --out white/prune/cat2_${i}
plink2 --bfile white/chr${i}_qc --extract white/category3_0.01_0.05.csv  --indep-pairwise 50 5 0.1 --out white/prune/cat3_${i}
plink2 --bfile white/chr${i}_qc --extract white/category4_0.05.csv --indep-pairwise 50 5 0.1 --out white/prune/cat4_${i}
done


for i in {1..22}; do
plink2 --bfile male/chr${i}_qc --extract male/category1_0.0001_0.001.csv --indep-pairwise 50 5 0.1 --out male/prune/cat1_${i}
plink2 --bfile male/chr${i}_qc --extract male/category2_0.001_0.01.csv --indep-pairwise 50 5 0.1 --out male/prune/cat2_${i}
plink2 --bfile male/chr${i}_qc --extract male/category3_0.01_0.05.csv  --indep-pairwise 50 5 0.1 --out male/prune/cat3_${i}
plink2 --bfile male/chr${i}_qc --extract male/category4_0.05.csv --indep-pairwise 50 5 0.1 --out male/prune/cat4_${i}
done



for i in {1..22}; do
plink2 --bfile female/chr${i}_qc --extract female/category1_0.0001_0.001.csv --indep-pairwise 50 5 0.1 --out female/prune/cat1_${i}
plink2 --bfile female/chr${i}_qc --extract female/category2_0.001_0.01.csv --indep-pairwise 50 5 0.1 --out female/prune/cat2_${i}
plink2 --bfile female/chr${i}_qc --extract female/category3_0.01_0.05.csv  --indep-pairwise 50 5 0.1 --out female/prune/cat3_${i}
plink2 --bfile female/chr${i}_qc --extract female/category4_0.05.csv --indep-pairwise 50 5 0.1 --out female/prune/cat4_${i}
done

 cat female/prune/cat1_*.prune.in > female/prune/cat1.prune.id
 cat female/prune/cat2_*.prune.in > female/prune/cat2.prune.id
 cat female/prune/cat3_*.prune.in > female/prune/cat3.prune.id
 cat female/prune/cat4_*.prune.in > female/prune/cat4.prune.id

 cat male/prune/cat1_*.prune.in > male/prune/cat1.prune.id
 cat male/prune/cat2_*.prune.in > male/prune/cat2.prune.id
 cat male/prune/cat3_*.prune.in > male/prune/cat3.prune.id
 cat male/prune/cat4_*.prune.in > male/prune/cat4.prune.id


 cat white/prune/cat1_*.prune.in > white/prune/cat1.prune.id
 cat white/prune/cat2_*.prune.in > white/prune/cat2.prune.id
 cat white/prune/cat3_*.prune.in > white/prune/cat3.prune.id
 cat white/prune/cat4_*.prune.in > white/prune/cat4.prune.id


 cat black/prune/cat1_*.prune.in > black/prune/cat1.prune.id
 cat black/prune/cat2_*.prune.in > black/prune/cat2.prune.id
 cat black/prune/cat3_*.prune.in > black/prune/cat3.prune.id
 cat black/prune/cat4_*.prune.in > black/prune/cat4.prune.id


## LD score calculation - scripts job submissions


### Categories based on quartiles as suggested in the paper  ## remove monomorphic variants before cuting to snps based on quartiles - code from gcta tutorial
## Male
require(data.table)
category1 = data.frame() 
category2 = data.frame() 
category3 = data.frame() 
category4 = data.frame() 

#./male/ld_sc/cat1_chr10_hqp.score.ld
#./female/ld_sc/cat1_chr10_hqp.score.ld
#./black/ld_sc/cat1_chr10_hqp.score.ld
#./white/ld_sc/cat1_chr10_hqp.score.ld

for(i in 1:22) {
lds_seg = read.table(paste0("/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/ancestry_specific/male/ld_sc/cat1_chr",i,"_hqp.score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
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

write.table(category1, "male/ld_sc/snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(category2, "male/ld_sc/snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(category3, "male/ld_sc/snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(category4, "male/ld_sc/snp_group4.txt", row.names=F, quote=F, col.names=F)



## Female
require(data.table)
category1 = data.frame() 
category2 = data.frame() 
category3 = data.frame() 
category4 = data.frame() 

for(i in 1:22) {
lds_seg = read.table(paste0("/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/ancestry_specific/female/ld_sc/cat1_chr",i,"_hqp.score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
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

write.table(category1, "female/ld_sc/snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(category2, "female/ld_sc/snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(category3, "female/ld_sc/snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(category4, "female/ld_sc/snp_group4.txt", row.names=F, quote=F, col.names=F



## Black
require(data.table)
category1 = data.frame() 
category2 = data.frame() 
category3 = data.frame() 
category4 = data.frame() 

for(i in 1:22) {
lds_seg = read.table(paste0("/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/ancestry_specific/black/ld_sc/cat1_chr",i,"_hqp.score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
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

write.table(category1, "black/ld_sc/snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(category2, "black/ld_sc/snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(category3, "black/ld_sc/snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(category4, "black/ld_sc/snp_group4.txt", row.names=F, quote=F, col.names=F)


## White
require(data.table)
category1 = data.frame() 
category2 = data.frame() 
category3 = data.frame() 
category4 = data.frame() 

for(i in 1:22) {
lds_seg = read.table(paste0("/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/ancestry_specific/white/ld_sc/cat1_chr",i,"_hqp.score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
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

write.table(category1, "white/ld_sc/snp_group1.txt", row.names=F, quote=F, col.names=F)
write.table(category2, "white/ld_sc/snp_group2.txt", row.names=F, quote=F, col.names=F)
write.table(category3, "white/ld_sc/snp_group3.txt", row.names=F, quote=F, col.names=F)
write.table(category4, "white/ld_sc/snp_group4.txt", row.names=F, quote=F, col.names=F)


### GRM generation for each category and ld score quartiles
for i in {1..22}
do
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat1_group1.txt --make-grm --make-grm-alg 1  --out black/GRM/cat1_chr${i}_hqp_q1 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat1_group2.txt --make-grm --make-grm-alg 1  --out black/GRM/cat1_chr${i}_hqp_q2 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat1_group3.txt --make-grm --make-grm-alg 1  --out black/GRM/cat1_chr${i}_hqp_q3 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat1_group4.txt --make-grm --make-grm-alg 1  --out black/GRM/cat1_chr${i}_hqp_q4 --thread-num 2
done



for i in {1..22}
do
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat2_group1.txt --make-grm --make-grm-alg 1  --out black/GRM/cat2_chr${i}_hqp_q1 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat2_group2.txt --make-grm --make-grm-alg 1  --out black/GRM/cat2_chr${i}_hqp_q2 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat2_group3.txt --make-grm --make-grm-alg 1  --out black/GRM/cat2_chr${i}_hqp_q3 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat2_group4.txt --make-grm --make-grm-alg 1  --out black/GRM/cat2_chr${i}_hqp_q4 --thread-num 2
done



for i in {1..22}
do
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat3_group1.txt --make-grm --make-grm-alg 1  --out black/GRM/cat3_chr${i}_hqp_q1 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat3_group2.txt --make-grm --make-grm-alg 1  --out black/GRM/cat3_chr${i}_hqp_q2 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat3_group3.txt --make-grm --make-grm-alg 1  --out black/GRM/cat3_chr${i}_hqp_q3 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat3_group4.txt --make-grm --make-grm-alg 1  --out black/GRM/cat3_chr${i}_hqp_q4 --thread-num 2
done




for i in {1..22}
do
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat4_group1.txt --make-grm --make-grm-alg 1  --out black/GRM/cat4_chr${i}_hqp_q1 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat4_group2.txt --make-grm --make-grm-alg 1  --out black/GRM/cat4_chr${i}_hqp_q2 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat4_group3.txt --make-grm --make-grm-alg 1  --out black/GRM/cat4_chr${i}_hqp_q3 --thread-num 2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile black/chr${i}_qc --extract black/ld_sc/snp_cat4_group4.txt --make-grm --make-grm-alg 1  --out black/GRM/cat4_chr${i}_hqp_q4 --thread-num 2
done

### merge list merged files
ls -l | grep cat1 | grep N.bin | grep q1 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat1_q1
ls -l | grep cat1 | grep N.bin | grep q2 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat1_q2
ls -l | grep cat1 | grep N.bin | grep q3 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat1_q3
ls -l | grep cat1 | grep N.bin | grep q4 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat1_q4


ls -l | grep cat2 | grep N.bin | grep q1 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat2_q1
ls -l | grep cat2 | grep N.bin | grep q2 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat2_q2
ls -l | grep cat2 | grep N.bin | grep q3 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat2_q3
ls -l | grep cat2 | grep N.bin | grep q4 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat2_q4

ls -l | grep cat3 | grep N.bin | grep q1 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat3_q1
ls -l | grep cat3 | grep N.bin | grep q2 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat3_q2
ls -l | grep cat3 | grep N.bin | grep q3 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat3_q3
ls -l | grep cat3 | grep N.bin | grep q4 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat3_q4


ls -l | grep cat4 | grep N.bin | grep q1 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat4_q1
ls -l | grep cat4 | grep N.bin | grep q2 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat4_q2
ls -l | grep cat4 | grep N.bin | grep q3 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat4_q3
ls -l | grep cat4 | grep N.bin | grep q4 | awk ' { print $9 } ' | sed 's|.grm.N.bin||g'  > merge_cat4_q4


### merge
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat1_q1  --make-grm --threads 10 --out cat1_hqp_q1
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat1_q2  --make-grm --threads 10 --out cat1_hqp_q2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat1_q3  --make-grm --threads 10 --out cat1_hqp_q3
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat1_q4  --make-grm --threads 10 --out cat1_hqp_q4

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat2_q1  --make-grm --threads 10 --out cat2_hqp_q1
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat2_q2  --make-grm --threads 10 --out cat2_hqp_q2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat2_q3  --make-grm --threads 10 --out cat2_hqp_q3
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat2_q4  --make-grm --threads 10 --out cat2_hqp_q4

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat3_q1  --make-grm --threads 10 --out cat3_hqp_q1
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat3_q2  --make-grm --threads 10 --out cat3_hqp_q2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat3_q3  --make-grm --threads 10 --out cat3_hqp_q3
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat3_q4  --make-grm --threads 10 --out cat3_hqp_q4

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat4_q1  --make-grm --threads 10 --out cat4_hqp_q1
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat4_q2  --make-grm --threads 10 --out cat4_hqp_q2
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat4_q3  --make-grm --threads 10 --out cat4_hqp_q3
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm merge_cat4_q4  --make-grm --threads 10 --out cat4_hqp_q4






### Subset to sparse GRM - cd to GRM folder
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat1_hqp_q1 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat1_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat1_hqp_q2 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat1_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat1_hqp_q3 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat1_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat1_hqp_q4 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat1_hqp_q4_0.05

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat2_hqp_q1 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat2_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat2_hqp_q2 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat2_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat2_hqp_q3 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat2_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat2_hqp_q4 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat2_hqp_q4_0.05

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat3_hqp_q1 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat3_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat3_hqp_q2 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat3_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat3_hqp_q3 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat3_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat3_hqp_q4 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat3_hqp_q4_0.05

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat4_hqp_q1 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat4_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat4_hqp_q2 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat4_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat4_hqp_q3 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat4_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm black/GRM/cat4_hqp_q4 --grm-cutoff 0.05 --make-grm --threads 10 --out black/GRM/cat4_hqp_q4_0.05



 
## PCs calculation for each category across each ancestry / sex group - cd to sparse GRM
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat1_hqp_q1_0.05 --pca 20 --threads 10 --out cat1_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat1_hqp_q2_0.05 --pca 20 --threads 10 --out cat1_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat1_hqp_q3_0.05 --pca 20 --threads 10 --out cat1_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat1_hqp_q4_0.05 --pca 20 --threads 10 --out cat1_hqp_q4_0.05


/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat2_hqp_q1_0.05 --pca 20 --threads 10 --out cat2_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat2_hqp_q2_0.05 --pca 20 --threads 10 --out cat2_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat2_hqp_q3_0.05 --pca 20 --threads 10 --out cat2_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat2_hqp_q4_0.05 --pca 20 --threads 10 --out cat2_hqp_q4_0.05


/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --grm cat3_hqp_q1_0.05 --pca 20 --threads 10 --out cat3_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat3_hqp_q2_0.05 --pca 20 --threads 10 --out cat3_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat3_hqp_q3_0.05 --pca 20 --threads 10 --out cat3_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat3_hqp_q4_0.05 --pca 20 --threads 10 --out cat3_hqp_q4_0.05


/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat4_hqp_q1_0.05 --pca 20 --threads 10 --out cat4_hqp_q1_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat4_hqp_q2_0.05 --pca 20 --threads 10 --out cat4_hqp_q2_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat4_hqp_q3_0.05 --pca 20 --threads 10 --out cat4_hqp_q3_0.05
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --grm cat4_hqp_q4_0.05 --pca 20 --threads 10 --out cat4_hqp_q4_0.05





