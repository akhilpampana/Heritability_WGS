########################################
# Script used to run GCTA/plink/bcftools for the analyses in the GREML-LDMS WGS 2022 paper
# It contains the most important commands for performing GREML-LDMS from WGS
# The Online Methods describe with greater details all the QC steps performed
# Additional R scripts are also available in this git
########################################

#### Variables
ncpu=12 
sample_list="/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/phenotypes/NTproBNP_14k.tsv" #List of samples with phenotypic information available

###### Create BED files and perform first filtering on SVM variants and other QC metrics


for i in  {1..22} ; do
    	bcftools view -m2 -M2  -Ou  -i'FILTER="PASS"' --threads "$ncpu" -S "$sample_list" --force-samples /data/project/Arora_lab/akhil/TOPMED/COMPLETE/GENOTYPES/freeze.10a.chr"$i".pass_only.gtonly.minDP10.vcf.gz | 
	bcftools  annotate --threads "$ncpu" -Ob -I +'%CHROM:%POS:%REF:%ALT' > /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/freeze10.14k.chr"$i".pass.bcf ;
	plink --bcf /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/heritability/freeze10.14k.chr"$i".pass.bcf --maf 0.0001 --allow-extra-chr --keep-allele-order --geno 0.05 --hwe 0.000001 --threads ${ncpu} -make-bed --out freeze10.14k.chr"$i".0.0001 ;
done

### Allele Frequency Calculations
for i in {1..22}; do
plink2 --bfile freeze10.14k.chr${i}.0.0001 --freq --out freeze10.14k.chr${i}
done


### Subset Variants to 4 bins [0.05],[0.01,0.05),[0.001,0.01),[0.0001,0.001)
#module load R => R
require(data.table)
category1 = data.frame() ## [0.0001,0.001)
category2 = data.frame() ## [0.001,0.01)
category3 = data.frame() ## [0.01,0.05)
category4 = data.frame() ## [0.05]
for(i in 1:22) {
	tmp = fread(paste0("freeze10.14k.chr",i,".afreq"))
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

save(category1,category2,category3,category4,file="allele_frequence_cutoff_4cat_09132022.rda")

write.table(category1,file="category1_0.0001_0.001.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category2,file="category2_0.001_0.01.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category3,file="category3_0.01_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)
write.table(category4,file="category4_0.05.csv",row.names=F,col.names=T,sep="\t",dec=".",quote=F)

## Category 1 ##merge not working 

for i in  {1..22} ; do 
plink2 --bfile plink_format/freeze10.14k.chr${i}.0.0001 --extract subset_for_h2_calc/category1_0.0001_0.001/category1_0.0001_0.001.csv --make-bed --out subset_for_h2_calc/category1_0.0001_0.001/category1_chr${i}  ; 
done

module load PLINK
cd subset_for_h2_calc/category1_0.0001_0.001/
for i in {1..22}; do
plink2 --bfile category1_chr${i} --set-all-var-ids @:#:'$r':'$a'  --new-id-max-allele-len 1000 --max-alleles 2 --make-bed --out cat1_chr${i}
done

ls -l | grep cat1 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile cat1_chr1 --merge-list merge --make-bed --out cat1_chr_all

### Prunning to generate high quality variants
for i in {1..22}; do
plink --bfile cat1_chr${i} --indep-pairwise 50 5 0.1 --out cat1_${i}
done

## subset to prunned variants
for i in {1..22}; do
plink --bfile cat1_chr${i} --extract cat1_${i}.prune.in --make-bed --out cat1_chr${i}_hqp
done


ls -l | grep hqp | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile cat1_chr1_hqp --merge-list merge --make-bed --out cat1_chr_all_hqp



## Category 2

for i in  {1..22} ; do 
plink2 --bfile plink_format/freeze10.14k.chr${i}.0.0001 --extract subset_for_h2_calc/category2_0.001_0.01/category2_0.001_0.01.csv --make-bed --out subset_for_h2_calc/category2_0.001_0.01/category2_chr${i}  ; 
done

module load PLINK
cd subset_for_h2_calc/category2_0.001_0.01/
for i in {1..22}; do
plink2 --bfile category2_chr${i} --set-all-var-ids @:#:'$r':'$a'  --new-id-max-allele-len 1000 --max-alleles 2 --make-bed --out cat2_chr${i}
done

ls -l | grep cat2 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile cat2_chr1 --merge-list merge --make-bed --out cat2_chr_all

### Prunning to generate high quality variants
plink --bfile cat2_chr_all --indep-pairwise 50 5 0.1 --out cat2

## subset to prunned variants
plink --bfile cat2_chr_all --extract cat2.prune.in --make-bed --out cat2_chr_all_hqp

### Category 3	
for i in  {1..22} ; do 
plink2 --bfile plink_format/freeze10.14k.chr${i}.0.0001 --extract subset_for_h2_calc/category3_0.01_0.05/category3_0.01_0.05.csv --make-bed --out subset_for_h2_calc/category3_0.01_0.05/category3_chr${i}  ; 
done

module load PLINK
cd subset_for_h2_calc/category3_0.01_0.05/
for i in {1..22}; do
plink2 --bfile category3_chr${i} --set-all-var-ids @:#:'$r':'$a'  --new-id-max-allele-len 1000 --max-alleles 2 --make-bed --out cat3_chr${i}
done

ls -l | grep cat3 | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile cat3_chr1 --merge-list merge --make-bed --out cat3_chr_all

### Prunning to generate high quality variants
plink --bfile cat3_chr_all --indep-pairwise 50 5 0.1 --out cat3

## subset to prunned variants
plink --bfile cat3_chr_all --extract cat3.prune.in --make-bed --out cat3_chr_all_hqp

### Category 4
for i in  {1..22} ; do 
plink2 --bfile plink_format/freeze10.14k.chr${i}.0.0001 --extract subset_for_h2_calc/category4_0.05/category4_0.05.csv --make-bed --out subset_for_h2_calc/category4_0.05/category4_chr${i}  ; 
done

cd subset_for_h2_calc/category4_0.05/
for i in {1..22}; do
plink2 --bfile category4_chr${i} --set-all-var-ids @:#:'$r':'$a'  --new-id-max-allele-len 1000 --max-alleles 2 --make-bed --out cat4_chr${i}
done

ls -l | grep bed | awk ' { print $9 } ' | sed 's|.bed||g' > merge
module load PLINK/1.90-foss-2016a
plink --bfile cat4_chr1 --merge-list merge --make-bed --out cat4_chr_all

### Prunning to generate high quality variants
plink --bfile cat4_chr_all --indep-pairwise 50 5 0.1 --out cat4


## subset to prunned variants
plink --bfile cat4_chr_all --extract cat4.prune.in --make-bed --out cat4_chr_all_hqp

### LD score calculation
/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/softwares/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat4_chr_all_hqp  

###### Merge BEDs ###### skip already done above
#${list_beds} contains the list of the autosomes (excepted chr 1) for merging

plink \
	--bfile ${BED_file_merged} \
	--merge-list ${list_beds} \
	--make-bed \
	--maf 0.0001 \
	--geno 0.05  \
	--hwe 0.000001 \
	--mind 0.05 \
	--out ${BED_file_merged_QC} \
	--threads ${ncpu}


###### PRS ###### Not using skip
#Simple PRS as a cohort QC check
#${scoreSNP} contains the effect sizes of the SNPs selected to construct the PRS

#plink \
#	--bfile ${BED_file_merged} \
#	--threads ${ncpu}  \
#	--score ${scoreSNP} \
#	--out P


###### Computing IBD ######
#Computing IBD segments. For computational reasons, it is recomended to prune/select common SNPs for faster runtime.
module load KING/2.1.2-foss-2016a

king -b cat4_chr_all_hqp.bed --ibdseg --prefix cat4_chr_all_hqp --cpus 4 --seglength 3
king -b ../category3_0.01_0.05/cat3_chr_all_hqp.bed --ibdseg --prefix ../category3_0.01_0.05/cat3_chr_all_hqp --cpus 4 --seglength 3

king -b ${BED_file_merged_QC}.bed  \
	--ibdseg \
	--prefix ${BED_file_merged_QC_out_name} \
	--cpus ${ncpu} \
	--seglength 3


###### GRM ######
#Construct a GRM, extracting a list of variants (see 01_LD_bins.R to create LD bins from a bed of a MAF range)
#GRM constructed by part (one part computed per script on an array job)

gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat3_chr_all_hqp --make-grm-alg 1 --thread-num 4 --out test
gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile cat4_chr_all_hqp --make-grm-alg 1 --thread-num 4 --out test



#Merge all GRM parts together

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


###### Create a file containing multiple GRMs in a directory (need full path) ######

for i in *.grm.bin ; do readlink -f "$i"  | cut -d'.' -f1-2 >>  ${mgrm_file_path}; done


###### Filter variants wihin a MAF range ######

plink \
	--bfile ${BED_file_merged_QC} \
	--maf 0.0001 \
	--max-maf 0.001 \
	--make-bed \
	--threads ${ncpu} \
	--out ${BED_file_MAF}


###### SNP pruning ######
#Select independant SNPs based on a different thresholds from a set of variants
#Indep-pairwise parameters can be modified depending if the pruning is performed on common or rare variants 

i={1..22}
plink \
	--bfile ${BED_file_merged_QC} \
	--chr "$i" \
	--extract ${list_variants_bin} \
	--indep-pairwise 50 5 0.1 \
	--out ${out_indep_var}_chr"$i" \
	--threads ${ncpu}


###### Compute PCs ######
#Compute PC from GRM (GCTA) or from BED (plink2)

GCTA \
	--grm ${GRM_out} \
	--threads ${ncpu} \
	--pca 20 \
	--out ${PCA_out}


plink2 \
	--bfile ${BED_file_merged_QC} \
	--extract ${list_variants_bin} \
	--pca 20 approx \
	--out ${PCA_out} \
	--thread-num ${ncpu}


###### REML ######
#Run unconstrained REML from multiple GRMs (for GREML-LDMS), fitting PCs as quantitative covariates
#REML-no-lrt doesn't calculate the reduced model

GCTA \
	--reml \
	--mgrm ${mgrm_file_path} \
	--reml-no-constrain \
	--pheno ${phenotype_file} \
	--out ${REML_output_file} \
	--thread-num ${ncpu} \
	--qcovar ${PCA_out} \
	--reml-no-lrt





