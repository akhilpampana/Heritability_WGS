#####################################################################################################################################################
#                                      Subset the datset to variants present in each loci for 41 tophits                                            #
#####################################################################################################################################################
#                            Create a file containing all variants which are in +/- 500kb region of sentinel snp                                    #  
#####################################################################################################################################################
#This extraction was done in R
# cd into the folder that containing sentinel snp with variants in +/- 500kb region
module load R
cd eqtl_colocquail/

#### Code to run in R
require(data.table)
files = list.files()
files = files[grep("500kb",files)]
df = data.frame()
for(file in files){
  tmp = fread(paste0(file))
  tmp = tmp[which(tmp$p.value < 0.05),]
  df = rbind(df,tmp)
}

df = df[!duplicated(df$SNPID),]
df$SNP = paste0(df$CHR,":",df$POS,":",df$Allele1,":",df$Allele2)
df1 = df[,c("SNPID","SNP","p.value")]
write.table(df1,file="Sentinel_snps_500kb_region_combined_02102023.tsv",row.names=F,col.names=F,sep="\t",dec=".",quote=F)

#####################################################################################################################################################
#                               Subset the Original genotypes to these variants and create a raw file using plink                                   #
#####################################################################################################################################################

module load PLINK/1.90-foss-2016a
for i in 1 4 8 12
do
plink --bfile ../heritability/plink_format/original/rsid/freeze10.14k.chr${i}.0.0001    --extract eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023_v1.csv --make-bed --out chr${i}_var_500kb
done

#cat eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023.csv | sed 's|,|\t|g' > eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023_v1.csv
#cat snplist.snplist | sort | uniq -d > duplicated_snps.snplist
#plink --bfile ../heritability/plink_format/original/rsid//freeze10.14k.chr${i}.0.0001 --exclude duplicated_snps.snplist --make-bed --out freeze10.14k.chr${i}.0.0001
#plink --bfile freeze10.14k.chr${i}.0.0001 --clump gwas_variants_12162022.csv --clump-kb 500 --clump-p1 5e-7 --clump-p2 0.05 --clump-r2 0.60 --out group3_chr${i} 


#####################################################################################################################################################
#                                                                 LD Matrix generation                                                              #
#####################################################################################################################################################
# Alternatively can use cor function in r based on .raw file generated in plink

for i in 1 4 8 12
do
plink --bfile chr${i}_var_500kb --r2 square --out chr${i}_var_500kb
done



#####################################################################################################################################################
#                               Use these correlation matrix for finemapping using  PAINTOR                                                         #
#####################################################################################################################################################

module load glibc/2.33-GCCcore-10.2.0

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/coloc/fine_mapping













