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
setwd("eqtl_colocquail/")
require(data.table)
files = list.files()
files = files[grep("500kb",files)]
#files = files[!(files %in% c("rs198389_500kb.z","Sentinel_snps_500kb_region_combined_02102023_v1.csv","Sentinel_snps_500kb_region_combined_02102023.csv"))]

df = data.frame()
for(file in files){
  var = gsub(".txt","",file)
  tmp = fread(paste0(file))
  tmp = tmp[which(tmp$p.value < 0.05),]
  tmp$flip = ifelse(tmp$AF_Allele2 > 0.50,1,0)
  tmp = tmp[,c(3,1,2,4,5,7,10,11,16)]
  colnames(tmp) = c("rsid","chromosome","position","allele1","allele2","maf","beta","se","flip")
  write.table(tmp,file=paste0("/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/coloc/",var,".z"),row.names=F,col.names=T,sep=" ",dec=".",quote=F)
}

#df = df[!duplicated(df$SNPID),]
#df$SNP = paste0(df$CHR,":",df$POS,":",df$Allele1,":",df$Allele2)
#df1 = df[,c("SNPID","SNP","p.value")]
#write.table(df1,file="Sentinel_snps_500kb_region_combined_02102023.tsv",row.names=F,col.names=F,sep="\t",dec=".",quote=F)

for(file in files){
  #var = gsub(".txt","",file)
  tmp = fread(paste0(file))
  #tmp = tmp[which(tmp$p.value < 0.05),]
  #tmp$flip = ifelse(tmp$AF_Allele2 > 0.50,1,0)
  
  tmp1 = tmp[which(tmp$flip == 1),]
  tmp1$Allele1 = tmp1$allele2
  tmp1$Allele2 = tmp1$allele1
  tmp1$final_maf = 1 - tmp1$maf
  tmp1$BETA = abs(tmp1$beta)
  
  tmp2 = tmp[which(tmp$flip == 0),]
  tmp2$Allele1 = tmp2$allele1
  tmp2$Allele2 = tmp2$allele2
  tmp2$final_maf = tmp2$maf
  tmp2$BETA = tmp2$beta
  
  tmp = rbind(tmp1,tmp2)
  tmp = tmp[,c(1,2,3,10,11,12,13,8,9)]
  colnames(tmp) = c("rsid","chromosome","position","allele1","allele2","maf","beta","se","flip")
  write.table(tmp,file=paste0(file,".txt"),row.names=F,col.names=T,sep=" ",dec=".",quote=F)
}



for(file in files){
  var = gsub(".txt","",file)
  tmp = fread(paste0(file))
  write.table(tmp,file=paste0(var,".z"),row.names=F,col.names=T,sep=" ",dec=".",quote=F)
}





#####################################################################################################################################################
#                               Subset the Original genotypes to these variants and create a raw file using plink                                   #
#####################################################################################################################################################

### combined subset
module load PLINK/1.90-foss-2016a
for i in 1 4 8 12
do
plink --bfile ../heritability/plink_format/original/rsid/freeze10.14k.chr${i}.0.0001    --extract eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023_v1.csv --make-bed --out chr${i}_var_500kb
done

### subset per loci
ls -l | grep z | awk ' { print $9 } ' | tr '\n' '\t' > rsid
ls -l | grep z | awk ' { print $9 } '  > rsid

for i in rs1009591_500kb.z       rs1023252_500kb.z       rs10550903_500kb.z      rs10689649_500kb.z      rs10858903_500kb.z      rs11105282_500kb.z      rs111478946_500kb.z     rs11555351_500kb.z      rs1208984_500kb.z       rs12402363_500kb.z      rs12411044_500kb.z      rs13107325_500kb.z      rs1318408_500kb.z       rs142005893_500kb.z     rs143149865_500kb.z     rs145485557_500kb.z     rs148333765_500kb.z     rs198383_500kb.z        rs198389_500kb.z        rs198392_500kb.z        rs202088699_500kb.z     rs2066462_500kb.z       rs2273286_500kb.z       rs28455075_500kb.z      rs34710782_500kb.z      rs34954501_500kb.z      rs3753584_500kb.z       rs41275462_500kb.z      rs41300100_500kb.z      rs4845875_500kb.z       rs4845876_500kb.z       rs4845881_500kb.z       rs4846063_500kb.z       rs55714388_500kb.z      rs61757273_500kb.z      rs72640280_500kb.z      rs72640281_500kb.z      rs7299091_500kb.z       rs7314459_500kb.z       rs79593079_500kb.z      rs9727993_500kb.z
do
plink --bfile /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/conditional_analysis/snps_for_conditioning/tophits_0.05/heritability/plink_format/original/rsid/freeze10.14k.chr1.0.0001    --extract ${i} --make-bed --out fine_mapping/${i} || plink --bfile ../heritability/plink_format/original/rsid/freeze10.14k.chr4.0.0001    --extract ${i} --make-bed --out fine_mapping/${i} || plink --bfile ../heritability/plink_format/original/rsid/freeze10.14k.chr8.0.0001    --extract ${i} --make-bed --out fine_mapping/${i} || plink --bfile ../heritability/plink_format/original/rsid/freeze10.14k.chr12.0.0001    --extract ${i} --make-bed --out fine_mapping/${i}
done

for i in rs1009591_500kb.z       rs1023252_500kb.z       rs10550903_500kb.z      rs10689649_500kb.z      rs10858903_500kb.z      rs11105282_500kb.z      rs111478946_500kb.z     rs11555351_500kb.z      rs1208984_500kb.z       rs12402363_500kb.z      rs12411044_500kb.z      rs13107325_500kb.z      rs1318408_500kb.z       rs142005893_500kb.z     rs143149865_500kb.z     rs145485557_500kb.z     rs148333765_500kb.z     rs198383_500kb.z        rs198389_500kb.z        rs198392_500kb.z        rs202088699_500kb.z     rs2066462_500kb.z       rs2273286_500kb.z       rs28455075_500kb.z      rs34710782_500kb.z      rs34954501_500kb.z      rs3753584_500kb.z       rs41275462_500kb.z      rs41300100_500kb.z      rs4845875_500kb.z       rs4845876_500kb.z       rs4845881_500kb.z       rs4846063_500kb.z       rs55714388_500kb.z      rs61757273_500kb.z      rs72640280_500kb.z      rs72640281_500kb.z      rs7299091_500kb.z       rs7314459_500kb.z       rs79593079_500kb.z      rs9727993_500kb.z
do
plink --bfile ${i} --recode A --out ${i}
done


### Final Version
for i in rs1009591_500kb.z       rs1023252_500kb.z       rs10550903_500kb.z      rs10689649_500kb.z      rs10858903_500kb.z      rs11105282_500kb.z      rs111478946_500kb.z     rs11555351_500kb.z      rs1208984_500kb.z       rs12402363_500kb.z      rs12411044_500kb.z      rs13107325_500kb.z      rs1318408_500kb.z       rs142005893_500kb.z     rs143149865_500kb.z     rs145485557_500kb.z     rs148333765_500kb.z     rs198383_500kb.z        rs198389_500kb.z        rs198392_500kb.z        rs202088699_500kb.z     rs2066462_500kb.z       rs2273286_500kb.z       rs28455075_500kb.z      rs34710782_500kb.z      rs34954501_500kb.z      rs3753584_500kb.z       rs41275462_500kb.z      rs41300100_500kb.z      rs4845875_500kb.z       rs4845876_500kb.z       rs4845881_500kb.z       rs4846063_500kb.z       rs55714388_500kb.z      rs61757273_500kb.z      rs72640280_500kb.z      rs72640281_500kb.z      rs7299091_500kb.z       rs7314459_500kb.z       rs79593079_500kb.z      rs9727993_500kb.z
do
plink --bfile /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/conditional_analysis/snps_for_conditioning/tophits_0.05/chr1.pass_only_all   --update-name  eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023_v1.csv 1 2 --extract fine_mapping/${i} --make-bed --out fine_mapping/${i} || plink --bfile /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/conditional_analysis/snps_for_conditioning/tophits_0.05/chr4.pass_only_all --update-name  eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023_v1.csv 1 2  --extract fine_mapping/${i} --make-bed --out fine_mapping/${i} || plink --bfile /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/conditional_analysis/snps_for_conditioning/tophits_0.05/chr8.pass_only_all --update-name  eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023_v1.csv 1 2  --extract fine_mapping/${i} --make-bed --out fine_mapping/${i} || plink --bfile /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/conditional_analysis/snps_for_conditioning/tophits_0.05/chr12.pass_only_all  --update-name  eqtl_colocquail/Sentinel_snps_500kb_region_combined_02102023_v1.csv 1 2   --extract fine_mapping/${i} --make-bed --out fine_mapping/${i}
done

#####################################################################################################################################################
#                                                                 LD Matrix generation                                                              #
#####################################################################################################################################################
# Alternatively can use cor function in r based on .raw file generated in plink

for i in rs1009591_500kb.z       rs1023252_500kb.z       rs10550903_500kb.z      rs10689649_500kb.z      rs10858903_500kb.z      rs11105282_500kb.z      rs111478946_500kb.z     rs11555351_500kb.z      rs1208984_500kb.z       rs12402363_500kb.z      rs12411044_500kb.z      rs13107325_500kb.z      rs1318408_500kb.z       rs142005893_500kb.z     rs143149865_500kb.z     rs145485557_500kb.z     rs148333765_500kb.z     rs198383_500kb.z        rs198389_500kb.z        rs198392_500kb.z        rs202088699_500kb.z     rs2066462_500kb.z       rs2273286_500kb.z       rs28455075_500kb.z      rs34710782_500kb.z      rs34954501_500kb.z      rs3753584_500kb.z       rs41275462_500kb.z      rs41300100_500kb.z      rs4845875_500kb.z       rs4845876_500kb.z       rs4845881_500kb.z       rs4846063_500kb.z       rs55714388_500kb.z      rs61757273_500kb.z      rs72640280_500kb.z      rs72640281_500kb.z      rs7299091_500kb.z       rs7314459_500kb.z       rs79593079_500kb.z      rs9727993_500kb.z
do
plink --bfile ${i} --r square spaces --out ${i}
done

### Correlation matrix based on R using ldmat function in hibayes ~ testing
install.packages("hibayes")

require(data.table)
library("hibayes")
files = list.files()
files = files[grep(".raw",files)]
for(file in files){
  
  ############# LD matrix generation
  data = as.data.frame(fread(file))
  colnames(data) = gsub("[/A-Z]+$","",colnames(data))
  colnames(data) = gsub("_","",colnames(data))
  data1 = data[,c(7:ncol(data))]
  data1 = as.matrix(data1)
  data2 = ldmat(data1)
  rownames(data2) = colnames(data1)
  colnames(data2) = colnames(data1)
  
  ############# Subset variants to rsids present in ldmatrix
  var = gsub(".raw",".z",file)
  data = fread(paste0(var))
  data = data[which(data$rsid %in% colnames(data2)),]
  
  }


### Validation of equal number of lines between ld file and bim file
for i in rs1009591_500kb.z       rs1023252_500kb.z       rs10550903_500kb.z      rs10689649_500kb.z      rs10858903_500kb.z      rs11105282_500kb.z      rs111478946_500kb.z     rs11555351_500kb.z      rs1208984_500kb.z       rs12402363_500kb.z      rs12411044_500kb.z      rs13107325_500kb.z      rs1318408_500kb.z       rs142005893_500kb.z     rs143149865_500kb.z     rs145485557_500kb.z     rs148333765_500kb.z     rs198383_500kb.z        rs198389_500kb.z        rs198392_500kb.z        rs202088699_500kb.z     rs2066462_500kb.z       rs2273286_500kb.z       rs28455075_500kb.z      rs34710782_500kb.z      rs34954501_500kb.z      rs3753584_500kb.z       rs41275462_500kb.z      rs41300100_500kb.z      rs4845875_500kb.z       rs4845876_500kb.z       rs4845881_500kb.z       rs4846063_500kb.z       rs55714388_500kb.z      rs61757273_500kb.z      rs72640280_500kb.z      rs72640281_500kb.z      rs7299091_500kb.z       rs7314459_500kb.z       rs79593079_500kb.z      rs9727993_500kb.z
do
if [ $(wc -l ${i}) -eq $w(wc -l ${i}.ld) ]; then echo "Warning: No Match!"; fi 
done

#####################################################################################################################################################
#                                               Fine Mapping using FINEMAP tools                                                                    #
#####################################################################################################################################################

/data/project/Arora_lab/akhil/SOFTWARES/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --in-files master1  --sss --flip-beta

#####################################################################################################################################################
#                               Coverting Variants position to HG19 to match with PAINTOR annotations                                               #
#####################################################################################################################################################
require(data.table)
files = list.files()
files = gsub("z.ld","z.z",files)
files - unique(files)
lifted = fread("../../../../gwas/lifted_grch38_grch37/Lifted_overall.bed",header=F)
for(file in files) {
  var = gsub("500kb.z.z","grch37",file)
  tmp = fread(paste0(file))
  tmp = merge(tmp,lifted,by.x=c("rsid"),by.y=c("V4"))
  tmp = tmp[,c("chromosome","V2","rsid","beta","se")]
  tmp$z = tmp$beta/tmp$se
  tmp = tmp[,c(1,2,3,6)]
  colnames(tmp) = c("CHR","POS","RSID","ZSCORE")
  write.table(tmp,file=paste0(var,".txt"),row.names=F,col.names=T,sep=" ",quote=F,dec=".")
}
#####################################################################################################################################################
#                               Use these correlation matrix for finemapping using  PAINTOR                                                         #
#####################################################################################################################################################
## Annotations for this dataset
module load glibc/2.33-GCCcore-10.2.0
sbatch --array=0-40 script_final_invitae_v1.sh






../../../../../../../../SOFTWARES/bgen.tgz/bin/bgenix

cat rs198389_500kb.txt | awk ' { print $3 "\t" $1 "\t" $2 "\t" $4 "\t" $5 "\t" $7 "\t" $10 "\t" $11 } ' >  rs198389_500kb.z

/data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/coloc/fine_mapping

/data/project/Arora_lab/akhil/SOFTWARES/PAINTOR_V3.0/PAINTOR

/data/project/Arora_lab/akhil/SOFTWARES/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 

PAINTOR_test/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 --in-files /data/project/Arora_lab/akhil/TOPMED/BNP/NTproBNP/NTproBNP_14k/gwas/coloc/fine_mapping/master_v1 --write-bcor --read-only-bgen

/data/project/Arora_lab/akhil/SOFTWARES/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --config --in-files master --rsids rs198389

/data/project/Arora_lab/akhil/ANNOTATIONS/Functional_Annotations


# Tissues Annotations removed
/data/project/Arora_lab/akhil/ANNOTATIONS/Functional_Annotations/RoadMap_Assayed_NarrowPeak/.E085-H3K36me3.narrowPeak.Fetal_Intestine_Small.yIuNle
/data/project/Arora_lab/akhil/ANNOTATIONS/Functional_Annotations/Hnisz_Cell2013_SuperEnhancer/mmc2.pdf
