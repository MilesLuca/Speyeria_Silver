# Speyeria_Silver
A collection of code used for our study detailing the genomic basis for the silvering polymorphism in Speyeria mormonia.

# Table of contents

# I. Figure 2 - GWAS, Fst, PCAs and genotype plots

### 1. Alignment and genotype calling

### 2. Alignment QC (depth measurments)

### 3. VCF filtering

### 4. Genome-wide association using GEMMA

### 5. Importing p-values and chromosome positions in R and plotting

### 6. Code for calculating Fst on morphs by populations

### 7. Plotting Fst and -log10 P from GWAS on chromosome region of interest in R

### 8. Principle component analysis for genome-wide and GWAS specific regions

### 9. Plot PCAs in R

### 10. Genotype plot for SNP illutsration at the GWAS interval

# II. Figure 4 - Haplotype based statistics, sweeps and introgression analyses (Fst, Dxy, Tajima's D, CLR, TWISST, Fd)

### 1. Vcf filtering and recoding

### 2. Fst recalculated on phased haplotypes

### 3. Dxy, Pi, and Tajima's calculated on phased haplotypes

### 4. Sweep analysis using SweeD

### 5. Topology Weighting by Iterative Subsampling (TWISST)

### 6. Proportion of introgression in sliding windows (fd)

### 7. Plot all haplotype based statistics in R


-----------------------------------------------------------------------------------------------------------------------------------

# I. Figure 2 - GWAS, Fst, PCAs and genotype plots

## 1. [Link Text](###Alignment and genotype calling)

The following code is used to generate genome indexes, align sequencing data, mark sequencing duplicates, call genotypes and produce vcf files for downstream analyses.

**First Index the genome for gatk, bwa and samtools:**

```
gatk CreateSequenceDictionary -R Smor.assembly.v1.0.unmasked.fasta
bwa index Smor.assembly.v1.0.unmasked.fasta
samtools faidx Smor.assembly.v1.0.unmasked.fasta
```

Alignemts and genotype calling are performed using a job array.

**Save file called “jobs” with following 75 IDs in (corresponding to each sample):**

```
nano jobs

Shyd_WA_05 Shyd_WA_06 Shyd_WA_08 Shyd_WA_09 Shyd_WA_10 Smor_WA_01 Smor_WA_02 Smor_WA_07 Smor_WA_11 Smor_WA_13 Smor_WA_14 Smor_WA_15 Smor_WA_16 Smor_WA_17 Smor_WA_18 Smor_WA_19 Smor_WA_20 Smor_WA_21 Smor_WA_22 Smor_WA_23 Smor_WA_24 Smor_WA_25 Smor_WA_26 Smor_WA_27 Smor_WA_28 Smor_WA_29 Smor_WA_30 Smor_WA_31 Smor_WA_32 Smor_WA_33 Smor_WA_34 Smor_WA_35 Smor_WA_37 Smor_WA_38 Smor_WA_39 Smor_WA_40 Smor_WA_41 Smor_WA_42 Smor_WA_43 Smor_WA_44 Smor_WA_45 Smor_WA_46 Smor_WA_47 Smor_WA_48 NVJ_51_Smor_B NVJ_52_Smor_S NVJ_53_Smor_B NVJ_54_Smor_B NVJ_55_Smor_S NVJ_56_Smor_S NVJ_57_Smor_B NVJ_58_Smor_B NVJ_59_Smor_B NVJ_60_Smor_S NVJ_61_Smor_B NVJ_62_Smor_S NVJ_63_Smor_S NVJ_64_Smor_B NVJ_65_Smor_S NVJ_66_Smor_S NVJ_67_Smor_S NVJ_68_Smor_B NVJ_69_Smor_B NVJ_70_Smor_B NVJ_71_Smor_B NVJ_72_Smor_S NVJ_73_Smor_B NVJ_74_Smor_B NVJ_75_Smor_S NVJ_76_Smor_B NVJ_77_Smor_S NVJ_78_Smor_S NVJ_79_Smor_S NVJ_80_Smor_S NVJ_81_Smor_S
```

**Make a directory to store the alignment statistics in**

```
mkdir alnstats
```

**Generate a script to run job as array (75 samples total)**

This script is run in two parts

The first is a script called: **makeGVCF.array.sh**

Part A (**makeGVCF.array.sh**): Call variants on each individual against the genome:
1. Aligns using bwa-mem (preferred aligner for genotyping)
2. Marks duplicate reads (using a gatk-instance of picard's MarkDuplicates)
3. Ensures that sensible names are assigned to each individual with "AddOrReplaceReadGroups" (this is important when merging the files in part B).
4. Produces the gvcf (genotype variant call file), one per individual. This file will score every site in the genome as either non-variant or variant.

The second is a script called: **genotypeGVCFs.sh**

Part B (**genotypeGVCFs.sh**): Merge gVCFs, call genotypes and produce the VCF, and then filter the VCF. 
1. Runs the gatk tool CombineGVCFs to merge all the gVCFs into one (very large) gvcf. This file includes a separate entry for every site in every individual so it is upsettingly large, but it runs quickly
2. Runs the gatk tool GenotypeGVCFs on the combined gVCF. This is the 'genotype calling' step. This program goes through each site in the genome, looks at the variants from the gVCF for each individual, and calls the genotype at each site. The output is one row per variant, one column per individual (ie, a VCF).
3-6, Filtering. Based on previosuly used filtering parameters optimised for Heliconius butterfly studies. These steps will create two final vcfs, one that contains just SNPs and one that contains just indels.

**Run part A**

```
makeGVCF.array.sh

#!/bin/bash
#SBATCH -J array_gvcf
#SBATCH -o array_gvcf.out
#SBATCH -e array_gvcf.err
#SBATCH -p defq
#SBATCH -n 32
#SBATCH -t 3-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lucalivraghi@gwu.edu
#SBATCH --array=0-74%10

module load samtools
module load bwa
module load gatk

names=($(cat jobs))

echo ${names[${SLURM_ARRAY_TASK_ID}]}

## align. Commands below grab each fastq file pair and align it to Smor genome. tee command outputs alignment stats
bwa mem -t 32 /lustre/groups/martinlab/Smor_GWAS/genome/Smor.assembly.v1.0.unmasked.fasta resequence_data/fastqs/${names[${SLURM_ARRAY_TASK_ID}]}.R1.fastq.gz resequence_data/fastqs/${names[${SLURM_ARRAY_TASK_ID}]}.R2.fastq.gz | \
samtools sort -@32 - | \
tee >(samtools flagstat - > alnstats/${names[${SLURM_ARRAY_TASK_ID}]}.stats.out) > ${names[${SLURM_ARRAY_TASK_ID}]}.bam

## mark the dups. Marks duplicates for future removal
gatk —java-options "-Xmx40G -XX:+UseParallelGC -XX:ParallelGCThreads=32" MarkDuplicates \
-I ${names[${SLURM_ARRAY_TASK_ID}]}.bam \
-M ${names[${SLURM_ARRAY_TASK_ID}]}.MarkDuplicates.out \
-O ${names[${SLURM_ARRAY_TASK_ID}]}.marked.bam

## rename read groups. Ensure that sensible names are assigned to each individual with "AddOrReplaceReadGroups" (this is important when merging the files)
gatk —java-options "-Xmx40G -XX:+UseParallelGC -XX:ParallelGCThreads=32" AddOrReplaceReadGroups \
-I ${names[${SLURM_ARRAY_TASK_ID}]}.marked.bam \
-O ${names[${SLURM_ARRAY_TASK_ID}]}.marked.RG.bam \
-RGID 4 \
-RGLB lib1 \
-RGPL illumina \
-RGPU unit1 \
-RGSM ${names[${SLURM_ARRAY_TASK_ID}]}

# make one gvcf per individual. This file will score every site in the genome as either non-variant or variant.
samtools index ${names[${SLURM_ARRAY_TASK_ID}]}.marked.RG.bam
gatk —java-options "-Xmx40G -XX:+UseParallelGC -XX:ParallelGCThreads=32" HaplotypeCaller \
--native-pair-hmm-threads 32 \
-R /lustre/groups/martinlab/Smor_GWAS/genome/Smor.assembly.v1.0.unmasked.fasta \
-I ${names[${SLURM_ARRAY_TASK_ID}]}.marked.RG.bam \
-O ${names[${SLURM_ARRAY_TASK_ID}]}.raw.g.vcf \
-ERC GVCF
```

**Run part B**

```
genotypeGVCFs.sh

#!/bin/bash
#SBATCH -J gatkGeGVCF
#SBATCH -o gatkGeGVCF.out
#SBATCH -e gatkGeGVCF.err
#SBATCH -p highThru
#SBATCH -n 8
#SBATCH -t 1-00:00:00
#SBATCH —mail-type=END
#SBATCH --mail-user=lucalivraghi@gwu.edu

module load samtools
module load gatk

#!/bin/bash
#SBATCH -J gatkGeGVCF
#SBATCH -o gatkGeGVCF.out
#SBATCH -e gatkGeGVCF.err
#SBATCH -p highThru
#SBATCH -n 8
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jjhanly@gwu.edu

module load samtools
module load gatk

### STEP 1 #### combine gVCFs into one file. This file includes a separate entry for every site in every individual so it is upsettingly large, but it runs quickly
gatk —java-options "-Xmx120G -XX:+UseParallelGC -XX:ParallelGCThreads=8" CombineGVCFs \
-R /lustre/groups/martinlab/Smor_GWAS/genome/Smor.assembly.v1.0.unmasked.fasta \
--variant Scor_WA_01.raw.g.vcf \
--variant Shyd_WA_05.raw.g.vcf \
--variant Shyd_WA_06.raw.g.vcf \
--variant Shyd_WA_08.raw.g.vcf \
--variant Shyd_WA_09.raw.g.vcf \
--variant Shyd_WA_10.raw.g.vcf \
--variant Smor_WA_01.raw.g.vcf \
--variant Smor_WA_02.raw.g.vcf \
--variant Smor_WA_07.raw.g.vcf \
--variant Smor_WA_11.raw.g.vcf \
--variant Smor_WA_13.raw.g.vcf \
--variant Smor_WA_14.raw.g.vcf \
--variant Smor_WA_15.raw.g.vcf \
--variant Smor_WA_16.raw.g.vcf \
--variant Smor_WA_17.raw.g.vcf \
--variant Smor_WA_18.raw.g.vcf \
--variant Smor_WA_19.raw.g.vcf \
--variant Smor_WA_20.raw.g.vcf \
--variant Smor_WA_21.raw.g.vcf \
--variant Smor_WA_22.raw.g.vcf \
--variant Smor_WA_23.raw.g.vcf \
--variant Smor_WA_24.raw.g.vcf \
--variant Smor_WA_25.raw.g.vcf \
--variant Smor_WA_26.raw.g.vcf \
--variant Smor_WA_27.raw.g.vcf \
--variant Smor_WA_28.raw.g.vcf \
--variant Smor_WA_29.raw.g.vcf \
--variant Smor_WA_30.raw.g.vcf \
--variant Smor_WA_31.raw.g.vcf \
--variant Smor_WA_32.raw.g.vcf \
--variant Smor_WA_33.raw.g.vcf \
--variant Smor_WA_34.raw.g.vcf \
--variant Smor_WA_35.raw.g.vcf \
--variant Smor_WA_37.raw.g.vcf \
--variant Smor_WA_38.raw.g.vcf \
--variant Smor_WA_39.raw.g.vcf \
--variant Smor_WA_40.raw.g.vcf \
--variant Smor_WA_41.raw.g.vcf \
--variant Smor_WA_42.raw.g.vcf \
--variant Smor_WA_43.raw.g.vcf \
--variant Smor_WA_44.raw.g.vcf \
--variant Smor_WA_45.raw.g.vcf \
--variant Smor_WA_46.raw.g.vcf \
--variant Smor_WA_47.raw.g.vcf \
--variant Smor_WA_48.raw.g.vcf \
--variant NVJ_51_Smor_B.raw.g.vcf \
--variant NVJ_52_Smor_S.raw.g.vcf \
--variant NVJ_53_Smor_B.raw.g.vcf \
--variant NVJ_54_Smor_B.raw.g.vcf \
--variant NVJ_55_Smor_S.raw.g.vcf \
--variant NVJ_56_Smor_S.raw.g.vcf \
--variant NVJ_57_Smor_B.raw.g.vcf \
--variant NVJ_58_Smor_B.raw.g.vcf \
--variant NVJ_59_Smor_B.raw.g.vcf \
--variant NVJ_60_Smor_S.raw.g.vcf \
--variant NVJ_61_Smor_B.raw.g.vcf \
--variant NVJ_62_Smor_S.raw.g.vcf \
--variant NVJ_63_Smor_S.raw.g.vcf \
--variant NVJ_64_Smor_B.raw.g.vcf \
--variant NVJ_65_Smor_S.raw.g.vcf \
--variant NVJ_66_Smor_S.raw.g.vcf \
--variant NVJ_67_Smor_S.raw.g.vcf \
--variant NVJ_68_Smor_B.raw.g.vcf \
--variant NVJ_69_Smor_B.raw.g.vcf \
--variant NVJ_70_Smor_B.raw.g.vcf \
--variant NVJ_71_Smor_B.raw.g.vcf \
--variant NVJ_72_Smor_S.raw.g.vcf \
--variant NVJ_73_Smor_B.raw.g.vcf \
--variant NVJ_74_Smor_B.raw.g.vcf \
--variant NVJ_75_Smor_S.raw.g.vcf \
--variant NVJ_76_Smor_B.raw.g.vcf \
--variant NVJ_77_Smor_S.raw.g.vcf \
--variant NVJ_78_Smor_S.raw.g.vcf \
--variant NVJ_79_Smor_S.raw.g.vcf \
--variant NVJ_80_Smor_S.raw.g.vcf \
--variant NVJ_81_Smor_S.raw.g.vcf \
-O Smor.g.vcf.gz

### STEP 2 #### Genotype that combined gVCF file. Run the gatk tool GenotypeGVCFs on the combined gVCF. This is the 'genotype calling' step. This program goes through each site in the genome, looks at the variants from the gVCF for each individual, and calls the genotype at each site. The output is one row per variant, one column per individual (ie, a VCF).
gatk —java-options "-Xmx115G -XX:+UseParallelGC -XX:ParallelGCThreads=8" GenotypeGVCFs \
-R /lustre/groups/martinlab/Smor_GWAS/genome/Smor.assembly.v1.0.unmasked.fasta \
-V Smor.g.vcf.gz \
-O Smor.vcf.gz

#### STEP 3 #### subset SNPs with SelectVariants
gatk —java-options "-Xmx120G -XX:+UseParallelGC -XX:ParallelGCThreads=8" SelectVariants \
-V Smor.vcf.gz \
-select-type SNP \
-O Smor.snps.vcf.gz

### STEP 4 #### subset Indels with SelectVariants
gatk —java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=64" SelectVariants \
-V Smor.vcf.gz \
-select-type INDEL \
-O Smor.indels.vcf.gz

#### STEP 5 #### Hard-filter SNPs on multiple expressions using VariantFiltration
gatk —java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=8" VariantFiltration \
-V Smor.snps.vcf.gz \
-filter "QD < 2.0" —filter-name "QD2" \
-filter "QUAL < 30.0" —filter-name "QUAL30" \
-filter "SOR > 3.0" —filter-name "SOR3" \
-filter "FS > 60.0" —filter-name "FS60" \
-filter "MQ < 40.0" —filter-name "MQ40" \
-filter "MQRankSum < -12.5" —filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" —filter-name "ReadPosRankSum-8" \
-O Smor.snps_filtered.vcf.gz

### STEP 6 #### Hard-filter Indels on multiple expressions using VariantFiltration
gatk —java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=64" VariantFiltration \
-V Smor.indels.vcf.gz \
-filter "QD < 2.0" —filter-name "QD2" \
-filter "QUAL < 30.0" —filter-name "QUAL30" \
-filter "FS > 200.0" —filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" —filter-name "ReadPosRankSum-20" \
-O Smor.indels_filtered.vcf.gz
```

## 2. Alignment QC (depth measurments)

**Checking alignments for overall depth using mosdepth**

```
check_depth.sh

#!/bin/bash
#SBATCH -J mosdepth
#SBATCH -o mosdepth.out
#SBATCH -e mosdepth.err
#SBATCH -p nano
#SBATCH -n 32
#SBATCH -t 30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lucalivraghi@gwu.edu

for i in *.bam; do
mosdepth --threads 32 --by 500 --fast-mode ./mosdepth/depth_${i} ${i}
done

```

Check the resulting mosdepth.summary.txt for coverage by chromosome, and average coverage.

**For a slightly more accurate (but way slower method) of calculating genome-wide coverage you can use Bedtools genomecov as follows below.**

**First generate a sizes.genome file with all chromosomes and their sizes:**

```
samtools faidx genome.fasta
cut -f1,2 genome.fa.fai > sizes.genome
```

**Then run bedtools for calculating coverage:**

```
bedtools_genomcov.sh

#!/bin/bash
#SBATCH -J bedtools
#SBATCH -o bedtools.out
#SBATCH -e bedtools.err
#SBATCH -p tiny
#SBATCH -n 32
#SBATCH -t 1:30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lucalivraghi@gwu.edu

for i in *.bam; do
bedtools genomecov -ibam ${i} -g sizes.genome -d > ${i}_coverage_output.txt
done
```

The -d argument specifies to calculate true coverage for each position in the genome.

**This generates pretty big files, which you can then awk to get an average for each position in the genome:**

```
awk '{sum += $3} END {print "Average coverage:", sum/NR}' SRR25297463_coverage_output.txt
```

Outputs from mosdepth and bedtools are very similar, so probably not worth running bedtools.

## 3. VCF filtering

From here on, we are using the SNP file generated by gatk, called "Smor.snps_filtered.vcf.gz"

**Check for missing data using vcftools:**

```
vcftools --gzvcf Smor.snps_filtered.vcf.gz --missing-indv --out data
```

There was no evidence of high missingness in any of the samples. 

**Filter with bcftools (removing minor allele freq of 2% and missing genotypes of 10%)**

```
filter_maf_miss.sh

#!/bin/bash
#SBATCH -J bcf_filt_maf_miss
#SBATCH -o bcf_filt_maf_miss.out
#SBATCH -e bcf_filt_maf_miss.err
#SBATCH -p nano
#SBATCH -n 32
#SBATCH -t 30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lucalivraghi@gwu.edu

module load bcftools

bcftools view --threads 32 -i 'MAF > 0.02' Smor.snps_filtered.vcf.gz -O z -o Smor.snps_filtered.maf.vcf.gz
bcftools view --threads 32 -i 'F_MISSING < 0.1' Smor.snps_filtered.maf.vcf.gz -O z -o Smor.snps_filtered.maf.geno.vcf.gz
```

**Next phase and impute missing genotypes using beagle**

```
phase_impute.sh

#!/bin/bash
#SBATCH -J beagle
#SBATCH -o beagle.out
#SBATCH -e beagle.err
#SBATCH -p tiny
#SBATCH -n 32
#SBATCH -t 4:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lucalivraghi@gwu.edu

java -Xmx50g -jar /CCAS/home/lucalivraghi/tools/beagle/beagle.05May22.33a.jar impute=TRUE window=1 overlap=0.1 nthreads=32 gt=Smor.snps_filtered.maf.geno.vcf.gz out=Smor.snps_filtered.maf.geno.beagled
```

In order to run the GWAS on only the correct population, I then removed the non-mormonia samples and split the vcfs into one containing only WA and one only NV samples. 
For illustration purposes, I'm only keeping the WA code here, but the process was performed the same way for NV samples.

**Purge non-variant sites and keep only bi-allelic sites**

```
bcf_filt_nonvar_biallele.sh

#!/bin/bash
#SBATCH -J bcf_filt
#SBATCH -o bcf_filt.out
#SBATCH -e bcf_filt.err
#SBATCH -p nano
#SBATCH -n 32
#SBATCH -t 30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lucalivraghi@gwu.edu

module load bcftools

#first remove non-biallelic sites
bcftools view --threads 32 -m2 -M2 -v snps Smor.snps_filtered.maf.geno.beagled.WA.vcf -O z -o Smor.snps_filtered.maf.geno.beagled.WA.biallele.vcf
#then remove non-variant sites
bcftools filter -e "MAC == 0" Smor.snps_filtered.maf.geno.beagled.WA.biallele.vcf -O z -o Smor.snps_filtered.maf.geno.beagled.WA.biallele.nonvar.vcf
```

## 4. Genome-wide association using GEMMA

**Generate plink files and run GWAS using GEMMA**

Phenotype file was created encoding each sample as either silver or unsilvered

```
nano PhenotypeFile.txt

IID	FID	Pheno
Smor_WA_01	Smor_WA_01	2
Smor_WA_02	Smor_WA_02	2
Smor_WA_07	Smor_WA_07	1
Smor_WA_11	Smor_WA_11	2
Smor_WA_13	Smor_WA_13	1
Smor_WA_14	Smor_WA_14	1
Smor_WA_15	Smor_WA_15	2
Smor_WA_16	Smor_WA_16	1
Smor_WA_17	Smor_WA_17	1
Smor_WA_18	Smor_WA_18	1
Smor_WA_19	Smor_WA_19	2
Smor_WA_20	Smor_WA_20	2
Smor_WA_21	Smor_WA_21	2
Smor_WA_22	Smor_WA_22	2
Smor_WA_23	Smor_WA_23	1
Smor_WA_24	Smor_WA_24	2
Smor_WA_25	Smor_WA_25	2
Smor_WA_26	Smor_WA_26	2
Smor_WA_27	Smor_WA_27	2
Smor_WA_28	Smor_WA_28	1
Smor_WA_29	Smor_WA_29	2
Smor_WA_30	Smor_WA_30	1
Smor_WA_31	Smor_WA_31	2
Smor_WA_32	Smor_WA_32	1
Smor_WA_33	Smor_WA_33	2
Smor_WA_34	Smor_WA_34	2
Smor_WA_35	Smor_WA_35	1
Smor_WA_37	Smor_WA_37	1
Smor_WA_38	Smor_WA_38	2
Smor_WA_39	Smor_WA_39	2
Smor_WA_40	Smor_WA_40	1
Smor_WA_41	Smor_WA_41	1
Smor_WA_42	Smor_WA_42	2
Smor_WA_43	Smor_WA_43	1
Smor_WA_44	Smor_WA_44	2
Smor_WA_45	Smor_WA_45	1
Smor_WA_46	Smor_WA_46	1
Smor_WA_47	Smor_WA_47	1
Smor_WA_48	Smor_WA_48	2

```

**Point plink to the correct files**

```
VCF=Smor.snps_filtered.maf.geno.beagled.WA.biallele.nonvar.vcf
PHENO=PhenotypeFile.txt
```

**Generate plink files**

```
plink --vcf ${VCF} --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --pheno $PHENO --make-bed --pca --out Smor_gemma_assoc
```

**Generate relatedness matrix**

```
gemma -bfile Smor_gemma_assoc -gk 1 -o Smor_relatedness.out > gemma.Smor_relatedness.out
```

**Run association, using relatedness matrix as covariate**

```
gemma -bfile Smor_gemma_assoc -lmm 1 -o Smor_silver.assoc.gemma -k Smor_relatedness.out.cXX.txt > Smor_silver.assoc.gemma.out
```

**Extract p-values**

```
awk '{print $1,$3,$12}' Smor_silver.assoc.gemma.assoc.txt > Smor_silver.assoc.gemma.assoc.txt.pvals
```

**Check minimum p-value**

```
awk 'min=="" || $3 < min {min=$3; minline=$0}; END{ print min, minline}' *pvals
```

**Optional: Reduce file size by filtering for min of 0.05**

```
awk '{ if($3 <= 0.05) { print }}' *pvals > Smor_silver.assoc.gemma.pval.noNA.0.05
```

## 5. Import p-values and chromosome positions in R and plot them

Functions to adjust scaffold positions

**Calculate chromosome coordinates:**

```
chrom.coords <- function(scafL,chromNames, gap = 2000000) {
  chromosome = vector()
  chromLengths  = vector()
  chromStarts = vector()
  chromEnds = vector()
  chromMid = vector()
  chrom = 1
  endLast = 0
  scafCurrent <- subset(scafL, chromosome == chromNames[1])
  chromosome[chrom] <- chrom
  chromLengths[chrom] <- sum(scafCurrent$length)
  chromStarts[chrom] <- endLast + 1
  chromEnds[chrom] <- endLast + chromLengths[chrom]
  chromMid[chrom] <- endLast + chromLengths[chrom]/2
  endLast = chromEnds[chrom]
  chrom = chrom + 1
  for (i in 2:length(chromNames)) {
    chromosome[chrom] <- chrom
    scafCurrent <- subset(scafL, chromosome == chromNames[i])
    chromLengths[chrom] <- sum(scafCurrent$length)
    chromStarts[chrom] <- endLast + gap + 1
    chromEnds[chrom] <- endLast + gap + chromLengths[chrom]
    chromMid[chrom] <- endLast + gap + chromLengths[chrom]/2
    endLast = chromEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(chromosome,chromLengths,chromStarts,chromEnds,chromMid))
  return(table)
}
```

**Calculate scaffold coordinates**

```
scaf.coords <- function(scafL,gap = 0) {
  scaffold = vector()
  scafStarts = vector()
  scafEnds = vector()
  chrom = 1
  endLast = 0
  for (e in 1:nrow(scafL)) {
    scaffold[chrom] <- levels(scafL$scaffold)[e]
    scafStarts[chrom] <- endLast + gap + 1
    scafEnds[chrom] <- endLast + gap + scafL$length[e]
    endLast = scafEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(scaffold,scafStarts,scafEnds))
  return(table)
}
```

**Specify PDF output with desired size**

```
pdf("GWAS_plot_full_WA.pdf", width = 8, height = 2)
```

Import your genome to create the object scaf_coords2 which will be merged with the data frame you want to plot. 
The scafL object below is just a table with the length of each of your scaffolds/chromosomes.

```
scafL<-read.table("Speyeria_scaffold_lengths_v1.txt",h=T)
chromNames <- c(1:30)
```

```
chrom_coords <- chrom.coords(scafL,chromNames)
scaf_coords2 <- merge(scafL,chrom_coords,by="chromosome", all.x=TRUE)
```

**Read in the file generated by GEMMA**

```
fileX <- read.table("Smor_silver.assoc.gemma.assoc.txt.pvals", h=T)
```

**Check file**

```
head(fileX)
```

**Filter for more signifanct snps to make plotting easier and omit NA values**

```
fileX <- na.omit(fileX) 
min(fileX$P)
nrow(fileX[fileX$P < 0.005,])
fileX <- fileX[fileX$P < 0.005,]
nrow(fileX)
```

**Make your life easier by having the CHROM / scaffold column called ?scaffold?, and your position called "pos". If the names don't match, it can be renamed like so:**

```
names(fileX)[names(fileX) == 'CHR'] <- 'scaffold'

fileX_merged <- merge(fileX, scaf_coords2,by = "scaffold", all.x=TRUE)
fileX_merged_pos <- cbind(fileX_merged, chromPos = fileX_merged$BP + fileX_merged$scafStart + fileX_merged$chromStarts-2)
```

Now here is the plotting code. 

top and bot set the y axis limits
in "begin" and "end", set the number of chromosmes. 

**Setup for plot: there is a variable to be plotted called -log10(comp[[i]]$P) , which is taking the column $P from fileX_merged_pos and doing a log10 transform on it.**

Your setup

```
comp=list(fileX_merged_pos)
names=c("fileXX_name")
col=c("darkgrey")

par(mfrow=c(1,1), mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 18
bot = 2.5
```

**Check how many chromosomes you have**

```
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[30]/1000000
```

**Plot stats**

```
for (i in 1:length(comp)){
  plot(0, pch = "", xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "", axes = FALSE)
  rect(chrom_coords$chromStarts/1000000, rep(bot, length(chrom_coords$chromStarts)),
       chrom_coords$chromEnds/1000000, rep(top, length(chrom_coords$chromStarts)),
       col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))
  
  #Plot your stats
  par(new = TRUE)
  
  #Logic for coloring points red if -log10(P) > 9.26
  point_colors <- ifelse(-log10(comp[[i]]$P) > 9.26, "red", adjustcolor(col[i], alpha=0.5))
  
  plot(comp[[i]]$chromPos/1000000, -log10(comp[[i]]$P), type = "p", pch = 19, 
       cex = 0.3, col = point_colors, xlim = c(begin,end), ylim = c(bot,top), 
       axes = FALSE, bty = "n", xlab = "", ylab = "", yaxt = "n", xaxt = "n")
  
  axis(2, cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("-Log10 P-val")), cex = 0.5, line = 0.8)
  
  #Add genome-wide significance line at y = 9.26 (Bonferroni-corrected for alpha 0.01)
  abline(h = 9.26, lty = 2, col = "red")
}

#Plot chromosome labels
axis(1, at = chrom_coords[,5][1:30]/1000000, labels = (1:30), lwd = 0, lwd.ticks = 0)

#Plot scale (adjust this if needed)
#segments(c(390, 390, 400), c(250.07, 250.0702, 250.0702), c(400, 390, 400), c(250.07, 250.0702, 250.0702), lwd = 1)
#text(395, 300.086, labels = "10Mb", cex = 1)
```

Close the PDF device

```
dev.off()
```

## 6 Code for calculating Fst on morphs by populations 

I used Simon Martin's scripts from his repo [genomics_general](https://github.com/simonhmartin/genomics_general) for Fst calculations in windows.
For Fst calculations, I went back to my original vcf file generated by gatk and kept nonvariant sites.
I then performed the following analyses on a vcf file filtered for chromosome 14 only, to speed up calculations.


**chromosome filtering and indexing using bcftools**

```
bcftools filter --threads 32 Smor.snps_filtered.vcf.gz -r Smor1400 -O z -o Smor.snps_filtered.Smor1400.vcf.gz
bcftools index --threads 40 Smor.snps_filtered.Smor1400.vcf.gz
```

**removing asetrisks from vcf file that cause problems downstream with bcftools**

```
bcftools filter --threads 40 -e 'ALT="*" || GT~"\*"' Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.vcf.gz -O z -o Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.vcf.gz
```

**phasing and imputing using beagle and indexing using bcftools**

```
java -Xmx50g -jar /CCAS/home/lucalivraghi/tools/beagle/beagle.05May22.33a.jar impute=FALSE window=1 overlap=0.1 nthreads=40 gt=Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.vcf.gz out=Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled
bcftools index --threads 40 Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.vcf.gz
```

**generating geno file format from vcf, required for Simon Martin's [genomics_general](https://github.com/simonhmartin/genomics_general) scripts**

```
python /CCAS/home/lucalivraghi/tools/genomics_general/VCF_processing/parseVCF.py \
--skipIndels -i \
Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.vcf.gz | gzip > Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.geno.gz
```

**Generate a popsfile that encodes each phenotype and population**

Popsfile for generating popgen statistics:

```
Scor_WA_01 Scor
Shyd_WA_05 Shyd
Shyd_WA_06 Shyd
Shyd_WA_08 Shyd
Shyd_WA_09 Shyd
Shyd_WA_10 Shyd
Smor_WA_01 Smor_WA_S
Smor_WA_02 Smor_WA_S
Smor_WA_07 Smor_WA_B
Smor_WA_11 Smor_WA_S
Smor_WA_13 Smor_WA_B
Smor_WA_14 Smor_WA_B
Smor_WA_15 Smor_WA_S
Smor_WA_16 Smor_WA_B
Smor_WA_17 Smor_WA_B
Smor_WA_18 Smor_WA_B
Smor_WA_19 Smor_WA_S
Smor_WA_20 Smor_WA_S
Smor_WA_21 Smor_WA_S
Smor_WA_22 Smor_WA_S
Smor_WA_23 Smor_WA_B
Smor_WA_24 Smor_WA
Smor_WA_25 Smor_WA_S
Smor_WA_26 Smor_WA_S
Smor_WA_27 Smor_WA_S
Smor_WA_28 Smor_WA_B
Smor_WA_29 Smor_WA_S
Smor_WA_30 Smor_WA_B
Smor_WA_31 Smor_WA_S
Smor_WA_32 Smor_WA_B
Smor_WA_33 Smor_WA_S
Smor_WA_34 Smor_WA_S
Smor_WA_35 Smor_WA_B
Smor_WA_37 Smor_WA_B
Smor_WA_38 Smor_WA_S
Smor_WA_39 Smor_WA_S
Smor_WA_40 Smor_WA_B
Smor_WA_41 Smor_WA_B
Smor_WA_42 Smor_WA_S
Smor_WA_43 Smor_WA_B
Smor_WA_44 Smor_WA_S
Smor_WA_45 Smor_WA_B
Smor_WA_46 Smor_WA_B
Smor_WA_47 Smor_WA_B
Smor_WA_48 Smor_WA_S
NVJ_51_Smor_B Smor_NVJ_B
NVJ_52_Smor_S Smor_NVJ_S
NVJ_53_Smor_B Smor_NVJ_B
NVJ_54_Smor_B Smor_NVJ_B
NVJ_55_Smor_S Smor_NVJ_S
NVJ_56_Smor_S Smor_NVJ_S
NVJ_57_Smor_B Smor_NVJ_B
NVJ_58_Smor_B Smor_NVJ_B
NVJ_59_Smor_B Smor_NVJ_B
NVJ_60_Smor_S Smor_NVJ_S
NVJ_61_Smor_B Smor_NVJ_B
NVJ_62_Smor_S Smor_NVJ_S
NVJ_63_Smor_S Smor_NVJ_S
NVJ_64_Smor_B Smor_NVJ_B
NVJ_65_Smor_S Smor_NVJ_S
NVJ_66_Smor_S Smor_NVJ_S
NVJ_67_Smor_S Smor_NVJ_S
NVJ_68_Smor_B Smor_NVJ_B
NVJ_69_Smor_B Smor_NVJ_B
NVJ_70_Smor_B Smor_NVJ_B
NVJ_71_Smor_B Smor_NVJ_B
NVJ_72_Smor_S Smor_NVJ_S
NVJ_73_Smor_B Smor_NVJ_B
NVJ_74_Smor_B Smor_NVJ_B
NVJ_75_Smor_S Smor_NVJ_S
NVJ_76_Smor_B Smor_NVJ_B
NVJ_77_Smor_S Smor_NVJ_S
NVJ_78_Smor_S Smor_NVJ_S
NVJ_79_Smor_S Smor_NVJ_S
NVJ_80_Smor_S Smor_NVJ_S
NVJ_81_Smor_S Smor_NVJ_S
```

**Run popgenWindows.py in 200bp windows with a 100bp slide**

```
/CCAS/home/lucalivraghi/tools/popgen_scripts/genomics_general/popgenWindows.py -w 200 -s 100 -g ../../geno_files/Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.geno.gz -o popgenstats_Smor_WA_S.vs.SmorWA_B_200w100s_nonphased.csv.gz -f phased -T 32 -p Smor_WA_S -p Smor_WA_B --popsFile ../../popsfiles/popsfileind.txt
/CCAS/home/lucalivraghi/tools/popgen_scripts/genomics_general/popgenWindows.py -w 200 -s 100 -g ../../geno_files/Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.geno.gz -o popgenstats_Smor_NVJ_S.vs.SmorNVJ_B_200w100s_nonphased.csv.gz -f phased -T 32 -p Smor_NVJ_S -p Smor_NVJ_B --popsFile ../../popsfiles/popsfileind.txt
```


## 7. Plot Fst and -log10 P from GWAS on chromosome region of interest in R

I then import the Fst data generated from Simon Martin's [genomics_general](https://github.com/simonhmartin/genomics_general) scripts and plot them alongside the GWAS P values.

```
#Code for plotting popgen stats at chromosome 14
library(ggplot2)
library(gridExtra)
library(tidyverse)

#Read the file contiaing the popgenstats for WA pops (Fst, Dxy, Pi)
file.popgen.stats <- read.csv("popgenstats_Smor_WA_S.vs.SmorWA_B_200w100s_nonphased.csv.gz", h=T)
head(file.popgen.stats)

#Firs read the GWAS pvalues from the exported table from GEMMA
file.chrm14 <- read.table("chrom14.WA.pvals", h=T)
head(file.chrm14)

#Make the plot for optix annotation (CDS and UTR)
p0 <- ggplot() + 
  geom_rect(aes(xmin = 1506454, xmax = 1506853, ymin = 0, ymax = 1), fill = "grey") + # No specific y-axis placement
  geom_rect(aes(xmin = 1492189, xmax = 1492995, ymin = 0, ymax = 1), fill = "grey") + 
  xlim(1480000, 1550000) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

#Make plot for GWAS on subset of chromosome 14

#log transform the pvalues
logp <- -log10(file.chrm14$P)

#Make the plot for the GWAS. 
p1 <- ggplot(data = file.chrm14) +
  geom_point(aes(x = BP, y = logp, color = logp > 9.26), pch = 19, size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  theme_classic() + 
  xlim(1480000, 1550000) +
  ylim(0, 16) +
  theme(
    legend.position = "none", 
    axis.title.y = element_text(size = 8), 
    axis.text.y = element_text(size = 6), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank()
  ) +
  ylab("-Log10 P-VALUE") + 
  xlab("Position (bp)") +
  geom_hline(yintercept = 9.26, linetype = "dashed", color = "red")

#Make the plot for FST between silver and buff INDIVIDUALS

#Calculate 99% quantile threshold for Fst_morm_WA.silver_morm_WA.buff

threshold <- quantile(file.popgen.stats$Fst_Smor_WA_S_Smor_WA_B, 0.99, na.rm = TRUE)

#Create a new column 'outlier' to flag points above the threshold
file.popgen.stats$outlier <- ifelse(file.popgen.stats$Fst_Smor_WA_S_Smor_WA_B > threshold, "outlier", "non-outlier")

p2 <- ggplot(data = file.popgen.stats, aes(x = mid, y = Fst_Smor_WA_S_Smor_WA_B)) +
  geom_point(aes(color = outlier), pch = 19, size = 1.5, alpha = 0.5) + # Map 'outlier' to color
  scale_color_manual(values = c("outlier" = "red", "non-outlier" = "grey")) + # Define colors for outliers and non-outliers
  geom_smooth(method = "loess", se = FALSE, color = "black", span = 0.05, size = 0.8) +
  theme_classic() + 
  xlim(1480000, 1550000) +
  ylim(0, 0.4) +
  theme(
    legend.position = "none", 
    axis.title.y = element_text(size = 8), 
    axis.text.y = element_text(size = 6), 
  ) +
  ylab("FST") + 
  xlab("Position (bp)")


p0 <- ggplotGrob(p0)
p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)

maxWidth = grid::unit.pmax(p0$widths[2:5], p1$widths[2:5], p2$widths[2:5])
p0$widths[2:5] <- as.list(maxWidth)
p1$widths[2:5] <- as.list(maxWidth)
p2$widths[2:5] <- as.list(maxWidth)

grid.arrange(p0, p1, p2, ncol = 1, heights = c(2.5, 10 ,10))
```

## 8. Principle component analysis for genome-wide and GWAS specific regions

Here I used plink to generate PCAs for both the genome-wide dataset and for the region spanning the GWAS signal.
These were run on the filtered set of vcfs produced at the end of step 4. 

**First prune SNPs for linkage using plink**

```
plink --vcf morm.NV.WA.filtered.phased.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out linkage_pruned
```

**Then run the PCA on the genomewide vcf containing both pops (WA and NV)**

```
plink --vcf morm.NV.WA.filtered.phased.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract linkage_pruned.prune.in \
--make-bed --pca --out PCA_WA_NV
```

**Now subset the vcf to the GWAS region, and re-run PCA for just this genomic interval**

```
bcftools filter --threads 32 morm.NV.WA.filtered.phased.vcf.gz -r Smor1400:1508454-1511684 -O z -o GWAS.morm.NV.WA.filtered.phased.vcf.gz

plink --vcf GWAS.morm.NV.WA.filtered.phased.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--make-bed --pca var-wts --out PCA_WA_NV_GWAS
```

## 9. Plot PCAs in R

library(tidyverse)

**Plot genome-wide PCA**

```
#Read in plink output

pca <- read_table("PCA_WA_NV.eigenvec", col_names = FALSE)
eigenval <- scan("PCA_WA_NV.eigenval")

# remove nuisance column
pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))
spp[grep("NVJ", pca$ind)] <- "Nevada"
spp[grep("WA", pca$ind)] <- "Washington"

# remake data.frame
pca <- as_tibble(data.frame(pca, spp))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot for eignevalues
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)


# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = spp)) + 
  geom_point(size = 5, alpha = 0.6) + 
  scale_shape_manual(values = c("Washington" = 17, "Nevada" = 16)) + 
  scale_colour_manual(values = c("red", "blue")) + 
  theme_classic()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  guides(color = "none", shape = "none")



dev.copy(pdf, 'genowidepca_triangle.pdf',width=12,height=7)
dev.off()


#Plot GWAS interval PCA#


#Read in plink output#
pca <- read_table("PCA_WA_NV_GWAS.eigenvec", col_names = FALSE)
eigenval <- scan("PCA_WA_NV_GWAS.eigenval")

# remove nuisance column
pca <- pca[,-1]


# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))
spp[grep("NVJ", pca$ind)] <- "Nevada"
spp[grep("WA", pca$ind)] <- "Washington"

# sort out the individual species and pops
# spp
pheno <- rep(NA, length(pca$ind))
pheno[grep("_S", pca$ind)] <- "Silver"
pheno[grep("_B", pca$ind)] <- "Buff_Het"

pheno[grep("NVJ_53_Smor_B.bam", pca$ind)] <- "Buff_Hom"
pheno[grep("NVJ_61_Smor_B.bam", pca$ind)] <- "Buff_Hom"
pheno[grep("NVJ_64_Smor_B.bam", pca$ind)] <- "Buff_Hom"
pheno[grep("NVJ_68_Smor_B.bam", pca$ind)] <- "Buff_Hom"
pheno[grep("NVJ_71_Smor_B.bam", pca$ind)] <- "Buff_Hom"
pheno[grep("NVJ_73_Smor_B.bam", pca$ind)] <- "Buff_Hom"


pheno[grep("Smor_WA_01", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_02", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_07", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_1", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_13", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_14", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_15", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_16", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_17", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_18", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_19", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_20", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_21", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_22", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_23", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_24", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_25", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_26", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_27", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_28", pca$ind)] <- "Buff_Hom"
pheno[grep("Smor_WA_29", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_30", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_31", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_32", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_33", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_34", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_35", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_37", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_38", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_39", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_40", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_41", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_42", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_43", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_44", pca$ind)] <- "Silver"
pheno[grep("Smor_WA_45", pca$ind)] <- "Buff_Hom"
pheno[grep("Smor_WA_46", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_47", pca$ind)] <- "Buff_Het"
pheno[grep("Smor_WA_48", pca$ind)] <- "Silver"


# remake data.frame
pca <- as_tibble(data.frame(pca, pheno, spp))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot for eignevalues
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)


# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = pheno, shape = spp)) + geom_point(size = 5)
b <- b + scale_colour_manual(values = c("red", "blue", "orange"))
b <- b + theme_classic()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  guides(shape = "none", color = "none")

dev.copy(pdf, 'PCA_GWAS_interval.pdf',width=12,height=7)
dev.off()
```

## 10. Genotype plot for SNP illutsration at the GWAS interval 

I jused Jim Whiting's genotype_plot code to plot each SNP along the GWAS interval as either reference, heterozygous for reference/alternative allele, or homozygous for alternative allele. 

**First, generate a popsfile encoding each sample**

```
nano popsfile.txt

Smor_WA_01      Silver_Smor_WA
Smor_WA_02      Silver_Smor_WA
Smor_WA_07      Buff_Smor_WA
Smor_WA_11      Silver_Smor_WA
Smor_WA_13      Buff_Smor_WA
Smor_WA_14      Buff_Smor_WA
Smor_WA_15      Silver_Smor_WA
Smor_WA_16      Buff_Smor_WA
Smor_WA_17      Buff_Smor_WA
Smor_WA_18      Buff_Smor_WA
Smor_WA_19      Silver_Smor_WA
Smor_WA_20      Silver_Smor_WA
Smor_WA_21      Silver_Smor_WA
Smor_WA_22      Silver_Smor_WA
Smor_WA_23      Buff_Smor_WA
Smor_WA_24      Silver_Smor_WA
Smor_WA_25      Silver_Smor_WA
Smor_WA_26      Silver_Smor_WA
Smor_WA_27      Silver_Smor_WA
Smor_WA_28      Buff_Smor_WA
Smor_WA_29      Silver_Smor_WA
Smor_WA_30      Buff_Smor_WA
Smor_WA_31      Silver_Smor_WA
Smor_WA_32      Buff_Smor_WA
Smor_WA_33      Silver_Smor_WA
Smor_WA_34      Silver_Smor_WA
Smor_WA_35      Buff_Smor_WA
Smor_WA_37      Buff_Smor_WA
Smor_WA_38      Silver_Smor_WA
Smor_WA_39      Silver_Smor_WA
Smor_WA_40      Buff_Smor_WA
Smor_WA_41      Buff_Smor_WA
Smor_WA_42      Silver_Smor_WA
Smor_WA_43      Buff_Smor_WA
Smor_WA_44      Silver_Smor_WA
Smor_WA_45      Buff_Smor_WA
Smor_WA_46      Buff_Smor_WA
Smor_WA_47      Buff_Smor_WA
Smor_WA_48      Silver_Smor_WA
NVJ_51_Smor_B   Buff_Smor_NVJ
NVJ_52_Smor_S   Silver_Smor_NVJ
NVJ_53_Smor_B   Buff_Smor_NVJ
NVJ_54_Smor_B   Buff_Smor_NVJ
NVJ_55_Smor_S   Silver_Smor_NVJ
NVJ_56_Smor_S   Silver_Smor_NVJ
NVJ_57_Smor_B   Buff_Smor_NVJ
NVJ_58_Smor_B   Buff_Smor_NVJ
NVJ_59_Smor_B   Buff_Smor_NVJ
NVJ_60_Smor_S   Silver_Smor_NVJ
NVJ_61_Smor_B   Buff_Smor_NVJ
NVJ_62_Smor_S   Silver_Smor_NVJ
NVJ_63_Smor_S   Silver_Smor_NVJ
NVJ_64_Smor_B   Buff_Smor_NVJ
NVJ_65_Smor_S   Silver_Smor_NVJ
NVJ_66_Smor_S   Silver_Smor_NVJ
NVJ_67_Smor_S   Silver_Smor_NVJ
NVJ_68_Smor_B   Buff_Smor_NVJ
NVJ_69_Smor_B   Buff_Smor_NVJ
NVJ_70_Smor_B   Buff_Smor_NVJ
NVJ_71_Smor_B   Buff_Smor_NVJ
NVJ_72_Smor_S   Silver_Smor_NVJ
NVJ_73_Smor_B   Buff_Smor_NVJ
NVJ_74_Smor_B   Buff_Smor_NVJ
NVJ_75_Smor_S   Silver_Smor_NVJ
NVJ_76_Smor_B   Buff_Smor_NVJ
NVJ_77_Smor_S   Silver_Smor_NVJ
NVJ_78_Smor_S   Silver_Smor_NVJ
NVJ_79_Smor_S   Silver_Smor_NVJ
NVJ_80_Smor_S   Silver_Smor_NVJ
NVJ_81_Smor_S   Silver_Smor_NVJ
```

**Then use the tool to plot the SNPs in R**

I used the same vcf that was used to generate the PCAs, focusing on the GWAS region.

```
library(GenotypePlot)

popmap <- read.table('popsfile.txt', h=F)

new_plot <- genotype_plot(vcf="morm.NV.WA.filtered.phased.vcf.gz",
                          chr="Smor1400",
                          start=1507000,
                          end=1510500,
                          popmap=popmap,
                          cluster=TRUE,
                          snp_label_size=1000,
                          invariant_filter = TRUE,
                          #polarise_genotypes='buff',
                          colour_scheme=c("#1babe2","#b0207c","#f89c30"))
combine_genotype_plot(new_plot)
dev.copy(pdf, 'optixPlot_GWAS_zoomin.pdf',width=20,height=10)
dev.off()
```

# II. Figure 4 - Haplotype based statistics, sweeps and introgression analyses (Fst, Dxy, Tajima's D, CLR, TWISST, Fd)

The following code is used to calculate popgen statistics in sliding windows on **phased haplotypes**, having assigned each haplotype group based on three diagnostic SNPs idenitifed in Figure 2.
For these statistics, a vcf file **including invariant SNPs** was used. 
Analyses were limited to chromosome 14 to speed up computation times.
Again, Simon Martin's [genomics_general](https://github.com/simonhmartin/genomics_general) scripts were used for calculations.

### 1. Vcf filtering and recoding

**The output vcf from gatk still contained some single "*" denoting single deletions, which were causing problems in the downstream analyses, so these are removed using bcftools:**

```
bcftools filter --threads 40 -e 'ALT="*" || GT~"\*"' Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.vcf.gz -O z -o Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.vcf.gz
```

**The vcf file is then phased using beagle and inexed using bcftools**

```
java -Xmx50g -jar /CCAS/home/lucalivraghi/tools/beagle/beagle.05May22.33a.jar impute=FALSE window=1 overlap=0.1 nthreads=40 gt=Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.vcf.gz out=Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled
bcftools index --threads 40 Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.vcf.gz
```

**A custom script is then used to separate each haplotype into a new entry in the vcf, recoding individual haplotypes as either "_A" or "_B" haplotypes.**

```
zcat Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.vcf.gz | \
awk '{
    if ($1 ~ /^##/) {
        print;
    } else if ($1=="#CHROM") {
        ORS="\t";
        for (i=1; i<=9; i++) print $i;
        for (i=10; i<NF; i++) {print $i"_A\t"$i"_B";}
        ORS="\n";
        print $NF"_A\t"$NF"_B";
    } else {
        ORS="\t";
        for (i=1; i<=9; i++) print $i;
        for (i=10; i<NF; i++) print substr($i, 1, 1)"\t"substr($i, 3, 1);
        ORS="\n";
        print substr($NF, 1, 1)"\t"substr($NF, 3, 1);
    }
}' > Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.haplo.vcf.gz
```

**Then recode this into geno format to run Simon Martin's [genomics_general](https://github.com/simonhmartin/genomics_general) scripts**

```
python /CCAS/home/lucalivraghi/tools/genomics_general/VCF_processing/parseVCF.py \
--skipIndels \
--ploidy 1 \
-i Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.haplo.vcf.gz | gzip > Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.geno.haplo.gz
```

**The resulting geno file was then inspected at each of the diagnostic SNP positions (Smor1400:1508598, Smor1400:1508581, Smor1400:1508599), and used to recode the popsfile**
**Buff alleles have "G,T,T" and silver alleles "T,C,A" at those positions.**

```
nano popsfilehaplo.txt

Scor_WA_01_A	coronis.silver
Scor_WA_01_B	coronis.silver
Shyd_WA_05_A	hydaspe.buff
Shyd_WA_05_B	hydaspe.buff
Shyd_WA_06_A	hydaspe.buff
Shyd_WA_06_B	hydaspe.buff
Shyd_WA_08_A	hydaspe.buff
Shyd_WA_08_B	hydaspe.buff
Shyd_WA_09_A	hydaspe.buff
Shyd_WA_09_B	hydaspe.buff
Shyd_WA_10_A	hydaspe.buff
Shyd_WA_10_B	hydaspe.buff
Smor_WA_01_A	morm_WA.silver
Smor_WA_01_B	morm_WA.silver
Smor_WA_02_A	morm_WA.silver
Smor_WA_02_B	morm_WA.silver
Smor_WA_07_A	morm_WA.buff
Smor_WA_07_B	morm_WA.silver
Smor_WA_11_A	morm_WA.silver
Smor_WA_11_B	morm_WA.silver
Smor_WA_13_A	morm_WA.silver
Smor_WA_13_B	morm_WA.buff
Smor_WA_14_A	morm_WA.buff
Smor_WA_14_B	morm_WA.silver
Smor_WA_15_A	morm_WA.silver
Smor_WA_15_B	morm_WA.silver
Smor_WA_16_A	morm_WA.silver
Smor_WA_16_B	morm_WA.buff
Smor_WA_17_A	morm_WA.silver
Smor_WA_17_B	morm_WA.buff
Smor_WA_18_A	morm_WA.silver
Smor_WA_18_B	morm_WA.buff
Smor_WA_19_A	morm_WA.silver
Smor_WA_19_B	morm_WA.silver
Smor_WA_20_A	morm_WA.silver
Smor_WA_20_B	morm_WA.silver
Smor_WA_21_A	morm_WA.silver
Smor_WA_21_B	morm_WA.silver
Smor_WA_22_A	morm_WA.silver
Smor_WA_22_B	morm_WA.silver
Smor_WA_23_A	morm_WA.buff
Smor_WA_23_B	morm_WA.silver
Smor_WA_24_A	morm_WA.silver
Smor_WA_24_B	morm_WA.silver
Smor_WA_25_A	morm_WA.silver
Smor_WA_25_B	morm_WA.silver
Smor_WA_26_A	morm_WA.silver
Smor_WA_26_B	morm_WA.silver
Smor_WA_27_A	morm_WA.silver
Smor_WA_27_B	morm_WA.silver
Smor_WA_28_A	morm_WA.buff
Smor_WA_28_B	morm_WA.buff
Smor_WA_29_A	morm_WA.silver
Smor_WA_29_B	morm_WA.silver
Smor_WA_30_A	morm_WA.silver
Smor_WA_30_B	morm_WA.buff
Smor_WA_31_A	morm_WA.silver
Smor_WA_31_B	morm_WA.silver
Smor_WA_32_A	morm_WA.silver
Smor_WA_32_B	morm_WA.buff
Smor_WA_33_A	morm_WA.silver
Smor_WA_33_B	morm_WA.silver
Smor_WA_34_A	morm_WA.silver
Smor_WA_34_B	morm_WA.silver
Smor_WA_35_A	morm_WA.silver
Smor_WA_35_B	morm_WA.buff
Smor_WA_37_A	morm_WA.silver
Smor_WA_37_B	morm_WA.buff
Smor_WA_38_A	morm_WA.silver
Smor_WA_38_B	morm_WA.silver
Smor_WA_39_A	morm_WA.silver
Smor_WA_39_B	morm_WA.silver
Smor_WA_40_A	morm_WA.silver
Smor_WA_40_B	morm_WA.buff
Smor_WA_41_A	morm_WA.buff
Smor_WA_41_B	morm_WA.silver
Smor_WA_42_A	morm_WA.silver
Smor_WA_42_B	morm_WA.silver
Smor_WA_43_A	morm_WA.buff
Smor_WA_43_B	morm_WA.silver
Smor_WA_44_A	morm_WA.silver
Smor_WA_44_B	morm_WA.silver
Smor_WA_45_A	morm_WA.buff
Smor_WA_45_B	morm_WA.buff
Smor_WA_46_A	morm_WA.buff
Smor_WA_46_B	morm_WA.silver
Smor_WA_47_A	morm_WA.silver
Smor_WA_47_B	morm_WA.buff
Smor_WA_48_A	morm_WA.silver
Smor_WA_48_B	morm_WA.silver
NVJ_51_Smor_B_A	morm_NVJ.buff
NVJ_51_Smor_B_B	morm_NVJ.silver
NVJ_52_Smor_S_A	morm_NVJ.silver
NVJ_52_Smor_S_B	morm_NVJ.silver
NVJ_53_Smor_B_A	morm_NVJ.buff
NVJ_53_Smor_B_B	morm_NVJ.buff
NVJ_54_Smor_B_A	morm_NVJ.buff
NVJ_54_Smor_B_B	morm_NVJ.silver
NVJ_55_Smor_S_A	morm_NVJ.silver
NVJ_55_Smor_S_B	morm_NVJ.silver
NVJ_56_Smor_S_A	morm_NVJ.silver
NVJ_56_Smor_S_B	morm_NVJ.silver
NVJ_57_Smor_B_A	morm_NVJ.buff
NVJ_57_Smor_B_B	morm_NVJ.silver
NVJ_58_Smor_B_A	morm_NVJ.buff
NVJ_58_Smor_B_B	morm_NVJ.silver
NVJ_59_Smor_B_A	morm_NVJ.buff
NVJ_59_Smor_B_B	morm_NVJ.silver
NVJ_60_Smor_S_A	morm_NVJ.silver
NVJ_60_Smor_S_B	morm_NVJ.silver
NVJ_61_Smor_B_A	morm_NVJ.buff
NVJ_61_Smor_B_B	morm_NVJ.buff
NVJ_62_Smor_S_A	morm_NVJ.silver
NVJ_62_Smor_S_B	morm_NVJ.silver
NVJ_63_Smor_S_A	morm_NVJ.silver
NVJ_63_Smor_S_B	morm_NVJ.silver
NVJ_64_Smor_B_A	morm_NVJ.buff
NVJ_64_Smor_B_B	morm_NVJ.buff
NVJ_65_Smor_S_A	morm_NVJ.silver
NVJ_65_Smor_S_B	morm_NVJ.silver
NVJ_66_Smor_S_A	morm_NVJ.silver
NVJ_66_Smor_S_B	morm_NVJ.silver
NVJ_67_Smor_S_A	morm_NVJ.silver
NVJ_67_Smor_S_B	morm_NVJ.silver
NVJ_68_Smor_B_A	morm_NVJ.buff
NVJ_68_Smor_B_B	morm_NVJ.buff
NVJ_69_Smor_B_A	morm_NVJ.silver
NVJ_69_Smor_B_B	morm_NVJ.buff
NVJ_70_Smor_B_A	morm_NVJ.buff
NVJ_70_Smor_B_B	morm_NVJ.silver
NVJ_71_Smor_B_A	morm_NVJ.buff
NVJ_71_Smor_B_B	morm_NVJ.buff
NVJ_72_Smor_S_A	morm_NVJ.silver
NVJ_72_Smor_S_B	morm_NVJ.silver
NVJ_73_Smor_B_A	morm_NVJ.buff
NVJ_73_Smor_B_B	morm_NVJ.buff
NVJ_74_Smor_B_A	morm_NVJ.silver
NVJ_74_Smor_B_B	morm_NVJ.buff
NVJ_75_Smor_S_A	morm_NVJ.silver
NVJ_75_Smor_S_B	morm_NVJ.silver
NVJ_76_Smor_B_A	morm_NVJ.buff
NVJ_76_Smor_B_B	morm_NVJ.silver
NVJ_77_Smor_S_A	morm_NVJ.silver
NVJ_77_Smor_S_B	morm_NVJ.silver
NVJ_78_Smor_S_A	morm_NVJ.silver
NVJ_78_Smor_S_B	morm_NVJ.silver
NVJ_79_Smor_S_A	morm_NVJ.silver
NVJ_79_Smor_S_B	morm_NVJ.silver
NVJ_80_Smor_S_A	morm_NVJ.silver
NVJ_80_Smor_S_B	morm_NVJ.silver
NVJ_81_Smor_S_A	morm_NVJ.silver
NVJ_81_Smor_S_B	morm_NVJ.silver
```

### 2. Fst recalculated on phased haplotypes

**Calculate Fst using the recoded individual haplotypes. These are calculated in 200bp windows with a 100bp slide**
**Statistics are calculated on all relevant population comparions, but code is illustrating WA population here**

```
/CCAS/home/lucalivraghi/tools/popgen_scripts/genomics_general/popgenWindows.py \
-w 200 -s 100 \
-g Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.geno.haplo.gz \
-o popgenstats_Smor_WA_S_haplo.vs.SmorWA_B_haplo_200w100s_phased.csv.gz -f haplo -T 32 -p morm_WA.silver -p morm_WA.buff --popsFile popsfilehaplo.txt
```

### 3. Dxy, Pi, and Tajima's calculated on phased haplotypes

**Calculate statistics using the recoded individual haplotypes as above. These are calculated in 200bp windows with a 100bp slide**

```
/CCAS/home/lucalivraghi/tools/popgen_scripts/genomics_general/popgenWindows.py \
-w 200 -s 100 \
-g Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.geno.haplo.gz \
-o popgenstats_Smor_NVJ_S_haplo.vs.SmorWA_B_haplo_200w_100s_phased.freq.csv.gz -f haplo -T 32 -p morm_WA.silver -p morm_WA.buff --popsFile popsfilehaplo.txt --analysis popFreq
```

### 4. Sweep analysis using SweeD

**The same phased and recoded haplotype was used to generate selective sweep statistics using the SweeD program**
**The vcf file is first split into single haplotype vcfs from each population, containing only unsilvered or silver haplotypes based on the diagnostic SNPs**

For example, using bcftools, the NV population vcf for unsilvered alleles is filtered like this:

```
bcftools view --threads 32 --samples NVJ_51_Smor_B_A,NVJ_53_Smor_B_A,NVJ_53_Smor_B_B,NVJ_54_Smor_B_A,NVJ_58_Smor_B_A,NVJ_59_Smor_B_A,NVJ_61_Smor_B_A,NVJ_61_Smor_B_B,NVJ_64_Smor_B_A,NVJ_64_Smor_B_B,NVJ_68_Smor_B_A,NVJ_68_Smor_B_B,NVJ_69_Smor_B_B,NVJ_70_Smor_B_A,NVJ_71_Smor_B_A,NVJ_71_Smor_B_B,NVJ_73_Smor_B_A,NVJ_73_Smor_B_B,NVJ_76_Smor_B_A Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.haplo.vcf.gz -O z -o Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.haplo.buff.NVJ.vcf.gz
```

**Next, vcftools is used to generate a site frequency spectrum file as input for SweeD**

Example on the NV vcf as above:

```
vcftools --counts2 --gzvcf Smor.snps_filtered.Smor1400.WA.NVJ.hyd.cor.removedasterisk.beagled.haplo.buff.NVJ.vcf.gz --stdout | awk 'NR<=1 {next} {print $2"\t"$6"\t"$4"\t0"}' > SF2.haplo.buff.NVJ.input
```

**SweeD is then run on this file using a grid size of 30kb, and including monomorphic sites**

```
SweeD -name SweeD.haplo.buff.mon.optix.NVJ -input SF2.haplo.buff.NVJ.input -grid 30000 -monomorphic
```

### 5. Topology Weighting by Iterative Subsampling (TWISST)

**The topology analysis was run on the WA population only, as this is the population for which we had relevant outgroups**
**Again, this is run using Simon Martin's [genomics_general](https://github.com/simonhmartin/genomics_general) scripts, by generating a geno file first, then recoding each haplotype based on diagnostic SNPs**

```
python /CCAS/home/lucalivraghi/tools/genomics_general/VCF_processing/parseVCF.py --skipIndels -i Chrm14.optix.filtered.WA.vcf.gz | gzip > Chrm14.optix.filtered.WA.geno.gz 

python /CCAS/home/lucalivraghi/tools/popgen_scripts/genomics_general/filterGenotypes.py \
--threads 40 \
--ploidy 2 \
-of bases \
-i Chrm14.optix.filtered.WA.geno.gz | gzip > Chrm14.optix.filtered.WA.geno.haplo.gz
```

**Then make the haplotype specific popsfile**

```
nano popsfilehaplo.txt

Scor_WA_01_A	coronis.silver
Scor_WA_01_B	coronis.silver
Shyd_WA_05_A	hydaspe.buff
Shyd_WA_05_B	hydaspe.buff
Shyd_WA_06_A	hydaspe.silver
Shyd_WA_06_B	hydaspe.buff
Shyd_WA_08_A	hydaspe.buff
Shyd_WA_08_B	hydaspe.buff
Shyd_WA_09_A	hydaspe.buff
Shyd_WA_09_B	hydaspe.buff
Shyd_WA_10_A	hydaspe.buff
Shyd_WA_10_B	hydaspe.silver
Smor_WA_01_A	morm_WA.silver
Smor_WA_01_B	morm_WA.silver
Smor_WA_02_A	morm_WA.silver
Smor_WA_02_B	morm_WA.silver
Smor_WA_07_A	morm_WA.buff
Smor_WA_07_B	morm_WA.silver
Smor_WA_11_A	morm_WA.silver
Smor_WA_11_B	morm_WA.silver
Smor_WA_13_A	morm_WA.silver
Smor_WA_13_B	morm_WA.buff
Smor_WA_14_A	morm_WA.silver
Smor_WA_14_B	morm_WA.buff
Smor_WA_15_A	morm_WA.silver
Smor_WA_15_B	morm_WA.silver
Smor_WA_16_A	morm_WA.silver
Smor_WA_16_B	morm_WA.buff
Smor_WA_17_A	morm_WA.buff
Smor_WA_17_B	morm_WA.silver
Smor_WA_18_A	morm_WA.silver
Smor_WA_18_B	morm_WA.buff
Smor_WA_19_A	morm_WA.silver
Smor_WA_19_B	morm_WA.silver
Smor_WA_20_A	morm_WA.silver
Smor_WA_20_B	morm_WA.silver
Smor_WA_21_A	morm_WA.silver
Smor_WA_21_B	morm_WA.silver
Smor_WA_22_A	morm_WA.silver
Smor_WA_22_B	morm_WA.silver
Smor_WA_23_A	morm_WA.buff
Smor_WA_23_B	morm_WA.silver
Smor_WA_24_A	morm_WA.silver
Smor_WA_24_B	morm_WA.silver
Smor_WA_25_A	morm_WA.silver
Smor_WA_25_B	morm_WA.silver
Smor_WA_26_A	morm_WA.silver
Smor_WA_26_B	morm_WA.silver
Smor_WA_27_A	morm_WA.silver
Smor_WA_27_B	morm_WA.silver
Smor_WA_28_A	morm_WA.buff
Smor_WA_28_B	morm_WA.buff
Smor_WA_29_A	morm_WA.silver
Smor_WA_29_B	morm_WA.silver
Smor_WA_30_A	morm_WA.silver
Smor_WA_30_B	morm_WA.buff
Smor_WA_31_A	morm_WA.silver
Smor_WA_31_B	morm_WA.silver
Smor_WA_32_A	morm_WA.silver
Smor_WA_32_B	morm_WA.buff
Smor_WA_33_A	morm_WA.silver
Smor_WA_33_B	morm_WA.silver
Smor_WA_34_A	morm_WA.silver
Smor_WA_34_B	morm_WA.silver
Smor_WA_35_A	morm_WA.silver
Smor_WA_35_B	morm_WA.buff
Smor_WA_37_A	morm_WA.silver
Smor_WA_37_B	morm_WA.buff
Smor_WA_38_A	morm_WA.silver
Smor_WA_38_B	morm_WA.silver
Smor_WA_39_A	morm_WA.silver
Smor_WA_39_B	morm_WA.silver
Smor_WA_40_A	morm_WA.buff
Smor_WA_40_B	morm_WA.silver
Smor_WA_41_A	morm_WA.silver
Smor_WA_41_B	morm_WA.buff
Smor_WA_42_A	morm_WA.silver
Smor_WA_42_B	morm_WA.silver
Smor_WA_43_A	morm_WA.buff
Smor_WA_43_B	morm_WA.silver
Smor_WA_44_A	morm_WA.silver
Smor_WA_44_B	morm_WA.silver
Smor_WA_45_A	morm_WA.buff
Smor_WA_45_B	morm_WA.buff
Smor_WA_46_A	morm_WA.silver
Smor_WA_46_B	morm_WA.buff
Smor_WA_47_A	morm_WA.buff
Smor_WA_47_B	morm_WA.silver
Smor_WA_48_A	morm_WA.silver
Smor_WA_48_B	morm_WA.silver
```

**Then generate the trees using 50bp windows**

```
module load phyML/3.3 

python /CCAS/home/lucalivraghi/tools/popgen_scripts/genomics_general/phylo/phyml_sliding_windows.py \
-g Chrm14.optix.filtered.WA.geno.haplo.gz \
--prefix Chrm14.optix.filtered.all.trees.haplo \
-T 32 \
--windType sites --model GTR -w 50 --optimise n
```

**The trees are then used as the inputs for the twisst analysis**

```
python /CCAS/home/lucalivraghi/tools/popgen_scripts/genomics_general/twisst/twisst.py -t Chrm14.optix.filtered.all.trees.haplo.trees.gz \
-w SmorWA.haplo.weights.tsv \
-g coronis.silver -g hydaspe.buff -g morm_WA.buff -g morm_WA.silver --groupsFile popsfilehaplo.txt \
--method fixed --iterations 500 --outgroup coronis.silver
```


### 6. Proportion of introgression in sliding windows (fd)

**The same vcf file that is used for the twisst analysis is also used for calculating fd in sliding windows.**
**Analysis is run in 1.5kb windows with a 150bp slide**

```
/CCAS/home/lucalivraghi/tools/genomics_general/ABBABABAwindows.py -g Chrm14.optix.filtered.WA.geno.haplo.gz \
-f haplo \
--ploidy 1 \
-o ABBABABBA.csv \
-w 1500 -m 10 -s 150  \
-P1 morm_WA.silver -P2 morm_WA.buff -P3 hydaspe.buff -O coronis.silver -T 25  \
--minData 0.5  \
--popsFile popsfilehaplo.txt 
```

### 7. Plot all haplotype based statistics in R

**These stats are all plotted in R. The code below shows is shown as representative of the WA pops plots**

#Code for plotting popgen stats at chromosome 14
library(ggplot2)
library(gridExtra)
library(tidyverse)

```
#Read the file contiaing the popgenstats for WA pops (Fst, Dxy, Pi)
file.popgen.stats <- read.csv("popgenstats_Smor_WA_S_haplo.vs.SmorWA_B_haplo_200w100s_phased.csv.gz", h=T)
head(file.popgen.stats)

file.popgen.stats.freq <- read.csv("popgenstats_Smor_WA_S_haplo.vs.SmorWA_B_haplo_200w_100s_phased.freq.csv.gz", h=T)
head(file.popgen.stats.freq)

SDB <- read.table("SweeD_Report.SweeD.haplo.buff.silver.mon.optix.WA", h=T)
head(SDB)

#Make the plot for optix annotation (CDS and UTR)


p0 <- ggplot() + 
  geom_rect(aes(xmin = 1506454, xmax = 1506853, ymin = 0, ymax = 1), fill = "grey") + # No specific y-axis placement
  geom_rect(aes(xmin = 1492189, xmax = 1492995, ymin = 0, ymax = 1), fill = "grey") + 
  xlim(1480000, 1550000) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )


#Make the plot for FST between silver and buff HAPLOTYPES

# Calculate 99% quantile threshold for Fst_morm_WA.silver_morm_WA.buff
threshold <- quantile(file.popgen.stats$Fst_morm_WA.silver_morm_WA.buff, 0.99, na.rm = TRUE)


# Create a new column 'outlier' to flag points above the threshold
file.popgen.stats$outlier <- ifelse(file.popgen.stats$Fst_morm_WA.silver_morm_WA.buff > threshold, "outlier", "non-outlier")


p1 <- ggplot(data = file.popgen.stats, aes(x = mid, y = Fst_morm_WA.silver_morm_WA.buff)) +
  geom_point(aes(color = outlier), pch = 19, size = 1, alpha = 0.5) + # Map 'outlier' to color
  scale_color_manual(values = c("outlier" = "red", "non-outlier" = "grey")) + # Define colors for outliers and non-outliers
  geom_smooth(method = "loess", se = FALSE, color = "black", span = 0.05, size = 0.8) +
  theme_classic() + 
  xlim(1480000, 1550000) +
  ylim(0, 1) +
  theme(
    legend.position = "none", 
    axis.title.y = element_text(size = 8), 
    axis.text.y = element_text(size = 6), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank()
  ) +
  ylab("FST") + 
  xlab("Position (bp)")

# Calculate 99% quantile threshold for dxy
threshold <- quantile(file.popgen.stats$dxy_morm_WA.silver_morm_WA.buff, 0.95, na.rm = TRUE)


# Create a new column 'outlier' to flag points above the threshold
file.popgen.stats$outlier <- ifelse(file.popgen.stats$dxy_morm_WA.silver_morm_WA.buff > threshold, "outlier", "non-outlier")


p2 <- ggplot(data = file.popgen.stats, aes(x = mid, y = dxy_morm_WA.silver_morm_WA.buff)) +
  geom_point(aes(color = outlier), pch = 19, size = 1, alpha = 0.5) + # Map 'outlier' to color
  scale_color_manual(values = c("outlier" = "red", "non-outlier" = "grey")) + # Define colors for outliers and non-outliers
  geom_smooth(method = "loess", se = FALSE, color = "black", span = 0.05, size = 0.8) +
  theme_classic() + 
  xlim(1480000, 1550000) +
  ylim(0, 0.26) +
  theme(
    legend.position = "none", 
    axis.title.y = element_text(size = 8), 
    axis.text.y = element_text(size = 6), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank()
  ) +
  ylab("Dxy") + 
  xlab("Position (bp)")

#Alternative plot as buff subtracted from silver Pi with outliers

# Generate the distribution
pi.buff.minus.pi.silver <- file.popgen.stats$pi_morm_WA.buff - file.popgen.stats$pi_morm_WA.silver

# Calculate 99th and 1th percentile thresholds
threshold_upper <- quantile(pi.buff.minus.pi.silver, 0.99, na.rm = TRUE)
threshold_lower <- quantile(pi.buff.minus.pi.silver, 0.01, na.rm = TRUE)

# Flag points as outliers if they are above the 95th percentile or below the 5th percentile
file.popgen.stats$outlier.delta.pi <- ifelse(pi.buff.minus.pi.silver > threshold_upper | 
                                               pi.buff.minus.pi.silver < threshold_lower, 
                                             "outlier", "non-outlier")


# Plot with outliers colored
p3 <- ggplot(data = file.popgen.stats) +
  geom_point(aes(x = mid, y = pi_morm_WA.buff - pi_morm_WA.silver, color = outlier.delta.pi), pch = 19, size = 1, alpha = 0.5) +
  geom_smooth(aes(x = mid, y = pi_morm_WA.buff - pi_morm_WA.silver), method = "loess", se = FALSE, color = "black", span = 0.05, size = 0.8) +
  theme_classic() + 
  scale_color_manual(values = c("outlier" = "red", "non-outlier" = "grey")) +  # Color outliers in red, others in black
  xlim(1480000, 1550000) +
  ylim(-0.15, 0.15) +
  theme(legend.position = "none", axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ylab("Δπ") + 
  xlab("Position (bp)")

#Make the plot for TajD for silver and buff haplotypes calculating outliers for their respective distributions

# Calculate 99th and 1st percentiles for TajD_morm_WA.silver
threshold_upper_silver <- quantile(file.popgen.stats.freq$TajD_morm_WA.silver, 0.99, na.rm = TRUE)
threshold_lower_silver <- quantile(file.popgen.stats.freq$TajD_morm_WA.silver, 0.01, na.rm = TRUE)

# Calculate 99th and 1st percentiles for TajD_morm_WA.buff
threshold_upper_buff <- quantile(file.popgen.stats.freq$TajD_morm_WA.buff, 0.99, na.rm = TRUE)
threshold_lower_buff <- quantile(file.popgen.stats.freq$TajD_morm_WA.buff, 0.01, na.rm = TRUE)

# Create separate flags for outliers in silver and buff
file.popgen.stats.freq$outlier_silver <- ifelse(
  file.popgen.stats.freq$TajD_morm_WA.silver > threshold_upper_silver | 
    file.popgen.stats.freq$TajD_morm_WA.silver < threshold_lower_silver, 
  "outlier", "non-outlier"
)

file.popgen.stats.freq$outlier_buff <- ifelse(
  file.popgen.stats.freq$TajD_morm_WA.buff > threshold_upper_buff | 
    file.popgen.stats.freq$TajD_morm_WA.buff < threshold_lower_buff, 
  "outlier", "non-outlier"
)

# Plot with full opacity for outliers and transparency for non-outliers, separately for silver and buff
p4 <- ggplot(data = file.popgen.stats.freq) +
  # Silver points
  geom_point(aes(x = mid, y = TajD_morm_WA.silver, shape = outlier_silver, alpha = outlier_silver), size = 1, color = "#1baae2") +
  geom_smooth(aes(x = mid, y = TajD_morm_WA.silver), method = "loess", se = FALSE, color = "#1baae2", span = 0.05, size = 0.8) +
  # Buff points
  geom_point(aes(x = mid, y = TajD_morm_WA.buff, shape = outlier_buff, alpha = outlier_buff), size = 1, color = "#f89b31") +
  geom_smooth(aes(x = mid, y = TajD_morm_WA.buff), method = "loess", se = FALSE, color = "#f89b31", span = 0.05, size = 0.8) +
  # Set shapes and alpha for outliers and non-outliers
  scale_shape_manual(values = c("outlier" = 19, "non-outlier" = 1)) +  # Full circle (19) for outliers, empty circle (1) for non-outliers
  scale_alpha_manual(values = c("outlier" = 1, "non-outlier" = 0.3)) + # Full opacity for outliers, transparency for non-outliers
  theme_classic() + 
  xlim(1480000, 1550000) +
  theme(legend.position = "none", axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ylab("TajD Silver / Buff") + 
  xlab("Position (bp)")

  #Make the plot for CLR for silver and buff haplotypes  

p5 <- ggplot(data = SDB) +
  geom_line(aes(x = Position.b, y = Likelihood.b), pch = 19, size = 1, color = "#f89b31", alpha = 0.7) +
  geom_line(aes(x = Position.s, y = Likelihood.s), pch = 19, size = 1, color = "#1baae2", alpha = 0.4) +
  theme_classic() + 
  xlim(1480000, 1550000) +
  ylim(0,100) +
  theme(legend.position = "none", axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.text.x = element_blank() ) +
  ylab("CLR Silver / Buff") + 
  xlab("Position (bp)")

############################## input files for twisst ######################################

source("plot_twisst.R")

#weights file with a column for each topology
weights_file <- "SmorWA.haplo.weights.tsv"

#coordinates file for each window
window_data_file <- "Chrm14.optix.filtered.all.trees.haplo.data.tsv"


################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)


############################## combined plots ##################################
# there are a functions available to plot both the weightings and the topologies

# make smooth weightings and plot those across chromosomes
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 23000, spacing = 1000)

# Load required libraries
library(reshape2)

# Extract weights from twisst_data_smooth
weights_df <- as.data.frame(twisst_data_smooth$weights[[1]])

# Add a region column to represent x-axis (for example, region number)
weights_df$region <- twisst_data_smooth$pos[[1]]

# Multiply topo1 and topo2 by -1 to invert them
weights_df$topo1 <- -weights_df$topo1
weights_df$topo2 <- -weights_df$topo2

# Melt the data frame to long format for ggplot
weights_long <- melt(weights_df, id.vars = "region", variable.name = "topology", value.name = "weight")

# Plot using ggplot

topo_colors <- c("topo1" = "blue", "topo2" = "cyan", "topo3" = "red")


p6 <- ggplot(weights_long, aes(x = region, y = weight, fill = topology)) +
  geom_area(position = "identity", alpha = 0.3) +
  scale_fill_manual(values = topo_colors) +
  labs(x = "Position (bp)", y = "Weight") +
  theme_classic() +
  xlim(1480000, 1550000) +
  theme(legend.position = "none", axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 8), axis.text.x = element_text(size = 6) )
  

##ABBA-BABA-fd###
file.popgen.stats <- read.csv("ABBABABBA.csv", h=T)
head(file.popgen.stats)

#fd is meaningless when D is negative, as it is designed to quantify the excess of ABBA over BABA only whgen an excess exists.
#We therefore convert all fd values to 0 at sites where D is negative.

file.popgen.stats$fd[file.popgen.stats$D < 0] <- 0


p7 <- ggplot(data = file.popgen.stats) +
  geom_point(aes(x = mid, y = fd), pch = 19, size = 1, color = "red", alpha = 0.4) +
  theme_classic() + 
  xlim(1480000, 1550000) +
  ylim(0, 1) +
  ylab("fd hyd to morm") +
  theme(legend.position = "none", axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.text.x = element_blank() ) +
  xlab("Position (bp)")



p0 <- ggplotGrob(p0)
p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)
p5 <- ggplotGrob(p5)
p6 <- ggplotGrob(p6)
p7 <- ggplotGrob(p7)


maxWidth = grid::unit.pmax(p0$widths[2:5], p1$widths[2:5], p2$widths[2:5], p3$widths[2:5], p4$widths[2:5], p5$widths[2:5], p6$widths[2:5], p7$widths[2:5])
p0$widths[2:5] <- as.list(maxWidth)
p1$widths[2:5] <- as.list(maxWidth)
p2$widths[2:5] <- as.list(maxWidth)
p3$widths[2:5] <- as.list(maxWidth)
p4$widths[2:5] <- as.list(maxWidth)
p5$widths[2:5] <- as.list(maxWidth)
p6$widths[2:5] <- as.list(maxWidth)
p7$widths[2:5] <- as.list(maxWidth)

grid.arrange(p0, p1, p2, p3, p4, p5, p7, p6, ncol = 1, heights = c(2.5, 10 ,10, 10, 10, 10, 10, 15))
```
