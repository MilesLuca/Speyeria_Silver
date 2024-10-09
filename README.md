# Speyeria_Silver
A collection of code used for our study detailing the genomic basis for the silvering polymorphism in Speyeria mormonia.

# 1. Quality control - fastqc

TBA

# 2. Alignment and genotype calling

The following code is used to generate genome indexes, align sequencing data, mark sequencing duplicates, call genotypes and produce vcf files for downstream analyses.

First Index the genome for gatk, bwa and samtools:

```
gatk CreateSequenceDictionary -R Smor.assembly.v1.0.unmasked.fasta
bwa index Smor.assembly.v1.0.unmasked.fasta
samtools faidx Smor.assembly.v1.0.unmasked.fasta
```

Alignemts and genotype calling are performed using a job array.
Prep a file for job array to grab sample IDs from
Save file called “jobs” with following 52 IDs in (corresponding to each sample):

```
nano jobs

Shyd_WA_05 Shyd_WA_06 Shyd_WA_08 Shyd_WA_09 Shyd_WA_10 Smor_WA_01 Smor_WA_02 Smor_WA_07 Smor_WA_11 Smor_WA_13 Smor_WA_14 Smor_WA_15 Smor_WA_16 Smor_WA_17 Smor_WA_18 Smor_WA_19 Smor_WA_20 Smor_WA_21 Smor_WA_22 Smor_WA_23 Smor_WA_24 Smor_WA_25 Smor_WA_26 Smor_WA_27 Smor_WA_28 Smor_WA_29 Smor_WA_30 Smor_WA_31 Smor_WA_32 Smor_WA_33 Smor_WA_34 Smor_WA_35 Smor_WA_37 Smor_WA_38 Smor_WA_39 Smor_WA_40 Smor_WA_41 Smor_WA_42 Smor_WA_43 Smor_WA_44 Smor_WA_45 Smor_WA_46 Smor_WA_47 Smor_WA_48 NVJ_51_Smor_B NVJ_52_Smor_S NVJ_53_Smor_B NVJ_54_Smor_B NVJ_55_Smor_S NVJ_56_Smor_S NVJ_57_Smor_B NVJ_58_Smor_B NVJ_59_Smor_B NVJ_60_Smor_S NVJ_61_Smor_B NVJ_62_Smor_S NVJ_63_Smor_S NVJ_64_Smor_B NVJ_65_Smor_S NVJ_66_Smor_S NVJ_67_Smor_S NVJ_68_Smor_B NVJ_69_Smor_B NVJ_70_Smor_B NVJ_71_Smor_B NVJ_72_Smor_S NVJ_73_Smor_B NVJ_74_Smor_B NVJ_75_Smor_S NVJ_76_Smor_B NVJ_77_Smor_S NVJ_78_Smor_S NVJ_79_Smor_S NVJ_80_Smor_S NVJ_81_Smor_S
```

Make a directory to store the alignment statistics in

```
mkdir alnstats
```

Generate a script to run job as array (75 samples total)
This script is run in two parts

The first is a script called: 

makeGVCF.array.sh

Part A (**makeGVCF.array.sh**): Call variants on each individual against the genome.
1. align using bwa-mem (preferred aligner for genotyping)
2. Mark duplicate reads (using a gatk-instance of picard's MarkDuplicates)
3. Ensure that sensible names are assigned to each individual with "AddOrReplaceReadGroups" (this is important when merging the files in part B).
4. produce the gvcf (genotype variant call file), one per individual. This file will score every site in the genome as either non-variant or variant.

The second is a script called: 

genotypeGVCFs.sh

Part B (**genotypeGVCFs.sh**): Merge gVCFs, call genotypes and produce the VCF, and then filter the VCF. 
1. Run the gatk tool CombineGVCFs to merge all the gVCFs into one (very large) gvcf. This file includes a separate entry for every site in every individual so it is upsettingly large, but it runs quickly
2. Run the gatk tool GenotypeGVCFs on the combined gVCF. This is the 'genotype calling' step. This program goes through each site in the genome, looks at the variants from the gVCF for each individual, and calls the genotype at each site. The output is one row per variant, one column per individual (ie, a VCF).
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
#SBATCH --array=0-51%10

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

# 3. Alignment QC (depth measurments)

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

For a slightly more accurate (but way slower method) of calculating genome-wide coverage you can use Bedtools genomecov as follows below.

First generate a sizes.genome file with all chromosomes and their sizes:

```
samtools faidx genome.fasta
cut -f1,2 genome.fa.fai > sizes.genome
```

Then run bedtools for calculating coverage:

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

This generates pretty big files, which you can then awk to get an average for each position in the genome:

```
awk '{sum += $3} END {print "Average coverage:", sum/NR}' SRR25297463_coverage_output.txt
```

Outputs from mosdepth and bedtools are very similar, so probably not worth running bedtools.

# 4. VCF filtering

From here on, we are using the SNP file generated by gatk, called "Smor.snps_filtered.vcf.gz"

Check for missing data using vcftools:

```
vcftools --gzvcf Smor.snps_filtered.vcf.gz --missing-indv --out data
```

No evidence of high missingness in any of the samples. 
Slightly higher in outgroups as expected:

INSERT IMAGE IF POSSIBLE

Filter with bcftools (removing minor allele freq of 2% and missing genotypes of 10%)

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

Next phase and impute missing genotypes using beagle

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

Purge non-variant sites and keep only bi-allelic sites

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

Generate plink files and run GWAS using GEMMA
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

Point plink to the correct files

```
VCF=Smor.snps_filtered.maf.geno.beagled.WA.biallele.nonvar.vcf
PHENO=PhenotypeFile.txt
```

Generate plink files

```
plink --vcf ${VCF} --double-id --allow-extra-chr --allow-no-sex --set-missing-var-ids @:# --pheno $PHENO --make-bed --pca --out Smor_gemma_assoc
```

Generate relatedness matrix

```
gemma -bfile Smor_gemma_assoc -gk 1 -o Smor_relatedness.out > gemma.Smor_relatedness.out
```
Run association, using relatedness matrix as covariate

```
gemma -bfile Smor_gemma_assoc -lmm 1 -o Smor_silver.assoc.gemma -k Smor_relatedness.out.cXX.txt > Smor_silver.assoc.gemma.out
```

Extract p-values

```
awk '{print $1,$3,$12}' Smor_silver.assoc.gemma.assoc.txt > Smor_silver.assoc.gemma.assoc.txt.pvals
```

Check minimum p-vale

```
awk 'min=="" || $3 < min {min=$3; minline=$0}; END{ print min, minline}' *pvals
```

Optional: Reduce file size by filtering for min of 0.05

```
awk '{ if($3 <= 0.05) { print }}' *pvals > Smor_silver.assoc.gemma.pval.noNA.0.05
```

5. Import p-values and chromosome positions in R and plot them

Functions to adjust scaffold positions

Calculate chromosome coordinates:

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

Calculate scaffold coordinates

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

Specify PDF output with desired size 

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

Read in the file generated by GEMMA

```
fileX <- read.table("Smor_silver.assoc.gemma.assoc.txt.pvals", h=T)
```

Check file

```
head(fileX)
```

Filter for more signifanct snps to make plotting easier and omit NA values

```
fileX <- na.omit(fileX) 
min(fileX$P)
nrow(fileX[fileX$P < 0.005,])
fileX <- fileX[fileX$P < 0.005,]
nrow(fileX)
```

Make your life easier by having the CHROM / scaffold column called ?scaffold?, and your position called "pos". If the names don't match, it can be renamed like so:

```
names(fileX)[names(fileX) == 'CHR'] <- 'scaffold'

fileX_merged <- merge(fileX, scaf_coords2,by = "scaffold", all.x=TRUE)
fileX_merged_pos <- cbind(fileX_merged, chromPos = fileX_merged$BP + fileX_merged$scafStart + fileX_merged$chromStarts-2)
```

#Now here is the plotting code. 

top and bot set the y axis limits
in "begin" and "end", set the number of chromosmes. 

plot your stats, there is a variable to be plotted called -log10(comp[[i]]$P) , which is taking the column $P from fileX_merged_pos and doing a log10 transform on it.

Your setup

```
comp=list(fileX_merged_pos)
names=c("fileXX_name")
col=c("darkgrey")

par(mfrow=c(1,1), mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 18
bot = 2.5
```

#Check how many chromosomes you have

```
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[30]/1000000
```

Plot background

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







