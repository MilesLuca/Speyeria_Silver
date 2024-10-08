# Speyeria_Silver
A collection of code used for our study detailing the genomic basis for the silvering polymorphism in Speyeria mormonia.

**1. Quality control - fastqc**
TBA

**2. Alignment and genotype calling**

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


