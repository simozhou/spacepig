Reference for parallelism:  [GATK PARALLEL](https://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster)

Previous steps:

- Sort

- Index

- Deduplicate

Deduplication is needed for an increased accuracy in variant calling, that will have to be performed later on. **But when do we want to perform deduplication?**

# Task #1

Easy peasy, example for Normal.bam

```bash
samtools stats Normal.bam | grep ^SN | cut -f 2- > bamstats.txt
```

# Task #2

Realignment command for both BAM files

##### Target creation  #2.1

This is an example for Normal.bam and realigner target creator

```bash
java –jar GenomeAnalysisTK.jar –T RealignerTargetCreator  –R human.fasta –I    Normal.bam    –known    indels.vcf –o realigner.intervals    
```

##### Target realigner #2.2

There's the possibility of using the smith waterman algorithm, but it's computationally expensive for what regards the realignment.

For what regards multitheading IndelRealigner has a scatter and gather approach.

```{bash}
java    
 –jar    
 GenomeAnalysisTK.jar    –T    IndelRealigner    –R    human.fasta    –I    original.bam    –known    indels.vcf    –targetIntervals    realigner.intervals     –o    realigned.bam    
```

For these two commands obviously the bam files become normal and tumor. The fasta file is the human genome file found on the annotation folder. Called **<mark>human_g1k_v37.fasta</mark>**

There are 

###### Task point #1

How to determine aligned reads?

- Grep for the aligned (@OC)

- Command not known

##### Recalibration #2.3

For the recalibration we need first a table.

For this command we can assign `-nct / --num_cpu_threads_per_data_thread`  controls the number of **CPU threads** allocated to each data thread.

Obvously we have to insert the whole genome

```bash
java    –jar    GenomeAnalysisTK.jar  –T     BaseRecalibrator    –R  human.fasta     –I  realigned.bam  –knownSites     dbsnp137.vcf       –knownSites gold.standard.indels.vcf –o    recal.table    
```

Now we have to recalibrate our bam file. PrintReads supports nct as the previous step.

```bash
java    –jar GenomeAnalysisTK.jar     –T    PrintReads    –R    human.fasta    –I     realigned.bam    –BQSR    recal.table    –o    recal.bam    
```

And to produce before and after plot. We have to use the tables in such code

```bash
java    –jar    GenomeAnalysisTK.jar    –T    AnalyzeCovariates    –R    human.fasta    –before    recal.table    –after    after_recal.table    –plots    recal_plots.pdf    
```

# Task #3

In order to perform the variant calling we can use

1.  GATK without the -L flag

2.  bcftools

### SNPs identification

#### bcftools code -->  directly to vcf

bcftools has also a **--threads** *INT*  option to define the threads

```bash
bcftools mpileup -Ou -a DP -f human_g1k_v37.fasta Sample.sorted.bam |
bcftools call -Ov -c -v > Sample.BCF.vcf
```

#### gatk code without the -L

The -L flag leads to a focus on a specific area

```bash
ava -jar GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R human_g1k_v37.fasta
-I Sample.sorted.bam -o Sample.GATK.vcf 
```

### Genotype identification ---> heterozygosity




