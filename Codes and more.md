Reference for parallelism:  [GATK PARALLEL](https://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster)

Previous steps:

- Correct wrong header format

- Sort

- Index

- Deduplicate

Deduplication is needed for an increased accuracy in variant calling, that will have to be performed later on. **But when do we want to perform deduplication?**



A main point that we have to implement inside our pipeline is the proper determination of the tumoral bam and the normal bam. Because in the task of the ancestry only the normal bam is needed.


To correct for wrong header format, samtools reheader can be used, flanked by a sed command. Example for Normal.bam

```bash
samtools reheader -c "sed '/^@RG/ s/\([A-Z][A-Z]\)ID/\1/g'" Normal.bam > Normal.RGcorrected.bam
```


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

How to determine realigned reads?

- Grep for the old CIGAR indicator ("OC:") in the additional fields of the realigned bam file

```{bash}
samtools view realigned.bam | grep -c '\sOC:'
```

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

The **-GF** field gives us information about the genotype. I don't think we should really care about.

I suppose we don't need to check for phased SNPs because they're relative to the haplotype.

This command, taken from [biostars](https://www.biostars.org/p/166260/) 

```bash
vcftools --gzvcf file.vcf.gz --extract-FORMAT-info GT | grep "0/1"
```

Is about regardless of the phasing. Given that we don't have paternal and maternal SNPs Maybe they're not supposed to be "|" separated.

<u>**Haplotype caller from GATK???**</u>

May be a possibility and it has multiple uses.



### Difference With clinvar file

Easy peasy with vcftools

```bash
vcftools --vcf Sample.BCF.recode.vcf --diff Sample.GATK.recode.vcf –diff-site

```

# Task #4

We can implement this task in this way:

1. We write the sample name into a file in this way 
   
   ```bash
   echo 'Normal.bam' > BAMs_List.txt
   ```

2. Use his script but **we have to install his package**
   
   ```bash
   Rscript RunEthSEQ.R
   ```

From which we can retrieve the report


