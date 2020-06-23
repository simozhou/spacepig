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

##### Target creation  #2.1.1

This is an example for Normal.bam and realigner target creator, on selected regions of the whole exome sequencing experiment

```bash
java –jar GenomeAnalysisTK.jar –T RealignerTargetCreator  –R human.fasta –I    Normal.bam    –known    indels.vcf –o realigner.intervals -L Haloplex_Regions.chr16.17.19.bed     
```

##### Target realigner #2.1.2

There's the possibility of using the smith waterman algorithm, but it's computationally expensive for what regards the realignment.

For what regards multitheading IndelRealigner has a scatter and gather approach.

```{bash}
java    
 –jar    
 GenomeAnalysisTK.jar    –T    IndelRealigner    –R    human.fasta    –I    original.bam    –known    indels.vcf    –targetIntervals    realigner.intervals     –o    realigned.bam    -L Haloplex_Regions.chr16.17.19.bed
```

For these two commands obviously the bam files become normal and tumor. The fasta file is the human genome file found on the annotation folder. Called **<mark>human_g1k_v37.fasta</mark>**

###### Quantify realinged reads #2.1.3

How to determine realigned reads?

- Grep for the old CIGAR indicator ("OC:") in the additional fields of the realigned bam file

```{bash}
samtools view realigned.bam | grep -c '\sOC:'
```

##### Recalibration #2.2

For the recalibration we need first a table.

For this command we can assign `-nct / --num_cpu_threads_per_data_thread`  controls the number of **CPU threads** allocated to each data thread.

Obvously we have to insert the whole genome

```bash
java    –jar    GenomeAnalysisTK.jar  –T     BaseRecalibrator    –R  human.fasta     –I  realigned.bam  –knownSites     dbsnp137.vcf       –knownSites gold.standard.indels.vcf –o    recal.table -L Haloplex_Regions.chr16.17.19.bed   
```

Now we have to recalibrate our bam file. PrintReads supports nct as the previous step.

```bash
java    –jar GenomeAnalysisTK.jar     –T    PrintReads    –R    human.fasta    –I     realigned.bam    –BQSR    recal.table    –o    recal.bam  -L Haloplex_Regions.chr16.17.19.bed
```

And to produce before and after plot. To do so we need to use again BaseRecalibrator to produce the after recalibration table.

```bash
java    –jar    GenomeAnalysisTK.jar  –T     BaseRecalibrator    –R  human.fasta     –I  realigned.bam  –knownSites     dbsnp137.vcf       –knownSites gold.standard.indels.vcf  -BQSR recal.table  –o    after_recal.table -L Haloplex_Regions.chr16.17.19.bed   
```

```bash
java    –jar    GenomeAnalysisTK.jar    –T    AnalyzeCovariates    –R    human.fasta    –before    recal.table    –after    after_recal.table    –plots    recal_plots.pdf    
```

# Task #3

### {{{[[[(((

In order to perform the variant calling we've been asked to use GATK UnifiedGenotyper

1.  GATK without the -L flag

2.  bcftools

### SNPs identification

#### bcftools code -->  directly to vcf

bcftools has also a **--threads** *INT*  option to define the threads

```bash
bcftools mpileup -Ou -a DP -f human_g1k_v37.fasta Sample.sorted.bam |
bcftools call -Ov -c -v > Sample.BCF.vcf
```
### )))]]]}}}
#### gatk code without the -L

In order to perform the variant calling we've been asked to use GATK UnifiedGenotyper

The -L flag leads to a focus on a specific area

```bash
java -jar GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R human_g1k_v37.fasta -I Sample.sorted.bam -o Sample.GATK.vcf -L Haloplex_Regions.chr16.17.19.bed
```

### Genotype identification ---> heterozygosity

The **-GF** field gives us information about the genotype. I don't think we should really care about.

I suppose we don't need to check for phased SNPs because they're relative to the haplotype.

Example for Sample.GATK.vcf

```bash
grep "0[/\|]1" Sample.GATK.vcf > Sample.GATK.Het.vcf
```

Is about regardless of the phasing. Given that we don't have paternal and maternal SNPs Maybe they're not supposed to be "|" separated.

<u>**Haplotype caller from GATK???**</u>

May be a possibility and it has multiple uses.



### Intersection With clinvar file

Easy peasy with vcftools

We may also apply a filter for quality to our results before intersecting to see how many variants we loose in the final result

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

# Task #5
To get raw CNV results from Normal.bam and Tumor.bam
```bash
samtools mpileup -q 1 -f ../../Annotations/human_g1k_v37.fasta Normal.bam Tumor.bam | java -jar ../../Tools/VarScan.v2.3.9.jar copynumber --output-file SCNA --mpileup 1
```

To adjust CNV data for GC content sequencing bias
```bash
java -jar ../../Tools/VarScan.v2.3.9.jar copyCaller SCNA.copynumber --output-file SCNA.copynumber.called
```

To retrieve and visualize CNV segments. Should modify the R script to work in the right directories internally or to 
accept them as input from shell.
```bash
Rscript CBS.R
```
