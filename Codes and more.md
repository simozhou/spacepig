Reference for parallelism:  [GATK PARALLEL](https://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster)

Previous steps:

- Correct wrong header format

- Sort

- Index

To correct for wrong header format, samtools reheader can be used, flanked by a sed command.

```bash
samtools reheader -c "sed '/^@RG/ s/\([A-Z][A-Z]\)ID/\1/g'" Normal.bam > Normal.RGcorrected.bam
samtools reheader -c "sed '/^@RG/ s/\([A-Z][A-Z]\)ID/\1/g'" Tumor.bam > Tumor.RGcorrected.bam
```


# Task #1


```bash
samtools stats Normal.RGcorrected.sorted.bam | grep ^SN | cut -f 2- > Normalbamstats.txt
samtools stats Tumor.RGcorrected.sorted.bam | grep ^SN | cut -f 2- > Tumorbamstats.txt
```

# Task #2

##### Target creation  #2.1.1

Command for realigner target creator, on selected regions of the whole exome sequencing experiment

```bash
java -jar ../../Tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../../Annotations/human_g1k_v37.fasta -I Normal.RGcorrected.sorted.bam -o Normal.realigner.intervals -L Haloplex_Regions.chr16.17.19.bed
java -jar ../../Tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../../Annotations/human_g1k_v37.fasta -I Tumor.RGcorrected.sorted.bam -o Tumor.realigner.intervals -L Haloplex_Regions.chr16.17.19.bed     
```

##### Target realigner #2.1.2

```bash
java -jar ../../Tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../../Annotations/human_g1k_v37.fasta -I Normal.RGcorrected.sorted.bam -targetIntervals Normal.realigner.intervals -o Normal.RGcorrected.sorted.realigned.bam -L Haloplex_Regions.chr16.17.19.bed
java -jar ../../Tools/GenomeAnalysisTK.jar -T IndelRealigner -R ../../Annotations/human_g1k_v37.fasta -I Tumor.RGcorrected.sorted.bam -targetIntervals Tumor.realigner.intervals -o Tumor.RGcorrected.sorted.realigned.bam -L Haloplex_Regions.chr16.17.19.bed
```

###### Quantify realinged reads #2.1.3

How to determine realigned reads?

- Grep for the old CIGAR indicator ("OC:") in the additional fields of the realigned bam file

```bash
samtools view Normal.RGcorrected.sorted.realigned.bam | grep -c '\sOC:' > Normal.realigned.reads.number
samtools view Tumor.RGcorrected.sorted.realigned.bam | grep -c '\sOC:' > Tumor.realigned.reads.number
```

##### Recalibration #2.2

For the recalibration we need first a table.

```bash
java -jar    ../../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R  ../../Annotations/human_g1k_v37.fasta -I Normal.RGcorrected.sorted.realigned.bam -knownSites ../../Annotations/hapmap_3.3.b37.vcf -o Normal.recal.table -L Haloplex_Regions.chr16.17.19.bed
java -jar    ../../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R  ../../Annotations/human_g1k_v37.fasta -I Tumor.RGcorrected.sorted.realigned.bam -knownSites ../../Annotations/hapmap_3.3.b37.vcf -o Tumor.recal.table -L Haloplex_Regions.chr16.17.19.bed
```

Now we have to recalibrate our bam file. PrintReads supports nct as the previous step.

```bash
java -jar    ../../Tools/GenomeAnalysisTK.jar -T PrintReads -R  ../../Annotations/human_g1k_v37.fasta -I Normal.RGcorrected.sorted.realigned.bam -BQSR Normal.recal.table -o Normal.RGcorrected.sorted.realigned.recal.bam -L Haloplex_Regions.chr16.17.19.bed
java -jar    ../../Tools/GenomeAnalysisTK.jar -T PrintReads -R  ../../Annotations/human_g1k_v37.fasta -I Tumor.RGcorrected.sorted.realigned.bam -BQSR Tumor.recal.table -o Tumor.RGcorrected.sorted.realigned.recal.bam -L Haloplex_Regions.chr16.17.19.bed
```

And to produce before and after plot. To do so we need to use again BaseRecalibrator with ```-BQSR``` option to produce the after recalibration table.

```bash
java -jar    ../../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R  ../../Annotations/human_g1k_v37.fasta -I Normal.RGcorrected.sorted.realigned.bam -knownSites ../../Annotations/hapmap_3.3.b37.vcf -BQSR Normal.recal.table -o Normal.after.recal.table -L Haloplex_Regions.chr16.17.19.bed
java -jar    ../../Tools/GenomeAnalysisTK.jar -T BaseRecalibrator -R  ../../Annotations/human_g1k_v37.fasta -I Tumor.RGcorrected.sorted.realigned.bam -knownSites ../../Annotations/hapmap_3.3.b37.vcf -BQSR Tumor.recal.table -o Tumor.after.recal.table -L Haloplex_Regions.chr16.17.19.bed
```
Eventually, before/after plots are generated (**Requires R<4.0.0**)

```bash
java -jar    ../../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R  ../../Annotations/human_g1k_v37.fasta -before Normal.recal.table -after Normal.after.recal.table -plots Normal.recal.plots.pdf
java -jar    ../../Tools/GenomeAnalysisTK.jar -T AnalyzeCovariates -R  ../../Annotations/human_g1k_v37.fasta -before Tumor.recal.table -after Tumor.after.recal.table -plots Tumor.recal.plots.pdf
```

# Task #3

In order to perform the variant calling we've been asked to use GATK UnifiedGenotyper

The -L flag leads to a focus on a specific area

```bash
java -jar ../../Tools/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../../Annotations/human_g1k_v37.fasta -I Normal.RGcorrected.sorted.realigned.recal.bam -o Normal.snps.vcf -L Haloplex_Regions.chr16.17.19.bed
```

### Genotype identification ---> heterozygosity

The **-GF** field gives us information about the genotype. ```bcftools filter``` allows us to extract heterozygous genotypes.

```bash
bcftools filter -i 'GT="0/1" | GT="0|1"' -o Normal.snps.het.vcf Normal.snps.vcf
```


### Intersection with clinvar file

Easy peasy with bcftools isec, but it needs BGZF files, so must bgzip both **Normal.snps.het.vcf** and **clinvar.vcf**
Plus, bcftools needs also an index for each BGZF file.
Plus, **clinvar.vcf** was provided without header, so we found one and we added it to the top of the file.

<u>We may also apply a filter for quality to our results before intersecting to see how many variants we loose in the final result</u>

```bash
bgzip Normal.snps.het.vcf
bgzip clinvar.vcf
bcftools index Normal.snps.het.vcf.gz
bcftools index clinvar.vcf.gz
```

After this preprocessing step, bcftools isec can be run

```bash
bcftools isec -p clinvar_intersection Normal.snps.het.vcf.gz clinvar.vcf.gz
```

Eventually, a BGZF, correctly named is produced from the intersection Normal.snps.het.vcf.gz --> clinvar.vcf.gz. Plus, an index file is produced.


```bash
cp clinvar_intersection/0002.vcf ./Normal.snps.het_clinvar_intersection.vcf
bgzip Normal.snps.het_clinvar_intersection.vcf
bcftools index Normal.snps.het_clinvar_intersection.vcf.gz
```

# Task #4

We can implement this task in this way:

1. We write the sample name into a file in this way 
   
   ```bash
   echo 'Normal.RGcorrected.sorted.realigned.recal.bam' > BAMs_List.txt
   ```

2. Use provided R script
   
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
The R script produces a table which can be visualized on IGV to isolate segments with heterozigous deletions.