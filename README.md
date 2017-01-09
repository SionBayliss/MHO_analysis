# MHO 001 MinION Analysis 

### Data Availability
The dataset supporting the conclusions of this article is available in the European Nucleotide Archive repository under project number PRJEB14152. Further supporting data is also available from the GigaScience GigaDB repository:

Bayliss SC, Hunt VL, Yokoyama M, Thorpe HA, Feil EJ. Supporting data for "The use of Oxford Nanopore native barcoding for complete genome assembly". GigaScience Database. 2016. http://dx.doi.org/10.5524/100269.


### Availability and requirements
*Project name: MHO_001 hybrid read assembly and analysis
*Project home page: https://github.com/SionBayliss/MHO_analysis
*Operating systems: Unix
*Programming language: R, perl
*Other requirements: Dependencies include Samtools (>=1.18), Trimmomatic, SPAdes v3.6.1, BWA (0.7.5a-r405), BioPerl, MAUVE, BLASR, prokka, Tablet/Artemis, (R Libraries: ggplot2, cowplot, genoPlotR)
*License: GNU GPL v3

# Workflow
For publication - Gigascience
  
### Set working directory (change as appropriate) 
```
DIR="/path/to/working/directory/"
cd $DIR
```

#### Download GitHub Directory
```
git clone http://githum.com/SionBayliss/MHO_analysis.git
```
### Download Reads
#### MinION Reads (fasta format)
ENA read accession: ERS1178418/ERR1424936
```
mkdir Reads
cd Reads
mkdir minION
cd minION
mkdir Raw
cd Raw
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA628/ERA628229/oxfordnanopore_native/MHO_001.tar.gz
tar xvzf MHO_001.tar.gz
mv MH0_001.pass.fasta reads_pass.fasta
mv MHO_001.fail.fasta reads_fail_raw.fasta 
rm MHO_001.tar.gz
```
#### Illumina Reads
ENA read accessions: ERS1180806 and ERS1180807
```
cd $DIR/Reads
mkdir Illumina; cd Illumina
mkdir Raw; cd Raw
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA630/ERA630327/fastq/2998-174_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA630/ERA630327/fastq/2998-174_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA630/ERA630327/fastq/2998-174_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA630/ERA630327/fastq/2998-174_2.fastq.gz
```

#### Set read paths
```
read1_1="$DIR/Reads/Illumina/Raw/2998-174_1.fastq.gz"
read2_1="$DIR/Reads/Illumina/Raw/2998-174_1._hiseq.fastq.gz"
read1_2=$DIR/Reads/Illumina/Raw/2998-174_2.fastq.gz
read2_2=$DIR/Reads/Illumina/Raw/2998-174_2.hiseq.fastq.gz
```
#### Combine Illumina reads for Miseq/Hiseq runs
Note: two Illumina sequencing runs were performed for the same library. The reads were combined before further processing. 
```
zcat $read1_1 $read2_1 | gzip - > $DIR/Reads/Illumina/Raw/illumina_1.fastq.gz
zcat $read1_2 $read2_2 | gzip - > $DIR/Reads/Illumina/Raw/illumina_2.fastq.gz
```

Assess quality using fastqc 

```
fastqc -t 2 $DIR/Reads/Illumina/Raw/illumina_1.fastq.gz $DIR/Reads/Illumina/Raw/illumina_2.fastq.gz
```
#### Note : Convert pass and fail reads to fasta format.
The following analysis assumes you have downloaded the fasta files from the ENA. The raw native MinION data is also available and can be converted to fasta format using PoRe, poRetools or conversion tool of your choice.

## QC, demultiplex and trim reads


#### Trim Illumina Reads using Trimmomatic.
```
cd $DIR/Reads/Illumina/
mkdir Trimmed
trimmomatic PE -threads 4 Raw/illumina_1.fastq.gz Raw/illumina_2.fastq.gz Trimmed/illumina_1.paired.fastq.gz Trimmed/illumina_1.unpaired.fastq.gz Trimmed/illumina_2.paired.fastq.gz Trimmed/illumina_2.unpaired.fastq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:20
```
NOTE - This wasn't a great library and a lot of reads are lost with even these quite permissive parameters. 
```
Results
Input Read Pairs: 713279 Both Surviving: 439480 (61.61%) Forward Only Surviving: 272446 (38.20%) Reverse Only Surviving: 691 (0.10%) Dropped: 662 (0.09%)
```

Set path to trimmed Illumina reads.
```
illumina_1="$DIR/Reads/Illumina/Trimmed/illumina_1.paired.fastq.gz"
illumina_2="$DIR/Reads/Illumina/Trimmed/illumina_2.paired.fastq.gz"
```


### Demultiplex MinION Reads 
Pass reads have previously been demultiplexed by Metrichore.

Note: Multiple diverse prokaryotic and eukaryotic sample libraries were included in this run. For this analysis we are only interested in the target barcode. To be conservative we remove other barcodes first. 

We will demultiplex the fail reads by finding the closest match to the barcode string (the threshold is less than 14 substitutions, insertions or deletions between the sequence and the barcode). Reads are also trimmed of adapter sequence if possible.

```
cd $DIR/Reads/minION/
perl $DIR/Scripts/FilterBarcodes.pl $DIR/Reads/minION/Raw/reads_fail_raw.fasta $DIR/Scripts/ONTBarCodes.fasta $DIR/Reads/minION/Demultiplexed_Fail/ 14

```

Check the summary file:

```
Barcode		Sequence	Match	Conflicted

No Barcode	NA	9501	0
NB1		GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTT	14	5
NB10	GGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCTTT	7	0
NB11	GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCTTT	1499	9  <---- This is the target barcode
NB12	GGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCTTT	904	3
NB2		GGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCTTT	23	7
NB3		GGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCTTT	1744	2
NB4		GGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCTTT	1087	6
NB5		GGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTT	591	18
NB6		GGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCTTT	1008	11
NB7		GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCTTT	2999	20
NB8		GGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCTTT	631	10
NB9		GGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCTTT	766	5
```

We will perform some analysis of these reads later on (the proportion of demultiplexed reads that are S. aureus). 

Concatenate filtered pass and fail reads. This will be our working file for the remainder of the workflow. 

``` 
cat $DIR/Reads/minION/Raw/reads_pass.fasta $DIR/Reads/minION/Demultiplexed_Fail/NB11.fasta > $DIR/Reads/minION/nanopore_combined.fasta
```
Set path to MinION reads.
```
nanopore_reads=$DIR/Reads/minION/nanopore_combined.fasta
```

#### Analyse Read Length Distribution


```
cd $DIR/Reads/minION 
awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' Demultiplexed_Fail/NB11.fasta > Read_length_fail.txt
awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' Raw/reads_pass.fasta > Read_length_pass.txt
```
Plot histograms and figures.
```
mkdir $DIR/Figures
Rscript $DIR/Scripts/PlotReadLHistograms_MinionReads.R Read_length_fail.txt Read_length_pass.txt $DIR/Figures
rm Rplots.pdf
```
```
Results:
----------------------
[1] "Summary File 1"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    250    5124    7502    7567    9829   23380 
[1] "Summary File 2"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    345    5562    7660    7746    9812   20590 
[1] "Combined File 2"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    250    5287    7577    7651    9825   23380 
[1] "Number of Reads File 1 = 1498"
[1] "Number of Reads File 2 = 1323"
[1] "Number of Reads Combined = 2821"
----------------------
```

### Assemble Genome

Assembly was performed using SPAdes v 3.61 in hybrid assembly mode  (i.e. with `-nanopore` option). Later versions of SPAdes produce slightly different output that may require modification of downstream scripts. To avoid confusion use v3.6.1.

Note: Set the appropriate number of threads for your desktop/cluster using the `-t` option.
```
cd $DIR
mkdir $DIR/Assembly
spades.py --pe1-1 $illumina_1 --pe1-2 $illumina_2 --nanopore $nanopore_reads --cov-cutoff 5 --careful -t 6 -k 21,33,55,77,99,127 -o $DIR/Assembly/
```
Summarise assembly
 
```
grep ">" $DIR/Assembly/scaffolds.fasta 
```
```
4 Contigs:
-----------------------
>NODE_1_length_2857339_cov_20.968_ID_575 <-- Chromosome 
>NODE_2_length_27751_cov_40.0464_ID_3259 <-- Plasmid 1
>NODE_3_length_3252_cov_3778.63_ID_5847 <-- Plasmid 2
>NODE_4_length_128_cov_1170_ID_6546 *** This is a long repeat of Cs - it will be removed
-----------------------
```
Separate out Chromosome/Plasmid_1/Plasmid_2
```
mkdir $DIR/Contigs
mkdir $DIR/Contigs/Raw
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $DIR/Assembly/scaffolds.fasta > $DIR/Contigs/Raw/all_contigs.fasta # Make one line multifasta
grep -A 1 '>NODE_1_' < $DIR/Contigs/Raw/all_contigs.fasta > $DIR/Contigs/Raw/Chromosome.fasta
grep -A 1 '>NODE_2_' < $DIR/Contigs/Raw/all_contigs.fasta > $DIR/Contigs/Raw/PlasmidA.fasta
grep -A 1 '>NODE_3_' < $DIR/Contigs/Raw/all_contigs.fasta > $DIR/Contigs/Raw/PlasmidB.fasta
```

BLAST contigs to find nearest reference genomes - Genbank and fasta files stored in References 
```
---------------------------------------------------
USA300_FPR3757 - Chromosome 
SAP046A - PlasmidA
SAP046B - PlasmidB
---------------------------------------------------
```
Note:  Reference sequences can be found in the $DIR/References/ directory.

#### Align contigs to reference sequences

Using MAUVE

```
cd $DIR/
mkdir $DIR/Alignments/
progressiveMauve --seed-weight=24 --output=$DIR/Alignments/AllRefsvsAllScaffolds $DIR/References/All_References.fasta $DIR/Assembly/scaffolds.fasta
mauve $DIR/Alignments/AllRefsvsAllScaffolds
```

### Circularise Contigs

Contigs are circularised using the redundant overlapping sequences at the end of the contigs (typically the length of the k-mer used to produce the contig). 
Sequence start sites are fixed at the origin of replication or, in the case of plasmids, to the beginning of the reference sequence.

```
perl $DIR/Scripts/CirculariseOnOverlaps.pl $DIR/Contigs/Raw/Chromosome.fasta $DIR/References/USA300_FPR3757.fasta 1000
perl $DIR/Scripts/CirculariseOnOverlaps.pl $DIR/Contigs/Raw/PlasmidA.fasta $DIR/References/SAP046A.fasta 1000
perl $DIR/Scripts/CirculariseOnOverlaps.pl $DIR/Contigs/Raw/PlasmidB.fasta $DIR/References/SAP046B.fasta 200
```
Move files to new folder.
```
mkdir $DIR/Contigs/Circularised
find $DIR/Contigs/Raw/ -type f -iname \*.trimmed -print0 | xargs -0 -I{} mv {} $DIR/Contigs/Circularised/
find $DIR/Contigs/Raw/ -type f -iname \*.clustalw -print0 | xargs -0 -I{} mv {} $DIR/Contigs/Circularised/
find $DIR/Contigs/Raw/ -type f -iname \*.circularised.fasta -print0 | xargs -0 -I{} mv {} $DIR/Contigs/Circularised/
cat $DIR/Contigs/Circularised/Chromosome.circularised.fasta $DIR/Contigs/Circularised/PlasmidA.circularised.fasta $DIR/Contigs/Circularised/PlasmidB.circularised.fasta > $DIR/Contigs/Circularised/All_Circularised.fasta
```
Check alignment in Mauve
```
progressiveMauve --seed-weight=24 --output=$DIR/Alignments/AllRefsvsAllScaffoldsCircularised $DIR/References/All_References.fasta $DIR/Contigs/Circularised/All_Circularised.fasta
mauve $DIR/Alignments/AllRefsvsAllScaffoldsCircularised
```
Note : Output SNP files and Gaps for further analysis.



#### Annotate genomes and references with prokka [Optional]

Annotate reference and analysis assemblies with prokka. 

MHO_001: 

```

# Rename contigs to be genbank + prokka compliant
mkdir $DIR/Annotation/
awk '/^>/{print ">Chromosome"; next}{print}' < $DIR/Contigs/Circularised/Chromosome.circularised.fasta > $DIR/Annotation/Chromosome.fasta
awk '/^>/{print ">PlasmidA"; next}{print}' < $DIR/Contigs/Circularised/PlasmidA.circularised.fasta > $DIR/Annotation/PlasmidA.fasta
awk '/^>/{print ">PlasmidB"; next}{print}' < $DIR/Contigs/Circularised/PlasmidB.circularised.fasta > $DIR/Annotation/PlasmidB.fasta
cat $DIR/Annotation/Chromosome.fasta $DIR/Annotation/PlasmidA.fasta $DIR/Annotation/PlasmidB.fasta > $DIR/Annotation/Genome.fasta
```
Run Prokka. 

```
prokka --outdir $DIR/Annotation/Genome/ --force --prefix "Genome" --locustag "Genome" --genus 'Staphylococcus' --species 'aureus' --kingdom 'Bacteria' --addgenes --cpus '6' --usegenus  $DIR/Annotation/Genome.fasta
prokka --outdir $DIR/Annotation/Chromosome/ --force --prefix "Chromosome" --locustag "Chromosome" --genus 'Staphylococcus' --species 'aureus' --kingdom 'Bacteria' --addgenes --cpus '6' $DIR/Annotation/Chromosome.fasta
prokka --outdir $DIR/Annotation/PlasmidA/ --force --prefix "PlasmidA" --locustag "PlasmidA" --genus "Staphylococcus" --species "aureus" --plasmid "USA300A" --kingdom Bacteria --cpus 6 $DIR/Annotation/PlasmidA.fasta
prokka --outdir $DIR/Annotation/PlasmidB/ --force --prefix "PlasmidB" --locustag "PlasmidB" --genus "Staphylococcus" --species "aureus" --plasmid "USA300B" --kingdom Bacteria --cpus 6 $DIR/Annotation/PlasmidB.fasta
```
Make reference chromosome and plasmids prokka complaint
```
awk '/^>/{print ">USA300_FPR3757"; next}{print}' < $DIR/References/USA300_FPR3757.fasta > $DIR/Annotation/USA300_FPR3757.fasta
awk '/^>/{print ">SAP046A"; next}{print}' < $DIR/References/SAP046A.fasta > $DIR/Annotation/SAP046A.fasta
awk '/^>/{print ">SAP046B"; next}{print}' < $DIR/References/SAP046B.fasta > $DIR/Annotation/SAP046B.fasta
cat $DIR/Annotation/USA300_FPR3757.fasta $DIR/Annotation/SAP046A.fasta $DIR/Annotation/SAP046B.fasta > $DIR/Annotation/USA300_genome.fasta

````

Run prokka

```
prokka --outdir $DIR/Annotation/USA300/ --force --prefix "USA300" --locustag "USA300" --genus "Staphylococcus" --species "aureus" --kingdom Bacteria --addgenes --usegenus --cpus 6 $DIR/Annotation/USA300_genome.fasta
prokka --outdir $DIR/Annotation/USA300_FPR3757/ --prefix "USA300_FPR3757" --force --locustag "USA300_FPR3757" --genus "Staphylococcus" --species "aureus" --kingdom Bacteria --addgenes --usegenus --cpus 6 $DIR/Annotation/USA300_FPR3757.fasta
prokka --outdir $DIR/Annotation/SAP046A/ --force --prefix "SAP046A" --locustag "SAP046A" --genus "Staphylococcus" --species "aureus" --kingdom Bacteria --addgenes --usegenus --cpus 6 $DIR/Annotation/SAP046A.fasta
prokka --outdir $DIR/Annotation/SAP046B/ --force --prefix "SAP046B" --locustag "SAP046B" --genus "Staphylococcus" --species "aureus" --kingdom Bacteria --addgenes --usegenus --cpus 6 $DIR/Annotation/SAP046B.fasta

```

Check alignment in Mauve.

```
progressiveMauve --seed-weight=24 --output=$DIR/Alignments/AllRefsvsAllAnnotated $DIR/Annotation/USA300/USA300.gbk $DIR/Annotation/Genome/Genome.gbk
mauve $DIR/Alignments/AllRefsvsAllAnnotated
```


### Read Coverage Stats

Map the long and short reads to the circularised genome for summary statistics.

Make mapping directory to work in.

```
mkdir $DIR/Mapping
cd $DIR/Mapping

```

Start by mapping Illumina reads. 

Copy and index files for mapping.

```
cp $DIR/Annotation/Genome.fasta $DIR/Mapping/Genome.fasta
bwa index $DIR/Mapping/Genome.fasta
samtools faidx $DIR/Mapping/Genome.fasta
```
Map the Illumina reads using BWA-mem
```
bwa mem -t 8 $DIR/Mapping/Genome.fasta $illumina_1 $illumina_2 -I 300,60,900,50 | samtools view -T $DIR/Mapping/Genome.fasta -bS - | samtools sort -T $DIR/Mapping/GenomeVsIllumina.bwa -o $DIR/Mapping/GenomeVsIllumina.bwa.bam -
```
Index the bam files for tablet/ACT or other BAM viewer.

```
samtools index $DIR/Mapping/GenomeVsIllumina.bwa.bam
```
Summarise the coverage statisticsfor short read data.
```
samtools depth $DIR/Mapping/GenomeVsIllumina.bwa.bam > $DIR/Mapping/GenomeVsIllumina.coverage 
echo "Illumina Reads"
grep "Chromosome" $DIR/Mapping/GenomeVsIllumina.coverage > $DIR/Mapping/ChromosomeVsIllumina.coverage 
echo "Chromosome"; awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }' < $DIR/Mapping/ChromosomeVsIllumina.coverage 
grep "PlasmidA" $DIR/Mapping/GenomeVsIllumina.coverage > $DIR/Mapping/PlasmidAVsIllumina.coverage 
echo "PlasmidA";  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }' < $DIR/Mapping/PlasmidAVsIllumina.coverage 
grep "PlasmidB" $DIR/Mapping/GenomeVsIllumina.coverage > $DIR/Mapping/PlasmidBVsIllumina.coverage 
echo "PlasmidB";  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }' < $DIR/Mapping/PlasmidBVsIllumina.coverage 

```
```

Results:
---------------
Illumina Reads

Chromosome
Average =  49.6106  Standard Deviation =  7.04348
PlasmidA
Average =  78.3502  Standard Deviation =  8.8514
PlasmidB
Average =  7302.04  Standard Deviation =  85.4383
---------------

```

Map Nanopore long reads reads. We will map long reads to each contig individually to avoid biases that maybe caused by mapping or high error rate reads to smaller sequences.

Copy files and index reads for mapping.

```
cp $DIR/Annotation/Chromosome.fasta $DIR/Mapping/Chromosome.fasta
bwa index $DIR/Mapping/Chromosome.fasta
samtools faidx $DIR/Mapping/Chromosome.fasta
cp $DIR/Annotation/PlasmidA.fasta $DIR/Mapping/PlasmidA.fasta
bwa index $DIR/Mapping/PlasmidA.fasta
samtools faidx $DIR/Mapping/PlasmidA.fasta
cp $DIR/Annotation/PlasmidB.fasta $DIR/Mapping/PlasmidB.fasta
bwa index $DIR/Mapping/PlasmidB.fasta
samtools faidx $DIR/Mapping/PlasmidB.fasta
```

Map reads using BWA-mem

```
bwa mem -t 8 -x ont2d $DIR/Mapping/Chromosome.fasta $nanopore_reads | samtools view -T $DIR/Mapping/Chromosome.fasta -bS - | samtools sort -T $DIR/Mapping/ChromosomeVsNanopore.bwa -o $DIR/Mapping/ChromosomeVsNanopore.bwa.bam -
bwa mem -t 8 -x ont2d $DIR/Mapping/PlasmidA.fasta $nanopore_reads | samtools view -T $DIR/Mapping/PlasmidA.fasta -bS - | samtools sort -T $DIR/Mapping/PlasmidAVsNanopore.bwa -o $DIR/Mapping/PlasmidAVsNanopore.bwa.bam -
bwa mem -t 8 -x ont2d $DIR/Mapping/PlasmidB.fasta $nanopore_reads | samtools view -T $DIR/Mapping/PlasmidB.fasta -bS - | samtools sort -T $DIR/Mapping/PlasmidBVsNanopore.bwa -o $DIR/Mapping/PlasmidBVsNanopore.bwa.bam -
```

Index the bam files for tablet/ACT or other BAM viewer.

```
samtools index $DIR/Mapping/ChromosomeVsNanopore.bwa.bam
samtools index $DIR/Mapping/PlasmidAVsNanopore.bwa.bam
samtools index $DIR/Mapping/PlasmidBVsNanopore.bwa.bam

```
Summarise the coverage statistics for long read data.

```
echo "Nanopore Reads"
samtools depth $DIR/Mapping/ChromosomeVsNanopore.bwa.bam > $DIR/Mapping/ChromosomeVsNanopore.coverage
echo "Chromosome"; awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }' < $DIR/Mapping/ChromosomeVsNanopore.coverage
samtools depth $DIR/Mapping/PlasmidAVsNanopore.bwa.bam > $DIR/Mapping/PlasmidAVsNanopore.coverage
echo "PlasmidA"; awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }' < $DIR/Mapping/PlasmidAVsNanopore.coverage
samtools depth $DIR/Mapping/PlasmidBVsNanopore.bwa.bam > $DIR/Mapping/PlasmidBVsNanopore.coverage
echo "PlasmidB"; awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }' < $DIR/Mapping/PlasmidBVsNanopore.coverage

```
```

Results:
---------------
Nanopore Reads

Chromosome
Average =  6.82736  Standard Deviation =  2.61292
PlasmidA
Average =  4.05767  Standard Deviation =  2.01433
PlasmidB
Average =  2.90293  Standard Deviation =  1.70347
---------------

```

Visualise the reads coverage with tablet.
```
# Illumina 
tablet $DIR/Mapping/GenomeVsIllumina.bwa.bam $DIR/Mapping/Genome.fasta

# Nanopore
tablet $DIR/Mapping/ChromosomeVsNanopore.bwa.bam $DIR/Mapping/Chromosome.fasta
tablet $DIR/Mapping/PlasmidAVsNanopore.bwa.bam $DIR/Mapping/PlasmidA.fasta
tablet $DIR/Mapping/PlasmidBVsNanopore.bwa.bam $DIR/Mapping/PlasmidB.fasta
```

Prepare Figure 2 for the paper. This script has a ridiculous amount of input files (17). 

```
# Align each contig to reference seperately.
progressiveMauve --output=$DIR/Alignments/ChromosomeVsUSA300 $DIR/Annotation/Chromosome.fasta $DIR/References/USA300_FPR3757.fasta
progressiveMauve --output=$DIR/Alignments/PlasmidAVsSAP046A $DIR/Annotation/PlasmidA.fasta $DIR/References/SAP046A.fasta
progressiveMauve --output=$DIR/Alignments/PlasmidBVsSAP046B  $DIR/Annotation/PlasmidB.fasta $DIR/References/SAP046B.fasta

# Plot Figure
Rscript $DIR/Scripts/PlotAlignment.R $DIR/Alignments/ChromosomeVsUSA300.backbone $DIR/Annotation/Chromosome/Chromosome.gbk $DIR/Annotation/USA300_FPR3757/USA300_FPR3757.gbk $DIR/Scripts/Figure_Annotation.tab $DIR/Mapping/ChromosomeVsIllumina.coverage $DIR/Mapping/ChromosomeVsNanopore.coverage $DIR/Alignments/PlasmidAVsSAP046A.backbone $DIR/Annotation/PlasmidA/PlasmidA.gbk $DIR/Annotation/SAP046A/SAP046A.gbk $DIR/Mapping/PlasmidAVsIllumina.coverage $DIR/Mapping/PlasmidAVsNanopore.coverage $DIR/Alignments/PlasmidBVsSAP046B.backbone $DIR/Annotation/PlasmidB/PlasmidB.gbk $DIR/Annotation/SAP046B/SAP046B.gbk $DIR/Mapping/PlasmidBVsIllumina.coverage $DIR/Mapping/PlasmidBVsNanopore.coverage $DIR/Figures/Figure2.pdf 

```


### Align the Nanopore reads to MHO_001 with BLASR

In order to assess the similarity of the nanopore reads to the final assembly in terms of percentage sequence similarity and alignment length we will use BLASR (designed for PacBio reads but useful for Nanopore reads). 
It provides an output analogous to BLAST which is useful for our needs.

Copy reads to working directory and tidy up read names for downstream analysis.

```
mkdir $DIR/BLASR
cd $DIR/BLASR
sed -e "s/\s.*//g" $DIR/Reads/minION/Demultiplexed_Fail/NB11.fasta > reads_fail.fasta
sed -e "s/\s.*//g" $DIR/Reads/minION/Raw/reads_pass.fasta > $DIR/BLASR/reads_pass.fasta
cp $DIR/Annotation/Genome.fasta $DIR/BLASR/Genome.fasta

```

Compare pass and fail reads separately to the final MHO_001 assembly. 

```

# Generate length of each read
awk '/^>/ {if (seqlen){print read,"\t",seqlen}; read=$0 ;seqlen=0;next; } { seqlen += length($0)}END{print read,"\t",seqlen}' $DIR/BLASR/reads_pass.fasta | sed -e "s/>//g" | sed -e "s/\s\t\s/\t/g" > $DIR/BLASR/pass.lengths.tab
awk '/^>/ {if (seqlen){print read,"\t",seqlen}; read=$0 ;seqlen=0;next; } { seqlen += length($0)}END{print read,"\t",seqlen}' $DIR/BLASR/reads_fail.fasta | sed -e "s/>//g" | sed -e "s/\s\t\s/\t/g" > $DIR/BLASR/fail.lengths.tab

# Run BLASR
headers="qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV"
echo $headers > $DIR/BLASR/pass.blasr; blasr $DIR/BLASR/reads_pass.fasta $DIR/BLASR/Genome.fasta -m 4 -useccs -bestn 1 -refineConcordantAlignments >> $DIR/BLASR/pass.blasr
echo $headers > $DIR/BLASR/fail.blasr; blasr $DIR/BLASR/reads_fail.fasta $DIR/BLASR/Genome.fasta -m 4 -useccs -bestn 1 -refineConcordantAlignments >> $DIR/BLASR/fail.blasr

# Plot figures for paper
Rscript $DIR/Scripts/BLASR_Analysis.R $DIR/BLASR/pass.blasr $DIR/BLASR/fail.blasr $DIR/BLASR/pass.lengths.tab $DIR/BLASR/fail.lengths.tab $DIR/BLASR/temp Pass Fail; mv $DIR/BLASR/temp.table.pdf $DIR/Figures/Table1.pdf; mv $DIR/BLASR/temp.figure.pdf $DIR/Figures/Figure1.pdf; rm Rplots.pdf

```

Compare reads which were demultiplexed to non-target samples and reads in which no barcode was detected. 

```

# Copy unfiltered reads and filtered non-NB11 reads to working dir.
cd $DIR/Reads/minION/Demultiplexed_Fail
cat NB10.fasta NB12.fasta NB1.fasta NB2.fasta NB3.fasta NB4.fasta NB5.fasta NB6.fasta NB7.fasta NB8.fasta NB9.fasta > $DIR/BLASR/combined_filtered.fasta
cd $DIR/BLASR
sed -e "s/\s.*//g" $DIR/BLASR/combined_filtered.fasta > $DIR/BLASR/temp.txt && mv $DIR/BLASR/temp.txt $DIR/BLASR/combined_filtered.fasta
sed -e "s/\s.*//g" $DIR/Reads/minION/Demultiplexed_Fail/unbarcoded.fasta > $DIR/BLASR/unbarcoded.fasta

# Generate length of each read
awk '/^>/ {if (seqlen){print read,"\t",seqlen}; read=$0 ;seqlen=0;next; } { seqlen += length($0)}END{print read,"\t",seqlen}' $DIR/BLASR/combined_filtered.fasta | sed -e "s/>//g" | sed -e "s/\s\t\s/\t/g" > $DIR/BLASR/filtered.lengths.tab
awk '/^>/ {if (seqlen){print read,"\t",seqlen}; read=$0 ;seqlen=0;next; } { seqlen += length($0)}END{print read,"\t",seqlen}' $DIR/BLASR/unbarcoded.fasta | sed -e "s/>//g" | sed -e "s/\s\t\s/\t/g" > $DIR/BLASR/unbarcoded.lengths.tab

# Run BLASR
headers="qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV"
echo $headers > $DIR/BLASR/filtered.blasr; blasr combined_filtered.fasta Genome.fasta -m 4 -useccs -bestn 1 -refineConcordantAlignments >> $DIR/BLASR/filtered.blasr
echo $headers > $DIR/BLASR/unbarcoded.blasr; blasr unbarcoded.fasta Genome.fasta -m 4 -useccs -bestn 1 -refineConcordantAlignments >> $DIR/BLASR/unbarcoded.blasr

# Plot figures for Supplementary
Rscript $DIR/Scripts/BLASR_Analysis.R $DIR/BLASR/filtered.blasr $DIR/BLASR/unbarcoded.blasr $DIR/BLASR/filtered.lengths.tab $DIR/BLASR/unbarcoded.lengths.tab $DIR/BLASR/temp Filtered Unbarcoded; mv $DIR/BLASR/temp.table.pdf $DIR/Figures/Supp_Table1.pdf; mv $DIR/BLASR/temp.figure.pdf $DIR/Figures/Supp_Figure1.pdf; rm Rplots.pdf

```

### Confirmation of MAUVE SNPs via mapping to USA300_FPR3757

The methodology used to produce the SNPs relative to USA300_FPR3757 have been detailed in the main manuscript. 

The pipeline is currently being prepared for publication (watch this space). 

The SNPs identified by the pipeline are available in VCF format as a part of the supplementary materials for the paper (Supplementary Table 2) or the GigaScience Repository (See start of README). 

Full workings for the paper are included in the supplementary spreadsheet.
