# MHO 001 MinION Analysis 

### Dependencies

#### Analysis 
*Samtools (>=1.18)
*Trimmomatic 
*SPAdes
*BWA (0.7.5a-r405)
*BioPerl
*Mauve

#### Plot Creation
*R
*ggplot
*cowplot
*genoPlotR

#### Optional 
*prokka
*pilon
*tablet/artemis

# Workflow
 For publication - TBA
  
### Set working directory (change as appropriate) 
```
DIR="/path/to/working/directory/"
cd $DIR
```

#### Download GitHub Directory
git clone http://githum.com/SionBayliss/MHO_analysis.git

### Download Reads
#### MinION Reads (fasta format)
ENA read accession: ERS1178418/ERR1424936
```
mkdir Reads && cd Reads
mkdir MinION && cd MinION
mkdir Raw && cd Raw
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA628/ERA628229/oxfordnanopore_native/MHO_001.tar.gz
tar xvzf MHO_001.tar.gz
```
#### Illumina Reads
ENA read accessions: ERS1180806 and ERS1180807
```
cd $DIR/Reads
mkdir Illumina && cd Illumina
mkdir Raw && cd Raw
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
Note: two Illumina runs were performed. The reads were combined before further processing. 
```
gunzip -c $read1_1 $read2_1 | gzip - > $DIR/Reads/Illumina/Raw/illumina_1.fastq.gz
gunzip -c $read1_2 $read2_2 | gzip - > $DIR/Reads/Illumina/Raw/illumina_2.fastq.gz
```

#### Note : Convert pass and fail reads to fasta format.
The following analysis assumes you have downloaded the fasta files from the ENA. The raw native MinIOn data is also available and can be converted to fasta format using PoRe, poRetools or conversion tool of your choice.

## QC, demultiplex and trim reads


#### Trim Illumina Reads using Trimmomatic #
```
cd $DIR/Reads/Illumina/
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

```
cd $DIR/Reads/minION/Raw
mkdir Demultiplexed_Stringency14
cd Demultiplexed_Stringency14
perl $DIR/Scripts/SplitBarcodes.pl ../reads_fail_raw.fasta 14    
```
Note :  SplitBarcodes.pl is not my script. It was produced by ONT and modified for this application.

```
Results - Stringency 14
GGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCT	918
GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCT	41
GGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCT	29
GGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCT	1029
GGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCT	695
GGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCT	2049
GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCT	3468
GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCT	1783 <-- This is our barcode
GGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCT	78
GGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCT	731
GGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCT	1240
GGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCT	1155
Sequences in:  21098
Sequences out: 13216
```

Rename file of demultiplexed sample.
``` 
cd $DIR/Reads/minION/
mkdir Trimmed
cp Raw/Demultiplexed_Stringency14/reads_fail_raw-GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCT.fasta Trimmed/reads_fail_str14.fasta
cp Raw/reads_pass.fasta Trimmed/reads_pass.fasta
```

Filter fail reads for any remaining unremoved barcodes/adapters.
```
cd $DIR/Reads/minION/Trimmed/
perl $DIR/Scripts/RemoveBarcodes_LongReads.pl reads_fail_str14.fasta $DIR/Scripts/ONTBarCodes.txt
```
```
Results

-----------------------
Seq 11 of 12 <-- Our barcode
Iteration 1 - 151 hits - 8616 bases removed from 1783 sequences
Iteration 2 - 8 hits - 481 bases removed from 1783 sequences
Iteration 3 - 0 hits - 0 bases removed from 1783 sequences
Total bases removed = 9097
Total barcodes removed = 159
-----------------------
```
Filter pass reads for remaining barcodes. 
```
perl $DIR/Scripts/RemoveBarcodes_LongReads.pl reads_pass.fasta $DIR/Scripts/ONTBarCodes.txt
```
```
-----------------------
Seq 11 of 12 Seq 11 of 12 <-- Our barcode
Iteration 1 - 391 hits - 22078 bases removed from 1324 sequences
Iteration 2 - 27 hits - 1653 bases removed from 1324 sequences
Iteration 3 - 0 hits - 0 bases removed from 1324 sequences
Total bases removed = 23731
Total barcodes removed = 418
-----------------------
```
Concatenate filtered pass and fail reads. This will be our working file for the remainder of the workflow. 
```
cat reads_fail_str14.fasta.barcodes_removed reads_pass.fasta.barcodes_removed > Reads_Combined_Str14.fasta # for Assembly
```

Set path to MinION reads.
```
nanopore_reads=$DIR/Reads/minION/Trimmed/Reads_Combined_Str14.fasta
```

#### Analyse Read Length Distribution


```
awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' reads_fail_str14.fasta.barcodes_removed > Read_length_fail_str14.txt
awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' reads_pass.fasta.barcodes_removed > Read_length_pass.txt
```
Plot histograms and figure for paper.
```
mkdir $DIR/Figures
Rscript $DIR/Scripts/PlotReadLHistograms_MinionReads.R Read_length_fail_str14.txt Read_length_pass.txt $DIR/Figures
```
```
Results:
----------------------
[1] "Summary File 1"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    237    5174    7568    7605    9867   23440 
[1] "Summary File 2"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    345    5562    7653    7728    9812   20590 
[1] "Combined File 2"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    237    5300    7588    7657    9837   23440 
[1] "Number of Reads File 1 = 1782"
[1] "Number of Reads File 2 = 1323"
[1] "Number of Reads Combined = 3105"
----------------------
```

### Assemble Genome

Assembly was performed using SPAdes in hybrid assembly mode  (i.e. with `-nanopore` option).

Note: Set the appropriate number of threads for your desktop/cluster using the `-t` option.
```
cd $DIR
mkdir $DIR/Assembly
spades.py --pe1-1 $illumina_1 --pe1-2 $illumina_2 --nanopore $nanopore_reads --cov-cutoff 5 --careful -t 8 -k 21,33,55,77,99,127 -o $DIR/Assembly/ 
```
 Summarise assembly
 
 ```
grep ">" $DIR/Assembly/scaffolds.fasta 
```
```
4 Contigs:
-----------------------
>NODE_1_length_2857339_cov_20.9061_ID_6763 <-- Chromosome 
>NODE_2_length_27751_cov_39.9639_ID_6765 <-- Plasmid 1
>NODE_3_length_3252_cov_3777.09_ID_6767 <-- Plasmid 2
>NODE_4_length_128_cov_1170_ID_6769 *** This is a long repeat of Cs - remove it
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

BLAST contigs to find nearest reference genomes - Genebank and fasta files stored in References 
```
---------------------------------------------------
USA300_FPR3757 - Chromosome 
SAP046A - PlasmidA
SAP046B - PlasmidB
---------------------------------------------------
```
Note:  Reference sequences can be found in the $DIR/References directory.

#### Align contigs to reference sequences

Using Mauve contig mover and visualisation GUI.
```
cd $DIR/
mkdir $DIR/References/Alignments/
progressiveMauve --seed-weight=24 --output=$DIR/References/Alignments/AllRefsvsAllScaffolds $DIR/References/All_References.fasta $DIR/Assembly/scaffolds.fasta
mauve $DIR/References/Alignments/AllRefsvsAllScaffolds
```

### Circularise Genome

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
progressiveMauve --seed-weight=24 --output=$DIR/References/Alignments/AllRefsvsAllScaffoldsCircularised $DIR/References/All_References.fasta $DIR/Contigs/Circularised/All_Circularised.fasta
mauve $DIR/References/Alignments/AllRefsvsAllScaffoldsCircularised
```
Note : Output SNP files and Gaps for further analysis.

#### Annotate genomes and references with prokka [Optional]

```
mkdir $DIR/Annotation/
```
Rename contigs to be genbank + prokka compliant.
```
awk '/^>/{print ">Chromosome"; next}{print}' < $DIR/Contigs/Circularised/Chromosome.circularised.fasta > $DIR/Annotation/Chromosome.fasta
awk '/^>/{print ">PlasmidA"; next}{print}' < $DIR/Contigs/Circularised/PlasmidA.circularised.fasta > $DIR/Annotation/PlasmidA.fasta
awk '/^>/{print ">PlasmidB"; next}{print}' < $DIR/Contigs/Circularised/PlasmidB.circularised.fasta > $DIR/Annotation/PlasmidB.fasta
cat $DIR/Annotation/Chromosome.fasta $DIR/Annotation/PlasmidA.fasta $DIR/Annotation/PlasmidB.fasta > $DIR/Annotation/Genome.fasta
```
Run Prokka. 

```
 prokka --outdir $DIR/Annotation/Genome/ --force --locustag "Genome" --genus 'Staphylococcus' --species 'aureus' --kingdom 'Bacteria' --addgenes --cpus '6' --addgenes --usegenus  $DIR/Annotation/Chromosome.fasta
prokka --outdir $DIR/Annotation/Chromosome/ --force --locustag "Chromosome" --genus 'Staphylococcus' --species 'aureus' --kingdom 'Bacteria' --addgenes --cpus '6' --addgenes $DIR/Annotation/Chromosome.fasta
prokka --outdir $DIR/Annotation/PlasmidA/ --force --locustag "PlasmidA" --genus "Staphylococcus" --species "aureus" --plasmid "USA300A" --kingdom Bacteria --addgenes --cpus 6 $DIR/Annotation/PlasmidA.fasta
prokka --outdir $DIR/Annotation/PlasmidB/ --force --locustag "PlasmidB" --genus "Staphylococcus" --species "aureus" --plasmid "USA300B" --kingdom Bacteria --addgenes --cpus 6 $DIR/Annotation/PlasmidB.fasta
```
Run prokka on the reference genome combined with USA300 plasmids.
```
awk '/^>/{print ">USA300_FPR3757"; next}{print}' < $DIR/References/USA300_FPR3757.fasta > $DIR/Annotation/USA300_FPR3757.fasta
awk '/^>/{print ">SAP046A"; next}{print}' < $DIR/References/SAP046A.fasta > $DIR/Annotation/SAP046A.fasta
awk '/^>/{print ">SAP046B"; next}{print}' < $DIR/References/SAP046B.fasta > $DIR/Annotation/SAP046B.fasta
cat $DIR/Annotation/USA300_FPR3757.fasta $DIR/Annotation/SAP046A.fasta $DIR/Annotation/SAP046B.fasta > $DIR/Annotation/USA300_genome.fasta
prokka --outdir $DIR/Annotation/USA300/ --force --locustag "USA300" --genus "Staphylococcus" --species "aureus" --kingdom Bacteria --addgenes --usegenus--cpus 6 $DIR/Annotation/USA300_genome.fasta
```
Check alignment in Mauve.
```
progressiveMauve --seed-weight=24 --output=$DIR/References/Alignments/AllRefsvsAllAnnotated $DIR/Annotation/USA300/USA300*.gbk $DIR/Annotation/Genome/Genome*.gbk
mauve $DIR/References/Alignments/AllRefsvsAllAnnotated
```



### Read Coverage Stats
Map the long and short reads to the circularised genome for summary statistics.

Make mapping directory to work in.
```
mkdir $DIR/Mapping
cd $DIR/Mapping
cp $DIR/Contigs/Annotation/All_Circularised.fasta $DIR/Contigs/Mapping/Genome.fasta
cp $DIR/Contigs/Annotation/Chromosome.circularised.fasta $DIR/Contigs/Mapping/Chromosome.fasta
cp $DIR/Contigs/Annotation/PlasmidA.circularised.fasta $DIR/Contigs/Mapping/PlasmidA.fasta
cp $DIR/Contigs/Annotation/PlasmidB.circularised.fasta $DIR/Contigs/Mapping/PlasmidB.fasta
```

Index
```
bwa index $DIR/Mapping/Genome.fasta
samtools faidx $DIR/Mapping/Genome.fasta
bwa index $DIR/Mapping/Chromosome.fasta
samtools faidx $DIR/Mapping/Chromosome.fasta
bwa index $DIR/Mapping/PlasmidA.fasta
samtools faidx $DIR/Mapping/PlasmidA.fasta
bwa index $DIR/Mapping/PlasmidB.fasta
samtools faidx $DIR/Mapping/PlasmidB.fasta
```
Map the Illumina reads.
```
bwa mem -t 8 $DIR/Mapping/Genome.fasta $illumina_1 $illumina_2 -I 300,60,900,50 | samtools view -T $DIR/Mapping/Genome.fasta -bS - | samtools sort -T $DIR/Mapping/GenomeVsIllumina.bwa -o $DIR/Mapping/GenomeVsIllumina.bwa.bam -
```

Map Nanopore reads. 
Note: The minION reads are mapped seperately to the chromosome and plasmids as mapping of the error prone long reads is problematic and spurious mapping to the longer chromosome rather than the plasmids seems likely. 
```
bwa mem -t 8 -x ont2d $DIR/Mapping/Chromosome.fasta $nanopore_reads | samtools view -T $DIR/Mapping/Chromosome.fasta -bS - | samtools sort -T $DIR/Mapping/ChromosomeVsNanopore.bwa -o $DIR/Mapping/ChromosomeVsNanopore.bwa.bam -
bwa mem -t 8 -x ont2d $DIR/Mapping/PlasmidA.fasta $nanopore_reads | samtools view -T $DIR/Mapping/PlasmidA.fasta -bS - | samtools sort -T $DIR/Mapping/PlasmidVsNanopore.bwa -o $DIR/Mapping/PlasmidAVsNanopore.bwa.bam -
bwa mem -t 8 -x ont2d $DIR/Mapping/PlasmidB.fasta $nanopore_reads | samtools view -T $DIR/Mapping/PlasmidB.fasta -bS - | samtools sort -T $DIR/Mapping/PlasmidBVsNanopore.bwa -o $DIR/Mapping/PlasmidBVsNanopore.bwa.bam -
```
 Index the bam files.
```
samtools index $DIR/Mapping/GenomeVsIllumina.bwa.bam
samtools index $DIR/Mapping/ChromosomeVsNanopore.bwa.bam
samtools index $DIR/Mapping/PlasmidAVsNanopore.bwa.bam
samtools index $DIR/Mapping/PlasmidBVsNanopore.bwa.bam
```
Summarise the coverage statistics
```
echo "ChromosomeVSIllumina"; samtools depth $DIR/Mapping/GenomeVsIllumina.bwa.bam | grep "NODE_1_" | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }'
echo "PlasmidAVSIllumina"; samtools depth $DIR/Mapping/GenomeVsIllumina.bwa.bam | grep "NODE_2_" |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }'
echo "PlasmidBVSIllumina"; samtools depth $DIR/Mapping/GenomeVsIllumina.bwa.bam | grep "NODE_3_" | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }'
echo "ChromosomeVSminION";samtools depth $DIR/Mapping/ChromosomeVsNanopore.bwa.bam | grep "NODE_1_" | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }'
echo "PlasmidAVSminION";samtools depth $DIR/Mapping/PlasmidAVsNanopore.bwa.bam | grep "NODE_2_" | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }'
echo "PlasmidBVSminION";samtools depth $DIR/Mapping/PlasmidBVsNanopore.bwa.bam | grep "NODE_3_" | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR," Standard Deviation = ",sqrt(sum/NR-(sum/NR**2)) }'
```

```
Results:
-----------------------

ChromosomeVSIllumina
Average =  49.6106  Standard Deviation =  7.04348
PlasmidAVSIllumina
Average =  78.5701  Standard Deviation =  8.86381
PlasmidAVSIllumina
Average =  7300.33  Standard Deviation =  85.4283
ChromosomeVSminION
Average =  7.44738  Standard Deviation =  2.72899
PlasmidAVSminION
Average =  4.50124  Standard Deviation =  2.12157
PlasmidBVSminION
Average =  2.89865  Standard Deviation =  1.70221

-----------------------
```

Visualise the reads coverage with tablet.
```
tablet $DIR/Mapping/GenomeVsIllumina.bwa.bam $DIR/Mapping/Genome.fasta
tablet $DIR/Mapping/ChromosomeVsNanopore.bwa.bam $DIR/Mapping/Chromosome.fasta
tablet $DIR/Mapping/PlasmidAVsNanopore.bwa.bam $DIR/Mapping/PlasmidA.fasta
tablet $DIR/Mapping/PlasmidBVsNanopore.bwa.bam $DIR/Mapping/PlasmidB.fasta
```

Note: Plasmid A has a region with very high coverage corresponding with an insertion sequence (of presumably high copy number in the genome). 

### Optional: Polish assembly with short reads 

Note : This seems to produce as many errors as it corrects - if fails to correct the regions between ribosomal genes and introduces SNPs relative to the reference genome. 

```
mkdir $DIR/Polished_Contigs
PILON_PATH="/opt/Pilon" ## Change for system
```
Run pilon for each contig seperately.
Warning : Pilon fixes some base positions BUT introduces many errors in highly repetitive elements (such as transposases). 
```
java -Xmx12288m -jar $PILON_PATH/pilon-1.16.jar --genome $DIR/Mapping/Chromosome.fasta --frags $DIR/Mapping/ChromosomeVsIllumina.bwa.bam --output $DIR/Polished_Contigs/Chromosome --changes --fix bases --threads 6
java -Xmx12288m -jar $PILON_PATH/pilon-1.16.jar --genome $DIR/Mapping/PlasmidA.fasta --frags $DIR/Mapping/PlasmiAVsIllumina.bwa.bam --output $DIR/Polished_Contigs/Plasmid_1 --changes --fix bases --threads 6
java -Xmx12288m -jar $PILON_PATH/pilon-1.16.jar --genome $DIR/Mapping/PlasmidA.fasta --frags $DIR/Mapping/PlasmidAVsIllumina.bwa.bam --output $DIR/Polished_Contigs/Plasmid_2 --changes --fix bases --threads 6
```
Convert to uppercase (lowercase bases are introduced by pilon).
```
awk '{ print toupper($0) }' $DIR/Polished_Contigs/Chromosome.fasta > $DIR/Polished_Contigs/Chromosome.uppercase.fasta
awk '{ print toupper($0) }' $DIR/Polished_Contigs/Plasmid_1.fasta > $DIR/Polished_Contigs/Plasmid_1.uppercase.fasta
awk '{ print toupper($0) }' $DIR/Polished_Contigs/Plasmid_2.fasta > $DIR/Polished_Contigs/Plasmid_2.uppercase.fasta
```
Visualise with Mauve .
```
cat $DIR/Polished_Contigs/Chromosome.uppercase.fasta $DIR/Polished_Contigs/Plasmid_1.uppercase.fasta $DIR/Polished_Contigs/Plasmid_2.uppercase.fasta > $DIR/Polished_Contigs/Genome_Polished.fasta
progressiveMauve --seed-weight=24 --output=$DIR/References/Alignments/AllRefsvsAllPolished $DIR/References/All_References.fasta $DIR/Polished_Contigs/Genome_Polished.fasta
mauve $DIR/References/Alignments/AllRefsvsAllPolished
```
