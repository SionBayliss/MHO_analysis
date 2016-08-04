#!/bin/perl 

# Circularises circular genomes where there is local overlap at the ends of contigs. 
# Accepts reference genome or DNA sequence to fix origin of replication. 
# Breaks experimental contigs, identifies overlaps and trims erroneous bases. 

# OUTPUTS:
# A) Overlap alignment (Clustalw)
# B) Trimmed contig - overlaps removed. 
# C) Origin Fixed

### NOTE: Manually check clustalw alignment to ensure it is sensible and cut positions are relevant ##

# BIOperl libs. 
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::Tools::Run::StandAloneBlastPlus;

## INPUTS ##

# File  
$CONTIG=shift;

# Contig or Gene to fix origin - The origin will be fixed on the first 200 bp. 
$FIX_START=shift;

# Overlap has been calculated using SPAdes hybrid assembly. 
# Overlap may need to me modified per experiment to get sensible results
$overlap=shift;
if($overlap eq ""){ $overlap=1000 }



## OUTPUT ##
# Find Path and Name
@suffixlist=(".fasta",".fa",".fas");
$name = fileparse($CONTIG,@suffixlist);
$dirname  = dirname($CONTIG);

# Circularised
$OUTPUT="$dirname/$name.circularised.fasta";

# Overlap alignment
$OUTPUT_ALIGN="$dirname/$name.clustalw";

## INPUT CHECK ##

# Open contig file.
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $CONTIG);
my $seq; my @seq_array; my $check=0;
while( $seq = $seq_in->next_seq() ) {
    $contig_info=$seq; ++$check;
}

# Check number of contigs is sane. 
if($check>1){ die "Too many contigs in file. Max = 1.\n"; }


## IDENTIFY OVERLAPS AND TRIM ##

# Contig id 
$sample=$contig_info->id;
$c_seq=$contig_info->seq;
$c_length=length($c_seq);
print "$sample is $c_length bp.\n";

# Take $overlap bp from each end to align.
$c1_end=substr($c_seq,0,$overlap);
$c2_end=substr($c_seq,$c_length-$overlap, $c_length);

# Align overlaps.
my @params = (quiet => 0, maxmb => '8000', maxiters => '10000'); # The parameters to be passed to MUSCLE
my $factory = Bio::Tools::Run::Alignment::Muscle->new(@params); # Factory

# Make seqIO objects
$seq_1 = Bio::Seq->new(-seq => $c1_end, -alphabet => 'dna' , -display_id=> "END1");
$seq_2 = Bio::Seq->new(-seq => $c2_end, -alphabet => 'dna' , -display_id=> "END2");
@seq_array=($seq_1, $seq_2); 

# Align overlaps using MUSCLE
my $aln = $factory->align(\@seq_array);

# Save alignment for visualisation
my $out = Bio::AlignIO->new(-file   => ">".$OUTPUT_ALIGN, -format => 'clustalw');
$out->write_aln($aln);

# Store Alignments
$aln_1=$aln->get_seq_by_id("END1");
$aln_2=$aln->get_seq_by_id("END2");

# FInd the number of matching overlapping bases. 
$match_line=$aln->match_line;
@matches = $match_line =~ /\*/g;
$match_l=length($match_line);
print scalar(@matches)," matching bases out of $overlap aligned bases\n";

# Find position of last matching aligned matching base.
$pos_count=0;
foreach(split(//,$match_line)){ ++$pos_count; if($_=~/\*/){ $cut=$pos_count} }

# Get reference position for the beginning of the contig.
$loc1=$aln_1->location_from_column($cut);
$cut_point=$loc1->start; # position to cut contig 1 - Removes unaligned and aligned i.e. redundant bases.
print "Trim beginning of contig to position $cut_point\n";
 
# Trim bases from beginning of contig. 
$contig_trimmed=substr($c_seq, $cut_point, $c_length-$cut_point);

# Trim bases from end of contig. (Produced alternative alignment if clustal shows non-identical alignment)
#$contig_trimmed=substr($c_seq, 0, $c_length-$cut_point);

# Save trimmed contig to use as  blast database. 
open TRIMMED, ">$dirname/$name.trimmed" or die "Could not write to $dirname/$name.trimmed";
print TRIMMED ">$sample\n$contig_trimmed\n";

## FIND ORIGIN IN CONTIG##

# Open $FIX_START 
my $seq_fix_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $FIX_START);
$seq_fix = $seq_fix_in->next_seq();
$origin=(substr($seq_fix->seq,0,200));

# Align origin to $CONTIG using BLAST.
my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(-db_data => "$dirname/$name.trimmed", -create => '1', -dbtype=>'nucl', -overwrite => 1);
$fac->make_db();
$seq_origin = Bio::Seq->new(-seq => $origin, -alphabet => 'dna' , -display_id=> "ORIGIN");
$result = $fac->blastn( -query=>$seq_origin, -db_data=>"$dirname/$name.trimmed", -method_args => ['-num_alignments' => 4]); 
$fac->cleanup();

# Store top hit.
while( my $hit = $result->next_hit ) {
	while( my $hsp = $hit->next_hsp ) {
		
		$strand=$hsp-> strand('hit');		
		$perc_id=$hsp->percent_identity;		
		$length=$hsp->length('total');		
		$start= $hsp->start('hit');
		
		#$end= $hsp->end('hit');
		print "Blast % id = $perc_id\tHit Pos = $start\tHit Length = $length\tStrand=$strand\n";
		last;
	}  
} 
print "Break Point = ",$start-1," bp\n";

# Break $CONTIG at origin. 
$c1=substr($contig_trimmed,0,$start-1);
$c2=substr($contig_trimmed,$start-1,length($contig_trimmed)-$start+1);

# Length of contig fragments. 
$c1_l=length($c1);
$c2_l=length($c2);
print "Contig Fragment 1 = $c1_l bp\nContig Fragment 2 = $c2_l bp\nTotal = ",$c1_l+$c2_l,"\n";

# Output is circularised contig with origin fixed on target sequence.
open OUTPUT , ">$OUTPUT" or die "OUTPUT did not open.\n";
print OUTPUT ">$sample\n$c2$c1\n";


exit;
