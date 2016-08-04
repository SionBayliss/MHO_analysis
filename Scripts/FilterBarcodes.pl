#!/bin/perl

# Filter fasta files that contain barcoded sequences using fuzz string matching.
# Arguement 1 = Query reads as fasta file. 
# Arguement 2 = Reference BARCODE sequences.

# Requires Bioperl
use Bio::SeqIO; 
use Text::LevenshteinXS qw(distance);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
 
# Input fasta 
my $query=$ARGV[0];  

# Barcodes
my $reference=$ARGV[1];

# Output directory
my $output=$ARGV[2];
mkdir $output;

# Max acceptable differences 
my $no_diff=$ARGV[3];

# Read barcode sequences as Seq Object. Identify number of barcodes.
# Create output files in output directory.
$seq_in  = Bio::SeqIO->new(-format => 'fasta', -file   => $reference);
my $total=0;
while( $seq = $seq_in->next_seq() ) {
		
	# Open file per barcode.
	$barcode_name=$seq->id;
	$output_files{$barcode_name} = Bio::SeqIO->new( -file => ">$output/$barcode_name.fasta",-format => "fasta");
	$output_files{"$barcode_name.conflict"} = Bio::SeqIO->new( -file => ">$output/$barcode_name.conflict.fasta",-format => "fasta");
	
	# Store barcode info
	++$total;
	$barcode_info{$barcode_name}=$seq;
	
} $output_files{"unbarcoded"}  = Bio::SeqIO->new( -file => ">$output/unbarcoded.fasta",-format => "fasta");

# Find total reads + store read headers and sequence in memory
my $reads_in  = Bio::SeqIO->new(-format => 'fasta', -file   => $query);
my $total_reads=0; my $removed=0;
while( $seq = $reads_in->next_seq() ) {

	++$total_reads;
	
	# If read is < 300 bp discard it.
	if(length($seq->seq())>=300){
		$reads{$seq->id}=$seq; # Sequence hash 
	}else{$removed++}
	
} print "$total barcodes\n$total_reads reads\n$removed reads failed length threshold\n";

# Loop through all barcode sequences
my $seq_in  = Bio::SeqIO->new(-format => 'fasta', -file   => $reference);
foreach ( sort keys %barcode_info ) {
        
	++$count;
	print "Barcode $count of $total\n";	
	
	# Barcode info
	$header=$barcode_info{$_}->id;
	$b_sequence=$barcode_info{$_}->seq; 
	$revcom_bc=reverse($b_sequence);
	$revcom_bc=~tr/ATGC/TACG/;
	
	# Length Barcode  
	$l_bc=length($b_sequence);
	
	# Loop through each read.
	$barcode=0; $read_pos=0;
	$temp_count=0;
	foreach(sort keys %reads){
		
		$read_id=$_;
		$sq=$reads{$read_id}->seq;
		
		# Look for barcode in first and last 150 bp + barcode length
		$l=length($sq);
				
		# Find Levenshtein distance between strings
		$contains=0;
		for my $pos(0..149){ 
		
			$first=substr($sq, $pos, $l_bc); # start of read
			$last=substr($sq, $l-($pos+$l_bc), $l_bc); # end of read

			for my $barcode($b_sequence,$revcom_bc){			
			
				# Compare start of read to barcode sequence
				my $distance = distance($barcode, $first);
				
				# Store barcode distance information for each barcode for each read
				if($distance<$no_diff){
					if( !$barcode_dist{$read_id}{$header} ){
						$barcode_dist{$read_id}{$header}=$distance;
						$barcode_pos{$read_id}{$header}=$pos+1;
					}
					elsif($barcode_dist{$read_id}{$header}>$distance){
						$barcode_dist{$read_id}{$header}=$distance;
						$barcode_pos{$read_id}{$header}=$pos+1;
					}
				}
								
				# Compare end of read to barcoded sequence. 			
				my $rc_distance = distance($barcode, $last);
				
				# Store barcode distance information for each barcode for each read
				if($rc_distance<$no_diff){
					if( !$barcode_dist{$read_id}{$header} ){
						$barcode_dist{$read_id}{$header}=$distance;
						$barcode_pos{$read_id}{$header}=$l-($pos+$l_bc);
					}
					elsif($barcode_dist{$read_id}{$header}>$distance){
						$barcode_dist{$read_id}{$header}=$distance;
						$barcode_pos{$read_id}{$header}=$l-($pos+$l_bc);
					}
				}			
			}
		}
	}
}

# Find barcoded reads to seperate files based upon their best match barcode.
# Unmatched barcodes are printed to a seperate file.
# Conficted barcodes (with same distance) are printed to file with a #CONFLICT tag.
foreach(keys %reads){
	$k=$_;
	
	# Print unbarcoded reads.
	if(!$barcode_dist{$k}){
		$output_files{"unbarcoded"}->write_seq($reads{$k});
		$unbarcoded++;
	}else{
		# Find minumum distance (best barcode match);
		$mi=min(values($barcode_dist{$k}));
		@matching_keys = grep { $barcode_dist{$k}{$_} eq $mi } keys $barcode_dist{$k};
		$no_matching=scalar(@matching_keys);	

		foreach( @matching_keys ){
			$k2=$_;
			if( $no_matching > 1 ){
				$output_files{"$k2.conflict"}->write_seq($reads{$k});
				$results{$k2}{"Conflict"}++;
			}else{
				# Trim sequence after unambigous barcodes from each read.
				$temp_seq=$reads{$k}->seq();
				if($barcode_pos{$k}{$k2}<=150){
					$trimmed=substr($temp_seq, $barcode_pos{$k}{$k2}, ( length($temp_seq) - ($barcode_pos{$k}{$k2}-1) ) );
				}else{
					$trimmed=substr($temp_seq, 0, $barcode_pos{$k}{$k2} );
				}
				$reads{$k}->seq($trimmed);
				
				# print to file
				$output_files{$k2}->write_seq($reads{$k});
				$results{$k2}{"Match"}++;
			}
		}
	}
}

# Summarise barcode matching
open SUMMARY, ">$output/summary.tab";
print SUMMARY "Barcode\tMatch\tConflicted\n";
print SUMMARY "No Match\t",$unbarcoded,"\t0\n";
foreach(sort keys %results ){
	if(!$results{$_}{"Conflict"}){
		print SUMMARY "$_\t",$barcode_info{$_}->seq,"\t", $results{$_}{"Match"},"\t0\n";
	}else{
		print SUMMARY "$_\t",$barcode_info{$_}->seq,"\t", $results{$_}{"Match"},"\t",$results{$_}{"Conflict"},"\n";
	}
}



exit
