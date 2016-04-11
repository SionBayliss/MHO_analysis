#!/bin/perl 

# Simple script for removing barcodes/adapters from minION long reads. 
# Arguement 1 = Query
# Arguement 2 = Reference DNA sequences.

# Requires Bioperl
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SearchIO; 

# Input 
my $query=$ARGV[0];  

# Sequence to filter.
my $reference=$ARGV[1];

# Output file
my $outfile="$query\.barcodes_removed" ;
my $seq_out = Bio::SeqIO->new(-file   => ">$outfile", -format => 'fasta');

# Create working file.
my $query_db = Bio::SeqIO->new(-format => 'fasta', -file   => $query);
while (my $query_entry = $query_db->next_seq) {
    $seq_out->write_seq($query_entry);
}

# Index database (output file)
my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(-db_data => $outfile, -create => '1', -dbtype=>'nucl', -overwrite => 1);
$fac->make_db();

# Read reference as Seq Object
my $seq_in  = Bio::SeqIO->new(-format => 'fasta', -file   => $reference);
my $total=0;
while( $seq = $seq_in->next_seq() ) {++$total} # find number of sequences;
my $seq_in  = Bio::SeqIO->new(-format => 'fasta', -file   => $reference);

# Loop through all reference sequences. 
my $count=0; my $total_removed=0; $total_barcodes=0;
my $seq; my @seq_array; my $header; my $sequence; my $bases_removed;my @position;
while( $seq = $seq_in->next_seq() ) {
        
	++$count;
	print "Seq $count of $total\n";	
		
	# Info 
	$header=$seq->id;
	$sequence=$seq->seq;	    
	
	# Variable that indicates whether all basrcodes have been removed. 
	$bar_check=0;
	
	$iteration=0;
	while($bar_check==0){
		
		++$iteration;	
		
		%blast_hash=();
		
		$removed_bases=0;
	
		# Blast against reference 
		my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(-db_data => $outfile, -create => '1', -dbtype=>'nucl', -overwrite => 1);
		$fac->make_db();
		$result_in = $fac->blastn( -query=>$seq, -db_data=>$outfile, -method_args => ['-num_threads'=>4, '-num_alignments' => 1000000 ]); 
		
		# Filter results.
		$total_coverage=0; $hits=0;
		 while( my $hit = $result_in->next_hit ) {
		   while( my $hsp = $hit->next_hsp ) {
		     if( $hsp->length('total') > (length($sequence)/2)) {
			if ( $hsp->percent_identity >= 50 ) {		  
			 		    
			    # Find position of hit
			    @position=$hsp->range('hit');
			    
			    # Hit name minus blast tag
			    $hit->name=~/lcl\|(.+)/;
			    $hit_name=$1;
			    
			    #Store info
			    $blast_hash{$hit_name}=join(":",@position);			    
			    
			    # Record Number of Hits
			    ++$hits; 
			    
			}
		     }
		   }  
		 }
		 
 		# Cleanup temporary files  
		$fac->cleanup();
		 
		# Remove barcode sequences from fasta file. 
		 
		# Identify hit query sequences- read trimmed and untrimmed sequences to output file. 
		my $query_db = Bio::SeqIO->new(-format => 'fasta', -file   => $outfile);
		my $seq_out = Bio::SeqIO->new(-file   => ">$outfile.temp", -format => 'fasta');
		$seq_check=0;		
		while( $q = $query_db->next_seq() ) {
			++$seq_check;			
			if(!$blast_hash{$q->id}){
				 $seq_out->write_seq($q);
			}else{				
				
				$seq_header=$q->id;
				$q_seq=$q->seq;
				
				$positions=$blast_hash{$q->id};
				@position=split(":",$positions);
				
				# Find nearest end of sequence.
				my $target = $position[0];
				my $l_seq=length($q_seq);
				my @vals = (0, $l_seq);
					     
				# Index of nearest end. 
				$idx=0;
				if( abs($target-$vals[0]) > abs($target-$vals[1]) ){  $idx=1  }
						  
				# Only capture bases outside of hit match.   
				if($idx==0){
					$output_sequence=substr($q->seq, $position[1], $l_seq-$position[1]);
					$removed_bases+=$position[1]; $total_bases+=$position[1];
				}
				if($idx==1){
					$output_sequence=substr($q->seq, 0, ($position[0]-1) );
					$removed_bases+=$l_seq-$position[0]; $total_bases+=$l_seq-$position[0];
				}
				
				# Print trimmed sequence to file.
				$q->seq($output_sequence);
				$seq_out->write_seq($q);
			}
		}
				
		# Rename to working file. 
		rename "$outfile.temp", $outfile;				
		
		# Iteration Summary. 
		print "Iteration $iteration - $hits hits - $removed_bases bases removed from $seq_check sequences\n";
		$total_barcodes+=$hits;
		
		if($hits==0){
			$bar_check=1;
		}
	}		
}

# Summary Stats.
print "Total bases removed = $total_bases\n";
print "Total barcodes removed = $total_barcodes\n";

exit