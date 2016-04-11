#!/usr/bin/env perl
#########
# -*- mode: cperl; tab-width: 8; indent-tabs-mode: nil; basic-offset: 2 -*-
# vim:ts=8:sw=2:et:sta:sts=2
#
# Copyright (c) 2014 Oxford Nanopore Technologies Ltd.
#
# Author:        dturner
# Last Modified: $Date$
# Id:            $Id$
# $HeadURL$
#
# Usage: ./split_barcodes <filename> <stringency> (default = 14, lower = more stringent)
#
use strict;
use warnings;
use Text::LevenshteinXS qw(distance);
use Bio::SeqIO;
use Carp;
use English qw(-no_match_vars);
use Readonly;
use Getopt::Long;

our $VERSION = '0.01';

Readonly::Scalar my $SEQ_LENGTH   => 200;
Readonly::Scalar my $MATCH_WITHIN => 150;
Readonly::Scalar my $STRINGENCY   => 14;

my $opts = {};
GetOptions($opts, qw(help));

if($opts->{help}) {
  print <<"EOT" or croak "Error printing: $ERRNO";
Usage: $PROGRAM_NAME <filename> <stringency> (default = 14, lower = more stringent)
EOT
  exit;
}

my $infile     = $ARGV[0];
my $stringency = $ARGV[1] || $STRINGENCY;
local $RS      = undef;
my @barcodes   = split /\s+/smx, <DATA>;
my ($basename) = $infile =~ m{([^\\/]+)$}smx;
$basename      =~ s{[.]\S+$}{}smx;

for my $bc (@barcodes) {
  unlink "$basename-$bc.fasta";
}

#########
# precompute a few things
#
my $bc_lengths  = {};
my $revcomps    = {};
my $out_handles = {};
for my $bc (@barcodes) {
  $bc_lengths->{$bc} = length $bc;

  my $revcomp_bc   = scalar reverse $bc;
  $revcomp_bc      =~ tr/[A,T,G,C]/[T,A,C,G]/;
  $revcomps->{$bc} = $revcomp_bc;

  $out_handles->{$bc} = Bio::SeqIO->new(
					-file   => ">>$basename-${bc}.fasta",
					-format => 'Fasta',
				       );
}

my $sequences_in = 0;
my $count        = {};
my $io_in        = Bio::SeqIO->new(
				   -file   => $infile,
				   -format => 'Fasta',
				  ) or croak "Could not open '$infile' $ERRNO\n";
while (my $rec = $io_in->next_seq()) {
  my $header  = $rec->id();
  my $seq     = $rec->seq();
  my $seq_len = length $seq;
  $sequences_in ++;

  if ($seq_len <= $SEQ_LENGTH) {
    #########
    # sequence too short to bother
    #
    next;
  }

  my $min_distance = $SEQ_LENGTH; # arbitrary high number of changes
  my $min_bc;

  #########
  # iterate over the barcodes we have
  #
  for my $bc (@barcodes) {
    my $bc_length  = $bc_lengths->{$bc};
    my $revcomp_bc = $revcomps->{$bc};

    #########
    # only try and match within the first $MATCH_WITHIN bases of the target sequence
    #
    for my $scan (0 .. $MATCH_WITHIN) {
      my $head = substr $seq, $scan, $bc_length;
      my $tail = substr $seq, $seq_len-$scan, $bc_length;

      for my $window ($head, $tail) {
	for my $barcode ($bc, $revcomp_bc) {
	  my $distance = distance($barcode, $window);

	  if($distance < $min_distance) { ## no critic (ProhibitDeepNests)
	    $min_distance = $distance;
	    $min_bc       = $bc;
	  }
	}
      }
    }
  }

  if ($min_distance < $stringency) {
    my $io_out = $out_handles->{$min_bc};
    $io_out->write_seq($rec);
    $count->{$min_bc}++;
  }
}

$io_in->close();

my $sequences_out = 0;
for my $key (sort keys %{$count}) {
  printf "%s\t%d\n", $key, $count->{$key} or croak qq[Error printing: $ERRNO];
  $sequences_out += $count->{$key};
}

printf "Sequences in:  %d\nSequences out: %d\n", $sequences_in, $sequences_out or croak qq[Error printing: $ERRNO];

1;

__DATA__
GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCT
GGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCT
GGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCT
GGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCT
GGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCT
GGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCT
GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCT
GGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCT
GGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCT
GGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCT
GGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCT
GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCT
