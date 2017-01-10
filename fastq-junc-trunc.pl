#!/usr/bin/perl

use Bio::SeqIO;
use Getopt::Std;
use vars qw( $opt_f $opt_j $opt_n $opt_e );
use strict;

my $fq_in = &init();
my $fq_out = Bio::SeqIO->new( '-format' => 'fastq' );
my ( $fq_seq, $first_jnc, $trunc_seq );
my $total_trunc = 0;
while( $fq_seq = $fq_in->next_seq() ) {
    $first_jnc = &find_first_jnc( $fq_seq, $opt_j );
    if ( $first_jnc >= 0 ) {
	$fq_out->write_seq( $fq_seq->trunc(1, ($first_jnc + $opt_n)) );
	$total_trunc++;
    }
    else {
	$fq_out->write_seq( $fq_seq );
    }
}

print STDERR ( "Total number of sequences truncated: $total_trunc\n" );

sub init {
    my $j_DEF = 'GATCGATC';
    my $n_DEF = 4;
    my $fq_in;
    getopts( 'f:j:n:e' );
    unless( -f $opt_f ) {
	print( "fastq-junc-trunc.pl -f <fastq filename>\n" );
	print( "                    -j <junction sequence; default = $j_DEF>\n" );
	print( "                    -n <number of junction bases to keep; default = $n_DEF>\n" );
	print( "                    -e <truncate at the end if part of the junction\n" );
	print( "                        sequence is seen>\n" );
	print( "This program takes as input a fastq format sequence file. It searches\n" );
	print( "each entry for the first occurance of the user-defined junction sequence.\n" );
	print( "If the junction sequence is found, it truncates the fastq entry at that\n" );
	print( "point, including the user-defined number of bases of the junction sequence.\n" );
	print( "It writes to STDOUT the truncated sequence. If the junction sequence is not\n" );
	print( "found, it writes to STDOUT the original sequence.\n" );
	print( "If the end of the sequence contains an incomplete junction sequence and the\n" );
	print( "-e option is specified, this is considered a junction sequence.\n" );
	exit( 0 );
    }
    ### Gzipped input?
    if ( $opt_f =~ /\.gz$/ ) {
	$fq_in = Bio::SeqIO->new( '-file' => "zcat $opt_f |",
				  '-format' => 'fastq' );
    }
    ### Bzipped input
    else {
	if ( $opt_f =~ /\.bz2$/ ) {
	    $fq_in = Bio::SeqIO->new( '-file' => "bzcat $opt_f |",
				      '-format' => 'fastq' );
	}
	else {
	    $fq_in = Bio::SeqIO->new( '-file' => $opt_f,
				      '-format' => 'fastq' );
	}
    }
    unless( defined( $opt_j ) ) {
	$opt_j = $j_DEF;
    }
    unless( defined( $opt_n ) ) {
	$opt_n = $n_DEF;
    }
    return $fq_in;
}

sub find_first_jnc {
    my $fq_seq = shift; # Bio::Seq [fastq] 
    my $jnc_seq_str = shift; # DNA sequence that occurs at the junction
    my $jnc_seq_len;
    my ( $pos, $jnc_seq_len, $trunc_jnc_seq, $seq_len );
    $pos = index( $fq_seq->seq, $jnc_seq_str );
    # If the full, correct junction sequence occurs, return the position of 
    # it's first occurance.
    if ( $pos != -1 ) {
	return $pos;
    }

    # The full, correct junction sequence did not occur. Let's check the end of 
    # the sequence to see if it ends in part of the junction sequence.
    $jnc_seq_len = length( $jnc_seq_str ) - 1;
    $seq_len = $fq_seq->length();
    while( $jnc_seq_len >= $opt_n ) {
	$trunc_jnc_seq = substr( $jnc_seq_str, $jnc_seq_str, $jnc_seq_len );
	if (rindex($fq_seq->seq, $trunc_jnc_seq) == ( $seq_len - $jnc_seq_len ) ) {
	    return ( $seq_len - $jnc_seq_len );
	}
	$jnc_seq_len--;
    }
    return -1;
}

