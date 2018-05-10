#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_f $opt_r $opt_s $opt_d $opt_m $opt_B );
use strict;

my ( $good_line, $f_line, $r_line );

&init();


### Open FILEHANDLES for forward and reverse sam files
### -f 64 -F 2048 == forward reads; no supplementary alignments
### -f 128 -F 2048 == reverse reads; no supplementary alignments
if ( -f $opt_B ) {
    open( FSAM, "samtools view -h -f 64 -F 2048 $opt_B |" );
    open( RSAM, "samtools view -h -f 128 -F 2048 $opt_B |" );
}
else {
    ### Open forward and reverse sam files
    open( FSAM, $opt_f ) or die( "$!: $opt_f\n" );
    open( RSAM, $opt_r ) or die( "$!: $opt_r\n" );
}

### Wind them through until past the header section
chomp( $f_line = <FSAM> );
while ( $f_line =~ /^@/ ) {
    chomp( $f_line = <FSAM> );
}

chomp( $r_line = <RSAM> );
while ( $r_line =~ /^@/ ) {
    chomp( $r_line = <RSAM> );
}

### Open output file handles
open( SAME_TABLE, ">$opt_s" ) or die( "$!: $opt_s\n" );
open( DIFF_TABLE, ">$opt_d" ) or die( "$!: $opt_d\n" );

$good_line = 1;
while( $good_line ) {
    $good_line = &process_lines( $f_line, $r_line );
    if ( ( $f_line = <FSAM> ) &&
	 ( $r_line = <RSAM> ) ) {
	$good_line = 1;
    }
    else {
	$good_line = 0;
    }
}

sub process_lines {
    my $f_line = shift;
    my $r_line = shift;
    my $innie;
    my ($f_scaf, $f_start, $f_end, $f_mq, $f_strand);
    my ($r_scaf, $r_start, $r_end, $r_mq, $r_strand);
    my ( @f_el, @r_el );
    @f_el = split( "\t", $f_line );
    @r_el = split( "\t", $r_line );

    ### Are these the same read? Check their IDs
    ### If not, return false
    unless ( $f_el[0] eq $f_el[0] ) {
	return 0; # We're done & stop
    }

    ### Do both reads pass the map-quality filter?
    unless ( ($f_el[4] >= $opt_m) &&
	     ($r_el[4] >= $opt_m) ) {
	return 1; # We're done, but keep on
    }

    ### All's good.
    ($f_scaf, $f_start, $f_end, $f_mq, $f_strand) = 
	&sam_line2table_data( @f_el );
    ($r_scaf, $r_start, $r_end, $r_mq, $r_strand) = 
	&sam_line2table_data( @r_el );

    ### Same scaffold or different?
    if ( $f_el[2] eq $r_el[2] ) {
	$innie = &is_it_innie( $f_start, $r_start, $f_strand );
	print SAME_TABLE ( join( "\t", ($f_el[0], 
					$f_scaf, $f_start, $f_end,
					$f_mq, $f_strand, 
					$r_scaf, $r_start, $r_end,
					$r_mq, $r_strand,
					$innie) ),
			   "\n" );
    }

    else {
	# Different scaffold
	print DIFF_TABLE ( join( "\t", ($f_el[0], 
					$f_scaf, $f_start, $f_end,
					$f_mq, $f_strand, 
					$r_scaf, $r_start, $r_end,
					$r_mq, $r_strand,
					'.') ),
			   "\n" );
    }
    return 1;
}

sub is_it_innie {
    my $f_start  = shift;
    my $r_start  = shift;
    my $f_strand = shift;

    if ( $f_start < $r_start ) {
	if ( $f_strand eq '+' ) {
	    return 1;
	}
	else {
	    return 0;
	}
    }
    else {
	if ( $f_strand eq '-' ) {
	    return 1;
	}
	else {
	    return 0;
	}
    }
}

sub sam_line2table_data {
    my @sam_elements = @_;
    my ( $scaf, $start, $end, $strand, $mq );
    $scaf   = $sam_elements[2];
    $start  = $sam_elements[3];
    $end    = $start + length( $sam_elements[9] );
    $mq     = $sam_elements[4];
    $strand = &bit_flag2strand( $sam_elements[1] );
    return ( $scaf, $start, $end, $mq, $strand );
}

sub bit_flag2strand {
    my $bit_flag = shift;
    my @bits = reverse(split( '', sprintf("%11b", ($bit_flag+2048) ) ));
    if ( $bits[4] ) {
	return '-';
    }
    else {
	return '+';
    }
}

sub init {
    my $m_DEF = 20;
    getopts( 'f:r:s:d:m:B:' );
    unless( -f $opt_B || (-f $opt_f && -f $opt_r) ) {
	print( "sesam2table.pl -f <forward sam file> -r <reverse sam file>\n" );
	print( "               -B <readname sorted BAM file>\n" );
	print( "               -s <output table for read pairs mapped to same scaffold>\n" );
	print( "               -d <output table for read pairs mapped to different scaffold>\n" );
	print( "               -m <map quality cutoff; default = $m_DEF>\n" );
	print( "Makes as output two files with tables that describe the map\n" );
	print( "positions of each forward / reverse read pair. One table\n" );
	print( "describes read pairs that both map to the same scaffold.\n" );
	print( "The other describes read pairs that map to different scaffolds.\n" );
	print( "The format of the tables are tab-delimited columns:\n" );
	print( "1. Read ID\n" );
	print( "2. Forward scaffold (or chromosome or conting, etc.)\n" );
	print( "3. Forward start coordinate\n" );
	print( "4. Forward end coordinate\n" );
	print( "5. Foward strand (+ or -)\n" );
	print( "6. Forward map-quality\n" );
	print( "7. Reverse scaffold (or chromosome or conting, etc.)\n" );
	print( "8. Reverse start coordinate\n" );
	print( "9. Reverse end coordinate\n" );
	print( "10. Foward strand (+ or -)\n" );
	print( "11. Reverse map-quality\n" );
	print( "12. Innie (1 = true, 0=false)\n" );
	exit( 0 );
    }
    unless( defined( $opt_m ) ) {
	$opt_m = $m_DEF;
    }
}
