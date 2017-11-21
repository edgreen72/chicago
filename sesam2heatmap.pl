#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_f $opt_r $opt_m $opt_b $opt_m $opt_M $opt_B );
use strict;

my ( $good_line, $f_line, $r_line, $bin, $i, $j, $chr );
my ( @sorted_chrs );
my ( %chr2len );
my ( %chr2bin );
my %heatmap;
my $HEATMAX;

&init();

### Open FILEHANDLES for forward and reverse sam files
if ( -f $opt_B ) {
    open( FSAM, "samtools view -h -f 64 $opt_B |" );
    open( RSAM, "samtools view -h -f 128 $opt_B |" );
}
else {
    if ( $opt_f =~ /.sam$/ ) {
	open( FSAM, $opt_f ) or die( "$!: $opt_f\n" );
    }
    elsif ( $opt_f =~ /.bam$/ ) {
	open( FSAM, "samtools view -h $opt_f |" );
    }
    else {
	print STDERR ("Don't know what kind of file $opt_f is\n" );
	exit( 1 );
    }
    if ( $opt_r =~ /.sam$/ ) {
	open( RSAM, $opt_r ) or die( "$!: $opt_r\n" );
    }
    elsif ( $opt_r =~ /.bam$/ ) {
	open( RSAM, "samtools view -h $opt_r |" );
    }
    else {
	print STDERR ("Don't know what kind of file $opt_r is\n" );
	exit( 1 );
    }
}

### Wind them through until past the header section
### populated %chr2len from header @SQ SN lines
### Don't include if it's < $opt_M
chomp( $f_line = <FSAM> );
while ( $f_line =~ /^@/ ) {
    if ( $f_line =~ /^\@SQ/ ) {
	&parse_SN_header( $f_line );
    }
    chomp( $f_line = <FSAM> );
}

chomp( $r_line = <RSAM> );
while ( $r_line =~ /^@/ ) {
    chomp( $r_line = <RSAM> );
}

### Sort chromosomes by length
@sorted_chrs = sort { $chr2len{$b} <=> $chr2len{$a} } (keys %chr2len);

### Populate %chr2bin with the first bin of each chromosome
$bin = 0;
foreach $chr ( @sorted_chrs ) {
    $chr2bin{ $chr } = $bin;
    $bin += int( $chr2len{ $chr } / $opt_b ) + 1;
}
$HEATMAX = $bin; # set the maximum bin; useful for output matrix

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

for( $i = 0; $i <= $HEATMAX; $i++ ) {
    for( $j = 0; $j <= $HEATMAX; $j++ ) {
	if ( defined( $heatmap{$i}->{$j} ) ) {
	    printf( "%d %d %d\n", $i, $j, $heatmap{$i}->{$j} );
	}
	else {
	    printf( "%d %d 0\n", $i, $j );
	}
    }
    print( "\n" );
}

foreach $chr ( @sorted_chrs ) {
    printf( "# %s %d\n", $chr, $chr2bin{ $chr } );
}

sub process_lines {
    my $f_line = shift;
    my $r_line = shift;
    my $innie;
    my ($f_scaf, $f_start, $f_end, $f_mq, $f_strand, $f_mid);
    my ($r_scaf, $r_start, $r_end, $r_mq, $r_strand, $r_mid);
    my ( $f_bin, $r_bin );
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

    ### Are both reads mapped to a sequence that passed
    ### length filter?
    unless( $chr2len{ $f_el[2] } &&
	    $chr2len{ $r_el[2] } ) {
	return 1;
    }
    
    ### All's good.
    ($f_scaf, $f_start, $f_end, $f_mq, $f_strand) = 
	&sam_line2table_data( @f_el );
    ($r_scaf, $r_start, $r_end, $r_mq, $r_strand) = 
	&sam_line2table_data( @r_el );

    $f_mid = ($f_start + $f_end) / 2;
    $r_mid = ($r_start + $r_end) / 2;

    $f_bin = $chr2bin{ $f_scaf } + int( $f_mid / $opt_b );
    $r_bin = $chr2bin{ $r_scaf } + int( $r_mid / $opt_b );

    if ( $f_bin == $r_bin ) { # same bin, just put it in once
	$heatmap{$f_bin}->{$r_bin}++;
    }
    else {
	$heatmap{$f_bin}->{$r_bin}++;
	$heatmap{$r_bin}->{$f_bin}++;
    }
    
    return 1;
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

sub parse_SN_header {
    my $sn_line = shift;
    my ( $chr, $len );

    ($chr, $len) = ( $sn_line =~ /SN:(\S+)\s+LN:(\d+)/);
    if ( $len >= $opt_M ) {
	$chr2len{ $chr } = $len;
    }
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
    my $b_DEF = 2000000;
    my $M_DEF = 100000;
    getopts( 'f:r:m:b:M:B:' );
    unless( -f $opt_B ||( -f $opt_f && -f $opt_r ) ) {
	print( "sesam2table.pl -f <forward sam/bam file>\n" );
	print( "               -r <reverse sam/bam file>\n" );
	print( "               -B <read name sorted bam file with both forward and reverse>\n" );
	print( "               -m <map quality cutoff; default = $m_DEF>\n" );
	print( "               -b <bin size; default = $b_DEF>\n" );
	print( "               -M <minimum length of chr/scaf/contig; default = $M_DEF>\n" );
	print( "Makes heatmap data file for gnuplot visualization.\n" );
	exit( 0 );
    }
    unless( defined( $opt_m ) ) {
	$opt_m = $m_DEF;
    }
    unless( defined( $opt_b ) ) {
	$opt_b = $b_DEF;
    }
    unless( defined( $opt_M ) ) {
	$opt_M = $M_DEF;
    }
}
