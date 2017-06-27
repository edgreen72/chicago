#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_s $opt_b );
use strict;

my ( $l, $dist );
my ( @el, @in, @out, @ss );
my $t_in  = 0; # total innies
my $t_out = 0; # total outies
my $t_ss  = 0; # total same-strand

&init();

### Open file-handle to the sesam2table.pl file of hits on
### the same chr/scaffold/etc
open( SAME, $opt_s ) or die( "$!: $_\n" );
while( chomp( $l = <SAME> ) ) {
    @el = split( "\t", $l );
    
    $dist = &find_dist( $el[2], $el[3], $el[7], $el[8] );
    # Innie, outie, or same-strand
    if ( $el[5] eq $el[10] ) { # Same strand
	$ss[int( $dist/$opt_b )]++;
	$t_ss++;
    }
    
    elsif ( $el[11] ) { # innie? check the last flag
	$in[int( $dist/$opt_b )]++;
	$t_in++;
    }
    
    else { # must be an outie (different strand, but not an innie)
	$out[int( $dist/$opt_b )]++;
	$t_out++;
    }
}
close( SAME );

### Write the output
printf( "### Fraction of short innies: %.2f\n", ($in[0] / ($t_ss + $t_in + $t_out)) );
printf( "### Short innies: %d\n", $in[0] );
printf( "### Total: %d\n", ($t_ss + $t_in + $t_out) );
print( "### Same-strand distance historgram for $opt_s.\n" );
print( "### Bin size $opt_b\n" );
print( "### Total same-strand hits: $t_ss\n" );
&print_array( \@ss, $t_ss );
print( "\n\n" );

print( "### Innies distance historgram for $opt_s.\n" );
print( "### Bin size $opt_b\n" );
print( "### Total innies hits: $t_in\n" );
&print_array( \@in, $t_in );
print( "\n\n" );

print( "### Outies distance historgram for $opt_s.\n" );
print( "### Bin size $opt_b\n" );
print( "### Total outies hits: $t_out\n" );
&print_array( \@out, $t_out );
print( "\n\n" );

sub init {
    my $b_DEF = 1000;
    getopts( 's:b:' );
    unless( -f $opt_s ) {
	print( "innies-outies-same-strand-lenhist.pl -s <SAME table>\n" );
	print( "                                     -b <bin size; default = $b_DEF>\n" );
	exit( 0 );
    }
    unless( defined( $opt_b ) ) {
	$opt_b = $b_DEF;
    }
}

### Takes the array of histogram data and prints it
sub print_array {
    my $ar_p = shift;
    my $tot  = shift;
    my $tot_seen = 0;
    my $bin;
    for( $bin = 0; $bin <= $#{ $ar_p }; $bin++ ) {
	$tot_seen += $ar_p->[$bin];
	printf( "%d %d %0.5f %.4f\n", 
		$bin, 
		$ar_p->[$bin], 
		$ar_p->[$bin]/$tot, 
		$tot_seen/$tot );
    }
}

### Takes the start and end coordinates for both read pairs
### and returns the distance between the outer coordinates
sub find_dist {
    my ( $s1, $e1, $s2, $e2 ) = @_;
    if ( $s1 < $s2 ) {
	return ( abs( $e2 - $s1 + 1 ) );
    }
    else {
	return ( abs( $e1 - $s2 + 1 ) );
    }
}
