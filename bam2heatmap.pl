#!/usr/local/bin/perl

use Getopt::Long;
use strict;

our( $outfile_prefix, # file path prefix for all output files
     $min_map_quality, # minimum MAPQ to use an alignment
     $min_scaffold_length, # minimum size to include a scaffold in heatmap
     $bin_size, # size of bins
     $scaffold_list_file, # path to file containing list of scaffolds to include
     $bamfile, # path to bamfile containing alignments
     @scaffolds_list ); # list of scaffolds to include

my ( $good_line, $f_line, $r_line, $bin, $i, $j, $chr );
my ( @sorted_chrs );
my ( %chr2len );
my ( %chr2bin );
my %heatmap;
my $HEATMAX;

&init(); # parse command line arguments

# if a scaffold list was provided, parse it into @scaffold_list
if ( defined $scaffold_list_file ) {
    open ( SCAFFOLDS, "<", $scaffold_list_file)
        or die "$scaffold_list_file is not readable.";
    while ( <SCAFFOLDS> ) {
        chomp;
        push @scaffolds_list, $_;
    }
}

# Use samtools view to split forward and reverse reads
if ( -e $bamfile ) {
    open( FSAM, "samtools view -h -f 64 $bamfile |" );
    open( RSAM, "samtools view -f 128 $bamfile |" );
} else {
    die "$bamfile is not readable."
}

# loop through forward read file's header, parsing the scaffold definition lines
# into %chr2len
chomp( $f_line = <FSAM> );
while ( $f_line =~ /^@/ ) {
    if ( $f_line =~ /^\@SQ/ ) {
        &parse_SN_header( $f_line );
    }
    chomp( $f_line = <FSAM> );
}
chomp( $r_line = <RSAM> );

# Sort chromosomes by length, unless a list of scaffolds was provided, in which
# case we'll just keep them in that order.
if ( @scaffolds_list ) {
    @sorted_chrs = @scaffolds_list
} else {
    @sorted_chrs = sort { $chr2len{$b} <=> $chr2len{$a} } (keys %chr2len);
}

### Populate %chr2bin with the first bin of each chromosome
$bin = 0;
foreach $chr ( @sorted_chrs ) {
    $chr2bin{ $chr } = $bin;
    $bin += int( $chr2len{ $chr } / $bin_size ) + 1;
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

my $heatmatrix_filename = $outfile_prefix . ".data";
my $gnuplot_filename = $outfile_prefix . ".gp";
my $png_filename = $outfile_prefix . ".png";

# print out the heat matrix
open ( HEATMATRIX, '>', $heatmatrix_filename );
for( $i = 0; $i <= $HEATMAX; $i++ ) {
    for( $j = 0; $j <= $HEATMAX; $j++ ) {
        if ( defined( $heatmap{$i}->{$j} ) ) {
            printf HEATMATRIX "%d %d %d\n", $i, $j, $heatmap{$i}->{$j};
        }
        else {
            printf HEATMATRIX "%d %d 0\n", $i, $j;
        }
    }
    print HEATMATRIX "\n";
}
close( HEATMATRIX );

# create the list of scaffolds and positions for the xtics/ytics in gnuplot
my @scaffold_strings = map { sprintf( "\"%s\" %d", $_, $chr2bin{ $_ }) } @sorted_chrs;
my $tics = join( ',', @scaffold_strings );
my $max = $chr2bin{ $sorted_chrs[$#sorted_chrs] } + 10;

# print out the gnuplot file
open( GNUPLOT, '>', $gnuplot_filename );
print GNUPLOT qq{set terminal png size 2048, 2048;
set output "$png_filename";
set palette rgbformula 34, 35, 36 negative;
unset colorbox;
set xtics rotate ($tics);
set ytics ($tics);
set grid front;
plot [0:$max][0:$max] '$heatmatrix_filename' using 1:2:(log(\$3)) notitle with image;};
close( GNUPLOT );

# finally, run gnuplot to get the output plot
system( "gnuplot", $gnuplot_filename );

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
    unless ( ($f_el[4] >= $min_map_quality) &&
             ($r_el[4] >= $min_map_quality) ) {
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

    $f_bin = $chr2bin{ $f_scaf } + int( $f_mid / $bin_size );
    $r_bin = $chr2bin{ $r_scaf } + int( $r_mid / $bin_size );

    $heatmap{$f_bin}->{$r_bin}++;
    $heatmap{$f_bin}->{$r_bin}++;

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

# parse a sam header line containing a scaffold definition
# Arguments: $sn_line -- string containing "@SQ SN:" header line
# Returns nothing, but adds scaffold to %chr2len if it belongs there.
sub parse_SN_header {
    my $sn_line = shift;
    my ( $chr, $len );

    ($chr, $len) = ( $sn_line =~ /SN:(\S+)\s+LN:(\d+)/);

    # if there is a list of scaffolds, only add this scaffold if it's part of
    # that list. Otherwise, just add it. Either way, though, only if it meets
    # the minimum length requirement
    if ( ! @scaffolds_list || grep( /^$chr$/, @scaffolds_list ) ) {
        if ( $len >= $min_scaffold_length ) {
            $chr2len{ $chr } = $len;
        }
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

# parses command line arguments
sub init {
    my $outfile_prefix_DEF = "heatmap";
    my $min_map_quality_DEF = 20;
    my $bin_size_DEF = 2000000;
    my $min_scaffold_length_DEF = 100000;
    (my $help_message =
    qq{bam2heatmap.pl <read name sorted bam file with both forward and reverse>
                   -p <outfile prefix; default = $outfile_prefix_DEF>
                   -m <map quality cutoff; default = $min_map_quality_DEF>
                   -b <bin size; default = $bin_size_DEF>
                   -M <minimum length of chr/scaf/contig; default = $min_scaffold_length_DEF>
                   -s <file with scaffolds to make heatmap of; one per line>
    Makes heatmap data file for gnuplot visualization.
    }) =~ s/^ {4}//gm;

    GetOptions('outfile-prefix|p' => \$outfile_prefix,
               'min-map-quality|m=i' => \$min_map_quality,
               'bin-size|b=i' => \$bin_size,
               'min-scaffold-length|M=i' => \$min_scaffold_length,
               'scaffold-list-file|s=s' => \$scaffold_list_file);
    # print help message if there's no argument left after parsing switches
    $bamfile = shift @ARGV or die $help_message;

    # set default values for arguments not specified on command line
    unless( defined( $outfile_prefix ) ) {
        $outfile_prefix = $outfile_prefix_DEF
    }
    unless( defined( $min_map_quality ) ) {
        $min_map_quality = $min_map_quality_DEF;
    }
    unless( defined( $bin_size ) ) {
        $bin_size = $bin_size_DEF;
    }
    unless( defined( $min_scaffold_length ) ) {
        $min_scaffold_length = $min_scaffold_length_DEF;
    }
}
