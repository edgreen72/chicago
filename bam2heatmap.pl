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
my ( %chr2len, # hash of scaffold name (string) to length (bp, int)
     %chr2bin, # hash of scaffold name to bin location (int)
     %chr2maxbin, # scaf name to last bin
     %chr2orientation); # hash of scaf name to orientation ("+"/"-")
my %heatmap;
my $HEATMAX;

&init(); # parse command line arguments

# if a scaffold list was provided, parse it into @scaffold_list
if ( defined $scaffold_list_file ) {
    open ( SCAFFOLDS, "<", $scaffold_list_file)
        or die "$scaffold_list_file is not readable.";
    while ( <SCAFFOLDS> ) {
        chomp;
        # if the scaffold name is preceded by an orientation, store this
        # orientation in addition to the scaffold name
        my ($scaffold_name, $orientation);
        if ( /^([\+-])(.+)$/ ) {
            $orientation = $1;
            $scaffold_name = $2;
        } else { # otherwise, the orientation is '+' by default
            $orientation = "+";
            $scaffold_name = $_;
        }
        push @scaffolds_list, $scaffold_name;
        $chr2orientation{ $scaffold_name } = $orientation;
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
    $chr2maxbin{ $chr } = $bin;
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
    my ($f_scaf, $f_start, $f_end, $f_mq, $f_strand, $f_mid);
    my ($r_scaf, $r_start, $r_end, $r_mq, $r_strand, $r_mid);
    my ( $f_bin, $r_bin );

    # make lists of sam fields in the forward and reverse reads
    my @forward_read_fields = split( "\t", $f_line );
    my @reverse_read_fields = split( "\t", $r_line );

    ### Are these the same read? Check their IDs
    ### If not, return false
    unless ( $forward_read_fields[0] eq $reverse_read_fields[0] ) {
        return 0; # We're done & stop
    }

    ### Do both reads pass the map-quality filter?
    unless ( ($forward_read_fields[4] >= $min_map_quality) &&
             ($reverse_read_fields[4] >= $min_map_quality) ) {
        return 1; # We're done, but keep on
    }

    ### Are both reads mapped to a sequence that passed
    ### length filter?
    unless( $chr2len{ $forward_read_fields[2] } &&
            $chr2len{ $reverse_read_fields[2] } ) {
        return 1;
    }

    # This read pair is fine, so we can add it to the proper bin

    # First, parse the sam lines into variables
    ($f_scaf, $f_start, $f_end, $f_mq, $f_strand) =
        &sam_line2table_data( @forward_read_fields );
    ($r_scaf, $r_start, $r_end, $r_mq, $r_strand) =
        &sam_line2table_data( @reverse_read_fields );

    # find the midpoints of the alignments (mean of start and end position)
    $f_mid = ($f_start + $f_end) / 2;
    $r_mid = ($r_start + $r_end) / 2;

    # determine which bins these alignments fall into
    $f_bin = &get_bin( $f_scaf, $f_mid );
    $r_bin = &get_bin( $r_scaf, $r_mid );

    # increment the heatmap matrix for this pair of bins, symmetrically
    $heatmap{$f_bin}->{$r_bin}++;
    $heatmap{$f_bin}->{$r_bin}++;

    return 1;
}

# get the bin number for a given location on a scaffold
# Arguments: scaffold -- name of scaffold, for calculating offset
# location -- base pair location of alignment midpoint on scaffold (int)
# Returns: the index of the proper bin
sub get_bin {
    my ( $scaffold, $location ) = @_;
    if ( $chr2orientation{ $scaffold } eq '+' ) {
        return $chr2bin{ $scaffold } + int( $location / $bin_size );
    } else {
        return $chr2maxbin{ $scaffold } - int( $location / $bin_size );
    }
}

# parse a list of sam fields into variables
# Arguments: an array of fields in a single line of a sam file, e.g., from
#            split ( "\t", $sam_file_line ).
# Returns: an array containing:
# - name of scaffold to which this read aligns
# - starting location of alignment on this scaffold, in bp
# - ending location of alignment, in bp
# - MAPQ of this alignment
# - strand of this alignment ('+' or '-')
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

# given a samflag (int), return '+' if it indicates an alignment on the forward
# strand or '-' on the reverse strand
sub bit_flag2strand {
    return ( ( $_ & 16  ) ? '-' : '+' );
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

    GetOptions('outfile-prefix|p=s' => \$outfile_prefix,
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
