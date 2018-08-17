#!/usr/bin/env perl

# Program: mapcount.pl v.2.2
# Count mapping hits against reference fasta (SAM file)
# Murray Cox (with some code contributions by Daniel Peterson)
# <murray.p.cox@gmail.com>
# March 2015 (code modified from September 2009 onward)

use strict;
use warnings;
use Getopt::Long;

# usage
my $usage = "$0 -f|fasta fasta_file -s|sam SAM_file -o|output filename [-t|trim] [-l|length 25] [-p|pair]\n
-f|fasta	fasta reference file
-s|sam		SAM file
-o|output	output filename
-t|trim		trim ID before a space
-l|length	minimum match length for counting
-p|pair		count paired matches as one match
\n";

# global values
my $fasta;
my $map;
my $trim;
my $outfile = "mapping_distribution.dat";
my $min_length = 0;
my $pair = 0;

# read user options
GetOptions(
	"f|fasta=s"   => \$fasta,
	"s|sam=s"     => \$map,
	"o|output=s"  => \$outfile,
	"t|trim"      => \$trim,
	"l|length=i"  => \$min_length,
	"p|pair"      => \$pair
);

# check for correct usage
if( !$fasta || !$map ){
	die $usage;
}

# open reference sequence file
unless( -e $fasta ){
	die "error: file $fasta does not exist\n";
}
open( FASTA, "<$fasta" ) or die "error: file $fasta cannot be opened: $!\n";

# check file format
my $return_value = &get_format(*FASTA);

# get and print sequences
my @returns = ();
my @ids = ();
until( eof(FASTA) ){
	
	@returns = &get_next(*FASTA, $return_value);
	
	my $id = $returns[0];
	
	if( $trim ){
		
		my @segments = split( /\s+/, $id );
		$id = $segments[0];
	}
	
	push( @ids, $id );	
}

# close reference sequence file
close FASTA or die "error: file $fasta cannot be closed: $!\n";

# open input SAM file
unless( -e $map ){
	die "error: file $map does not exist\n";
}
open( MAP, "<$map" ) or die "error: cannot open $map for reading: $!\n";

# open output counts file
if( -e $outfile ){
	die "error: file $outfile already exists\n";
}
open( OUT, ">$outfile" ) or die "error: cannot open $outfile for writing: $!\n";

# make hash
my %counter = ();
foreach my $entry ( @ids ){
	$counter{$entry} = 0;
}

my $chaff = 0;
my $double_hit = 0;
my $unmapped = 0;
my $previous_match = "";

# step through map output, skipping first line
while( <MAP> ){
	
	my $line = $_;
	my @line_array = split( /\s+/, $line );
	
	if( $line !~ m/^@/){
		
		my $match_name = $line_array[2]; 
		
		if( $match_name eq "*" ){
			$unmapped++;
			next;
		}
		
		my @bit_flag = &getBits($line_array[1], 11);
		my $first_of_pair = $bit_flag[6];
		my $second_of_pair = $bit_flag[7];
		
		if( ( $first_of_pair == $second_of_pair ) && $pair ){
			die "error: bitflags not set correctly at $line_array[0]\n";
		}
	
		if( !$pair ){
			$first_of_pair = 1;
		}
	
		if( $trim ){
			
			my @segments = split( /\s+/, $match_name );
			$match_name = $segments[0];
		}
		
		my $match_length = 0;
		my @match_array = split ( //, $line_array[5] );
		my $current_sum = 0;
 
		foreach( @match_array ){
			if( $_ =~ /\d/ ){
				$current_sum = ($current_sum*10) + $_;
			}elsif( $_ eq "M" ){
				$match_length += $current_sum;
				$current_sum = 0;
			}else{
				$current_sum = 0;
			}
		}
	
		if( $first_of_pair && ( $match_length >= $min_length ) ){
			$counter{$match_name}++;
			$previous_match = $match_name;
			
		}elsif( $second_of_pair && ( $match_length >= $min_length ) && ( $match_name ne $previous_match ) ){
			$counter{$match_name}++;
	
		}elsif( $first_of_pair ){
			$chaff++;
		
		}elsif( $second_of_pair && ( $match_length >= $min_length ) ){
			$double_hit++;	
		
		}else{
			$chaff++;
		}
	
		if( $second_of_pair ){
			$previous_match = "";
		}
	}
}

# print counts to file
print OUT "id\tcount\n";
foreach my $value ( @ids ){
	print OUT $value, "\t", $counter{$value}, "\n";
}

print STDOUT "Unmapped reads:\t$unmapped\n";

if( $min_length ){
	print STDOUT "Reads failing length cutoff:\t$chaff\n";
}

if( $pair ){
	my $true_doubles = $double_hit;
	print STDOUT "Paired end hits:\t$true_doubles\n";
}

# check if anything left over
my @difference = ();
my %count = ();
foreach my $element (@ids, keys %counter) { $count{$element}++ }
foreach my $element (keys %count) {
	if( $count{$element} <= 1 ){
		push( @difference, $element );
	}
}

# print difference to screen
if( @difference ){
	print STDOUT "\nSequence names not matching fasta file...\n\n";
	print STDOUT "id\tcount\n";
	foreach my $entry ( @difference ){
		print STDOUT $entry, "\t", $counter{$entry}, "\n";
	}
}

# close remaining files and terminate
close MAP or die "error: cannot close $map: $!\n";
close OUT or die "error: cannot close $outfile: $!\n";
exit 0 or die "error: $0 ended abnormally: $!\n";


# get bit information
sub getBits {
        my $byte;
        my $num_bits;
        ( $byte, $num_bits ) = @_;
        my @bits;
        for( my $i = 0; $i < $num_bits; $i++){
	        push(@bits, ( $byte & 2**$i ) );
	}
	return @bits;
}

# determine FASTA/FASTQ format
# usage: &get_format(*FILEHANDLE)
# return: "fasta" or "fastq"
#         0 on failure
sub get_format(*){
	
	# set function variables
	local *FILEHANDLE = shift;
	
	# retrieve file position
	my $position = tell FILEHANDLE;
	
	# retrieve first line
	seek(FILEHANDLE, 0, 0);
	my $first_line = <FILEHANDLE>;
	
	# retrieve first character
	my $first_character = substr($first_line, 0, 1);
	
	# reset filehandle
	seek(FILEHANDLE, $position, 0);
	
	# return format
	if( $first_character eq ">" ){
		return "fasta";
	}elsif( $first_character eq "@" ){
		return "fastq";
	}else{
		return 0;
	}
}

# retrieve next FASTA/FASTQ entry
# usage: &get_next(*FILEHANDLE, $format)
# return: @array
#         $array[0] = header, $array[1] = sequence, $array[2] = quality
#         0 on failure
sub get_next(*$){
	
	# set function variables
	local *FILEHANDLE = shift;
	my $format = shift;
	
	my @data = ("", "", "");
	
	if( $format eq "fasta" ){
		
		# retrieve first line
		my $first_line = <FILEHANDLE>;
		chomp($first_line);
		
		# retrieve file position
		my $position = tell FILEHANDLE;
		
		# retrieve header
		$data[0] = substr($first_line, 1, length($first_line)-1);
		
		# retrieve sequence
		until( eof(FILEHANDLE) ){
			
			# retrieve line
			my $line = <FILEHANDLE>;
			chomp($line);
			
			# retrieve first character
			my $first_character = substr($line, 0, 1);
			
			# step through multiline fasta
			if( $first_character eq ">" ){
				seek(FILEHANDLE, $position, 0);
				last;
			}else{
				$data[1] .= $line;
				$position = tell FILEHANDLE;
			}
		}
	}elsif( $format eq "fastq" ){
		
		# retrieve first line
		my $first_line = <FILEHANDLE>;
		chomp($first_line);
		
		# retrieve header
		$data[0] = substr($first_line, 1, length($first_line)-1);
		
		# retrieve sequence
		my $second_line = <FILEHANDLE>;
		chomp($second_line);
		$data[1] = $second_line;

		# ignore next header line
		my $third_line = <FILEHANDLE>;

		# retrieve quality
		my $fourth_line = <FILEHANDLE>;
		chomp($fourth_line);
		$data[2] = $fourth_line;
	}else{
		return 0;
	}
	
	# return data
	return @data;
}

# print FASTA/FASTQ entry
# usage: &print_data(*FILEHANDLE, $format, $fasta_wrap, @data)
# $fasta_wrap gives number of characters to wrap, or 0 for no wrapping
# return: 0 on failure
sub print_data(*$$@){

	# set function variables
	local *FILEHANDLE = shift;
	my $format = shift;
	my $fasta_wrap = shift;
	
	my $header = shift;
	my $sequence = shift;
	my $quality = shift;
	
	if( $format eq "fasta" ){
		
		print FILEHANDLE ">", $header, "\n";
		
		if( $fasta_wrap == 0 ){
			print FILEHANDLE $sequence, "\n";
		}else{
			for(my $i = 0; $i < length $sequence; $i += $fasta_wrap){
				print FILEHANDLE substr($sequence, $i, $fasta_wrap), "\n";
			}
		}
	}elsif( $format eq "fastq" ){
		
		if( !$quality ){
			die "error in print_data(): sequences not in fastq format\n";
		}
		
		print FILEHANDLE "@", $header, "\n";
		print FILEHANDLE $sequence, "\n";
		print FILEHANDLE "+", $header, "\n";
		print FILEHANDLE $quality, "\n";
	}else{
		return 0;
	}

	return 1;
}
