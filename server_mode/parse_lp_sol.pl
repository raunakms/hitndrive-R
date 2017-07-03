#######################################################################
### THIS SCRIPT PARSES THE CPLEX SOL FILE TO EXTRACT REQUIRED DATA ###
#######################################################################

### FIRST ARGUMENT : PATH TO SOL FILES
### SECOND ARGUMENT : PATH TO PARSED FILES
### THIRD ARGUMENT : BATCH NAME

#!/usr/bin/perl

use strict;
use warnings;

use File::Temp;

my ($solution_dir, $output_dir, $batch_name) = @ARGV;

my @files = <$solution_dir/*>;

foreach my $sol_file(@files) {
	my $name = $batch_name."_";
	my ($sample_name) = $sol_file =~ /($name\d+[^.sol])/;
	
	my $out_file = $output_dir."/".$sample_name.".txt";
	
	my $tmp = File::Temp->new();
	my $tmp_filename = $tmp->filename;

	system("grep 'variable name' $sol_file > $tmp_filename");

	open(FILE, $tmp_filename) or die "ERROR: CANNOT OPEN THE TEMPORARY FILE";
	open(OUT_FILE, ">$out_file") or die "out_file error";

	my $sol_count = 1;
	my $const_count = 1;
	
	#print OUT_FILE "Solution", "\t", "Type", "\t", "Variable", "Value";
		
	while(<FILE>){
		chomp;
		my $line = $_;
		my ($var) = $line =~ /<variable name="([^"]+)/;
		my ($val) = $line =~ /value="([^"]+)/;

		if($var =~ m/@/){
			my @fields = split("@", $var);
			print OUT_FILE $sol_count, "\t", "NODE", "\t", $fields[1], "\t", $val, "\n";
		}else{
			my @fields1 = split(";", $var);
			if($fields1[0] eq "Y"){
				print OUT_FILE $sol_count, "\t", "EDGE", "\t", $var, "\t", $val, "\n";
			}

			if($fields1[0] eq "Z"){
				print OUT_FILE $sol_count, "\t", "ZETA", "\t", $var, "\t", $val, "\n";
			}
		}
		$const_count++;
	}	
	close FILE;
	close OUT_FILE;

	print "FILE GENERATED :",$out_file, "\n";
}

