#!/bin/perl

##Leandro Nascimento Lemos
##Victor Pylro
##Luiz Roesch
##10-04-2014 versao0.2
#QiimeToUparse
#How to use this script?
#perl bmp-otuName.pl -i something.fasta -o something.uparse.fasta
#!/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $inputfile;
my $outputfile;

Getopt::Long::Configure('bundling');
GetOptions(
    'i|input_file=s'         => \$inputfile,
    'o|output_file_prefix=s' => \$outputfile
);

if (!defined($inputfile)) {
    die("Usage: otuName.pl -i <input file> -o <output file>\n");
}

if (!defined($outputfile)) {
    die("Usage: otuName.pl -i <input file> -o <output file>\n");
}

$/ = ">";

open(FASTA, $inputfile) or die "Cannot open the fasta file...";
open(SAIDA, ">" . $outputfile) or die "!";

$_ = <FASTA>;  # Skip the first empty record
my %used_random_numbers;  # Hash to store used random numbers
my $range = 1000000000;

while (my $inputfile = <FASTA>) {
    chomp($inputfile);
    my ($header, $sequence) = split("\n", $inputfile, 2);
    $sequence =~ s/\n//g;  # Remove any newline characters from the sequence
    
    my $random_number;
    
    do {
        $random_number = int(rand($range));
    } while (exists $used_random_numbers{$random_number});
    
    $used_random_numbers{$random_number} = 1;  # Mark the random number as used
    
    print SAIDA ">OTU$random_number\n$sequence\n";
}

close(SAIDA);
close(FASTA);
