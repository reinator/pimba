#!/bin/perl

##Leandro Nascimento Lemos
##Victor Pylro
##Luiz Roesch
##10-04-2014 versao0.2
#QiimeToUparse
#How to use this script?
#perl bmp-otuName.pl -i something.fasta -o something.uparse.fasta
use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;


my $inputfile;
my $outputfile;


Getopt::Long::Configure ('bundling');
GetOptions ('i|input_file=s' => \$inputfile,
	    'o|output_file_prefix=s' => \$outputfile);

if(!defined($inputfile)) {
    die ("Usage: otuName.pl -i <input file> -o <output file>\n");
}

if(!defined($outputfile)) {
    die ("Usage: otuName.pl -i <input file> -o <output file>\n");
}


$/ = ">"; 

open(FASTA, $inputfile) or die "Cannot open the fasta file...";
my $f = $inputfile;
$f =~ s/\..*$//g;
open(SAIDA, ">" . $outputfile) or die "!";
$_ = <FASTA>;


	    while(my $inputfile = <FASTA>){
	chomp($inputfile);
	(my $header, my $sequence) = split("\n", $inputfile);
	(my $nome, my $campo02, my $campo3, my $campo4, my $barcode01, my $barcode02, my$campo6) = split(" ", $header);
		my $novo_nome = $nome;
	
	$novo_nome =~ s/\_.*$//g;
	
	my $range = 1000000000;
	my $random_number = int(rand($range));
	
	print SAIDA ">OTU$random_number\n$sequence\n";

}
close(SAIDA);
close(FASTA);
