#!/usr/bin/perl

# Seperates a fastq file of paired reads (for example, output from matchAdapters) into
# a first file containing the 5' reads and a second file containing the 3' reads

unless(@ARGV) {
	print "-i <file>\tInput file (combined PE PEAT data)\n";
	print "-o1 <file>\tFirst read output file.\n";
	print "-o2 <file>\tSecond read output file.\n";
	
	exit;
}

my $inputFile = "input.error";
my $outputFile1 = "output1.error";
my $outputFile2 = "output2.error";

for(my $i=0; $i<(@ARGV); $i++) {
	if($ARGV[$i] eq "-i") {
		$inputFile = $ARGV[$i+1];
	} elsif($ARGV[$i] eq "-o1") {
		$outputFile1 = $ARGV[$i+1];
	} elsif($ARGV[$i] eq "-o2") {
		$outputFile2 = $ARGV[$i+1];
	}
}

open(I, $inputFile) || die("Unable to open $inputFile for reading.\n");
open(O1, "> $outputFile1") || die("Unable to open $outputFile1 for writing.\n");
open(O2, "> $outputFile2") || die("Unable to open $outputFile2 for writing.\n");

my $write = 1;
my $counter = 0;

while($line1=<I>) {
	$line2 = <I>;
	$line3 = <I>;
	$line4 = <I>;
	$line5 = <I>;
	$line6 = <I>;
	$line7 = <I>;
	$line8 = <I>;

	chomp($line1);
	chomp($line2);
	chomp($line3);
	chomp($line4);
	chomp($line5);
	chomp($line6);
	chomp($line7);
	chomp($line8);

	if($line1 =~ m/F$/) {
		print O1 $line1 . "\n";
		print O1 $line2 . "\n";
		print O1 $line3 . "\n";
		print O1 $line4 . "\n";
		print O2 $line5 . "\n";
		print O2 $line6 . "\n";
		print O2 $line7 . "\n";
		print O2 $line8 . "\n";
	} elsif($line5 =~ m/F$/) {
		print O2 $line1 . "\n";
		print O2 $line2 . "\n";
		print O2 $line3 . "\n";
		print O2 $line4 . "\n";
		print O1 $line5 . "\n";
		print O1 $line6 . "\n";
		print O1 $line7 . "\n";
		print O1 $line8 . "\n";
	}
}
close(O2);
close(O1);
close(I);
