#!/usr/bin/perl -w
#
# select.pl
# Version of December 21, 2017.  Volker Brendel

use strict;
use Getopt::Std;


my $USAGE="\nUsage: $0 [-p] [-i] [-o outfile] -r readfile\


** This script does selects reads or read pairs from a FASTQ input
   file based on the condition that the read (or first read for read
   pair input) starts with a G (indicating capped transcript end in
   CAGE experiments).
   \n";


my %args;
getopts('pio:r:', \%args);

my ($readfile,$outputfile,$outputfile1,$outputfile2);
my $isPaired = 0;
my $outInterleafed = 0;
my ($FHINF,$FHOUTF);

if (defined($args{p})) {
  $isPaired = 1;
}
if (defined($args{i})) {
  $outInterleafed = 1;
}
if (!defined($args{r})) {
  print "\nPlease specify a read input file.\n\n";
  die $USAGE;
}
else {
  $readfile = $args{r};
 if (! -e $readfile) {
   print "\nRead file $readfile does not exist.\n\n";
   die $USAGE;
 }
}
if (!defined($args{o})) {
  $outputfile = $readfile . "_select";
  print "\nOutput file base set to:\t$outputfile\n\n";
}
else {
  $outputfile = $args{o};
  print "\nOutput file base set to:\t$outputfile\n\n";
}


if ($isPaired == 0) {
  $outInterleafed = 0;
}



open FHINF,  "<$readfile" || die ("Cannot open file: $readfile"); 
if ($outInterleafed) {
  open FHOUTF,  ">$outputfile" || die ("Cannot open file: $outputfile"); 
}
else {
  if ($isPaired == 0) {
    $outputfile1 = $outputfile . ".fq";
    open FHOUTF1,  ">$outputfile1" || die ("Cannot open file: $outputfile1"); 
  }
  else {
    $outputfile1 = $outputfile . "_R1.fq";
    open FHOUTF1,  ">$outputfile1" || die ("Cannot open file: $outputfile1"); 
    $outputfile2 = $outputfile . "_R2.fq";
    open FHOUTF2,  ">$outputfile2" || die ("Cannot open file: $outputfile2"); 
  }
}



my $line1     = "";
my $line2     = "";
my $line3     = "";
my $line4     = "";
my $line5     = "";
my $line6     = "";
my $line7     = "";
my $line8     = "";

my $readstotal = 0;
my $readsaccepted = 0;
my $readsrejected = 0;

while (defined($line1=<FHINF>)) {
  if ($line1 =~ /^@/) {
    defined($line2=<FHINF>) || die ("trouble");
    defined($line3=<FHINF>) || die ("trouble");
    defined($line4=<FHINF>) || die ("trouble");
    if ($isPaired) {
      defined($line5=<FHINF>) || die ("trouble");
      defined($line6=<FHINF>) || die ("trouble");
      defined($line7=<FHINF>) || die ("trouble");
      defined($line8=<FHINF>) || die ("trouble");
    }
    $readstotal++;

    if ($line2 =~ /^G/) {
      $readsaccepted++;
      if ($outInterleafed) {
        print FHOUTF $line1;
        print FHOUTF $line2;
        print FHOUTF $line3;
        print FHOUTF $line4;
        if ($isPaired) {
          print FHOUTF $line5;
          print FHOUTF $line6;
          print FHOUTF $line7;
          print FHOUTF $line8;
        }
      }
      else {
        print FHOUTF1 $line1;
        print FHOUTF1 $line2;
        print FHOUTF1 $line3;
        print FHOUTF1 $line4;
        if ($isPaired) {
          print FHOUTF2 $line5;
          print FHOUTF2 $line6;
          print FHOUTF2 $line7;
          print FHOUTF2 $line8;
        }
      }
    }
    else {
      $readsrejected++;
    }
  }
  else {
    print FHOUTF "should not happen";
  }
}

printf ("\n\nTotal reads checked: %10d\n\n", $readstotal);
printf ("Reads accepted     : %10d (%6.2f%%)\n", $readsaccepted,  100*$readsaccepted/$readstotal);
printf ("Reads rejected     : %10d (%6.2f%%)\n", $readsrejected,  100*$readsrejected/$readstotal);
