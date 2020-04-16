#!/usr/bin/perl -w
#
# scrubSAMfile.pl
# Version of January 15, 2018.  Volker Brendel

use strict;
use Getopt::Std;


my $USAGE="\nUsage: $0 [-p] [-i min_intronlength] [-I min_longintron]\n
           [-t max_TLEN]  [-o outfile] -s SAMfile\


** This script does selects reads or read pairs from a SAM file by the
   following criteria:
     - reads must be paired
     - the mappings should not allow presumed introns longer than min_longintron
     - the SAM TLEN value should not exceed max_TLEN
   \n";


my %args;
getopts('pi:I:t:o:s:', \%args);

my ($SAMfile,$outputfile,$outputfile1,$outputfile2,$outputfile3,
             $outputfile4,$outputfile5);
my $isPaired = 0;
my $min_longintron   = 1000;
my $min_intronlength = 50;
my $max_TLEN         =  500;
my ($FHINF,$FHOUTF);

if (defined($args{p})) {
  $isPaired = 1;
}
if (!defined($args{s})) {
  print "\nPlease specify a SAM input file (sorted by read name of format id#0:T/F).\n\n";
  die $USAGE;
}
else {
  $SAMfile = $args{s};
  if (! -e $SAMfile) {
    print "\nSAM file $SAMfile does not exist.\n\n";
    die $USAGE;
  }
}

if (defined($args{i})) {
  $min_intronlength = $args{i};
}
print "\nMinimal 'intron' size (for special reporting) set to :\t$min_intronlength";
if (defined($args{I})) {
  $min_longintron = $args{I};
}
print "\nLong 'intron' threshold set to :\t$min_longintron";
if (defined($args{t})) {
  $max_TLEN = $args{t};
}
print "\nMaximal TLEN (column 9 value) set to :\t$max_TLEN\n";

if (!defined($args{o})) {
  $outputfile = $SAMfile;
  print "\nOutput file base set to:\t$outputfile\n";
}
else {
  $outputfile = $args{o};
  print "\nOutput file base set to:\t$outputfile\n";
}



open FHINF,  "<$SAMfile" || die ("Cannot open file: $SAMfile"); 
  if ($isPaired == 0) {
    $outputfile1 = $outputfile . ".fq";
    open FHOUTF1,  ">$outputfile1" || die ("Cannot open file: $outputfile1"); 
  }
  else {
    $outputfile1 = "Singlets-" .$outputfile;
    open FHOUTF1,  ">$outputfile1" || die ("Cannot open file: $outputfile1"); 
    $outputfile2 = "LongIntrons-" . $outputfile;
    open FHOUTF2,  ">$outputfile2" || die ("Cannot open file: $outputfile2"); 
    $outputfile3 = "LongTLEN-" . $outputfile;
    open FHOUTF3,  ">$outputfile3" || die ("Cannot open file: $outputfile3"); 
    $outputfile4 = "ShortIntrons-" . $outputfile;
    open FHOUTF4,  ">$outputfile4" || die ("Cannot open file: $outputfile4"); 
    $outputfile5 = "Scrubbed-" . $outputfile;
    open FHOUTF5,  ">$outputfile5" || die ("Cannot open file: $outputfile5"); 
  }



my $line1     = "";
my $line2     = "";
my $tmpline   = "";
my @read1vars;
my @read2vars;
my $read1name = "";
my $read2name = "";
my $read1nbrN;
my $read2nbrN;
my $read1TLEN;
my $read2TLEN;

my $skip = 0;
my $readstotal = 0;
my $singletreads = 0;
my $readslongintron = 0;
my $readslongTLEN = 0;
my $readsshortintron = 0;
my $readsaccepted = 0;
my $readsrejected = 0;

while (defined($line1=<FHINF>)) {
  if ($line1 =~ /^@/) {
    next;
  }
  $readstotal++;
  if ($isPaired) {
    if ($skip) {
      $tmpline = $line1; $line1 = $line2; $line2 = $tmpline;
    } else {
      defined($line2=<FHINF>) || die ("trouble");
      $readstotal++;
    }
    @read1vars = split('\t',$line1);
    @read2vars = split('\t',$line2);
    ($read1name) = $read1vars[0] =~ m/([^#]*)#.*/;
    ($read2name) = $read2vars[0] =~ m/([^#]*)#.*/;
    if ($read1name ne $read2name) {
      print FHOUTF1 $line1;
      $singletreads++;
      $skip = 1;
    } else {
      $skip = 0;
      ($read1nbrN) = $read1vars[5] =~ m/\d*M(\d*)N.*/;
      ($read2nbrN) = $read2vars[5] =~ m/\d*M(\d*)N.*/;
      $read1TLEN = $read1vars[8];
      $read2TLEN = $read2vars[8];
      if ((defined($read1nbrN) &&  $read1nbrN >= $min_longintron)  || 
          (defined($read2nbrN) &&  $read2nbrN >= $min_longintron)  
         ) {
            print FHOUTF2 $line1;
            print FHOUTF2 $line2;
            $readslongintron++;
	    next;
      }
      if (abs($read1TLEN) > $max_TLEN || 
          abs($read2TLEN) > $max_TLEN
         ) {
            print FHOUTF3 $line1;
            print FHOUTF3 $line2;
            $readslongTLEN++;
	    next;
      }
      if ((defined($read1nbrN) &&  $read1nbrN >= $min_intronlength)  || 
          (defined($read2nbrN) &&  $read2nbrN >= $min_intronlength)  
         ) {
            print FHOUTF4 $line1;
            print FHOUTF4 $line2;
            $readsshortintron++;
      }
      print FHOUTF5 $line1;
      print FHOUTF5 $line2;
      $readsaccepted++;
    }
  }
}

$readsrejected = $singletreads + 2*($readslongintron + $readslongTLEN);
printf ("\n\nTotal reads checked: %10d\n\n", $readstotal);
if ($isPaired) {
  printf ("Accepted read pairs            : %10d (%6.2f%%)\n", $readsaccepted,  200*$readsaccepted/$readstotal);
  printf ("Rejected reads                 : %10d (%6.2f%%)\n", $readsrejected,  100*$readsrejected/$readstotal);

  printf ("\n\n  Rejection categories:\n\n");
  printf ("  Singlet reads                : %10d (%6.2f%%)\n", $singletreads,  100*$singletreads/$readstotal);
  printf ("  Read pairs with long introns : %10d (%6.2f%%)\n", $readslongintron,  200*$readslongintron/$readstotal);
  printf ("  Read pairs with long TLEN    : %10d (%6.2f%%)\n", $readslongTLEN,  200*$readslongTLEN/$readstotal);

  printf ("\n\nAccepted reads with special features:\n\n");
  printf ("  Read pairs with short introns: %10d (%6.2f%%)\n", $readsshortintron,  200*$readsshortintron/$readstotal);
}
