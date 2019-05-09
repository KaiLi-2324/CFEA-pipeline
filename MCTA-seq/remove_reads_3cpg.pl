#!/usr/bin/perl
use strict;
use warnings;


# this script is designed to remove reads with more than 3 CpGs
my ($line1,$line2,$line3,$line4,$count);
my $usage = "perl $0 <IN> <OUT>
<IN>:fastq.gz
<OUT>:fastq.gz";
open IN,"zcat $ARGV[0]|" or die $!;
open OUT,"| gzip > $ARGV[1]" or die $!;
while (<IN>)
{
	$count=0;
	chomp;
	$line1 = $_;
	chomp($line2=<IN>);
	chomp($line3=<IN>);
	chomp($line4=<IN>);
	while ($line2 =~m/CA/g)
	{
		$count++;
	}
	while ($line2 =~m/CT/g)
        {
                $count++;
        }
	while ($line2 =~m/CC/g)
        {
                $count++;
        }
	next if ($count >= 3);
	print OUT "$line1\n$line2\n$line3\n$line4\n";
}
close IN;
