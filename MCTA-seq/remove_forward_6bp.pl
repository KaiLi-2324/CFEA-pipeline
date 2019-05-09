#! /usrt/bin/perl -w
use strict;
use warnings FATAL => 'all';


# this script is designed to remove forward xbp, backward zbp, remain more than ybp
# for read1, remove the forward 6bp and backward 12 bp, and for read2, remove forward
# 12bp and backward 6bp

my $usage = "perl $0 <IN1> <OUT> x y z";
die $usage unless @ARGV==5;
open IN,"zcat $ARGV[0] |" or die $!;
open OUT,"| gzip > $ARGV[1]" or die $!;
while(<IN>)
{
	chomp;
	my $line1=$_;
	chomp(my $line2=<IN>);
	chomp(my $line3=<IN>);
	chomp(my $line4=<IN>);
	my $len_2=length($line2);
	if($len_2 >=$ARGV[2]+$ARGV[4]+$ARGV[3])
	{
		my $second=substr($line2,$ARGV[2],$len_2-$ARGV[2]-$ARGV[4]);
		my $line4_1=substr($line4,$ARGV[2],$len_2-$ARGV[2]-$ARGV[4]);
		my $len=length($second);
		#my $third=substr($second,0,$len_2-24);
		#my $line4_2=substr($line4_1,0,$len_2-24);
		print OUT "$line1\n$second\n$line3\n$line4_1\n";
	}
	else
	{
		next;
	}
}
close IN;
close OUT;
