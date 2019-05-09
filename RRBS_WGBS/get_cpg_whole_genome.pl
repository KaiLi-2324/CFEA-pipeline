#!/usr/bin/perl
use warnings;
use strict;

# this script is written to extract CpG positions from whole genome fasta file
# the input of this script is a whole genome fastq file, foe example, hg19.fa
# you can simply run the script like perl getCG-1.0.0.pl hg19.fa CpG.hg19.txt
# __author__ == wangxinyu

my $marker = "CG";										#marker

die "Thin program help you get whole genome CpGs position.
(perl) getCG genome.fa geneome_CpG.txt
(perl) getCG hg19/ genome_CpG.txt
	(require .fa files in hg19/)\n" unless(@ARGV);
my ( $infile , $outfile ) = @ARGV;
my ($in_fh,$out_fh);

$marker = uc($marker);									#To capital letter.

if(-e $infile and -d $infile){
	print "Input is a directory.\n";
	
	open $out_fh,">",$outfile or die "$outfile open ERROR!\n";
	
	my @files = glob("$infile/*.fa");
	foreach my $chr_file(@files){
		open $in_fh,"<",$chr_file or die "$chr_file open ERROR!\n";
		my ($chr , $chr_string) = ( "" , "" );
		
		#new chrom start	
		$chr = <$in_fh>;
		chomp $chr;
		$chr =~ s/>//;
		print $chr." is running...\n";
		
		while(my $line = <$in_fh>){
			chomp $line;
			$chr_string .= uc($line);
		}
		print_marker_position($chr,$chr_string,$marker,$out_fh);
	}
}
elsif(-e $infile){
	print "Input is a file.\n";
	
	open $in_fh,"<",$infile or die "$infile open ERROR!\n";
	open $out_fh,">",$outfile or die "$outfile open ERROR!\n";
	
	my ($chr , $chr_string) = ( "" , "" );
	while(my $line = <$in_fh>){
		chomp $line;
		if($line =~ />/){
			$line =~ s/>//;
			print $line."  is running...\n";
			if($chr){
				print_marker_position($chr,$chr_string,$marker,$out_fh);
			}
			
			#new chrom start
			$chr_string = "";
			$chr = $line;
		}
		else{
			$chr_string .= uc($line);
		}
	}
	#the last chrom
	print_marker_position($chr,$chr_string,$marker,$out_fh);
	
	close $in_fh;
	close $out_fh;
}
else{
	print "INPUT ERROR!\nCheck you input:\n$infile\n";
}

#print_marker_position($chr,$chr_string,$marker,$out_fh)
sub print_marker_position{
	my ($chr,$chr_string,$marker,$out_fh) = @_;
	
	my $position = index($chr_string,$marker) + 1;
	while($position ne 0){
		print $out_fh "$chr\t$position\n";
		$position = index($chr_string,$marker,$position) + 1;
	}
}
