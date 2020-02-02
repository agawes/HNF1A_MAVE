#!/usr/bin/perl
use strict;

my $RG_info="RG_info";

open(IN,'<',$RG_info) or die;
while(my $line = <IN>){
	chomp $line;
	my @line = split("\t", $line);
	my $name=$1 if $line[2] =~ m/SM="(.+)"/;
	system "cat raw_fastq/$line[0]\_1.fastq.gz >> concat_fastq/$name.1.fastq.gz";
	system "cat raw_fastq/$line[0]\_2.fastq.gz >> concat_fastq/$name.2.fastq.gz";

}
close IN;
