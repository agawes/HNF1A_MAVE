#!/usr/bin/perl
#use strict;
use File::Basename;
use Getopt::Long;


my $MAPQ=30;
my $baseQ=20;

GetOptions("sam=s" => \$sam, "MAPQ=i" => \$MAPQ, "baseQ=i" => \$baseQ);

my $base= basename($sam, ".sam");

my %AA= ( "ATT" => "I", "ATC" => "I", "ATA" => "I",
	"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L", "TTA" => "L", "TTG" => "L",
	"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
	"TTT" => "F", "TTC" => "F",
	"ATG" => "M",
	"TGT" => "C", "TGC" => "C",
	"GCT" => "A", "GCC" => "A", "GCA" => "A","GCG" => "A",
	"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G",
	"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
	"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
	"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "AGT" => "S", "AGC" => "S",
	"TAT" => "Y", "TAC" => "Y",
	"TGG" => "W",
	"CAA" => "Q", "CAG" => "Q",
	"AAT" => "N", "AAC" => "N",
	"CAT" => "H", "CAC" => "H",
	"GAA" => "E", "GAG" => "E",
	"GAT" => "D", "GAC" => "D",
	"AAA" => "K", "AAG" => "K",
	"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R", "AGA" => "R", "AGG" => "R",
	"TAA" => "*", "TAG" => "*", "TGA" => "*"
);

### read the HNF1A reference sequence
open(IN,'<','/well/mccarthy/production/genetic-screens/resources/HNF1A_whole_construct.fna') or die;
my $fline= <IN>;
my $HNF1A_seq=<IN>;
chomp $HNF1A_seq;
close IN;

my @codon_pos;
my %ref_seq;
## i - pos within HNF1A seq; 113 - 1st codon; 2003 - codon 631
for (my $i=113; $i<=2003; $i+=3){
	push(@codon_pos, $i);
	$ref_seq{$i} = substr($HNF1A_seq, $i-1,3);
	#print "$i\t",$ref_seq{$i},"\n";
}


my %detected_variants; ### save all non-ref alleles here; 
					### key: Pos_WT_ALT 
					### value: supporting read id
					### then can be sorted on key; and allele count by count of unique read id's (in this way we don't count twice the same variant
					### if appears in both reads from the pair

open(IN,'<',$sam) or die;
my $id='';
LINE:while(my $line=<IN>){ 
	if ($line !~ /^@/){
	chomp $line;
	my @line=split("\t",$line);
	if ($line[4]>=$MAPQ){ 

#	print $line,"\n";	
	my $read_id = $line[0];
#	print $read_id,"\n";
	my $start = $line[3];
	my $cigar = $line[5];
	my $seq = $line[9];
	my $qual = $line[10];
	
	if ($cigar == "*"){next LINE;}
	# hash with bases called at each reference position; and QUAL
	#print "##### $start\t$cigar ######\n";

	my %hash;
	if ($cigar =~m/[MSID]/){ #next LINE;

		(my @cigar) = $cigar=~ m/(\d+[IDSM])/g;
		my $pos_in_read=1;
		my $pos_in_ref = $start;
		foreach my $c(@cigar){
			my ($count, $cig)= $c=~ m/(\d+)([IDSM])/;
			#print "####### $count\t$cig #######\n";
			if ($cig eq "M"){
				for (my $k=$pos_in_read; $k<($pos_in_read+$count); $k++){
					$hash{$pos_in_ref}{'seq'}=substr($seq,$k-1,1);
					#print "CHECK: ",substr($seq,$k-1,1),"\t",$hash{$pos_in_ref}{'seq'},"\n";
					$hash{$pos_in_ref}{'qual'}=ord(substr($qual,$k-1,1))-33;
					#print $k,"\t",$pos_in_ref,"\t*",$hash{$pos_in_ref}{'seq'},"*\t",substr($HNF1A_seq,$pos_in_ref-1,1),"\n";
					$pos_in_ref++;
				}
				$pos_in_read=$pos_in_read+$count;	
				#print "# pos in read: $pos_in_read\tpos in ref: $pos_in_ref\n";		
			}
			if ($cig eq "I"){
				for (my $k=$pos_in_read; $k<($pos_in_read+$count); $k++){
					#print $k,"\tI_",$pos_in_ref,"\t","I","\t","I","\n";
				}
				$pos_in_read=$pos_in_read+$count;					
				#print "# pos in read: $pos_in_read\tpos in ref: $pos_in_ref\n";		
				
			}						
			if ($cig eq "D"){ ## do nothing to hash; this is not in the read
				$hash{$pos_in_ref}{'seq'}="*";
				$hash{$pos_in_ref}{'qual'}=10;
				#print $pos_in_read,"\t",$pos_in_ref,"\tD\n";	
				$pos_in_ref=$pos_in_ref+$count;	
				#print "# pos in read: $pos_in_read\tpos in ref: $pos_in_ref\n";		

			}
			if ($cig eq "S"){ ## don't save these positions in hash; this is soft-clipped
				#$pos_in_ref=$pos_in_ref+$count;	## if S at the beginning - start already marks the position after the soft-clip
													## if S at the end - it doesn't matter, we don't do anything else to these bases
				$pos_in_read=$pos_in_read+$count;	
				#print "# pos in read: $pos_in_read\tpos in ref: $pos_in_ref\n";		
			}
		}
	} else {
				print "Here\n";

		for (my $k=1; $k<=length($seq); $k++){$hash{$k}=$k+$start-1;}
	}
=pod
	foreach my $k (sort {$a <=> $b} keys %hash){
		#print $k,"\t",$hash{$k}{'seq'},"\t",$hash{$k}{'qual'},"\n";
	}
=cut

	### now go through the read and call full codons:
	
	## first find start position of first full codon:
	my $codon_start=0;
	foreach my $s(@codon_pos){
		if ($s<=$start){$codon_start=$s;}
		if ($start < 113){$codon_start=113;}
	}
#	print "first codon: $codon_start\t", $ref_seq{$codon_start},"\n";
	for (my $i = $codon_start; $i < ($codon_start+length($seq)); $i+=3){
		#print $i,"\n";
		### check if $i not beyond the coding sequence:
		if ($i <=2003){
		if (defined $hash{$i} && defined $hash{$i+1} && defined $hash{$i+2}){
			#print $hash{$i}{'qual'},"\t",$hash{$i+1}{'qual'},"\t",$hash{$i+2}{'qual'},"\n";
			if ($hash{$i}{'qual'}>=$baseQ && $hash{$i+1}{'qual'}>=$baseQ && $hash{$i+2}{'qual'}>=$baseQ){
				my $codon =  $hash{$i}{'seq'}.$hash{$i+1}{'seq'}.$hash{$i+2}{'seq'};
#				if ($ref_seq{$i} ne $codon){
					#print "$i\t",$ref_seq{$i},"\t", $codon,"\t",$AA{$codon},"\n";
					if (exists $detected_variants{$i."_".$ref_seq{$i}."_".$codon}){
						$detected_variants{$i."_".$ref_seq{$i}."_".$codon}.=",".$read_id;
					} else {$detected_variants{$i."_".$ref_seq{$i}."_".$codon}=$read_id;}
#				}
			} #else {print "$i\t",$ref_seq{$i},"\t","low_QUAL\n";}
		} #else {print "$i\t",$ref_seq{$i},"\t","not_in_seq\n";}
		} # else {print "outside of HNF1A coding seq\n";}
	}
}}}	

#########  print final results  ###############
#########  compare with the TWIST design ######
#########  and print out the variants2reads.log - this can then be used to see if any variants co-occur more frequently than expected ######

### read the TWIST design into mem
my $design = '/well/mccarthy/production/genetic-screens/resources/Genewiz_results.faulty_library.tab';
my %genewiz;
open(IN,'<',$design) or die;
my $fline=<IN>;
while (my $line = <IN>){
	chomp $line;
	my @line = split("\t",$line);
	my $ref_pos = $line[0]+ 85;

	if ($line[5] eq "NA"){
		$genewiz{$ref_pos."_".$line[1]."_".$line[2]}{'design'} = 0;
	} else {$genewiz{$ref_pos."_".$line[1]."_".$line[2]}{'design'} = 1; }
	$genewiz{$ref_pos."_".$line[1]."_".$line[2]}{'genewiz'} = $line[6];
}
close IN;

__END__

## if for some reason we would like to know which sequencing reads contain the specific codon calls,
##  this would generate a file with this info:

open(OUT,'>',"$base.variants2reads.log") or die;
foreach my $var (sort {$a <=> $b} keys %detected_variants){
	print $var,"\t";
	my @var=split("_",$var);
	## ref pos, AA pos,
	print $var[0],"\t",($var[0]-112+2)/3,"\t",$var[1],"\t",$AA{$var[1]},"\t",$var[2],"\t",$AA{$var[2]},"\t";
	my @reads = split(",",$detected_variants{$var});
	my %seen=();
	my @unique = grep { ! $seen{$_} ++ } @reads;
	print scalar(@unique),"\t";	
	## check if in design:
	if (exists $genewiz{$var}{'design'}){print $genewiz{$var}{'design'},"\t",$genewiz{$var}{'genewiz'},"\n";}
	else {print "0\t0\n";}
	### for the non-ref vars 
	if ($var[1] ne $var[2]){
		print OUT $var,"\t",join(",",@unique),"\n";
	}
}
close OUT;

