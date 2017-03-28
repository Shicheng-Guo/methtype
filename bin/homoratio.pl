#!/usr/bin/perl -w
use strict;
use Cwd;

# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2017-03-23

my %sam;
my %hapinfo;
my %loc;
my %hap;
my %cpg;

chdir "C:\\Users\\shicheng\\Downloads\\test";
# die &USAGE if scalar @ARGV<1;
# my $input=shift @ARGV;
# my $output=shift @ARGV;
my $input="input.txt";
my $output="output.txt";

open F1,$input;
while(<F1>){
	my ($sam,$file)=split/\t/;
	open F2,$file;
	while(<F2>){
	chomp;
    next if /^\s+$/;
	my ($loc,$haptype,$count,$cpg)=split /\s+/;
	my $T_number = () = $haptype =~ /T/gi;	
		$hapinfo{$loc}{$cpg}{$haptype}{$sam}+=$count;
		$loc{$loc}=$loc;
		$cpg{$cpg}=$cpg;
		$hap{$haptype}=$haptype;
		$sam{$sam}=$sam;
	}
	close F2;
}
close F1;

#  
my @sam=sort keys %sam;
my $header=join("\t",@sam);
print "GenomeInterval\tCpGLocList\tHaplotype\t$header\n";
foreach my $loc(sort keys %hapinfo){
	my %Hapinfo;   # reformed hapinfo(hapsplit and hapmerge)
	foreach my $cpg(sort keys %{$hapinfo{$loc}}){
	foreach my $hap(sort keys %{$hapinfo{$loc}{$cpg}}){
		my @cpg=split/,/,$cpg;
		my $i=0;
		while($i<($#cpg)){
			my $j=2;
			while($j<=$#cpg-$i+1){
			foreach my $sam(sort keys %sam){
			my $string_hap=substr($hap,$i,$j);
			# print "$string_hap\t";
			my $string_cpg=join(",",@cpg[$i..($i+$j-1)]);
			$hapinfo{$loc}{$string_cpg}{$string_hap}{$sam}=$hapinfo{$loc}{$cpg}{$hap}{$sam};
			if($hapinfo{$loc}{$cpg}{$hap}{$sam}){
			$Hapinfo{$loc}{$string_cpg}{$string_hap}{$sam}+=$hapinfo{$loc}{$string_cpg}{$string_hap}{$sam};
			}
			}
			$j++;
			}
			$i++;
		}
		#print "\n";
	}
}

#### calculate Cn/(Cn+Tn+1) for each continous CpG sites(n=2,3,4,...) 
	open OUT,">$output";
	my %homoratio;
	foreach my $loc(sort keys %Hapinfo){
		foreach my $cpg(sort keys %{$Hapinfo{$loc}}){
			print "$loc\t$cpg";
			foreach my $sam(sort keys %sam){
				my $C=0;
				my $T=0;
				my $N=0;
				foreach my $hap(sort keys %{$Hapinfo{$loc}{$cpg}}){
				$Hapinfo{$loc}{$cpg}{$hap}{$sam}=0 if (! defined $Hapinfo{$loc}{$cpg}{$hap}{$sam});
				my $num=length($hap);
    			if($hap eq "C"x $num){
    				$C=$C+$Hapinfo{$loc}{$cpg}{$hap}{$sam};
    			}elsif($hap eq "T"x $num){
    				$T=$T+$Hapinfo{$loc}{$cpg}{$hap}{$sam};
    			}else{
    				$N=$N+$Hapinfo{$loc}{$cpg}{$hap}{$sam};
    			}
				}
				my $ratio1=sprintf("%.2f",$C/($T+$C));
				my $ratio2=sprintf("%.2f",(1-$N/($T+$C+$N)));
				print "\t$ratio1\t$ratio2";
			}
			print "\n";
		}
	}
}	
close OUT;



sub USAGE{
print "Usage: perl $0 M-Primary-Hapinfo M-Plasma-Hapinfo ExcludeDataBase-Hapinfo OutputFileName\n";
print "Status Screening for (Shared HMH:High Methylation Haplotype by Tissue and Plam) in Other Tissues\n";
print "Move All Hapinfo Files into One Folder and assign Matched Tissue and Plasma\n";

print '
Format For: M-Primary-Hapinfo and M-Plasma-Hapinfo: 
chr1:100231312-100231328	CCT		1	100231316,100231324,100231328
chr1:100231312-100231328	CC		1	100231324,100231328
chr1:100231312-100231328	CCCC	4	100231312,100231316,100231324,100231328
';
print '
Format For: ExcludeDataBase-Hapinfo: Hapinfo File List to Excluded: 
STL001BL-01.hapInfo.txt
STL001FT-01.hapInfo.txt
STL001GA-01.hapInfo.txt
';
}

