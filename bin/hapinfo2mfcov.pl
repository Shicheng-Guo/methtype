#!/usr/bin/perl -w
use strict;
use Cwd;
chdir getcwd;

# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2017-04-5

my %sam;
my %hapinfo;
my %loc;
my %hap;
my %cpg;

 die &USAGE if scalar @ARGV<3;
 my $input=shift @ARGV;
 my $bed=shift @ARGV;
 my $output=shift @ARGV;
# my $input="input.txt";
# my $output="output.txt";

open F,$bed;
my $window=50;
my %hash;
my %HapinfoRegion;
my %Hapinfocov;
my %Hapinfomf;
while(<F>){
        chomp;
        next if /^\s+$/;
        my ($chr,$start,$end,$id,$gene,$block)=split/\s+/;
        my $bin1=int($start/$window);
        my $bin2=int($end/$window);
        foreach my $i($bin1..$bin2){
        $hash{$chr}{$i-1}="$chr:$start-$end";
        $hash{$chr}{$i}="$chr:$start-$end";
        $hash{$chr}{$i+1}="$chr:$start-$end";
		# print "$i\n";
	}
}
close F;

open F1,$input || die "cannot open $input\n";
print "\nStart to reading haplotype files....\n";
while(<F1>){
    chomp;
    next if /^\s*$/;
    my ($sam,$file)=split/\t/;
    print "\t$file reading completed....\n";
	open F2,$file || die "cannot open $file\n";
	while(<F2>){
		chomp;
		my ($loc,$haptype,$count,$cpg)=split /\s+/;                      # load and read haplotype files
		my ($chr,undef,undef)=split/[:-]/,$loc;
		my @tmp=split/,/,$cpg;
		my $bin1=int($tmp[0]/$window);
		my $bin2=int($tmp[$#tmp]/$window);
		my $inner=0;
		while(! $inner and $bin1<=$bin2){
			$inner=1 if defined $hash{$chr}{$bin1};
			$bin1++;
		}
    if($inner){ 	                                                    #  filter out haplotype within interest bed regions
		$hapinfo{$loc}{$cpg}{$haptype}{$sam}+=$count;
		$loc{$loc}=$loc;
		$cpg{$cpg}=$cpg;
		$hap{$haptype}=$haptype;
		$sam{$sam}=$sam;
		# print "$sam\t$loc\t$cpg\t$hapinfo{$loc}{$cpg}{$haptype}{$sam}\n";
		}

	}
	close F2;
}
close F1;

# print  
open OUT1,">$output.region.mf";
open OUT2,">$output.region.cov";

my @sam=sort keys %sam;
my $header=join("\t",@sam);
print OUT1 "Chr\tStart\tEnd\t$header\n";
print OUT2 "Chr\tStart\tEnd\t$header\n";

print "Start to Performance ...\n";
my @loc=keys %loc;
my $loci=0;
my $part=(int($#loc/10+1));
# print "$#loc\tpart=$part\n";

foreach my $loc(sort keys %hapinfo){
	## Step 1. performance haplotype split and recombination
    $loci++;	
	my ($chr,undef,undef)=split/[:-]/,$loc;
	my %Hapinfo;
	my %Hapinfocov;
    foreach my $cpg(sort keys %{$hapinfo{$loc}}){
	foreach my $hap(sort keys %{$hapinfo{$loc}{$cpg}}){
		my @cpg=split/,/,$cpg;
		my $i=0;		
		while($i<=($#cpg)){
			my $string_hap=substr($hap,$i,1);
			my $string_cpg=$cpg[$i];
			my $bin=int($cpg[$i]/$window);
			# print "$cpg[$i]\t$bin\n";
			foreach my $sam(sort keys %sam){
			if( defined $hash{$chr}{$bin} and defined $hapinfo{$loc}{$cpg}{$hap}{$sam}){			
			$Hapinfocov{$hash{$chr}{$bin}}{$string_cpg}{$string_hap}{$sam}+=$hapinfo{$loc}{$cpg}{$hap}{$sam};
			}
			}
			$i++;
			}
		#print "\n";
		}
	}
	
	## Step 2. summary and report homozygote haplotype freqeunce.
	# print "test step2\n";
	foreach my $loc(sort keys %Hapinfocov){
        #print "$loc\n";
        my ($chr,$start,$end)=split/:|-/,$loc;
		print OUT1 "$chr\t$start\t$end";
		print OUT2 "$chr\t$start\t$end";
		foreach my $sam(sort keys %sam){
		my $C=0;
		my $T=0;
		my $mf;
		my $coverage;
		foreach my $cpg(sort keys %{$Hapinfocov{$loc}}){
			foreach my $hap(sort keys %{$Hapinfocov{$loc}{$cpg}}){
               # next if length($hap)<1;
				if (defined $Hapinfocov{$loc}{$cpg}{$hap}{$sam}){
					if($hap eq "C"){
						$C=$C+$Hapinfocov{$loc}{$cpg}{$hap}{$sam};
						 print "$sam\t$loc\t$cpg\t$hap\t$Hapinfocov{$loc}{$cpg}{$hap}{$sam}\t";
					}elsif($hap eq "T"){
						$T=$T+$Hapinfocov{$loc}{$cpg}{$hap}{$sam};
						 print "$sam\t$loc\t$cpg\t$hap\t$Hapinfocov{$loc}{$cpg}{$hap}{$sam}\t";
					}
				}
				 print "\n";
			}
			
		}
			print "$C\t$T\n";
			$coverage= ($T+$C)>0 ? $C+$T : 0;
			$mf= ($T+$C)>0 ? sprintf("%.3f",$C/($T+$C+0.001)) : "NA";
		#	print  "\t$mf";
		#	print  "\t$coverage";
			print OUT1 "\t$mf";
			print OUT2 "\t$coverage";
		}
		#print  "\n";
		#print  "\n";
		print OUT1 "\n";
		print OUT2 "\n";
	}
	
	# print "$loci\n";
	my $status=sprintf("\t........%2d%%(%d/%d) completed!",100*$loci/($#loc+1),$loci,$#loc+1);
	print "$status\n" if $loci % $part eq 0;
}

print "\n";	
close OUT1;
close OUT2;


sub USAGE{
print "Usage: perl $0 samplesheet.txt interest.bed prefix\n";
print "Methylation Haplotype and Counts in Each File will be Collected with Mode of Matrix\n";

print '
Format For: samplesheet file: 
Indx01  Indx01.sortc.bam.hapInfo.txt
Indx02  Indx02.sortc.bam.hapInfo.txt
';
print '
Format For: Interest genomic interval file: 
chr10   76532564        76532591        Col18a1 Exon
chr10   76532243        76532564        Col18a1 Intron
';
}





