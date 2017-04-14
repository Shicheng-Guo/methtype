#!/usr/bin/perl -w
# Hapinfo to methylation haplotype load (MHL)
# Run the script to the Hapinfo directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2017-04-13

use strict;
use Cwd;
die &USAGE if @ARGV <1;
my %mch_load_matrix;
my %probe_HMH_samples;
my %hap_count_matrix;
my $hapinfList=shift @ARGV;
open FF,$hapinfList;
chomp(my @hapInfo_files=<FF>);
close FF;

my @sample_list;
foreach my $hapInfo_file(sort @hapInfo_files){
        my @line=split /\//,$hapInfo_file;
	my $sample_name = $line[$#line];
	$sample_name =~ s/.hapInfo.txt//;
	push(@sample_list, $sample_name);
	open(INFILE, "$hapInfo_file") || die("Error in opening $hapInfo_file!");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\s+/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<1);
		$hap_count_matrix{$probeID}->{$sample_name}->{$hapString}=$fields[2];
	}
	close(INFILE);
}

my @unmethylated_haps= ("T"x1,"T"x2,"T"x3,"T"x4,"T"x5,"T"x6,"T"x7,"T"x8,"T"x9,"T"x10,"T"x11,"T"x12,"T"x13,"T"x14);
my @methylated_haps  = ("C"x1,"C"x2,"C"x3,"C"x4,"C"x5,"C"x6,"C"x7,"C"x8,"C"x9,"C"x10,"C"x11,"C"x12,"C"x13,"C"x14);

foreach my $probeID (keys(%hap_count_matrix)){
	foreach my $sample_name (keys(%{$hap_count_matrix{$probeID}})){
		my %k_mer_counts;
		my $mc_hap_load=0;
		foreach my $hapString (keys(%{$hap_count_matrix{$probeID}->{$sample_name}})){
			for(my $word_size = 1; $word_size<=length($hapString); $word_size++){
				next if($word_size>9);
				for(my $i=0; $i<=length($hapString)-$word_size; $i++){
					my $sub_hapString = substr($hapString,$i,$word_size);
					#next if($sub_hapString =~ /[NAG]/i);
					$k_mer_counts{$word_size}->{$sub_hapString}+=$hap_count_matrix{$probeID}->{$sample_name}->{$hapString};					
				}
			}
		}
		my $norm_factor=0;
		foreach my $word_size (sort keys(%k_mer_counts)){
			$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]});
			$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]});
			my $total_count=0;
			foreach my $allele (sort keys(%{$k_mer_counts{$word_size}})){
				$total_count+=$k_mer_counts{$word_size}->{$allele};
			}
			next if($total_count<1);
			my $mh_fraction = $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}/$total_count;
			my $weight = $word_size;
            #print "$word_size\t$methylated_haps[$word_size-1]\t$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}\t$total_count\t$weight\n";
			$mc_hap_load += $weight*$mh_fraction;
			$norm_factor+=$weight;
			#print "($weight:$mh_fraction:$mc_hap_load)\t";
		}
			#print "\n";
		    next if(!$norm_factor);
            #print "$mc_hap_load\t$norm_factor\n";
		    $mc_hap_load/=$norm_factor;
		    $mch_load_matrix{$probeID}->{$sample_name}=$mc_hap_load;
	}
}


print "Probe_id\t", join("\t", sort @sample_list), "\n";
foreach my $probeID (sort keys(%mch_load_matrix)){
	print "$probeID";
	foreach my $sample_name(sort @sample_list){
		$mch_load_matrix{$probeID}->{$sample_name}="NA" if(! defined($mch_load_matrix{$probeID}->{$sample_name}));
		print "\t", $mch_load_matrix{$probeID}->{$sample_name};
	}
	print "\n";
}

sub USAGE{
print "\nperl $0 Hapinfo_File_list > Ouput.txt\n";
print "Just use: ls *hapInfo.txt > Hapinfo_File_list to Get Hapinfo_File_list\n";
}
