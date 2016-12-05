#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;

my $hapinfo=shift @ARGV;
open F,$hapinfo;

my %mf;
my %key;

while(<F>){
my ($id,$haptype,$num,$pos)=split /\s+/;
my ($chr,undef,undef)=split/[:-]/,$id;
my @haptype=split//,$haptype;
my @pos=split/,/,$pos;
# Accumulate
foreach my $i(0..$#pos){
my $key="$chr:$pos[$i]";
$key{$key}=$key;
$mf{$key}{$haptype[$i]} += $num;
}
}

foreach my $key(sort keys %key){
        my($chr,$start)=split /:/,$key;
        my $end=$start+1;
        $mf{$key}{"T"}=0 if ! defined $mf{$key}{"T"};
        $mf{$key}{"C"}=0 if ! defined $mf{$key}{"C"};
        my $nme=$mf{$key}{"T"};
        my $me=$mf{$key}{"C"};
        my $sum=$me+$nme;
        my $ml=sprintf("%.3f",$me/($sum+1));
        print "$chr\t$start\t$end\t$ml\t$me\t$sum\n";
}
