#!/usr/bin/perl -w


use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

use Motif;

my $prog = basename ($0);

if (@ARGV != 1)
{
	print "$prog <in.mat>\n" if @ARGV != 1;
	exit (1);
}

my ($in) = @ARGV;

my $matrices = readMotifFile ($in);

print join ("\t", 'name', 'consensus', 'information_content', 'crosslink_counts', 'num_sites', 'score'), "\n";
foreach my $m (@$matrices)
{
        #print Dumper ($m), "\n";
    my $name = $m->{'AC'};
	my $desc = $m->{'DE'};
	$desc =~/N=(\d+)\, Consensus=(.*?)\, Score=(.*?)$/;
	my $nSites = $1;
	my $consensus = $2;
	my $score = $3;

	my $ic = $m->{'IC'};
	my $xl = $m->{'XL'};

	print join ("\t", $name, $consensus, join(";", @$ic), join (";", @$xl), $nSites, $score), "\n";
	#Carp::croak Dumper ($m), "\n";
}


