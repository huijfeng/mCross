#!/usr/bin/perl -w


=head1 NAME


=head1 AUTHOR

Chaolin Zhang (cz2294@columbia.edu)
Created on Nov 22, 2017

=cut



use strict;
use Getopt::Long;
use Carp;
use File::Basename;
use Data::Dumper;


use MyConfig;
use Common 1.02;
use Bed;
use Motif;



my $prog = basename ($0);
my $ext = 10;
my $word = "";
my $padding = 0;
my $mismatch = 0;

my $cache = MyConfig::getDefaultCache ($prog);

my $verbose = 0;

GetOptions ("l=i"=>\$ext,
		"w=s"=>\$word,
		"p:i"=>\$padding,
		"m:i"=>\$mismatch,
		"c:s"=>\$cache,
		"v"=>\$verbose);


if (@ARGV != 2)
{
	print "Get crosslink frequency in motif ...\n";
	print "Usage: $prog [options] <seq_file> <out_file>\n";
	print " -l       [int]   : sequence extension around crosslink site ($ext)\n";
	print " -w       [string]: word to search for\n";
	print " -p       [int]   : number of nucleotide to be padded on both sides ($padding)\n";
	print " -m       [int]   : number of mismatches allowed in the core motif ($mismatch)\n";
	print " -c             [string]: cache dir ($cache)\n";
	print " -v                     : verbose\n";
	exit (1);
}


my ($inFastaFile, $outFile) = @ARGV;


print "loading top words ...\n" if $verbose;


my %crosslinkedBase;

system ("mkdir $cache");
 
print "search motifs ...\n" if $verbose;

my $motifLen = length ($word);

my $tmpSiteFile = "$cache/site.bed";


my $cmd = "PatternMatch -c $word -m $mismatch $inFastaFile > $tmpSiteFile";
print "$cmd\n" if $verbose;

my $ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;



my %crosslinkCount;
for (my $i = 0; $i < $motifLen + $padding * 2; $i++)
{
	$crosslinkCount{$i - $padding} = 0;
}

print "get crosslink profile ...\n" if $verbose;
my $fin;
open ($fin, "<$tmpSiteFile") || Carp::croak "cannot open file $tmpSiteFile to read\n";
my $n = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	my $s = lineToBed ($line);
	
	my $start = $s->{'chromStart'};

	my $crosslinkPos = $ext - $start; #crosslink position relative to motif site
	next if $crosslinkPos < -$padding || $crosslinkPos >= $motifLen + $padding;

	$crosslinkCount{$crosslinkPos}++;
	$n++;
}
close ($fin);


print "write to output...\n" if $verbose;

my @pos = sort {$a <=> $b} keys %crosslinkCount;

my $fout;
open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
print $fout join ("\t", "#pos", "base", "count", "freq"), "\n";
foreach my $p (@pos)
{
	my $base = $p < 0 || $p >= $motifLen ? 'N' : substr($word, $p, 1);

	print $fout join ("\t", sprintf("%02d", $p), $base, $crosslinkCount{$p}, $crosslinkCount{$p} / $n), "\n";
}
close ($fout);

system ("rm -rf $cache");

