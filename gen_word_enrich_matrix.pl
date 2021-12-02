#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use Carp;
use Data::Dumper;
use Quantas;

my $prog = basename ($0);
my $verbose = 0;
my $log2 = 0;
my $count = 0;

my $base = "";
my $suffix = "";
#my $method = "mean";  #sum

GetOptions (
	"base:s"=>\$base,
	"suffix:s"=>\$suffix,
	"v|verbose"=>\$verbose
);

if (@ARGV != 2)
{
	print "generate word enrichment matrix\n";
	print "Usage $prog [options] <in.conf> <out.txt>\n";
	print " -base         [string] : base dir of input data\n";
	print " -suffix       [string] : add suffix to the file names\n";
	print " -v                     : verbose\n";
	exit (1);
}

my ($configFile, $outFile) = @ARGV;

print "loading configuration file from $configFile ...\n" if $verbose;
Carp::croak "contig file $configFile does not exist\n" unless -f $configFile;
my $groups = readExprConfigFile ($configFile, $base, $suffix);

print "done.\n" if $verbose;

print "loading data of individual samples ...\n" if $verbose;

my @sampleData;

my @groupNames = sort {$groups->{$a}->{"id"} <=> $groups->{$b}->{"id"}} keys %$groups;

my $iter = 0;
foreach my $gName (@groupNames)
{
	my $samples = $groups->{$gName}->{"samples"};
	my $s = $samples->[0];	#only one sample per group

	print "$iter: group=$gName, sample=$s\n" if $verbose;
	
	my $inputFile = $s;
	$inputFile = "$base/$inputFile" if $base ne '';
	
	my $sdata = readWordEnrichDataFile ($inputFile);

	$sampleData[$iter]{'name'} = $gName;
	$sampleData[$iter]{'data'} = $sdata;	

	$iter++;
}

print "$iter samples loaded.\n" if $verbose;

print "generate output matrix ...\n" if $verbose;


my %wordHash;
foreach my $s (@sampleData)
{
	my $sdata = $s->{'data'};
	map {$wordHash{$_} = 1} keys %$sdata;
}


my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

print $fout join ("\t", "word", @groupNames), "\n";

foreach my $w (sort keys %wordHash)
{
	my @out;
	foreach my $s (@sampleData)
	{
		my $sdata = $s->{'data'};
		my $z = exists $sdata->{$w} ? $sdata->{$w} : 'NA';
		push @out, $z;
	}

	print $fout join ("\t", $w, @out), "\n";
}

close ($fout);


sub readWordEnrichDataFile
{
	my ($inputFile, $verbose) = @_;
	my $fin;
	open ($fin, "<$inputFile") || Carp::croak "cannot open file $inputFile to read\n";
	<$fin>; #header;
	my %data;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my @cols = split (/\t/, $line);
		my $w = $cols[0];
		my $z = $cols[5];
		$data{$w} = $z;	
	}
	close ($fin);	
	return \%data;	
}


