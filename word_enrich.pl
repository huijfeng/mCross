#!/usr/bin/perl -w

use strict;

#use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use Carp;

use Data::Dumper;

use Common;
use Sequence 1.02;
use MyConfig;

my $verbose = 0;
my $prog = basename ($0);

my $wordSize = 6;
my $inclusive = 0;
my $testMethod = "binom"; #"hyperg"

my $noMask = 0;

my $cache = getDefaultCache ($prog);



GetOptions ('w|word-size:i'=>\$wordSize,
		'inc'=>\$inclusive,
		'test:s'=>\$testMethod,
		'no-mask'=>\$noMask,
		'v|verbose'=>\$verbose);


if (@ARGV != 3)
{
	print "test word enrichment in fg vs. bg sequences\n";
	print "Usage: $prog [options] <fg.fa> <bg.fa> <out.txt>\n";
	print " -w    [int]   : word size ($wordSize)\n";
	print " -inc          : inclusive (fg is included in bg)\n";
	print " -test [string]: test method ([binom]|hyperg)\n"; 
	print " --no-mask     : do not exclude masked words for testing\n";
	print " -c    [dir]   : cache dir ($cache)\n";
	print " -v            : verose\n";
	exit (1);
}


my ($fgFastaFile, $bgFastaFile, $outFile) = @ARGV;

#system ("mkdir $cache");

print "count words in fg file $fgFastaFile ...\n" if $verbose;
my $fgWordHash = wordcount ($fgFastaFile, $wordSize, {fasta=>1, noMask=>$noMask});

print "count words in bg file $bgFastaFile ...\n" if $verbose;
my $bgWordHash = wordcount ($bgFastaFile, $wordSize, {fasta=>1, noMask=>$noMask});

print "combine word lists from fg and bg ...\n" if $verbose;
my %wordHash;

map {$wordHash{$_}->{'fg'} = $fgWordHash->{$_}} keys %$fgWordHash;
map {$wordHash{$_}->{'bg'} = $bgWordHash->{$_}} keys %$bgWordHash;

map {
	$wordHash{$_}->{'fg'} = exists $wordHash{$_}->{'fg'}? $wordHash{$_}->{'fg'} : 0; 
	$wordHash{$_}->{'bg'} = exists $wordHash{$_}->{'bg'}? $wordHash{$_}->{'bg'} : 0;
} keys %wordHash;

#print Dumper (\%wordHash), "\n";

my ($fgSum, $bgSum) = 0;
map {$fgSum += $wordHash{$_}->{'fg'};
	 $bgSum += $wordHash{$_}->{'bg'};
} keys %wordHash;

$bgSum += $fgSum unless $inclusive;

print "number of words: $fgSum (fg), $bgSum (fg+bg)\n" if $verbose;

print "testing word enrichment ...\n" if $verbose;

my $iter = 0;

my @cleanWords;
foreach my $w (sort keys %wordHash)
{
	next if $w=~/[^ACGTU]/g;

	my $fg = $wordHash{$w}->{'fg'};
	my $bg = $wordHash{$w}->{'bg'};
	$bg += $fg unless $inclusive;
	
	#k - fgW
	#m - fgW+bgW
	#n - fgSum
	#N - fgSum+bgSum
	
	my $log2FC = ($fgSum == 0 || $bg-$fg == 0 || $fg == 0 || $bgSum - $fgSum == 0) ? 'NA' : log ($fg * ($bgSum - $fgSum)/ $fgSum / ($bg-$fg)) / log(2);
	
	print "$iter: $w ($fg of $bg in fg, log2fc=$log2FC) ...\n" if $verbose;

	my $pval = 1;
	my $r = $fgSum / $bgSum;
	my $zscore = ($fg - $bg * $r) / sqrt ($bg * $r * (1-$r)); 

	if ($testMethod eq 'hyperg')
	{
		$pval = hypergeoTest ($bg, $fg, $bgSum, $fgSum);
	}
	elsif ($testMethod eq 'binom')
	{
		$pval = binomTest ($fg, $bg, $r);
	}
	$wordHash{$w} = {fg=>$fg, bg=>$bg, log2FC=>$log2FC, zscore=>$zscore, pval=>$pval};
	push @cleanWords, $w;

	#print $fout join ("\t", $w, $fg, $bg, $log2FC, $pval), "\n";
	$iter++;
}

my $fout;

open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
print $fout join ("\t", 'word', 'di_entropy', 'fg', 'fg+bg', 'log2FC', 'zscore', 'pval'), "\n";
foreach my $w (sort {$wordHash{$b}->{'zscore'} <=> $wordHash{$a}->{'zscore'}} @cleanWords)
{
	my $wh = $wordHash{$w};
	print $fout join ("\t", $w, seqEntropy ($w, 2), $wh->{'fg'}, $wh->{'bg'}, $wh->{'log2FC'}, $wh->{'zscore'}, $wh->{'pval'}), "\n";
}
close ($fout);


#this function depends on emboss, so not used anymore
sub enumerateWords
{
    my ($fastaFile, $wordSize, $cache) = @_;
    my $tmpFile = "$cache/word.count." . rand();
    my $cmd = "wordcount -sequence $fastaFile -wordsize $wordSize -outfile $tmpFile >& /dev/null";
    my $ret = system ($cmd);
    Carp::croak "$cmd failed\n" unless $ret == 0;

    my $fin;
    open ($fin, "<$tmpFile") || Carp::croak "cannot open file $tmpFile to read\n";

    my %wordHash;
    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line =~/^\s*$/;
        my ($w, $c) = split (/\s+/, $line);
        $w=uc($w);

		if ($noMask == 0)
		{
        	next if $w=~/[^ACGTU]/g;
		}
        $wordHash{$w} = $c;
    }
    close ($fin);

    unlink $tmpFile;
    return \%wordHash;
}


sub countKnownMotifs
{
	my ($fastaFile, $motifListFile) = @_;

	my %wordHash;
	my $cmd = "PatternMatch -l $motifListFile";

	my $fin;
	open ($fin, "$cmd |") ||Carp::croak "CMD=$cmd failed\n";
	
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		my @cols = split ("\t", $line);
		my $name = $cols[3];
		my ($consensus, $site) = split (":", $name);

		$wordHash{$consensus}++;	
	}
	close ($fin);
	return \%wordHash;
}

