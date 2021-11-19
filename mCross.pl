#!/usr/bin/perl -w


=head1 NAME

mCross.pl - de novo motif discovery by precisely registering crosslink sites

=head1 AUTHOR

Chaolin Zhang (cz2294@columbia.edu)
Created on May 22, 2016

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
my $progDir = dirname ($0);
my $ext = 10;
my $padding = 0;
my $mismatch = 1;
my $topWordFile = "";
my $bgFastaFile = "";

my $maxN = 20;
my $maxIter = 100;
my $pseudoCount = 1e-6;

my $clusterSeedMotifs = 0;
my $crosslinkModel = 1; #1=simple, 2=nucleotide-specific
my $scoreMethod = "log"; # or 'sqrt'

my $prefix = "RBP";

my $singleOutputFile = 0;

my $cache = MyConfig::getDefaultCache ($prog);

my $verbose = 0;

GetOptions ("l=i"=>\$ext,
		"seed:s"=>\$topWordFile,
		"bg:s"=>\$bgFastaFile,
		"p:i"=>\$padding,
		"m:i"=>\$mismatch,
		"N:i"=>\$maxN,
		"cluster-seeds"=>\$clusterSeedMotifs,
		"xl-model:i"=>\$crosslinkModel,
		"score-method:s"=>\$scoreMethod,
		"prefix:s"=>\$prefix,
		"single-out-file"=>\$singleOutputFile,
		"c:s"=>\$cache,
		"v"=>\$verbose);


if (@ARGV != 2)
{
	print "Motif discovery anchored by crosslink sites ...\n";
	print "Usage: $prog [options] <seq_file> <out_file or out_file_stem>\n";
	print " -l       [int]   : sequence extension around crosslink site ($ext)\n";
	print " --seed   [file]  : top_nmer_file\n";
	print " --bg     [file]  : if top_nmer not provided, fg and bg file are used to get the list\n";
	print " -p       [int]   : pad the seed motif on both sides ($padding)\n";
	print " -m       [int]   : number of mismatches allowed in the core motif ($mismatch)\n";
	print " -N       [int]   : max number of seed words to search ($maxN)\n";
	print " --cluster-seeds  : cluster seed word\n";
	print " --xl-model     [int]   : crosslink model (1=simple(default), 2=nucleotide-specific)\n";
	print " --score-method [string]: [log]|sqrt\n";
	print " --prefix       [string]: prefix of the motif name ($prefix)\n";
	print " --single-output-file   : write all motifs to a single file\n";
	print " -c             [string]: cache dir ($cache)\n";
	print " -v                     : verbose\n";
	exit (1);
}


my ($inFastaFile, $outFile) = @ARGV;


Carp::croak "either top_nmer or bg fasta file has to be provided\n" unless (-f $topWordFile) || (-f $bgFastaFile);


system ("mkdir $cache");


my %alphabet2index = (
	'A'=>0,
	'C'=>1,
	'G'=>2,
	'T'=>3);




my @seedMotifs;

print "loading top words ...\n" if $verbose;


if (-f $topWordFile)
{
	my $fin; 
	open ($fin, "<$topWordFile") || Carp::croak "cannot open file $topWordFile to read\n";
	my $iter = 0;
	while (my $line = <$fin>)
	{
		if ($maxN > 0)
		{
			last if $iter >= $maxN;
		}
	
		chomp $line;

		$line =~s/\"//g;

		next if $line=~/^\s*$/;
		my @cols = split (/\s/, $line);

		next if @cols != 2;
		my $w = $cols[0];
		$w = "N"x$padding . $w . "N"x$padding if $padding > 0;

		push @seedMotifs, $w;
		$iter++;
	}
	close ($fin);
}
else
{
	my $fg = "$cache/fg.nomask.fa";
	my $bg = "$cache/bg.nomask.fa";
	my $tmpOutFile = "$cache/w7.txt";

	my $cmd = "perl $progDir/hardMask.pl --no-N --no-desc $inFastaFile $fg";
	print "$cmd\n" if $verbose;
	my $ret = system ($cmd);
	Carp::croak "cmd=$cmd failed:$?\n" if $ret != 0;
	
	$cmd = "perl $progDir/hardMask.pl --no-N --no-desc $bgFastaFile $bg";
	print "$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed:$?\n" if $ret != 0;

	$cmd = "perl $progDir/word_enrich.pl -w 7 -test binom -v $fg $bg $tmpOutFile";
	print "$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "cmd=$cmd failed:$?\n" if $ret != 0;

	my $fin;

	open ($fin, "<$tmpOutFile") || Carp::croak "cannot open file $tmpOutFile to read\n";
	<$fin>; #headline
	my %wordHash;
	while (my $line = <$fin>)
	{

		chomp $line;
		next if $line =~/^\s*$/;

		my @cols = split ("\t", $line);
		$wordHash{$cols[0]} = {log2fc=> $cols[4], z=>$cols[5], p=>$cols[6]}; #zscore
	}
	close ($fin);
	#Carp::croak Dumper (\%wordHash), "\n";	

	my $iter = 0;
	foreach my $w (sort {$wordHash{$b}->{'z'} <=> $wordHash{$a}->{'z'}} keys %wordHash)
	{
		if ($maxN > 0)
		{
			last if $iter >= $maxN;
		}

		print join("\t", $w, $wordHash{$w}->{'log2fc'},
			$wordHash{$w}->{'z'},
			$wordHash{$w}->{'p'}), "\n" if $verbose;

		$w = "N"x$padding . $w . "N"x$padding if $padding > 0;
		push @seedMotifs, $w;
		$iter++;	
	}
}


my $n = @seedMotifs;
print "$n top words loaded\n" if $verbose;

exit (0) if $n < 1;


my @seedMotifClusters;
if ($clusterSeedMotifs)
{
	print "clustering top words ...\n" if $verbose;
	my $clusters = clusterWords (\@seedMotifs, 2);

	my $iter = 0;

	foreach my $c (@$clusters)
	{
		my @words = map {$seedMotifs[$_]} @$c;
		print "$iter: ", join (", ", @words), "\n" if $verbose;
	
		push @seedMotifClusters, \@words;
		$iter++;
	}
}
else
{
	map {push @seedMotifClusters, [$_]} @seedMotifs;
}


my %crosslinkedBase;

if ($crosslinkModel == 2)
{
	my $seqIO = Bio::SeqIO->new (-file => $inFastaFile, -format=> 'Fasta');
	while (my $seq = $seqIO->next_seq())
	{
		my $id = $seq->id;
		my $base = substr ($seq->seq, $ext, 1);
		$crosslinkedBase{$id} = $base;
	}
}


 
print "refine motifs ...\n" if $verbose;

my $fout;

my $outFile2 = $singleOutputFile ? $outFile : "$outFile.00.mat";

open ($fout, ">$outFile2") || Carp::croak "cannot open file $outFile2 to write\n";
for (my $i = 0; $i < @seedMotifClusters; $i++)
{
	my $consensus = $seedMotifClusters[$i];

	print "$i : ", join (",", @$consensus), "...\n" if $verbose;

	my $motif = refineMotif ($consensus, $crosslinkModel, $mismatch, $inFastaFile, $ext, $maxIter, $cache);
	my $motifName = "$prefix.$i";
	$motif->{'name'} = $motifName;
	
	#output the motif
	my $old = select($fout);

	#print "\n##MOTIF $i\n";

	printMCrossMotif ($motif);
	select ($old);	

	if ($singleOutputFile ==0 && $i < @seedMotifClusters - 1)
	{
		close ($fout);
		$outFile2 = $outFile . "." . sprintf ("%02d", $i+1) . ".mat";
		open ($fout, ">$outFile2") || Carp::croak "cannot open file $outFile2 to write\n";
	}

	#exit (1);
}
close ($fout);

system ("rm -rf $cache");


sub refineMotif
{
	my ($consensus, $crosslinkModel, $mismatch, $seqFile, $ext, $maxIter, $cache) = @_;

	my $motifLen = length ($consensus->[0]);
	

	print "search sites ...\n" if $verbose;

	my $tmpSiteFile = "$cache/site.bed";
	my $consensusFile = "$cache/consensus.txt";

	my $fout;

	#we limit the degeneracy here
	#need to be fixed later
	
	if (@$consensus > 16)
	{
		$mismatch = 0;
	}
	elsif (@$consensus > 4)
	{
		$mismatch = min($mismatch, 1);
	}
	print "mismatch allowed=$mismatch\n" if $verbose;	

	open ($fout, ">$consensusFile") || Carp::croak "cannot open file $consensusFile to write\n";
	foreach my $c (@$consensus)
	{
		print $fout join ("\t", $c, $mismatch), "\n";
	}
	close ($fout);	

	my $cmd = "PatternMatch -l $consensusFile $seqFile > $tmpSiteFile";
	print "$cmd\n" if $verbose;

	my $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed:$?\n" if $ret != 0;
		
	my $fin;
	open ($fin, "<$tmpSiteFile") || Carp::croak "cannot open file $tmpSiteFile to read\n";
	my @sites;
	my %siteHash;
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my $s = lineToBed ($line);
		$s->{'inc'} = 1;
		my $siteKey = $s->{'chrom'} . "." . $s->{'chromStart'};
		push @sites, $s unless exists $siteHash{$siteKey};	#make sure we do not add the same site multiple times
		$siteHash{$siteKey} = 1;
	}
	close ($fin);
	
	my $n = @sites;
	print "$n sites matched\n" if $verbose;

	##
	print "build motif model ...\n" if $verbose;

	my $motif = buildMotifModel (\@sites, $crosslinkModel, $ext);
	$motif->{'consensus'} = $consensus;

	printSimpleMotif ($motif);

	##
	print "refine motif model ...\n" if $verbose;
	my $iter;
	for ($iter = 0; $iter < $maxIter; $iter++)
	{

		#my $s = scoreMotif ($motif);
		#print "score=$s\n";
	
		#$motif->{'score'} = $s;

		my $ntoggle = 0;
		foreach my $s (@sites)
		{
			$ntoggle += toggleSite ($motif, $s, $crosslinkModel, $scoreMethod);
		}

		print "$ntoggle sites changed in the motif model\n" if $verbose;

		last if $ntoggle < 1;

		#rebuild the motif
		$motif = buildMotifModel (\@sites, $crosslinkModel, $ext);

		$motif->{'consensus'} = $consensus;

		print "\niter=$iter\n";
		printSimpleMotif ($motif);
	}

	print "motif not converged\n" if $iter == $maxIter &&$verbose;


	#printMotif ($motif);
	#Carp::croak Dumper ($model), "\n";

	return $motif;
}


#if inclusion status of a a site need to be changed
#return 1 if yes, 0 if no

sub toggleSite
{
	my ($motif, $site, $crosslinkModel, $scoreMethod) = @_;
	my $removeSite = $site->{'inc'} == 1 ? 1 : 0;

	my $newMotif = copyMotif ($motif);

	my $seqId = $site->{'chrom'};
	my $start = $site->{'chromStart'};
	my @cols = split (":", $site->{'name'});
	my $siteSeq = uc($cols[1]);

	#Carp::croak "siteSeq=$siteSeq\n";

	return 0 if $siteSeq =~/[^ACGT]/g;
	
	my @b = split(//, $siteSeq);
	#Carp::croak Dumper(\@b), "\n";

	my $motifCountMatrix = $newMotif->{'count'};
	my $crosslinkCount = $newMotif->{'croslink'};
	my $motifLen = @b;

	if ($removeSite)
	{
		map {$motifCountMatrix->[$_][$alphabet2index{$b[$_]}]--} (0..($motifLen-1));
	
		my $crosslinkPos = $ext - $start; #crosslink position relative to motif site
		if ($crosslinkModel == 1)
		{
			$crosslinkCount->{$crosslinkPos}--;
		}
		elsif ($crosslinkModel == 2)
		{
			my $base = $crosslinkedBase{$seqId};
			$crosslinkCount->{$crosslinkPos}[$alphabet2index{$base}]--;		
		}
		$newMotif->{'n'}--;
	}
	else
	{
		map {$motifCountMatrix->[$_][$alphabet2index{$b[$_]}]++} (0..($motifLen-1));
	
		my $crosslinkPos = $ext - $start; #crosslink position relative to motif site
		if ($crosslinkModel == 1)
        {
            $crosslinkCount->{$crosslinkPos}++;
        }
        elsif ($crosslinkModel == 2)
        {
            my $base = $crosslinkedBase{$seqId};
            $crosslinkCount->{$crosslinkPos}[$alphabet2index{$base}]++;
        }
		$newMotif->{'n'}++;
	}
	
	my $score = $motif->{'score'};
	my $newScore = scoreMotif ($newMotif, $crosslinkModel, $scoreMethod);
	
	my $toggle = 0;
	if ($score < $newScore)
	{
		$toggle = 1;
		#my $n = $motif->{'n'};
		#my $n2 = $newMotif->{'n'};
		#print "seq=$siteSeq, score = $score ($n sites), newScore=$newScore ($n2 sites)\n";	
		$site->{'inc'} = 1 - $site->{'inc'};
	}
	return $toggle;
}



sub copyMotif
{
	my $motif = $_[0];
	my %motifNew = %$motif;
	
	#copy crosslink data
	my %crosslinkCount;

	if ($crosslinkModel == 1)
	{
	 	%crosslinkCount = %{$motif->{'crosslink'}};
	}
	elsif ($crosslinkModel == 2)
	{
		foreach my $p (keys %{$motif->{'crosslink'}})
		{
			my @c = @{$motif->{'crosslink'}{$p}};
			$crosslinkCount{$p} = \@c;
		}
	}
	$motifNew{'crosslink'} = \%crosslinkCount;

	#copy count matrix
	my @motifCountMatrix;
	for (my $i = 0; $i < @{$motif->{'count'}}; $i++)
	{
		my @c = @{$motif->{'count'}[$i]};
		$motifCountMatrix[$i] = \@c;
	}
	$motifNew{'count'} = \@motifCountMatrix;
	return \%motifNew;
}


sub scoreMotif
{
	my ($motif, $crosslinkModel, $scoreMethod) = @_;
	my $n = $motif->{'n'};
	
	my $motifCountMatrix=$motif->{'count'};
    my $crosslinkCount=$motif->{'crosslink'};
	
	my $s1 = 0;
	for (my $i = 0; $i < @$motifCountMatrix; $i++)
	{
		my $prob = count2prob ($motifCountMatrix->[$i], $pseudoCount);
		map {$s1+= $_ * log($_/0.25)} @$prob;
	}

	my $s2 = 0;
	
	if ($crosslinkModel == 1)
	{
		my $V = keys %$crosslinkCount;
		my @count = values %$crosslinkCount;
		my $prob = count2prob (\@count, $pseudoCount);
	
		map {$s2+= $_ * log($_ * $V)} @$prob;
	}
	elsif ($crosslinkModel == 2)
	{
		my $V = keys %$crosslinkCount;
		my $totalCount = 0;
		
		foreach my $i (keys %$crosslinkCount)
		{
			for (my $j = 0; $i < 3; $j++)
			{
				$totalCount += $crosslinkCount->{$i}[$j];
			}
		}
		
		Carp::croak "total site count = 0\n" if $totalCount == 0;

		foreach my $i (keys %$crosslinkCount)
        {
            for (my $j = 0; $i < 3; $j++)
            {
                my $p = $crosslinkCount->{$i}[$j]/$totalCount;
				$p = $pseudoCount if $p < $pseudoCount;
				$s2 += $p * log ($p * $V * 4);
            }
        }
	}

	#print "s1=$s1, s2=$s2\n";
	#return $s1 * $n;

	my $scaling = 0;
	$scaling = $scoreMethod eq 'log' ? log($n) : sqrt ($n) if $n > 0;

	return ($s1+$s2) * $scaling;
}


#print motif in transfac format
sub printMCrossMotif
{
	my $motif = $_[0];
	
	my $motifCountMatrix=$motif->{'count'};
	my $crosslinkCount=$motif->{'crosslink'};

	print join ("\t", "AC", $motif->{'name'}), "\n";
	print "XX\n";
	print join ("\t", "TY", "Motif"), "\n";
	print "XX\n";
	print join ("\t", "ID", $motif->{'name'}), "\n";
	print "XX\n";
	print join ("\t", "NA", $motif->{'name'}), "\n";
	print "XX\n";
	
	print "DE", "\tN=", $motif->{'n'}, ", Consensus=", join(",",@{$motif->{'consensus'}}), ", Score=", $motif->{'score'}, "\n";

	my @consensus = split (//, $motif->{'consensus'}[0]);
	#this could be multiple seed words
	#need to get the real consensus later

	print join ("\t", "P0", qw(A C G T)), "\n";

	my @infoContent;
	for (my $i = 0; $i < @$motifCountMatrix; $i++)
	{
		my $n = sum ($motifCountMatrix->[$i]);
		my $info = 2 - entropy ($motifCountMatrix->[$i]) - 1.5 / $n / log(2);
		$info = 0 if $info < 0;

		#https://en.wikipedia.org/wiki/Sequence_logo
		push @infoContent, $info;
		print join ("\t", sprintf("%02d", $i+1), @{$motifCountMatrix->[$i]}, $consensus[$i]), "\n";
	}

	print "XX\n";
	print "IC\t", join (" ", map {sprintf("%.3f", $_)} @infoContent), "\n";
	print "XX\n";

	#print "#Crosslink=\n";
	print join ("\t", "XL"), "\n"; 
	my @pos = sort {$a <=> $b} keys %$crosslinkCount;

    if ($crosslinkModel == 1)
    {
        foreach my $p (@pos)
        {
			print join ("\t", sprintf("%02d", $p), $crosslinkCount->{$p}), "\n";
        }
    }
    elsif ($crosslinkModel == 2)
    {
        print join ("\t", "POS", qw(A C G T)), "\n";
        foreach my $p (@pos)
        {
            print join ("\t", sprintf("%02d", $p), @{$crosslinkCount->{$p}}), "\n";
        }
    }

	print "XX\n";
	print "//\n";
}

sub trimMCrossMotif
{
	#not finished yet

	my ($motif, $threshold) = @_;
	my $motifCountMatrix=$motif->{'count'};

	my ($left, $right) = (0, 0);
	for (my $i = 0; $i < @$motifCountMatrix; $i++)
    {
		
    }	
}

sub printSimpleMotif
{
	my $motif = $_[0];
	
	my $motifCountMatrix=$motif->{'count'};
	my $crosslinkCount=$motif->{'crosslink'};

	print "N=", $motif->{'n'}, ", Consensus=", join (",", @{$motif->{'consensus'}}), ", Score=", $motif->{'score'}, "\n";
	print "#CountMatrix=\n";
	print join ("\t", "POS", qw(A C G T)), "\n";
	for (my $i = 0; $i < @$motifCountMatrix; $i++)
	{
		print join ("\t", $i, @{$motifCountMatrix->[$i]}), "\n";
	}

	print "#Crosslink=\n";

	my @pos = sort {$a <=> $b} keys %$crosslinkCount;
	
	if ($crosslinkModel == 1)
	{
		foreach my $p (@pos)
		{
			print join ("\t", $p, $crosslinkCount->{$p}), "\n";
		}
	}
	elsif ($crosslinkModel == 2)
	{
		print join ("\t", "POS", qw(A C G T)), "\n";
    	foreach my $p (@pos)
    	{
        	print join ("\t", $p, @{$crosslinkCount->{$p}}), "\n";
   		}
	}
}



sub buildMotifModel
{
	my ($sites, $crosslinkModel, $ext) = @_;
	my $motifLen = $sites->[0]->{'chromEnd'} - $sites->[0]->{'chromStart'} + 1;
	
	my @motifCountMatrix;
	my %crosslinkCount;
	
	#init
	for (my $i = 0; $i < $motifLen; $i++)
	{
		for (my $j = 0; $j < 4; $j++)
		{
			$motifCountMatrix[$i][$j] = 0;
		}
	}

	for (my $i = 0; $i < $ext * 2 + 1 - $motifLen + 1; $i++)
	{
		if ($crosslinkModel == 1)
		{
			$crosslinkCount{$ext-$i} = 0;
		}
		elsif ($crosslinkModel == 2)
		{
			for (my $j = 0; $j < 4; $j++)
			{
				$crosslinkCount{$ext-$i}[$j] = 0;
			}
		}
		else
		{
			Carp::croak "crosslink model unknown\n";
		}
	}

	my $n = 0;

	foreach my $s (@$sites)
	{
		my $seqId = $s->{'chrom'};
		my $start = $s->{'chromStart'};
		my @cols = split (":", $s->{'name'});
		my $siteSeq = uc($cols[1]);

		#Carp::croak "siteSeq=$siteSeq\n";

		next if $siteSeq =~/[^ACGT]/g || $s->{'inc'} != 1;		

		my @b = split(//, $siteSeq);
		#Carp::croak Dumper(\@b), "\n";

		#print "siteSeq = $siteSeq\n";
		map {$motifCountMatrix[$_][$alphabet2index{$b[$_]}]++} (0..($motifLen-1));
		
		my $crosslinkPos = $ext - $start; #crosslink position relative to motif site
		if ($crosslinkModel == 1)
		{
			$crosslinkCount{$crosslinkPos}++;
		}
		elsif ($crosslinkModel == 2)
		{
			my $base = $crosslinkedBase{$seqId};
			$crosslinkCount{$crosslinkPos}[$alphabet2index{$base}]++;
		}
		$n++;
	}

	my %motif = (n=> $n, count=>\@motifCountMatrix, crosslink=>\%crosslinkCount);
	my $score = scoreMotif (\%motif, $crosslinkModel, $scoreMethod);
	$motif{'score'} = $score;
	return \%motif;
}


=head2
	my $clusters = cluserWords (\@words);
	
=cut
sub clusterWords
{
	my ($words, $maxDiff) = @_;
	my @adjMatrix;  #adjacency matrix, 1 if dist < $maxDiff
	my $ignoreCase = 1;
	
	for (my $i = 0; $i < @$words; $i++)
	{
		for (my $j = $i+1; $j < @$words; $j++)
		{
			my $dist = countMismatch ($words->[$i], $words->[$j], $ignoreCase, $maxDiff);
			$adjMatrix[$i][$j] = $dist < $maxDiff ? 1 : 0;
			$adjMatrix[$j][$i] = $adjMatrix[$i][$j];
		}
	}
	my $clusters = matrix2clusters (\@adjMatrix);

	return $clusters;
}



