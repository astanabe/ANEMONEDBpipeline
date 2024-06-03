use strict;
use Cwd 'getcwd';
use Time::Piece;
use utf8;
use open ':encoding(utf8)';
use open ':std';

my $filehandleinput1;
my $filehandleoutput1;

my $outputfolder = $ARGV[-1];
my $solutionvoltable = 'solutionvoltable.tsv';
my $watervoltable = 'watervoltable.tsv';
my $stdconctable = 'stdconctable.tsv';
my $negativesamplelist = 'negativesamplelist.txt';
my $replicatelist = 'replicatelist.txt';

my %solutionvol;
my %watervol;
my %stdconc;
my %negativesample;
my %parentsample;
my %daughtersample;
{
	if (-e $stdconctable) {
		my $lineno = 1;
		my @label;
		open($filehandleinput1, "< $stdconctable") or die(__LINE__ . ": Cannot open \"$stdconctable\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			my @row = split(/\t/, $_);
			if ($lineno == 1 && scalar(@row) > 2 && $_ !~ /\t(?:\d+|\d+\.\d+)\t/ && $_ !~ /\t(?:\d+|\d+\.\d+)$/) {
				@label = @row;
			}
			elsif ($lineno > 1 && @label && scalar(@row) >= scalar(@label)) {
				for (my $i = 1; $i < scalar(@label); $i ++) {
					if ($row[$i] =~ /^(?:\d+|\d+\.\d+)$/) {
						$stdconc{$row[0]}{$label[$i]} = $row[$i];
					}
					else {
						&errorMessage(__LINE__, "\"$stdconctable\" is invalid.");
					}
				}
			}
			elsif (!@label && scalar(@row) == 2 && $row[1] =~ /^(?:\d+|\d+\.\d+)$/) {
				$stdconc{$row[0]} = $row[1];
			}
			elsif (!@label && scalar(@row) == 3 && $row[2] =~ /^(?:\d+|\d+\.\d+)$/) {
				$stdconc{$row[0]}{$row[1]} = $row[2];
			}
			elsif (@row) {
				&errorMessage(__LINE__, "\"$stdconctable\" is invalid.");
			}
			$lineno ++;
		}
		close($filehandleinput1);
	}
	if (-e $solutionvoltable) {
		open($filehandleinput1, "< $solutionvoltable") or die(__LINE__ . ": Cannot open \"$solutionvoltable\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			my @row = split(/\t/, $_);
			if (scalar(@row) >= 2) {
				for (my $i = 1; $i < scalar(@row); $i ++) {
					if ($row[$i] =~ /^(?:\d+|\d+\.\d+)$/) {
						$solutionvol{$row[0]} += $row[$i];
					}
				}
			}
		}
		close($filehandleinput1);
	}
	if (-e $watervoltable) {
		open($filehandleinput1, "< $watervoltable") or die(__LINE__ . ": Cannot open \"$watervoltable\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			my @row = split(/\t/, $_);
			if (scalar(@row) >= 2) {
				for (my $i = 1; $i < scalar(@row); $i ++) {
					if ($row[$i] =~ /^(?:\d+|\d+\.\d+)$/) {
						$watervol{$row[0]} += $row[$i];
					}
				}
			}
		}
		close($filehandleinput1);
	}
	if (-e $negativesamplelist) {
		open($filehandleinput1, "< $negativesamplelist") or die(__LINE__ . ": Cannot open \"$negativesamplelist\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			if ($_) {
				$negativesample{$_} = 1;
			}
		}
		close($filehandleinput1);
	}
	if (-e $replicatelist) {
		open($filehandleinput1, "< $replicatelist") or die(__LINE__ . ": Cannot open \"$replicatelist\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			my @temp = split(/\t/, $_);
			for (my $i = 0; $i < scalar(@temp); $i ++) {
				$parentsample{$temp[$i]} = $temp[0];
				push(@{$daughtersample{$temp[0]}}, $temp[$i]);
			}
		}
		close($filehandleinput1);
	}
}
my @samplenames;
{
	if (-e 'community_qc3nn_target.tsv') {
		my %samplenames;
		open($filehandleinput1, "< community_qc3nn_target.tsv") or die(__LINE__ . ": Cannot open \"community_qc3nn_target.tsv\".\n");
		my $lineno = 1;
		while (<$filehandleinput1>) {
			if ($lineno > 1) {
				s/\r?\n?$//;
				my @row = split(/\t/, $_);
				$samplenames{$row[0]} = 1;
			}
			$lineno ++;
		}
		close($filehandleinput1);
		@samplenames = sort({$a cmp $b} keys(%samplenames));
	}
	else {
		die(__LINE__ . "\"community_qc3nn_target.tsv\" does not exist.\n");
	}
}
my $root = getcwd();
my $locus;
my $team;
my $project;
my $run;
{
	my @path = split('/', $root);
	$team = $path[-3];
	$project = $path[-2];
	$run = $path[-1];
}

foreach my $locus (&readSeq('forwardprimer.fasta')) {
	if (!-e $outputfolder) {
		mkdir($outputfolder);
	}
	if (-e $outputfolder && !-d $outputfolder) {
		die(__LINE__ . ": \"$outputfolder\" is not a directory.\n");
	}
	foreach my $samplename (@samplenames) {
		if (!-e "$outputfolder/$locus") {
			if (!mkdir("$outputfolder/$locus")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$locus\".\n");
			}
		}
		if (!-e "$outputfolder/$locus/$team") {
			if (!mkdir("$outputfolder/$locus/$team")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team\".\n");
			}
		}
		if (!-e "$outputfolder/$locus/$team/$project") {
			if (!mkdir("$outputfolder/$locus/$team/$project")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team/$project\".\n");
			}
		}
		if (!-e "$outputfolder/$locus/$team/$project/$project$run") {
			if (!mkdir("$outputfolder/$locus/$team/$project/$project$run")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team/$project/$project$run\".\n");
			}
		}
		if (!-e "$outputfolder/$locus/$team/$project/$project$run/$samplename") {
			if (!mkdir("$outputfolder/$locus/$team/$project/$project$run/$samplename")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team/$project/$project$run/$samplename\".\n");
			}
			foreach my $strand ('.forward', '.reverse', '') {
				if ($daughtersample{$samplename}) {
					foreach my $daughtersample (@{$daughtersample{$samplename}}) {
						if (-e "$locus/demultiplexed/$daughtersample$strand.fastq.xz") {
							system("cp $locus/demultiplexed/$samplename$strand.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$samplename/$samplename$strand.fastq.xz");
						}
					}
				}
				elsif (-e "$locus/demultiplexed/$samplename$strand.fastq.xz") {
					system("cp $locus/demultiplexed/$samplename$strand.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$samplename/$samplename$strand.fastq.xz");
				}
			}
			if (-e "$locus/$samplename/phylum.png") {
				system("cp $locus/$samplename/phylum.png $outputfolder/$locus/$team/$project/$project$run/$samplename/wordcloud.phylum.png");
			}
			if (-e "$locus/$samplename/class.png") {
				system("cp $locus/$samplename/class.png $outputfolder/$locus/$team/$project/$project$run/$samplename/wordcloud.class.png");
			}
			if (-e "$locus/$samplename/order.png") {
				system("cp $locus/$samplename/order.png $outputfolder/$locus/$team/$project/$project$run/$samplename/wordcloud.order.png");
			}
			if (-e "$locus/$samplename/family.png") {
				system("cp $locus/$samplename/family.png $outputfolder/$locus/$team/$project/$project$run/$samplename/wordcloud.family.png");
			}
			if (-e "$locus/$samplename/genus.png") {
				system("cp $locus/$samplename/genus.png $outputfolder/$locus/$team/$project/$project$run/$samplename/wordcloud.genus.png");
			}
			if (-e "$locus/$samplename/species.png") {
				system("cp $locus/$samplename/species.png $outputfolder/$locus/$team/$project/$project$run/$samplename/wordcloud.species.png");
			}
			foreach my $communityfile ('community_qc_target.tsv', 'community_qc3nn_target.tsv', 'community_qc_nontarget.tsv', 'community_qc3nn_nontarget.tsv', 'community_standard.tsv') {
				my $tempcommunityfile = $communityfile;
				$tempcommunityfile =~ s/_target/_fishes/;
				if (-e "$locus/$communityfile") {
					if ((-s "$locus/$communityfile") > 27) {
						system("head -n 1 $locus/$communityfile > $outputfolder/$locus/$team/$project/$project$run/$samplename/$communityfile");
						system("grep -P '^$samplename\\t' $locus/$communityfile >> $outputfolder/$locus/$team/$project/$project$run/$samplename/$communityfile");
					}
				}
				elsif (-e "$locus/$tempcommunityfile") {
					if ((-s "$locus/$tempcommunityfile") > 27) {
						system("head -n 1 $locus/$tempcommunityfile > $outputfolder/$locus/$team/$project/$project$run/$samplename/$communityfile");
						system("grep -P '^$samplename\\t' $locus/$tempcommunityfile >> $outputfolder/$locus/$team/$project/$project$run/$samplename/$communityfile");
					}
				}
				if (-e "$outputfolder/$locus/$team/$project/$project$run/$samplename/$communityfile") {
					system("xz -9e $outputfolder/$locus/$team/$project/$project$run/$samplename/$communityfile");
				}
			}
			if (exists($negativesample{$samplename})) {
				open($filehandleoutput1, "> $outputfolder/$locus/$team/$project/$project$run/$samplename/README.txt") or die(__LINE__ . ": Cannot open \"$outputfolder/$locus/$team/$project/$project$run/$samplename/README.txt\".\n");
				print($filehandleoutput1 "This sample is NEGATIVE CONTROL and/or ERRONEOUS SAMPLE. DO NOT USE for normal analysis.");
				close($filehandleoutput1);
			}
			if ($samplename =~ /__[ACGT]+\+[ACGT]+__/) {
				open($filehandleoutput1, "> $outputfolder/$locus/$team/$project/$project$run/$samplename/README.txt") or die(__LINE__ . ": Cannot open \"$outputfolder/$locus/$team/$project/$project$run/$samplename/README.txt\".\n");
				print($filehandleoutput1 "This sample may be caused by index-hopping. DO NOT USE for normal analysis.");
				close($filehandleoutput1);
			}
		}
		else {
			print(STDERR __LINE__ . ": \"$samplename\" ($samplename) was multiply used.\n");
		}
	}
}

sub readSeq {
	my $seqfile = shift(@_);
	my @list;
	$filehandleinput1 = &readFile($seqfile);
	while (<$filehandleinput1>) {
		s/\r?\n?$//;
		if (/^> *(.+)/) {
			my $seqname = $1;
			$seqname =~ s/;+size=\d+;*//g;
			push(@list, $seqname);
		}
	}
	close($filehandleinput1);
	return(@list);
}
