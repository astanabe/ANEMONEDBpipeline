use strict;
use Cwd 'getcwd';
use Time::Piece;
use utf8;
use open ':encoding(utf8)';
use open ':std';

my $filehandleinput1;
my $filehandleoutput1;

my @inputfiles;
for (my $i = 0; $i < (scalar(@ARGV) - 1); $i ++) {
	$inputfiles[$i] = $ARGV[$i];
}
my $outputfolder = $ARGV[-1];
my $solutionvoltable = 'solutionvoltable.tsv';
my $watervoltable = 'watervoltable.tsv';
my $stdconctable = 'stdconctable.tsv';
my $blanklist = 'blanklist.txt';
my $negativesamplelist = 'negativesamplelist.txt';
my $replicatelist = 'replicatelist.txt';

my %solutionvol;
my %watervol;
my %stdconc;
my %blank;
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
				my $tempsamplename = shift(@row);
				@{$solutionvol{$tempsamplename}} = @row;
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
				my $tempsamplename = shift(@row);
				@{$watervol{$tempsamplename}} = @row;
			}
		}
		close($filehandleinput1);
	}
	if (-e $blanklist) {
		open($filehandleinput1, "< $blanklist") or die(__LINE__ . ": Cannot open \"$blanklist\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			if ($_) {
				$blank{$_} = 1;
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
my $team;
my $project;
my $run;
my $locus;
{
	my @path = split('/', $root);
	$team = $path[-4];
	$project = $path[-3];
	$run = $path[-2];
	$locus = $path[-1];
}

my %table;
foreach my $tempinputfile (@inputfiles) {
	my @label;
	my $lineno = 1;
	open(IN, "< $tempinputfile") or die(__LINE__ . ": Cannot open \"$tempinputfile\".\n");
	while (<IN>) {
		my $line = $_;
		$line =~ s/\r?\n?$//;
		my @row = split(/\t/, $line);
		if ($lineno == 1) {
			@label = @row;
		}
		elsif (@label) {
			for (my $i = 0; $i < scalar(@row); $i ++) {
				if ($row[$i] && $row[$i] !~ /^\s*$/) {
					$table{$line}{$label[$i]} = $row[$i];
				}
			}
		}
		else {
			print(STDERR "line $lineno\n");
			die(__LINE__ . ": Unknown error.\n");
		}
		$lineno ++;
	}
	close(IN);
}

my %samplename2organization;
foreach my $line (keys(%table)) {
	if (exists($table{$line}{'samplename'}) && exists($table{$line}{'organization'})) {
		$samplename2organization{$table{$line}{'samplename'}} = $table{$line}{'organization'};
	}
}

{
	if (!-e $outputfolder) {
		mkdir($outputfolder);
	}
	if (-e $outputfolder && !-d $outputfolder) {
		die(__LINE__ . ": \"$outputfolder\" is not a directory.\n");
	}
	foreach my $samplename (@samplenames) {
		if (!exists($samplename2organization{$samplename})) {
			$samplename2organization{$samplename} = 'OTHERS';
		}
		if (!-e "$outputfolder/$samplename2organization{$samplename}") {
			if (!mkdir("$outputfolder/$samplename2organization{$samplename}")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$samplename2organization{$samplename}\".\n");
			}
		}
		if (!-e "$outputfolder/$samplename2organization{$samplename}/$locus") {
			if (!mkdir("$outputfolder/$samplename2organization{$samplename}/$locus")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$samplename2organization{$samplename}/$locus\".\n");
			}
		}
		if (!-e "$outputfolder/$samplename2organization{$samplename}/$locus/$team") {
			if (!mkdir("$outputfolder/$samplename2organization{$samplename}/$locus/$team")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$samplename2organization{$samplename}/$locus/$team\".\n");
			}
		}
		if (!-e "$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project") {
			if (!mkdir("$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project\".\n");
			}
		}
		if (!-e "$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run") {
			if (!mkdir("$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run\".\n");
			}
		}
		if (!-e "$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename") {
			if (!mkdir("$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename")) {
				die(__LINE__ . ": Cannot make \"$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename\".\n");
			}
			foreach my $strand ('.forward', '.reverse', '') {
				if ($daughtersample{$samplename}) {
					foreach my $daughtersample (@{$daughtersample{$samplename}}) {
						if (-e "demultiplexed/$daughtersample$strand.fastq.xz") {
							system("cp demultiplexed/$samplename$strand.fastq.xz $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$samplename$strand.fastq.xz");
						}
					}
				}
				elsif (-e "demultiplexed/$samplename$strand.fastq.xz") {
					system("cp demultiplexed/$samplename$strand.fastq.xz $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$samplename$strand.fastq.xz");
				}
			}
			if (-e "$samplename/phylum.png") {
				system("cp $samplename/phylum.png $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/wordcloud.phylum.png");
			}
			if (-e "$samplename/class.png") {
				system("cp $samplename/class.png $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/wordcloud.class.png");
			}
			if (-e "$samplename/order.png") {
				system("cp $samplename/order.png $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/wordcloud.order.png");
			}
			if (-e "$samplename/family.png") {
				system("cp $samplename/family.png $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/wordcloud.family.png");
			}
			if (-e "$samplename/genus.png") {
				system("cp $samplename/genus.png $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/wordcloud.genus.png");
			}
			if (-e "$samplename/species.png") {
				system("cp $samplename/species.png $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/wordcloud.species.png");
			}
			foreach my $communityfile ('community_qc_target.tsv', 'community_qc3nn_target.tsv', 'community_qc_nontarget.tsv', 'community_qc3nn_nontarget.tsv', 'community_standard.tsv') {
				my $tempcommunityfile = $communityfile;
				$tempcommunityfile =~ s/_target/_fishes/;
				if (-e "$communityfile") {
					if ((-s "$communityfile") > 27) {
						system("head -n 1 $communityfile > $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$communityfile");
						system("grep -P '^$samplename\\t' $communityfile >> $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$communityfile");
					}
				}
				elsif (-e "$tempcommunityfile") {
					if ((-s "$tempcommunityfile") > 27) {
						system("head -n 1 $tempcommunityfile > $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$communityfile");
						system("grep -P '^$samplename\\t' $tempcommunityfile >> $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$communityfile");
					}
				}
				if (-e "$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$communityfile") {
					system("xz -9e $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/$communityfile");
				}
			}
			if (exists($blank{$samplename})) {
				open($filehandleoutput1, "> $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/README.txt") or die(__LINE__ . ": Cannot open \"$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/README.txt\".\n");
				print($filehandleoutput1 "This sample is FIELD BLANK. DO NOT USE for normal analysis.");
				close($filehandleoutput1);
			}
			elsif (exists($negativesample{$samplename})) {
				open($filehandleoutput1, "> $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/README.txt") or die(__LINE__ . ": Cannot open \"$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/README.txt\".\n");
				print($filehandleoutput1 "This sample is NEGATIVE CONTROL and/or ERRONEOUS SAMPLE. DO NOT USE for normal analysis.");
				close($filehandleoutput1);
			}
			if ($samplename =~ /__[ACGT]+\+[ACGT]+__/) {
				open($filehandleoutput1, "> $outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/README.txt") or die(__LINE__ . ": Cannot open \"$outputfolder/$samplename2organization{$samplename}/$locus/$team/$project/$project$run/$samplename/README.txt\".\n");
				print($filehandleoutput1 "This sample may be caused by index-hopping. DO NOT USE for normal analysis.");
				close($filehandleoutput1);
			}
		}
		else {
			print(STDERR __LINE__ . ": \"$samplename\" ($samplename) was multiply used.\n");
		}
	}
}
