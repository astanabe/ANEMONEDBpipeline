use strict;
use Cwd 'getcwd';
use Time::Piece;
use utf8;
use open ':encoding(utf8)';
use open ':std';

my $filehandleinput1;
my $filehandleoutput1;

my %table;
my %sterivex2samplename;
my @inputfiles;
for (my $i = 0; $i < (scalar(@ARGV) - 2); $i ++) {
	$inputfiles[$i] = $ARGV[$i];
}
my $inputfile = $ARGV[-2];
my $outputfolder = $ARGV[-1];
my $solutionvoltable = 'solutionvoltable.tsv';
my $watervoltable = 'watervoltable.tsv';
my $stdconctable = 'stdconctable.tsv';
my $negativesamplelist = 'negativesamplelist.txt';

my @sampleinfo = ('samp_name', 'project_name', 'correspondence', 'worldmesh', 'lat_lon', 'samp_taxon_id', 'samp_collec_device', 'samp_mat_process', 'collection_date', 'collection_date_local', 'collection_date_utc', 'temp', 'salinity', 'samp_size', 'size_frac_low', 'num_filter', 'filter_prod_name', 'biocide_used', 'samp_store_temp', 'samp_store_sol');
my %sampletable;

my @experimentinfo = ('samp_vol_we_dna_ext', 'dna_extract_vol', 'pool_dna_extracts', 'nucl_acid_ext', 'nucl_acid_amp', 'lib_layout', 'target_gene', 'pcr_primers', 'pcr_primers2', 'pcr_cond', 'pcr_cond2', 'seq_meth', 'sop');
my %experimenttable;

my %solutionvol;
my %watervol;
my %stdconc;
my %negativesample;
{
	if (-e $stdconctable) {
		my $lineno = 1;
		my @label;
		open($filehandleinput1, "< $stdconctable") or die(__LINE__ . ": Cannot open \"$stdconctable\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			my @row = split(/\t/, $_);
			if ($lineno ==1 && scalar(@row) > 2 && $_ !~ /\t(?:\d+|\d+\.\d+)\t/ && $_ !~ /\t(?:\d+|\d+\.\d+)$/) {
				@label = @row;
			}
			elsif ($lineno > 1 && @label && scalar(@row) == scalar(@label)) {
				for (my $i = 1; $i < scalar(@label); $i ++) {
					if ($row[$i] =~ /^(?:\d+|\d+\.\d+)$/) {
						$stdconc{$row[0]}{$label[$i]} = $row[$i];
					}
					else {
						die(__LINE__ . ": \"$stdconctable\" is invalid.\n");
					}
				}
			}
			elsif (scalar(@row) == 2 && $row[1] =~ /^(?:\d+|\d+\.\d+)$/) {
				$stdconc{$row[0]} = $row[1];
			}
			elsif (scalar(@row) == 3 && $row[2] =~ /^(?:\d+|\d+\.\d+)$/) {
				$stdconc{$row[0]}{$row[1]} = $row[2];
			}
			elsif (@row) {
				die(__LINE__ . ": \"$stdconctable\" is invalid.\n");
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
			if (scalar(@row) == 2 && $row[1] =~ /^(?:\d+|\d+\.\d+)$/) {
				$solutionvol{$row[0]} = $row[1];
			}
		}
		close($filehandleinput1);
	}
	if (-e $watervoltable) {
		open($filehandleinput1, "< $watervoltable") or die(__LINE__ . ": Cannot open \"$watervoltable\".\n");
		while (<$filehandleinput1>) {
			s/\r?\n?$//;
			my @row = split(/\t/, $_);
			if (scalar(@row) == 2 && $row[1] =~ /^(?:\d+|\d+\.\d+)$/) {
				$watervol{$row[0]} = $row[1];
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
}
my $root = getcwd();
my $locus;
my $team;
my $project;
my $run;
{
	my @path = split('/', $root);
	$locus = $path[-4];
	$team = $path[-3];
	$project = $path[-2];
	$run = $path[-1];
}

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
				if ($row[$i] && $row[$i] !~ /^\s+$/ && $label[$i] =~ /sterivex/) {
					$row[$i] =~ s/[\-ー－―‐]0*//g;
					$table{$line}{$label[$i]} = lc($row[$i]);
				}
				elsif ($row[$i]) {
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

foreach my $line (keys(%table)) {
	if (exists($table{$line}{'samplename'})) {
		my @samplenames;
		if (exists($table{$line}{'sterivex1'}) && exists($table{$line}{'watervol1'}) && exists($table{$line}{'sterivex2'}) && exists($table{$line}{'watervol2'}) && exists($table{$line}{'sterivexNC'}) && exists($table{$line}{'watervolNC'})) {
			$sterivex2samplename{$table{$line}{'sterivex1'}} = $table{$line}{'samplename'};
			$sterivex2samplename{$table{$line}{'sterivex2'}} = $table{$line}{'samplename'};
			$sterivex2samplename{$table{$line}{'sterivexNC'}} = $table{$line}{'samplename'} . '-NC';
			push(@samplenames, $table{$line}{'samplename'});
			push(@samplenames, $table{$line}{'samplename'} . '-NC');
		}
		elsif (exists($table{$line}{'sterivex1'}) && exists($table{$line}{'watervol1'}) && exists($table{$line}{'sterivex2'}) && exists($table{$line}{'watervol2'})) {
			$sterivex2samplename{$table{$line}{'sterivex1'}} = $table{$line}{'samplename'};
			$sterivex2samplename{$table{$line}{'sterivex2'}} = $table{$line}{'samplename'};
			push(@samplenames, $table{$line}{'samplename'});
		}
		elsif (exists($table{$line}{'sterivex1'}) && exists($table{$line}{'watervol1'}) && exists($table{$line}{'sterivexNC'}) && exists($table{$line}{'watervolNC'})) {
			$sterivex2samplename{$table{$line}{'sterivex1'}} = $table{$line}{'samplename'};
			$sterivex2samplename{$table{$line}{'sterivexNC'}} = $table{$line}{'samplename'} . '-NC';
			push(@samplenames, $table{$line}{'samplename'});
			push(@samplenames, $table{$line}{'samplename'} . '-NC');
		}
		elsif (exists($table{$line}{'sterivex2'}) && exists($table{$line}{'watervol2'}) && exists($table{$line}{'sterivexNC'}) && exists($table{$line}{'watervolNC'})) {
			$sterivex2samplename{$table{$line}{'sterivex2'}} = $table{$line}{'samplename'};
			$sterivex2samplename{$table{$line}{'sterivexNC'}} = $table{$line}{'samplename'} . '-NC';
			push(@samplenames, $table{$line}{'samplename'});
			push(@samplenames, $table{$line}{'samplename'} . '-NC');
		}
		elsif (exists($table{$line}{'sterivex1'}) && exists($table{$line}{'watervol1'})) {
			$sterivex2samplename{$table{$line}{'sterivex1'}} = $table{$line}{'samplename'};
			push(@samplenames, $table{$line}{'samplename'});
		}
		elsif (exists($table{$line}{'sterivex2'}) && exists($table{$line}{'watervol2'})) {
			$sterivex2samplename{$table{$line}{'sterivex2'}} = $table{$line}{'samplename'};
			push(@samplenames, $table{$line}{'samplename'});
		}
		elsif (exists($table{$line}{'sterivexNC'}) && exists($table{$line}{'watervolNC'})) {
			if ($table{$line}{'samplename'} =~ /-NC$/) {
				$sterivex2samplename{$table{$line}{'sterivexNC'}} = $table{$line}{'samplename'};
				push(@samplenames, $table{$line}{'samplename'});
			}
			else {
				$sterivex2samplename{$table{$line}{'sterivexNC'}} = $table{$line}{'samplename'} . '-NC';
				push(@samplenames, $table{$line}{'samplename'} . '-NC');
			}
		}
		foreach my $samplename (@samplenames) {
			foreach my $sampleinfo (@sampleinfo) {
				if ($sampleinfo eq 'samp_name') {
					$sampletable{$samplename}{$sampleinfo} = $samplename;
				}
				elsif ($sampleinfo eq 'project_name') {
					$sampletable{$samplename}{$sampleinfo} = "$team $project";
				}
				elsif ($sampleinfo eq 'correspondence' && exists($table{$line}{'correspondence'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'correspondence'};
				}
				elsif ($sampleinfo eq 'correspondence' && exists($table{$line}{'collector1'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'collector1'};
				}
				elsif ($sampleinfo eq 'samp_taxon_id' && $samplename =~ /-NC$/) {
					$sampletable{$samplename}{$sampleinfo} = 'blank sample';
				}
				elsif ($sampleinfo eq 'samp_taxon_id') {
					$sampletable{$samplename}{$sampleinfo} = 'metagenome';
				}
				elsif ($sampleinfo eq 'worldmesh' && $samplename !~ /-NC$/ && exists($table{$line}{'lat'}) && exists($table{$line}{'lon'})) {
					$sampletable{$samplename}{$sampleinfo} = &meshcode2($table{$line}{'lat'}, $table{$line}{'lon'});
				}
				elsif ($sampleinfo eq 'lat_lon' && $samplename !~ /-NC$/ && exists($table{$line}{'lat'}) && exists($table{$line}{'lon'})) {
					my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode2latlon(&meshcode2($table{$line}{'lat'}, $table{$line}{'lon'}));
					$sampletable{$samplename}{$sampleinfo} = "$latC $lonC";
				}
				elsif ($sampleinfo eq 'samp_collec_device' && $samplename !~ /-NC$/ && exists($table{$line}{'locality'}) && $table{$line}{'locality'} =~ /-(?:Bottom|Intermediate)/) {
					$sampletable{$samplename}{$sampleinfo} = 'bottle';
				}
				elsif ($sampleinfo eq 'samp_collec_device' && $samplename !~ /-NC$/) {
					$sampletable{$samplename}{$sampleinfo} = 'bucket';
				}
				elsif ($sampleinfo eq 'samp_mat_process' && $samplename =~ /GFF/) {
					$sampletable{$samplename}{$sampleinfo} = 'filtering water, storing in freezer';
				}
				elsif ($sampleinfo eq 'samp_mat_process') {
					$sampletable{$samplename}{$sampleinfo} = 'filtering water, adding RNAlater to the filter, storing in freezer';
				}
				elsif ($sampleinfo eq 'collection_date' && exists($table{$line}{'date'}) && exists($table{$line}{'time'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'date'} . 'T' . $table{$line}{'time'} . ':00+09:00';
				}
				elsif ($sampleinfo eq 'collection_date_local' && exists($table{$line}{'date'}) && exists($table{$line}{'time'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'date'} . 'T' . $table{$line}{'time'} . ':00+09:00';
				}
				elsif ($sampleinfo eq 'collection_date_utc' && exists($table{$line}{'date'}) && exists($table{$line}{'time'})) {
					my $t = Time::Piece->strptime($table{$line}{'date'} . 'T' . $table{$line}{'time'} . '+0900', '%Y-%m-%dT%H:%M%z');
					$sampletable{$samplename}{$sampleinfo} = $t->strftime('%Y-%m-%dT%H:%M:00Z');
				}
				elsif ($sampleinfo eq 'temp' && exists($table{$line}{'watertemp'}) && $samplename !~ /-NC$/) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'watertemp'};
				}
				elsif ($sampleinfo eq 'salinity' && exists($table{$line}{'salinity'}) && $samplename !~ /-NC$/) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'salinity'};
				}
				elsif ($sampleinfo eq 'samp_store_temp') {
					$sampletable{$samplename}{$sampleinfo} = -20;
				}
				elsif ($sampleinfo eq 'samp_store_sol' && $samplename !~ /GFF/) {
					$sampletable{$samplename}{$sampleinfo} = 'RNAlater';
				}
				elsif ($sampleinfo eq 'samp_size' && $samplename =~ /-NC$/ && exists($table{$line}{'watervolNC'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'watervolNC'};
				}
				elsif ($sampleinfo eq 'samp_size' && exists($table{$line}{'watervol1'}) && exists($table{$line}{'watervol2'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'watervol1'} + $table{$line}{'watervol2'};
				}
				elsif ($sampleinfo eq 'samp_size' && exists($table{$line}{'watervol1'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'watervol1'};
				}
				elsif ($sampleinfo eq 'samp_size' && exists($table{$line}{'watervol2'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'watervol2'};
				}
				elsif ($sampleinfo eq 'size_frac_low' && $samplename =~ /GFF/) {
					$sampletable{$samplename}{$sampleinfo} = 0.7;
				}
				elsif ($sampleinfo eq 'size_frac_low') {
					$sampletable{$samplename}{$sampleinfo} = 0.45;
				}
				elsif ($sampleinfo eq 'num_filter' && exists($table{$line}{'nfilter'})) {
					$sampletable{$samplename}{$sampleinfo} = $table{$line}{'nfilter'};
				}
				elsif ($sampleinfo eq 'num_filter' && $samplename =~ /-NC$/) {
					$sampletable{$samplename}{$sampleinfo} = 1;
				}
				elsif ($sampleinfo eq 'num_filter') {
					if (exists($table{$line}{'sterivex1'}) && exists($table{$line}{'sterivex2'})) {
						$sampletable{$samplename}{$sampleinfo} = 2;
					}
					else {
						$sampletable{$samplename}{$sampleinfo} = 1;
					}
				}
				elsif ($sampleinfo eq 'filter_prod_name' && $samplename =~ /GFF/) {
					$sampletable{$samplename}{$sampleinfo} = 'GF/F';
				}
				elsif ($sampleinfo eq 'filter_prod_name') {
					$sampletable{$samplename}{$sampleinfo} = 'Sterivex-HV';
				}
				elsif ($sampleinfo eq 'biocide_used' && exists($table{$line}{'BAC'}) && $table{$line}{'salinity'} eq 'yes') {
					$sampletable{$samplename}{$sampleinfo} = 'Osban S';
				}
			}
			foreach my $experimentinfo (@experimentinfo) {
				if ($experimentinfo eq 'samp_vol_we_dna_ext' && $samplename =~ /-NC$/ && exists($table{$line}{'watervolNC'})) {
					$experimenttable{$samplename}{$experimentinfo} = $table{$line}{'watervolNC'};
				}
				elsif ($experimentinfo eq 'samp_vol_we_dna_ext' && exists($table{$line}{'watervol1'}) && exists($table{$line}{'watervol2'})) {
					$experimenttable{$samplename}{$experimentinfo} = $table{$line}{'watervol1'} + $table{$line}{'watervol2'};
				}
				elsif ($experimentinfo eq 'samp_vol_we_dna_ext' && exists($table{$line}{'watervol1'})) {
					$experimenttable{$samplename}{$experimentinfo} = $table{$line}{'watervol1'};
				}
				elsif ($experimentinfo eq 'samp_vol_we_dna_ext' && exists($table{$line}{'watervol2'})) {
					$experimenttable{$samplename}{$experimentinfo} = $table{$line}{'watervol2'};
				}
				elsif ($experimentinfo eq 'pool_dna_extracts' && $samplename =~ /-NC$/) {
					$experimenttable{$samplename}{$experimentinfo} = 1;
				}
				elsif ($experimentinfo eq 'pool_dna_extracts') {
					if (exists($table{$line}{'sterivex1'}) && exists($table{$line}{'sterivex2'})) {
						$experimenttable{$samplename}{$experimentinfo} = 2;
					}
					else {
						$experimenttable{$samplename}{$experimentinfo} = 1;
					}
				}
				elsif ($experimentinfo eq 'dna_extract_vol' && $samplename =~ /-NC$/) {
					$experimenttable{$samplename}{$experimentinfo} = 200;
				}
				elsif ($experimentinfo eq 'dna_extract_vol') {
					if (exists($table{$line}{'sterivex1'}) && exists($table{$line}{'sterivex2'})) {
						$experimenttable{$samplename}{$experimentinfo} = 400;
					}
					else {
						$experimenttable{$samplename}{$experimentinfo} = 200;
					}
				}
				elsif ($experimentinfo eq 'nucl_acid_ext') {
					$experimenttable{$samplename}{$experimentinfo} = 'https://ednasociety.org/manual/';
				}
				elsif ($experimentinfo eq 'nucl_acid_amp') {
					$experimenttable{$samplename}{$experimentinfo} = 'https://ednasociety.org/manual/';
				}
				elsif ($experimentinfo eq 'lib_layout') {
					$experimenttable{$samplename}{$experimentinfo} = 'paired';
				}
				elsif ($experimentinfo eq 'target_gene') {
					$experimenttable{$samplename}{$experimentinfo} = '12S rRNA';
				}
				elsif ($experimentinfo eq 'pcr_primers') {
					if ($project =~ /^(?:2017|2018)/) {
						$experimenttable{$samplename}{$experimentinfo} = 'FWD:ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNGTTGGTAAATCTCGTGCCAGC;REV:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNCATAGTGGGGTATCTAATCCTAGTTTG;FWD:ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNGTCGGTAAAACTCGTGCCAGC;REV:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNCATAGTGGGGTATCTAATCCCAGTTTG;';
					}
					else {
						$experimenttable{$samplename}{$experimentinfo} = 'FWD:ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNRGTTGGTAAATCTCGTGCCAGC;REV:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNGCATAGTGGGGTATCTAATCCTAGTTTG;FWD:ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNGTCGGTAAAACTCGTGCCAGC;REV:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNCATAGTGGGGTATCTAATCCCAGTTTG;FWD:ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNGCCGGTAAAACTCGTGCCAGC;REV:GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNCATAGGAGGGTGTCTAATCCCCGTTTG;';
					}
				}
				elsif ($experimentinfo eq 'pcr_primers2') {
					$experimenttable{$samplename}{$experimentinfo} = 'FWD:AATGATACGGCGACCACCGAGATCTACACXXXXXXXXACACTCTTTCCCTACACGACGCTCTTCCGATCT;REV:CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT;';
				}
				elsif ($experimentinfo eq 'pcr_cond') {
					$experimenttable{$samplename}{$experimentinfo} = 'initial denaturation:95_3;denaturation:98_0.33;annealing:65_0.25;elongation:72_0.25;final elongation:72_5;35';
				}
				elsif ($experimentinfo eq 'pcr_cond2') {
					$experimenttable{$samplename}{$experimentinfo} = 'initial denaturation:95_3;denaturation:98_0.33;annealing:65_0.25;elongation:72_0.25;final elongation:72_5;10';
				}
				elsif ($experimentinfo eq 'seq_meth') {
					if ($project =~ /^(?:2019|2020|2021)/ && $run =~ /^RUN01$/) {
						$experimenttable{$samplename}{$experimentinfo} = 'NextSeq 500';
					}
					else {
						$experimenttable{$samplename}{$experimentinfo} = 'MiSeq';
					}
				}
				elsif ($experimentinfo eq 'sop') {
					$experimenttable{$samplename}{$experimentinfo} = 'https://anemone.bio/';
				}
			}
		}
	}
}

{
	if (!-e $outputfolder) {
		mkdir($outputfolder);
	}
	if (-e $outputfolder && !-d $outputfolder) {
		die(__LINE__ . ": \"$outputfolder\" is not a directory.\n");
	}
	open(IN, "< $inputfile") or die(__LINE__ . ": Cannot open \"$inputfile\".\n");
	while (<IN>) {
		my $line = $_;
		$line =~ s/\r?\n?$//;
		if ($line =~ /^\S+$/) {
			my $analysisname = $line;
			my @temp = split(/__/, lc($analysisname));
			$temp[1] =~ s/[\-ー－―‐]0*//g;
			if (exists($sterivex2samplename{$temp[1]})) {
				my $newsamplename = $sterivex2samplename{$temp[1]};
				if (!-e "$outputfolder/$locus") {
					if (!mkdir("$outputfolder/$locus")) {
						die(__LINE__ . ": Cannot make \"$outputfolder/$locus\".");
					}
				}
				if (!-e "$outputfolder/$locus/$team") {
					if (!mkdir("$outputfolder/$locus/$team")) {
						die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team\".");
					}
				}
				if (!-e "$outputfolder/$locus/$team/$project") {
					if (!mkdir("$outputfolder/$locus/$team/$project")) {
						die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team/$project\".");
					}
				}
				if (!-e "$outputfolder/$locus/$team/$project/$project$run") {
					if (!mkdir("$outputfolder/$locus/$team/$project/$project$run")) {
						die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team/$project/$project$run\".");
					}
				}
				if (!-e "$outputfolder/$locus/$team/$project/$project$run/$newsamplename") {
					if (!mkdir("$outputfolder/$locus/$team/$project/$project$run/$newsamplename")) {
						die(__LINE__ . ": Cannot make \"$outputfolder/$locus/$team/$project/$project$run/$newsamplename\".");
					}
					if (-e "demultiplexed/$analysisname.forward.fastq.xz") {
						system("cp demultiplexed/$analysisname.forward.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.forward.fastq.xz");
					}
					else {
						my $temp = $analysisname;
						$temp =~ s/__$locus$//;
						system("cp demultiplexed/$temp" . "re_1__$locus.forward.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.r1.forward.fastq.xz");
						system("cp demultiplexed/$temp" . "re_2__$locus.forward.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.r2.forward.fastq.xz");
						system("cp demultiplexed/$temp" . "re_3__$locus.forward.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.r3.forward.fastq.xz");
					}
					if (-e "demultiplexed/$analysisname.reverse.fastq.xz") {
						system("cp demultiplexed/$analysisname.reverse.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.reverse.fastq.xz");
					}
					else {
						my $temp = $analysisname;
						$temp =~ s/__$locus$//;
						system("cp demultiplexed/$temp" . "re_1__$locus.reverse.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.r1.reverse.fastq.xz");
						system("cp demultiplexed/$temp" . "re_2__$locus.reverse.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.r2.reverse.fastq.xz");
						system("cp demultiplexed/$temp" . "re_3__$locus.reverse.fastq.xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$newsamplename.r3.reverse.fastq.xz");
					}
					if (-e "$analysisname/phylum.png") {
						system("cp $analysisname/phylum.png $outputfolder/$locus/$team/$project/$project$run/$newsamplename/wordcloud.phylum.png");
					}
					if (-e "$analysisname/class.png") {
						system("cp $analysisname/class.png $outputfolder/$locus/$team/$project/$project$run/$newsamplename/wordcloud.class.png");
					}
					if (-e "$analysisname/order.png") {
						system("cp $analysisname/order.png $outputfolder/$locus/$team/$project/$project$run/$newsamplename/wordcloud.order.png");
					}
					if (-e "$analysisname/family.png") {
						system("cp $analysisname/family.png $outputfolder/$locus/$team/$project/$project$run/$newsamplename/wordcloud.family.png");
					}
					if (-e "$analysisname/genus.png") {
						system("cp $analysisname/genus.png $outputfolder/$locus/$team/$project/$project$run/$newsamplename/wordcloud.genus.png");
					}
					if (-e "$analysisname/species.png") {
						system("cp $analysisname/species.png $outputfolder/$locus/$team/$project/$project$run/$newsamplename/wordcloud.species.png");
					}
					foreach my $communityfile ('community_qc_target.tsv', 'community_qc3nn_target.tsv', 'community_standard.tsv') {
						my $tempcommunityfile = $communityfile;
						$tempcommunityfile =~ s/_target/_fishes/;
						if (-e $communityfile) {
							if ((-s $communityfile) > 27) {
								system("head -n 1 $communityfile > $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$communityfile");
								system("grep -P '^$analysisname\\t' $communityfile | perl -npe 's/^$analysisname\\t/$newsamplename\\t/' >> $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$communityfile");
								system("xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$communityfile");
							}
						}
						elsif (-e $tempcommunityfile) {
							if ((-s $tempcommunityfile) > 27) {
								system("head -n 1 $tempcommunityfile > $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$communityfile");
								system("grep -P '^$analysisname\\t' $tempcommunityfile | perl -npe 's/^$analysisname\\t/$newsamplename\\t/' >> $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$communityfile");
								system("xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/$communityfile");
							}
						}
					}
					open($filehandleoutput1, "> $outputfolder/$locus/$team/$project/$project$run/$newsamplename/sample.tsv") or die(__LINE__ . ": Cannot open \"$outputfolder/$locus/$team/$project/$project$run/$newsamplename/sample.tsv\".\n");
					print($filehandleoutput1 "samplename\tkey\tvalue\n");
					foreach my $sampleinfo (@sampleinfo) {
						if (exists($sampletable{$newsamplename}{$sampleinfo}) && $sampletable{$newsamplename}{$sampleinfo}) {
							print($filehandleoutput1 $newsamplename . "\t$sampleinfo\t" . $sampletable{$newsamplename}{$sampleinfo} . "\n");
						}
					}
					close($filehandleoutput1);
					system("xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/sample.tsv");
					open($filehandleoutput1, "> $outputfolder/$locus/$team/$project/$project$run/$newsamplename/experiment.tsv") or die(__LINE__ . ": Cannot open \"$outputfolder/$locus/$team/$project/$project$run/$newsamplename/experiment.tsv\".\n");
					print($filehandleoutput1 "samplename\tkey\tvalue\n");
					foreach my $experimentinfo (@experimentinfo) {
						if (exists($experimenttable{$newsamplename}{$experimentinfo}) && $experimenttable{$newsamplename}{$experimentinfo}) {
							print($filehandleoutput1 $newsamplename . "\t$experimentinfo\t" . $experimenttable{$newsamplename}{$experimentinfo} . "\n");
						}
					}
					if (exists($stdconc{$analysisname})) {
						print($filehandleoutput1 "$newsamplename\tpcr_standard\tMiFish_STD_01:GTCGGTAAAACTCGTGCCAGCCACCGCGGTTATACGACAGGCCCAAGTTGAACGCAGTCGGCGTAAAGAGTGGTTAAAAGGGTGAGTTGCTAAAGCCGAAACCCTGCTCCGCTGTTATACGTTACGTAAAACGAGAATATCGCTTACGAAAGTAGCTTTAAATAGGATCCAGGAACCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATG;MiFish_STD_02:GTCGGTAAAACTCGTGCCAGCCACCGCGGTTATACGACAGGCCCAAGTTGATCTTGAACGGCGTAAAGAGTGGTTAGATTTCCCTACTGCTAAAGCCGAAGCACCGCCGTGCTGTTATACGTGCCTCACAGTGAGAAGGGCGAAAACGAAAGTAGCTTTATTTCCGTCACGCGAACCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATG;MiFish_STD_03:GTCGGTAAAACTCGTGCCAGCCACCGCGGTTATACGACAGGCCCAAGTTGAAGCGACGCGGCGTAAAGAGTGGTTATCACTAGAAGATCCTAAAGCCGAAGCACTTCGCCGCTGTTATACGCCCACCGCTATCGGAATCACCCCAACGAAAGTAGCTTTATGCGTTAAATCGGAACCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATG;MiFish_STD_04-2:GTCGGTAAAACTCGTGCCAGCCACCGCGGTTATACGACAGGCCCAAGTTGAGATCCCACGGCGTAAAGAGTGGTTAGAACCGGAAACGCGTAAAGCCGAAGAACATCAGTGCTGTTATACGCATTCGATTAGGTGAATTTCAGTAACGAAAGTAGCTTTACCGATCAAATCCGAACCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATG;\n");
						print($filehandleoutput1 "$newsamplename\tpcr_standard_conc\t");
						foreach my $standardname (sort(keys(%{$stdconc{$analysisname}}))) {
							my $temp = $standardname;
							$temp =~ s/_\d+copies//;
							print($filehandleoutput1 "$temp:" . $stdconc{$analysisname}{$standardname} . ";");
						}
						print($filehandleoutput1 "\n");
					}
					close($filehandleoutput1);
					system("xz $outputfolder/$locus/$team/$project/$project$run/$newsamplename/experiment.tsv");
					if (exists($negativesample{$analysisname})) {
						open($filehandleoutput1, "> $outputfolder/$locus/$team/$project/$project$run/$newsamplename/README.txt") or die(__LINE__ . ": Cannot open \"$outputfolder/$locus/$team/$project/$project$run/$newsamplename/README.txt\".\n");
						if ($newsamplename =~ /-NC$/) {
							print($filehandleoutput1 "This sample is NEGATIVE CONTROL. DO NOT USE for normal analysis.");
						}
						else {
							print($filehandleoutput1 "This sample is ERRONEOUS. DO NOT USE for normal analysis.");
						}
						close($filehandleoutput1);
					}
				}
				else {
					print(STDERR __LINE__ . ": \"$newsamplename\" ($analysisname) was multiply used.\n");
				}
			}
		}
		else {
			die(__LINE__ . ": Unknown error.\n");
		}
	}
	close(IN);
}

sub meshcode2 {
	my $lat = shift(@_);
	my $lon = shift(@_);
	my $o = 0;
	if ($lat < 0) {
		$o += 4;
	}
	if ($lon <= -100) {
		$o += 3;
	}
	elsif ($lon < 0) {
		$o += 2;
	}
	elsif ($lon >= 100) {
		$o ++;
	}
	my $z = $o % 2;
	my $y = (($o - $z) / 2) % 2;
	my $x = ($o - ($y * 2) - $z) / 4;
	my $a = ($lat * (1 - ($x * 2)) * 60) / 40;
	my $b = ($lon * (1 - ($y * 2))) - ($z * 100);
	my $p = int($a);
	my $u = int($b);
	my $q = int((($a - $p) * 40) / 5);
	my $v = int((($b - $u) * 60) / 7.5);
	return(sprintf("%d%03d%02d%d%d", ($o + 1), $p, $u, $q, $v));
}

sub meshcode2latlon {
	my $meshcode2 = shift(@_);
	my $code0 = substr($meshcode2, 0, 1);
	$code0 --;
	my $code12 = substr($meshcode2, 1, 3);
	if ($code12 =~ /^00/) {
		$code12 = substr($meshcode2, 3, 1);
	}
	elsif ($code12 =~ /^0/) {
		$code12 = substr($meshcode2, 2, 2);
	}
	else {
		$code12 = substr($meshcode2, 1, 3);
	}
	my $code34;
	if (substr($meshcode2, 4, 1) eq '0') {
		$code34 = substr($meshcode2, 5, 1);
	}
	else {
		$code34 = substr($meshcode2, 4, 2);
	}
	my $latwidth = (2 / 3) / 8;
	my $lonwidth = 1 / 8;
	my $code5 = substr($meshcode2, 6, 1);
	my $code6 = substr($meshcode2, 7, 1);
	my $z = $code0 % 2;
	my $y = (($code0 - $z) / 2) % 2;
	my $x = ($code0 - (2 * $y) - $z) / 4;
	my $latN = ($code12 * 2) / 3;
	$latN += ((($code5 - $x + 1) * 2) / 3) / 8;
	$latN *= (1 - (2 * $x));
	my $lonW = $code34 + (100 * $z);
	$lonW += ($code6 + $y) / 8;
	$lonW *= 1 - (2 * $y);
	my $dlat = (2 / 3) / 8;
	my $dlon = 1 / 8;
	my $latS = sprintf("%.8f", ($latN - $dlat));
	my $latC = sprintf("%.8f", ($latN - ($dlat / 2)));
	my $lonE = sprintf("%.8f", ($lonW + $dlon));
	my $lonC = sprintf("%.8f", ($lonW + ($dlon / 2)));
	$latN = sprintf("%.8f", $latN);
	$lonW = sprintf("%.8f", $lonW);
	return($latN, $latC, $latS, $lonW, $lonC, $lonE);
}
