use strict;
use Time::Piece;
use utf8;
use open ':encoding(utf8)';
use open ':std';

my %replacetable;
my @inputfiles;
for (my $i = 0; $i < (scalar(@ARGV) - 2); $i ++) {
	$inputfiles[$i] = $ARGV[$i];
}
my $inputfile = $ARGV[-2];
my $outputfile = $ARGV[-1];

if ($inputfile eq $outputfile) {
	die(__LINE__ . ": Cannot set the same file name to output file as that of input file.\n");
}

foreach my $tempinputfile (@inputfiles) {
	open(IN, "< $tempinputfile") or die(__LINE__ . ": Cannot open \"$tempinputfile\".\n");
	while (<IN>) {
		my $line = $_;
		$line =~ s/\r?\n?$//;
		my @row = split(/\t/, $line);
		if (scalar(@row) == 2) {
			$replacetable{$row[0]} = $row[1];
		}
	}
	close(IN);
}

{
	my @label;
	my $lineno = 1;
	open(IN, "< $inputfile") or die(__LINE__ . ": Cannot open \"$inputfile\".\n");
	open(OUT, "> $outputfile") or die(__LINE__ . ": Cannot open \"$outputfile\".\n");
	while (<IN>) {
		my $line = $_;
		$line =~ s/\r?\n?$//;
		$line =~ tr/　！”＃＄％＆’（）＊＋，－．／０-９：；＜＝＞？＠Ａ-Ｚ［￥］＾＿｀ａ-ｚ｛｜｝/ -}/;
		$line =~ tr/ーｰ/--/;
		my @row = split(/\t/, $line);
		if ($lineno == 1) {
			@label = @row;
			print(OUT join("\t", @label) . "\n");
		}
		elsif (@label) {
			my %table;
			for (my $i = 0; $i < scalar(@row); $i ++) {
				if ($row[$i]) {
					$table{$label[$i]} = $row[$i];
				}
			}
			if (exists($table{'locality'}) && exists($replacetable{$table{'locality'}}) && exists($table{'date'}) && exists($table{'time'}) && (exists($table{'sterivex1'}) || exists($table{'sterivex2'}) || exists($table{'sterivexNC'}))) {
				$table{'locality'} = $replacetable{$table{'locality'}};
				my $t = Time::Piece->strptime($table{'date'} . 'T' . $table{'time'} . '+0900', '%Y-%m-%dT%H:%M%z');
				$table{'samplename'} = $t->strftime('%Y%m%dT%H%M-') . $table{'locality'};
				if (exists($table{'correspondence'}) && exists($replacetable{$table{'correspondence'}})) {
					$table{'correspondence'} = $replacetable{$table{'correspondence'}};
				}
				if (exists($table{'collector1'}) && exists($replacetable{$table{'collector1'}})) {
					$table{'collector1'} = $replacetable{$table{'collector1'}};
				}
				if (exists($table{'collector2'}) && exists($replacetable{$table{'collector2'}})) {
					$table{'collector2'} = $replacetable{$table{'collector2'}};
				}
				if (exists($table{'collector3'}) && exists($replacetable{$table{'collector3'}})) {
					$table{'collector3'} = $replacetable{$table{'collector3'}};
				}
				if (exists($table{'collector4'}) && exists($replacetable{$table{'collector4'}})) {
					$table{'collector4'} = $replacetable{$table{'collector4'}};
				}
				if (exists($table{'BAC'}) && $table{'BAC'} eq '使用') {
					$table{'BAC'} = 'yes';
				}
				if (exists($table{'BAC'}) && $table{'BAC'} eq '不使用') {
					$table{'BAC'} = 'no';
				}
				if (exists($table{'lat'}) && ($table{'lat'} < 20 || $table{'lat'} > 50) && $table{'lat'} ne 'na') {
					print("Latitude is weired at line" . __LINE__ . ".\n$line\n");
				}
				if (exists($table{'lon'}) && ($table{'lon'} < 120 || $table{'lon'} > 150) && $table{'lon'} ne 'na') {
					print("Longitude is weired at line" . __LINE__ . ".\n$line\n");
				}
				if (exists($table{'watertemp'}) && ($table{'watertemp'} < -5 || $table{'watertemp'} > 50) && $table{'watertemp'} ne 'na') {
					print("Water temperature is weired at line" . __LINE__ . ".\n$line\n");
				}
				if (exists($table{'salinity'}) && $table{'salinity'} > 10 && $table{'salinity'} < 50) {
					$table{'salinity'} = $table{'salinity'} / 10;
				}
				if (exists($table{'salinity'}) && ($table{'salinity'} < 0 || $table{'salinity'} > 5) && $table{'salinity'} ne 'na') {
					print("Salinity is weired at line" . __LINE__ . ".\n$line\n");
				}
				for (my $i = 0; $i < scalar(@label); $i ++) {
					if (exists($table{$label[$i]}) && $table{$label[$i]} ne 'na') {
						print(OUT $table{$label[$i]});
					}
					if ($i + 1 == scalar(@label)) {
						print(OUT "\n");
					}
					else {
						print(OUT "\t");
					}
				}
			}
		}
		else {
			print(STDERR "line $lineno\n");
			die(__LINE__ . ": Unknown error.\n");
		}
		$lineno ++;
	}
	close(OUT);
	close(IN);
}
