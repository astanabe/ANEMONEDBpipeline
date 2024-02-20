use strict;

my $inputfileA = $ARGV[0];
my $inputfileB = $ARGV[1];
my $outputfile = $ARGV[2];

my $label;

my %tableA;
my %tableB;
my @temp;

{
	my $lineno = 1;
	open(IN, "< $inputfileA") or die(__LINE__ . "Cannot open file.");
	while (<IN>) {
		s/\r?\n?$//;
		my @row;
		if (/^(.+)\t(\S+)$/) {
			@row = ($1, $2);
		}
		if ($lineno == 1) {
			$label = $row[0];
		}
		else {
			$tableA{$row[0]} = $row[1];
			push(@temp, $row[0]);
		}
		$lineno ++;
	}
	close(IN);
}

{
	my $lineno = 1;
	open(IN, "< $inputfileB") or die(__LINE__ . "Cannot open file.");
	while (<IN>) {
		s/\r?\n?$//;
		my @row;
		if (/^(.+)\t(\S+)$/) {
			@row = ($1, $2);
		}
		if ($lineno == 1) {
			if ($row[0] ne $label) {
				die(__LINE__ . "The label is not matched.");
			}
		}
		else {
			$tableB{$row[0]} = $row[1];
		}
		$lineno ++;
	}
	close(IN);
}

{
	open(OUT, "> $outputfile") or die(__LINE__ . "Cannot write file.");
	print(OUT "$label\tnreads\tncopiesperml\n");
	foreach my $temp (@temp) {
		print(OUT $temp);
		if (exists($tableA{$temp}) && $tableA{$temp} > 0) {
			print(OUT "\t" . $tableA{$temp});
		}
		else {
			die();
		}
		if (exists($tableB{$temp}) && $tableB{$temp} > 0) {
			print(OUT "\t" . $tableB{$temp} . "\n");
		}
		else {
			print(OUT "\tNA\n");
		}
	}
	close(OUT);
}
