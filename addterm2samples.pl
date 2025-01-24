#use strict;
use File::Spec;
my $devnull = File::Spec->devnull();

my $numthreads = 32;

my $filehandleinput1;
my $filehandleoutput1;

my $inputfolder = $ARGV[0];
if (!-d $inputfolder) {
	die(__LINE__);
}
my $url = $ARGV[1];

# store taxonomic hierarchy
my %parents;
my %havedaughters;
my %ranks;
my %rankname2taxa;
my %yearmonth;
my %taxonid;
my %yearmonthid;
my %taxon2slug;
# store number of samples of taxon-worldmesh, taxon-month-worldmesh, taxon-year-worldmesh, taxon-year-month-worldmesh
my %nsamples;
my %totalnsamples;
# for duplication check
my %combinedtaxonnames;
my %taxonnames;

my @loci;
foreach my $temp (glob("$inputfolder/*")) {
	if (-d $temp) {
		$temp =~ /([^\/]+$)/;
		push(@loci, $1);
	}
}
foreach my $locus (@loci) {
	my @teams;
	foreach my $temp (glob("$inputfolder/$locus/*")) {
		if (-d $temp) {
			$temp =~ /([^\/]+$)/;
			push(@teams, $1);
		}
	}
	foreach my $team (@teams) {
		my @projects;
		foreach my $temp (glob("$inputfolder/$locus/$team/*")) {
			if (-d $temp) {
				$temp =~ /([^\/]+$)/;
				push(@projects, $1);
			}
		}
		foreach my $project (@projects) {
			my @runs;
			foreach my $temp (glob("$inputfolder/$locus/$team/$project/*")) {
				if (-d $temp) {
					$temp =~ /([^\/]+$)/;
					push(@runs, $1);
				}
			}
			foreach my $run (@runs) {
				my @samples;
				foreach my $temp (glob("$inputfolder/$locus/$team/$project/$run/*")) {
					if (-d $temp) {
						$temp =~ /([^\/]+$)/;
						push(@samples, $1);
					}
				}
				foreach my $sample (@samples) {
					# get year, month, and worldmesh
					my $year;
					my $month;
					my $worldmesh;
					$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/sample.tsv.xz");
					{
						my $lineno = 1;
						while (<$filehandleinput1>) {
							s/\r?\n?$//;
							my @row = split(/\t/, $_);
							if ($lineno > 1 && $row[1] eq 'collection_date_utc' && $row[2]) {
								$row[2] =~ /^(\d\d\d\d)\-(\d\d)/;
								($year, $month) = ($1, $2);
								if (!exists($yearmonth{$year}{$month})) {
									$yearmonth{$year}{$month} = 1;
								}
							}
							elsif ($lineno > 1 && $row[1] eq 'worldmesh' && $row[2]) {
								$worldmesh = $row[2];
							}
							$lineno ++;
						}
					}
					close($filehandleinput1);
					$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/community_qc3nn_target.tsv.xz");
					{
						my $lineno = 1;
						my %rankname;
						my %dedup;
						while (<$filehandleinput1>) {
							s/\r?\n?$//;
							my @row = split(/\t/, $_);
							if ($lineno == 1) {
								for (my $i = 0; $i < scalar(@row); $i ++) {
									if ($row[$i] eq 'superkingdom' || $row[$i] eq 'kingdom' || $row[$i] eq 'phylum' || $row[$i] eq 'class' || $row[$i] eq 'order' || $row[$i] eq 'family' || $row[$i] eq 'genus' || $row[$i] eq 'species') {
										$rankname{$i} = $row[$i];
									}
								}
							}
							else {
								my %combinedtaxonname;
								{
									my $combinedtaxonname;
									for (my $i = 0; $i < scalar(@row); $i ++) {
										if (exists($rankname{$i}) && $row[$i] !~ /^unidentified /) {
											if ($combinedtaxonname) {
												$combinedtaxonname .= ";$row[$i]";
											}
											else {
												$combinedtaxonname = $row[$i];
											}
											$combinedtaxonname{$i} = $combinedtaxonname;
										}
									}
								}
								my $parent;
								for (my $i = 0; $i < scalar(@row); $i ++) {
									if (exists($rankname{$i}) && $row[$i] !~ /^unidentified /) {
										my $num = 1;
										while (!exists($combinedtaxonnames{$combinedtaxonname{$i}}) && exists($taxonnames{$row[$i] . '-' . $num})) {
											$num ++;
										}
										while (exists($combinedtaxonnames{$combinedtaxonname{$i}}) && exists($taxonnames{$row[$i] . '-' . $num}) && exists($parents{$row[$i] . '-' . $num}) && $parents{$row[$i] . '-' . $num} ne $parent) {
											$num ++;
										}
										my $taxonname = $row[$i] . '-' . $num;
										$combinedtaxonnames{$combinedtaxonname{$i}} = 1;
										$taxonnames{$taxonname} = 1;
										$dedup{$taxonname} = 1;
										if ($parent) {
											if (!exists($parents{$taxonname})) {
												$parents{$taxonname} = $parent;
											}
											elsif ($parents{$taxonname} ne $parent) {
												print(STDERR "Parent \"$parent\" of \"$taxonname\" is not matched to \"$parents{$taxonname}\".\n");
											}
											$havedaughters{$parent} = 1;
										}
										if (!exists($ranks{$taxonname})) {
											$ranks{$taxonname} = $rankname{$i};
										}
										elsif ($ranks{$taxonname} ne $rankname{$i}) {
											print(STDERR "Rank \"$rankname{$i}\" of \"$taxonname\" is not matched to \"$ranks{$taxonname}\".\n");
										}
										$rankname2taxa{$rankname{$i}}{$taxonname} = 1;
										$parent = $taxonname;
									}
								}
							}
							$lineno ++;
						}
						foreach my $taxonname (keys(%dedup)) {
							$nsamples{$taxonname}{$worldmesh} ++;
							$nsamples{$taxonname}{$month}{$worldmesh} ++;
							$nsamples{$taxonname}{$year}{$worldmesh} ++;
							$nsamples{$taxonname}{$year}{$month}{$worldmesh} ++;
						}
					}
					close($filehandleinput1);
				}
			}
		}
	}
}
# make taxonomy taxon
{
	my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name,slug,parent taxon`;
	foreach my $rankname ('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') {
		foreach my $taxonname (sort(keys(%{$rankname2taxa{$rankname}}))) {
			my $nonum= $taxonname;
			$nonum =~ s/-\d+$//;
			my $slug = lc($taxonname);
			$slug =~ s/[^a-z0-9\-]/-/g;
			$slug =~ s/-{2,}/-/g;
			$taxon2slug{$taxonname} = $slug;
			#print(STDERR "$rankname : $taxonname : $nonum : $slug\n");
			if ($slug =~ /[^a-z0-9\-]/) {
				print(STDERR "\"$slug\" is erroneous slug.\n");
			}
			if ($terms !~ /\d+,\"?$nonum\"?,$slug,\d+\n/) {
				if (exists($parents{$taxonname}) && $parents{$taxonname} && exists($taxonid{$parents{$taxonname}}) && $taxonid{$parents{$taxonname}}) {
					$taxonid{$taxonname} = `wp term create --url='$url/' --porcelain --slug='$slug' --parent=$taxonid{$parents{$taxonname}} taxon '$nonum'`;
				}
				elsif ($rankname eq 'superkingdom') {
					$taxonid{$taxonname} = `wp term create --url='$url/' --porcelain --slug='$slug' taxon '$nonum'`;
				}
				elsif (exists($parents{$taxonname})) {
					print(STDERR "There is no taxid of \"$parents{$taxonname}\" which is parent of \"$taxonname\".\n");
				}
				else {
					print(STDERR "There is no parent of \"$taxonname\".\n");
				}
				if (exists($taxonid{$taxonname})) {
					$taxonid{$taxonname} =~ s/\r?\n?$//;
				}
			}
			elsif ($terms =~ /(\d+),\"?$nonum\"?,$slug,(\d+)\n/) {
				$taxonid{$taxonname} = $1;
				my $parenttaxid = $2;
				if ($rankname ne 'superkingdom') {
					if (!exists($parents{$taxonname})) {
						print(STDERR "There is no parent of \"$taxonname\".\n");
					}
					elsif (!exists($taxonid{$parents{$taxonname}})) {
						print(STDERR "There is no taxid of \"$parents{$taxonname}\" which is parent of \"$taxonname\".\n");
					}
					elsif ($parenttaxid ne $taxonid{$parents{$taxonname}}) {
						print(STDERR "The parent taxid $parenttaxid of \"$taxonname\" is not matched to $taxonid{$parents{$taxonname}}.\n");
					}
				}
			}
		}
	}
}
# make taxonomy yearmonth
{
	my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name yearmonth`;
	foreach my $year (sort(keys(%yearmonth))) {
		if ($terms !~ /\d+,$year\n/) {
			$yearmonthid{$year} = `wp term create --url='$url/' --porcelain yearmonth $year`;
			$yearmonthid{$year} =~ s/\r?\n?$//;
		}
		elsif ($terms =~ /(\d+),$year\n/) {
			$yearmonthid{$year} = $1;
		}
		foreach my $month (sort(keys(%{$yearmonth{$year}}))) {
			if ($terms !~ /\d+,$year-$month\n/) {
				if (exists($yearmonthid{$year}) && $yearmonthid{$year}) {
					$yearmonthid{"$year-$month"} = `wp term create --url='$url/' --porcelain --parent=$yearmonthid{$year} yearmonth $year-$month`;
					$yearmonthid{"$year-$month"} =~ s/\r?\n?$//;
				}
			}
			elsif ($terms =~ /(\d+),$year-$month\n/) {
				$yearmonthid{"$year-$month"} = $1;
			}
		}
	}
}
{
	my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name yearmonth`;
	foreach my $year (sort(keys(%yearmonth))) {
		foreach my $month (sort(keys(%{$yearmonth{$year}}))) {
			if (!exists($yearmonthid{$month})) {
				if ($terms !~ /\d+,$month\n/) {
					$yearmonthid{$month} = `wp term create --url='$url/' --porcelain yearmonth $month`;
					$yearmonthid{$month} =~ s/\r?\n?$//;
				}
				elsif ($terms =~ /(\d+),$month\n/) {
					$yearmonthid{$month} = $1;
				}
			}
		}
	}
}
{
	my $child = 0;
	$| = 1;
	$? = 0;
	foreach my $locus (@loci) {
		my @teams;
		foreach my $temp (glob("$inputfolder/$locus/*")) {
			if (-d $temp) {
				$temp =~ /([^\/]+$)/;
				push(@teams, $1);
			}
		}
		my $locusid;
		{
			my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name project`;
			if ($terms !~ /\d+,$locus\n/) {
				$locusid = `wp term create --url='$url/' --porcelain project $locus`;
				$locusid =~ s/\r?\n?$//;
			}
			elsif ($terms =~ /(\d+),$locus\n/) {
				$locusid = $1;
			}
		}
		foreach my $team (@teams) {
			my @projects;
			foreach my $temp (glob("$inputfolder/$locus/$team/*")) {
				if (-d $temp) {
					$temp =~ /([^\/]+$)/;
					push(@projects, $1);
				}
			}
			my $teamid;
			if ($locusid) {
				my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name project`;
				if ($terms !~ /\d+,$team\n/) {
					$teamid = `wp term create --url='$url/' --parent=$locusid --porcelain project $team`;
					$teamid =~ s/\r?\n?$//;
				}
				elsif ($terms =~ /(\d+),$team\n/) {
					$teamid = $1;
				}
			}
			foreach my $project (@projects) {
				my @runs;
				foreach my $temp (glob("$inputfolder/$locus/$team/$project/*")) {
					if (-d $temp) {
						$temp =~ /([^\/]+$)/;
						push(@runs, $1);
					}
				}
				my $projectid;
				if ($teamid) {
					my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name project`;
					if ($terms !~ /\d+,$project\n/) {
						$projectid = `wp term create --url='$url/' --parent=$teamid --porcelain project $project`;
						$projectid =~ s/\r?\n?$//;
					}
					elsif ($terms =~ /(\d+),$project\n/) {
						$projectid = $1;
					}
				}
				foreach my $run (@runs) {
					my @samples;
					foreach my $temp (glob("$inputfolder/$locus/$team/$project/$run/*")) {
						if (-d $temp) {
							$temp =~ /([^\/]+$)/;
							push(@samples, $1);
						}
					}
					my $runid;
					if ($projectid) {
						my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name project`;
						if ($terms !~ /\d+,$run\n/) {
							$runid = `wp term create --url='$url/' --parent=$projectid --porcelain project $run`;
							$runid =~ s/\r?\n?$//;
						}
						elsif ($terms =~ /(\d+),$run\n/) {
							$runid = $1;
						}
					}
					foreach my $sample (@samples) {
						if (my $pid = fork()) {
							$child ++;
							if ($child == $numthreads) {
								if (wait == -1) {
									$child = 0;
								} else {
									$child --;
								}
							}
							if ($?) {
								&errorMessage(__LINE__);
							}
							next;
						}
						else {
							# correspondence, lat_lon, worldmesh, collection_date_local, collection_date_utc, temp, salinity, samp_taxon_id, filter_prod_name, size_frac_low, samp_size, num_filter
							my %sampletable;
							my $year;
							my $month;
							$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/sample.tsv.xz");
							{
								my $lineno = 1;
								while (<$filehandleinput1>) {
									s/\r?\n?$//;
									my @row = split(/\t/, $_);
									if ($lineno > 1) {
										$sampletable{$row[1]} = $row[2];
									}
									if ($lineno > 1 && $row[1] eq 'collection_date_utc' && $row[2]) {
										$row[2] =~ /^(\d\d\d\d)\-(\d\d)/;
										($year, $month) = ($1, $2);
									}
									$lineno ++;
								}
							}
							close($filehandleinput1);
							my $meshcode2id;
							if (exists($sampletable{'worldmesh'}) && $sampletable{'worldmesh'} > 0) {
								my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name meshcode2`;
								if ($terms !~ /\d+,$sampletable{'worldmesh'}\n/) {
									$meshcode2id = `wp term create --url='$url/' --porcelain meshcode2 $sampletable{'worldmesh'}`;
									$meshcode2id =~ s/\r?\n?$//;
								}
								elsif ($terms =~ /(\d+),$sampletable{'worldmesh'}\n/) {
									$meshcode2id = $1;
								}
							}
							{
								my $samples = `wp post list --url='$url/' --post_type='sample' --format=csv --fields=ID,post_title`;
								my $sampleid;
								if ($samples =~ /(\d+),$sample\n/) {
									$sampleid = $1;
								}
								# add meshcode2
								if ($sampleid && exists($sampletable{'worldmesh'})) {
									print(STDERR "Adding meshcode2 \"$sampletable{'worldmesh'}\" to sample $sampleid...\n");
									#my $terms = `wp post term list --url='$url/' --format=csv --fields=term_id,name $sampleid meshcode2`;
									#if ($terms =~ /(\d+),(\d+)\n/ && ($1 != $meshcode2id || $2 != $sampletable{'worldmesh'})) {
										system("wp post term remove --url='$url/' --all $sampleid meshcode2");
									#	$terms = '';
									#}
									if (exists($sampletable{'worldmesh'}) && $sampletable{'worldmesh'}) {
										#if ($terms !~ /$meshcode2id,$sampletable{'worldmesh'}\n/) {
											system("wp post term add --url='$url/' $sampleid meshcode2 $sampletable{'worldmesh'}");
										#}
									}
								}
								# add project
								if ($sampleid && $run && $runid) {
									print(STDERR "Adding project \"$run\" to sample $sampleid...\n");
									#my $terms = `wp post term list --url='$url/' --format=csv --fields=term_id,name $sampleid project`;
									#if ($terms =~ /(\d+),(\S+)\n/ && ($1 != $runid || $2 ne $run)) {
										system("wp post term remove --url='$url/' --all $sampleid project");
									#	$terms = '';
									#}
									#if ($terms !~ /$runid,$run\n/) {
										system("wp post term add --url='$url/' $sampleid project $run");
									#}
								}
								# add year, month, and year-month
								if ($sampleid && $month && $year) {
									print(STDERR "Adding month $month and year $year to sample $sampleid...\n");
									#my $terms = `wp post term list --url='$url/' --format=csv --fields=term_id,name $sampleid yearmonth`;
									system("wp post term remove --url='$url/' --all $sampleid yearmonth");
									#if ($terms !~ /$yearmonthid{$year},$year\n/) {
									#	system("wp post term add --url='$url/' $sampleid yearmonth $year");
									#}
									#if ($terms !~ /$yearmonthid{$month},$month\n/) {
										system("wp post term add --url='$url/' $sampleid yearmonth $month");
									#}
									my $yearmonth = "$year-$month";
									#if ($terms !~ /$yearmonthid{$yearmonth},$yearmonth\n/) {
										system("wp post term add --url='$url/' $sampleid yearmonth $yearmonth");
									#}
								}
								#$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/community_qc3nn_target.tsv.xz");
								#while (<$filehandleinput1>) {
								#	if (/([ACGT]{20,})/) {
								#		my $sequence = $1;
								#		system("wp post meta add --url='$url/' $sampleid sequence $sequence");
								#	}
								#}
								#close($filehandleinput1);
								print(STDERR "Reading taxonomy profile of sample $sampleid...\n");
								my @dedup;
								$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/community_qc3nn_target.tsv.xz");
								{
									my %dedup;
									my $lineno = 1;
									my %rankname;
									while (<$filehandleinput1>) {
										s/\r?\n?$//;
										my @row = split(/\t/, $_);
										if ($lineno == 1) {
											for (my $i = 0; $i < scalar(@row); $i ++) {
												if ($row[$i] eq 'superkingdom' || $row[$i] eq 'kingdom' || $row[$i] eq 'phylum' || $row[$i] eq 'class' || $row[$i] eq 'order' || $row[$i] eq 'family' || $row[$i] eq 'genus' || $row[$i] eq 'species') {
													$rankname{$i} = $row[$i];
												}
											}
										}
										else {
											my $line;
											for (my $i = 0; $i < scalar(@row); $i ++) {
												if (exists($rankname{$i}) && $row[$i] !~ /^unidentified /) {
													if ($line) {
														$line .= ";$row[$i]";
													}
													else {
														$line = $row[$i];
													}
												}
											}
											$dedup{$line} = 1;
										}
										$lineno ++;
									}
									my @temp = sort(keys(%dedup));
									if (scalar(@temp) > 1) {
										for (my $i = 0; $i < scalar(@temp) - 1; $i ++) {
											if ($temp[($i + 1)] !~ /^$temp[$i]/) {
												push(@dedup, $temp[$i]);
											}
										}
										push(@dedup, $temp[-1]);
									}
									elsif (scalar(@temp) == 1) {
										push(@dedup, $temp[0]);
									}
								}
								close($filehandleinput1);
								{
									my $terms = `wp post term list --url='$url/' --format=csv --fields=term_id,slug $sampleid taxon`;
									foreach (@dedup) {
										my @row = split(/;/, $_);
										my %combinedtaxonname;
										{
											my $combinedtaxonname;
											for (my $i = 0; $i < scalar(@row); $i ++) {
												if ($combinedtaxonname) {
													$combinedtaxonname .= ";$row[$i]";
												}
												else {
													$combinedtaxonname = $row[$i];
												}
												$combinedtaxonname{$i} = $combinedtaxonname;
											}
										}
										my $parent;
										for (my $i = 0; $i < scalar(@row); $i ++) {
											my $num = 1;
											while (!exists($combinedtaxonnames{$combinedtaxonname{$i}}) && exists($taxonnames{$row[$i] . '-' . $num})) {
												$num ++;
											}
											while (exists($combinedtaxonnames{$combinedtaxonname{$i}}) && exists($taxonnames{$row[$i] . '-' . $num}) && exists($parents{$row[$i] . '-' . $num}) && $parents{$row[$i] . '-' . $num} ne $parent) {
												$num ++;
											}
											$parent = $row[$i] . '-' . $num;
										}
										if ($parent) {
											if ($terms !~ /$taxonid{$parent},$taxon2slug{$parent}\n/) {
												print(STDERR "Adding taxonomy \"$parent\" to sample $sampleid...\n");
												system("wp post term add --url='$url/' $sampleid taxon $taxon2slug{$parent}");
											}
										}
									}
								}
							}
							exit;
						}
					}
				}
			}
		}
	}
}
# join
while (wait != -1) {
	if ($?) {
		&errorMessage(__LINE__, 'Cannot run BLAST search correctly.');
	}
}


sub readFile {
	my $filehandle;
	my $filename = shift(@_);
	if ($filename =~ /\.gz$/i) {
		unless (open($filehandle, "pigz -p 8 -dc $filename 2> $devnull |")) {
			die(__LINE__ . ": Cannot open \"$filename\".\n");
		}
	}
	elsif ($filename =~ /\.bz2$/i) {
		unless (open($filehandle, "lbzip2 -n 8 -dc $filename 2> $devnull |")) {
			die(__LINE__ . ": Cannot open \"$filename\".\n");
		}
	}
	elsif ($filename =~ /\.xz$/i) {
		unless (open($filehandle, "xz -T 8 -dc $filename 2> $devnull |")) {
			die(__LINE__ . ": Cannot open \"$filename\".\n");
		}
	}
	else {
		unless (open($filehandle, "< $filename")) {
			die(__LINE__ . ": Cannot open \"$filename\".\n");
		}
	}
	return($filehandle);
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
