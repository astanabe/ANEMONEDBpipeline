#use strict;
use File::Spec;
my $devnull = File::Spec->devnull();

my $numthreads = 16;

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
my %worldmesh;
my %meshcode2id;
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
								$worldmesh{$worldmesh} = 1;
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
# make taxonomy meshcode2
{
	my $terms = `wp term list --url='$url/' --format=csv --fields=term_id,name meshcode2`;
	foreach my $worldmesh (sort(keys(%worldmesh))) {
		if ($terms !~ /\d+,$worldmesh\n/) {
			$meshcode2id{$worldmesh} = `wp term create --url='$url/' --porcelain meshcode2 $worldmesh`;
			$meshcode2id{$worldmesh} =~ s/\r?\n?$//;
		}
		elsif ($terms =~ /(\d+),$worldmesh\n/) {
			$meshcode2id{$worldmesh} = $1;
		}
	}
}
# Make draft pages
{
	my $child = 0;
	$| = 1;
	$? = 0;
	my $maps = `wp post list --url='$url/' --post_type='map' --format=csv --fields=ID,post_name`;
	foreach my $rankname ('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') {
		foreach my $taxonname (sort(keys(%{$rankname2taxa{$rankname}}))) {
			my $nonum= $taxonname;
			$nonum =~ s/-\d+$//;
			my $slug = lc($taxonname);
			$slug =~ s/[^a-z0-9\-]/-/g;
			$slug =~ s/-{2,}/-/g;
			if ($maps !~ /\d+,\"?$slug\"?\n/) {
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
					system("wp post create --url='$url/' --post_author=1 --post_type='map' --comment_status='closed' --ping_status='closed' --post_status='draft' --post_title='$nonum' --post_name='$slug' --post_parent='' --post_excerpt='' --post_content='temp'");
					exit;
				}
			}
		}
	}
}
# join
while (wait != -1) {
	if ($?) {
		&errorMessage(__LINE__, 'Cannot wp run correctly.');
	}
}
# Make maps for each taxon
my %mapid;
{
	my $child = 0;
	$| = 1;
	$? = 0;
	my $maps = `wp post list --url='$url/' --post_type='map' --format=csv --fields=ID,post_name`;
	foreach my $rankname ('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') {
		foreach my $taxonname (sort(keys(%{$rankname2taxa{$rankname}}))) {
			my $nonum= $taxonname;
			$nonum =~ s/-\d+$//;
			my $slug = lc($taxonname);
			$slug =~ s/[^a-z0-9\-]/-/g;
			$slug =~ s/-{2,}/-/g;
			my $nsamples = 0;
			open($filehandleoutput1, "> $slug.html") or die(__LINE__ . ": Cannot open \"$slug.html\"\n");
			print($filehandleoutput1 "<h2>Browse all metabarcoding samples including $nonum on the earth</h2>[leaflet-map fitbounds][zoomhomemap][fullscreen]");
			foreach my $tempvalue (sort({$a <=> $b} keys(%{$nsamples{$taxonname}}))) {
				# case of worldmesh
				if (length($tempvalue) == 8) {
					my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode2latlon($tempvalue);
					if ($nsamples{$taxonname}{$tempvalue} > 160) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FC8E50\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					elsif ($nsamples{$taxonname}{$tempvalue} > 80) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF7D93\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					elsif ($nsamples{$taxonname}{$tempvalue} > 40) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF72C7\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					elsif ($nsamples{$taxonname}{$tempvalue} > 20) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF74F1\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					elsif ($nsamples{$taxonname}{$tempvalue} > 10) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#DF86FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					elsif ($nsamples{$taxonname}{$tempvalue} > 5) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#9B9FFF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					elsif ($nsamples{$taxonname}{$tempvalue} > 1) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					elsif ($nsamples{$taxonname}{$tempvalue} == 1) {
						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?meshcode2=$tempvalue\">$nsamples{$taxonname}{$tempvalue} sample is available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					}
					$nsamples += $nsamples{$taxonname}{$tempvalue};
				}
			}
			if (exists($havedaughters{$taxonname}) && $havedaughters{$taxonname} || $rankname eq 'species') {
				print($filehandleoutput1 "<h2>List of submaps</h2>[subpages post_type=\"map\" depth=\"1\" post_status=\"publish\" sort_column=\"post_title\" sort_order=\"asc\"]");
			}
			close($filehandleoutput1);
			if ($maps =~ /(\d+),\"?$slug\"?\n/) {
				$mapid{$taxonname} = $1;
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
					system("wp post update --url='$url/' --post_type='map' --post_author=1 --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$nonum ($nsamples)' --post_name='$slug' --post_parent='$mapid{$parents{$taxonname}}' --post_excerpt='' $mapid{$taxonname} $slug.html");
					unlink("$slug.html");
					exit;
				}
			}
		}
	}
}
# join
while (wait != -1) {
	if ($?) {
		&errorMessage(__LINE__, 'Cannot wp run correctly.');
	}
}
# Make draft pages
{
	my $child = 0;
	$| = 1;
	$? = 0;
	my $maps = `wp post list --url='$url/' --post_type='map' --format=csv --fields=ID,post_name`;
	#foreach my $rankname ('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') {
	foreach my $rankname ('genus', 'species') {
		foreach my $taxonname (sort(keys(%{$rankname2taxa{$rankname}}))) {
			my $nonum= $taxonname;
			$nonum =~ s/-\d+$//;
			my $slug = lc($taxonname);
			$slug =~ s/[^a-z0-9\-]/-/g;
			$slug =~ s/-{2,}/-/g;
			foreach my $tempvalue (sort({$a <=> $b} keys(%{$nsamples{$taxonname}}))) {
				# case of month or year
				#if (length($tempvalue) == 2 || length($tempvalue) == 4) {
				if (length($tempvalue) == 4) {
					if ($maps !~ /\d+,\"?$tempvalue-$slug\"?\n/) {
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
							system("wp post create --url='$url/' --post_author=1 --post_type='map' --comment_status='closed' --ping_status='closed' --post_status='draft' --post_title='$tempvalue - $nonum' --post_name='$tempvalue-$slug' --post_parent='' --post_excerpt='' --post_content='temp'");
							exit;
						}
					}
					#if (length($tempvalue) == 4) {
					#	foreach my $month (sort(keys(%{$nsamples{$taxonname}{$tempvalue}}))) {
					#		if (length($month) == 2) {
					#			if ($maps !~ /\d+,\"?$tempvalue-$month-$slug\"?\n/) {
					#				if (my $pid = fork()) {
					#					$child ++;
					#					if ($child == $numthreads) {
					#						if (wait == -1) {
					#							$child = 0;
					#						} else {
					#							$child --;
					#						}
					#					}
					#					if ($?) {
					#						&errorMessage(__LINE__);
					#					}
					#					next;
					#				}
					#				else {
					#					system("wp post create --url='$url/' --post_author=1 --post_type='map' --comment_status='closed' --ping_status='closed' --post_status='draft' --post_title='$tempvalue - $month - $nonum' --post_name='$tempvalue-$month-$slug' --post_parent='' --post_excerpt='' --post_content='temp'");
					#					exit;
					#				}
					#			}
					#		}
					#	}
					#}
				}
			}
		}
	}
}
# join
while (wait != -1) {
	if ($?) {
		&errorMessage(__LINE__, 'Cannot wp run correctly.');
	}
}
# Make maps for each year, month or year-month
{
	my $child = 0;
	$| = 1;
	$? = 0;
	my $maps = `wp post list --url='$url/' --post_type='map' --format=csv --fields=ID,post_name`;
	#foreach my $rankname ('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') {
	foreach my $rankname ('genus', 'species') {
		foreach my $taxonname (sort(keys(%{$rankname2taxa{$rankname}}))) {
			my $nonum= $taxonname;
			$nonum =~ s/-\d+$//;
			my $slug = lc($taxonname);
			$slug =~ s/[^a-z0-9\-]/-/g;
			$slug =~ s/-{2,}/-/g;
			foreach my $tempvalue (sort({$a <=> $b} keys(%{$nsamples{$taxonname}}))) {
				# case of month or year
				#if (length($tempvalue) == 2 || length($tempvalue) == 4) {
				if (length($tempvalue) == 4) {
					open($filehandleoutput1, "> $tempvalue-$slug.html") or die(__LINE__ . ": Cannot open \"$slug-$tempvalue.html\"\n");
					#if (length($tempvalue) == 2) {
					#	print($filehandleoutput1 "<h2>Browse all metabarcoding samples of $tempvalue-th month across all years including $nonum on the earth</h2>[leaflet-map fitbounds][zoomhomemap][fullscreen]");
					#}
					#elsif (length($tempvalue) == 4) {
						print($filehandleoutput1 "<h2>Browse all metabarcoding samples of $tempvalue including $nonum on the earth</h2>[leaflet-map fitbounds][zoomhomemap][fullscreen]");
					#}
					my $tempmapid;
					{
						my $nsamples = 0;
						foreach my $worldmesh (sort(keys(%{$nsamples{$taxonname}{$tempvalue}}))) {
							if (length($worldmesh) == 8) {
								my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode2latlon($worldmesh);
								if ($nsamples{$taxonname}{$tempvalue}{$worldmesh} > 160) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FC8E50\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								elsif ($nsamples{$taxonname}{$tempvalue}{$worldmesh} > 80) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF7D93\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								elsif ($nsamples{$taxonname}{$tempvalue}{$worldmesh} > 40) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF72C7\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								elsif ($nsamples{$taxonname}{$tempvalue}{$worldmesh} > 20) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF74F1\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								elsif ($nsamples{$taxonname}{$tempvalue}{$worldmesh} > 10) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#DF86FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								elsif ($nsamples{$taxonname}{$tempvalue}{$worldmesh} > 5) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#9B9FFF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								elsif ($nsamples{$taxonname}{$tempvalue}{$worldmesh} > 1) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								elsif ($nsamples{$taxonname}{$tempvalue}{$worldmesh} == 1) {
									print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$worldmesh} sample is available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
								}
								$nsamples += $nsamples{$taxonname}{$tempvalue}{$worldmesh};
							}
						}
						if (length($tempvalue) == 2 || length($tempvalue) == 4) {
							print($filehandleoutput1 "<h2>List of sibling maps</h2>[siblings post_type=\"map\" depth=\"1\" post_status=\"publish\" sort_column=\"post_title\" sort_order=\"asc\"]");
						}
						#if (length($tempvalue) == 4) {
						#	print($filehandleoutput1 "<h2>List of submaps</h2>[subpages post_type=\"map\" depth=\"1\" post_status=\"publish\" sort_column=\"post_title\" sort_order=\"asc\"]");
						#}
						close($filehandleoutput1);
						if ($maps =~ /(\d+),\"?$tempvalue-$slug\"?\n/) {
							$tempmapid = $1;
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
								system("wp post update --url='$url/' --post_type='map' --post_author=1 --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$tempvalue - $nonum ($nsamples)' --post_name='$tempvalue-$slug' --post_parent='$mapid{$taxonname}' --post_excerpt='' $tempmapid $tempvalue-$slug.html");
								unlink("$tempvalue-$slug.html");
								exit;
							}
						}
					}
					#if (length($tempvalue) == 4) {
					#	foreach my $month (sort(keys(%{$nsamples{$taxonname}{$tempvalue}}))) {
					#		if (length($month) == 2) {
					#			my $nsamples = 0;
					#			open($filehandleoutput1, "> $tempvalue-$month-$slug.html") or die(__LINE__ . ": Cannot open \"$tempvalue-$month-$slug.html\"\n");
					#			print($filehandleoutput1 "<h2>Browse all metabarcoding samples of $month-th month, $tempvalue including $nonum on the earth</h2>[leaflet-map fitbounds][zoomhomemap][fullscreen]");
					#			foreach my $worldmesh (sort(keys(%{$nsamples{$taxonname}{$tempvalue}{$month}}))) {
					#				if (length($worldmesh) == 8) {
					#					my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode2latlon($worldmesh);
					#					if ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} > 160) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FC8E50\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					elsif ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} > 80) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF7D93\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					elsif ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} > 40) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF72C7\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					elsif ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} > 20) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF74F1\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					elsif ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} > 10) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#DF86FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					elsif ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} > 5) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#9B9FFF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					elsif ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} > 1) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					elsif ($nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} == 1) {
					#						print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.3\"]<a href=\"/taxon/$slug/?yearmonth=$tempvalue-$month\&meshcode2=$worldmesh\">$nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh} sample is available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
					#					}
					#					$nsamples += $nsamples{$taxonname}{$tempvalue}{$month}{$worldmesh};
					#				}
					#			}
					#			print($filehandleoutput1 "<h2>List of sibling maps</h2>[siblings post_type=\"map\" depth=\"1\" post_status=\"publish\" sort_column=\"post_title\" sort_order=\"asc\"]");
					#			close($filehandleoutput1);
					#			if ($maps =~ /(\d+),\"?$tempvalue-$month-$slug\"?\n/) {
					#				my $tempmapid2 = $1;
					#				if (my $pid = fork()) {
					#					$child ++;
					#					if ($child == $numthreads) {
					#						if (wait == -1) {
					#							$child = 0;
					#						} else {
					#							$child --;
					#						}
					#					}
					#					if ($?) {
					#						&errorMessage(__LINE__);
					#					}
					#					next;
					#				}
					#				else {
					#					system("wp post update --url='$url/' --post_type='map' --post_author=1 --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$tempvalue - $month - $nonum ($nsamples)' --post_name='$tempvalue-$month-$slug' --post_parent='$tempmapid' --post_excerpt='' $tempmapid2 $tempvalue-$month-$slug.html");
					#					unlink("$tempvalue-$month-$slug.html");
					#					exit;
					#				}
					#			}
					#		}
					#	}
					#}
				}
			}
		}
	}
}
# join
while (wait != -1) {
	if ($?) {
		&errorMessage(__LINE__, 'Cannot wp run correctly.');
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
