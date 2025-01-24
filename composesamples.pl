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
	my $terms = `wp term list --format=csv --fields=term_id,name,slug,parent taxon`;
	foreach my $rankname ('superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') {
		foreach my $taxonname (sort(keys(%{$rankname2taxa{$rankname}}))) {
			my $nonum = $taxonname;
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
					$taxonid{$taxonname} = `wp term create --porcelain --slug='$slug' --parent=$taxonid{$parents{$taxonname}} taxon '$nonum'`;
				}
				elsif ($rankname eq 'superkingdom') {
					$taxonid{$taxonname} = `wp term create --porcelain --slug='$slug' taxon '$nonum'`;
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
	my $terms = `wp term list --format=csv --fields=term_id,name yearmonth`;
	foreach my $year (sort(keys(%yearmonth))) {
		if ($terms !~ /\d+,$year\n/) {
			$yearmonthid{$year} = `wp term create --porcelain yearmonth $year`;
			$yearmonthid{$year} =~ s/\r?\n?$//;
		}
		elsif ($terms =~ /(\d+),$year\n/) {
			$yearmonthid{$year} = $1;
		}
		foreach my $month (sort(keys(%{$yearmonth{$year}}))) {
			if ($terms !~ /\d+,$year-$month\n/) {
				if (exists($yearmonthid{$year}) && $yearmonthid{$year}) {
					$yearmonthid{"$year-$month"} = `wp term create --porcelain --parent=$yearmonthid{$year} yearmonth $year-$month`;
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
	my $terms = `wp term list --format=csv --fields=term_id,name yearmonth`;
	foreach my $year (sort(keys(%yearmonth))) {
		foreach my $month (sort(keys(%{$yearmonth{$year}}))) {
			if (!exists($yearmonthid{$month})) {
				if ($terms !~ /\d+,$month\n/) {
					$yearmonthid{$month} = `wp term create --porcelain yearmonth $month`;
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
	my $terms = `wp term list --format=csv --fields=term_id,name meshcode2`;
	foreach my $worldmesh (sort(keys(%worldmesh))) {
		if ($terms !~ /\d+,$worldmesh\n/) {
			$meshcode2id{$worldmesh} = `wp term create --porcelain meshcode2 $worldmesh`;
			$meshcode2id{$worldmesh} =~ s/\r?\n?$//;
		}
		elsif ($terms =~ /(\d+),$worldmesh\n/) {
			$meshcode2id{$worldmesh} = $1;
		}
	}
}
# make draft pages
{
	my $child = 0;
	$| = 1;
	$? = 0;
	my $samples = `wp post list --post_type='sample' --format=csv --fields=ID,post_title`;
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
			my $terms = `wp term list --format=csv --fields=term_id,name project`;
			if ($terms !~ /\d+,$locus\n/) {
				$locusid = `wp term create --porcelain project $locus`;
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
				my $terms = `wp term list --format=csv --fields=term_id,name project`;
				if ($terms !~ /\d+,$team\n/) {
					$teamid = `wp term create --parent=$locusid --porcelain project $team`;
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
					my $terms = `wp term list --format=csv --fields=term_id,name project`;
					if ($terms !~ /\d+,$project\n/) {
						$projectid = `wp term create --parent=$teamid --porcelain project $project`;
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
						my $terms = `wp term list --format=csv --fields=term_id,name project`;
						if ($terms !~ /\d+,$run\n/) {
							$runid = `wp term create --parent=$projectid --porcelain project $run`;
							$runid =~ s/\r?\n?$//;
						}
						elsif ($terms =~ /(\d+),$run\n/) {
							$runid = $1;
						}
					}
					foreach my $sample (@samples) {
						# make sample draft page
						if ($samples !~ /\d+,\"?$sample\"?\n/) {
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
								system("wp post create --post_author=1 --post_type='sample' --comment_status='closed' --ping_status='closed' --post_status='draft' --post_title='$sample' --post_excerpt='' --post_content='temp'");
								exit;
							}
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
		&errorMessage(__LINE__, 'Cannot wp run correctly.');
	}
}
# update and publish sample pages
{
	my $child = 0;
	$| = 1;
	$? = 0;
	my $samples = `wp post list --post_type='sample' --format=csv --fields=ID,post_title`;
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
			my $terms = `wp term list --format=csv --fields=term_id,name project`;
			if ($terms !~ /\d+,$locus\n/) {
				$locusid = `wp term create --porcelain project $locus`;
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
				my $terms = `wp term list --format=csv --fields=term_id,name project`;
				if ($terms !~ /\d+,$team\n/) {
					$teamid = `wp term create --parent=$locusid --porcelain project $team`;
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
					my $terms = `wp term list --format=csv --fields=term_id,name project`;
					if ($terms !~ /\d+,$project\n/) {
						$projectid = `wp term create --parent=$teamid --porcelain project $project`;
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
						my $terms = `wp term list --format=csv --fields=term_id,name project`;
						if ($terms !~ /\d+,$run\n/) {
							$runid = `wp term create --parent=$projectid --porcelain project $run`;
							$runid =~ s/\r?\n?$//;
						}
						elsif ($terms =~ /(\d+),$run\n/) {
							$runid = $1;
						}
					}
					foreach my $sample (@samples) {
						if ($samples =~ /(\d+),\"?$sample\"?\n/) {
							my $sampleid = $1;
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
								# target_gene, pcr_primers, pcr_standard, pcr_standard_conc
								my %experimenttable;
								$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/experiment.tsv.xz");
								{
									my $lineno = 1;
									my @label;
									while (<$filehandleinput1>) {
										s/\r?\n?$//;
										my @row = split(/\t/, $_);
										if ($lineno > 1) {
											$experimenttable{$row[1]} = $row[2];
										}
										$lineno ++;
									}
								}
								close($filehandleinput1);
								my $wordcloudphylum;
								my $wordcloudclass;
								my $wordcloudorder;
								my $wordcloudfamily;
								my $wordcloudgenus;
								my $wordcloudspecies;
								if (-e "$inputfolder/$locus/$team/$project/$run/$sample/wordcloud.phylum.png") {
									$wordcloudphylum = "/dist/$locus/$team/$project/$run/$sample/wordcloud.phylum.png";
								}
								if (-e "$inputfolder/$locus/$team/$project/$run/$sample/wordcloud.class.png") {
									$wordcloudclass = "/dist/$locus/$team/$project/$run/$sample/wordcloud.class.png";
								}
								if (-e "$inputfolder/$locus/$team/$project/$run/$sample/wordcloud.order.png") {
									$wordcloudorder = "/dist/$locus/$team/$project/$run/$sample/wordcloud.order.png";
								}
								if (-e "$inputfolder/$locus/$team/$project/$run/$sample/wordcloud.family.png") {
									$wordcloudfamily = "/dist/$locus/$team/$project/$run/$sample/wordcloud.family.png";
								}
								if (-e "$inputfolder/$locus/$team/$project/$run/$sample/wordcloud.genus.png") {
									$wordcloudgenus = "/dist/$locus/$team/$project/$run/$sample/wordcloud.genus.png";
								}
								if (-e "$inputfolder/$locus/$team/$project/$run/$sample/wordcloud.species.png") {
									$wordcloudspecies = "/dist/$locus/$team/$project/$run/$sample/wordcloud.species.png";
								}
								open($filehandleoutput1, "> $inputfolder/$locus/$team/$project/$run/$sample/wordpress.html") or die(__LINE__ . ": Cannot open \"$inputfolder/$locus/$team/$project/$run/$sample/wordpress.html\"\n");
								my $teampath = lc("/project/$locus/$team/");
								my $projectpath = lc("/project/$locus/$team/$project/");
								my $runpath = lc("/project/$locus/$team/$project/$run/");
								my $locuspath = lc("/project/$locus/");
								print($filehandleoutput1 <<"_END");
<h2>Data file download URL</h2>
<ul>
<li><a href=\"/dist/$locus/$team/$project/$run/$sample/\"><?php echo esc_url( home_url() ); ?>/dist/$locus/$team/$project/$run/$sample/</a></li>
</ul>
<h2>Sample information</h2>
<ul>
<li>Sample name: $sample</li>
<li>Ambrella project or team: <a href=\"$teampath\">$team</a></li>
<li>Project: <a href=\"$projectpath\">$project</a></li>
<li>Sequencing run: <a href=\"$runpath\">$run</a></li>
<li>Correnpondence: $sampletable{'correspondence'}</li>
<li>Collection date (UTC): $sampletable{'collection_date_utc'}</li>
<li>Collection date (Local): $sampletable{'collection_date_local'}</li>
<li>Worldmesh: <a href=\"$locuspath?meshcode2=$sampletable{'worldmesh'}/\">$sampletable{'worldmesh'}</a></li>
<li>Latitude and longitude (Center of mesh): $sampletable{'lat_lon'}</li>
<li>Taxonomy: $sampletable{'samp_taxon_id'}</li>
<li>Temperature: $sampletable{'temp'} degree Celsius</li>
<li>Salinity: $sampletable{'salinity'} percent</li>
<li>Filter: $sampletable{'filter_prod_name'}</li>
<li>Filter pore size: $sampletable{'size_frac_low'} micrometer</li>
<li>Filtered water volume: $sampletable{'samp_size'} milliliter</li>
<li>Number of filters used: $sampletable{'num_filter'}</li>
</ul>
<h2>Experiment information</h2>
<ul>
<li>Target locus: <a href=\"$locuspath\">$locus</a> ($experimenttable{'target_gene'})</li>
_END
								print($filehandleoutput1 "<li>PCR primer sequences: <ul>\n");
								foreach my $pcr_primer (split(/;/, $experimenttable{'pcr_primers'})) {
									print($filehandleoutput1 "<li>$pcr_primer</li>\n");
								}
								print($filehandleoutput1 "</ul></li>\n");
								print($filehandleoutput1 "<li>PCR internal standard sequences: <ul>\n");
								foreach my $pcr_standard (split(/;/, $experimenttable{'pcr_standard'})) {
									print($filehandleoutput1 "<li>$pcr_standard</li>\n");
								}
								print($filehandleoutput1 "</ul></li>\n");
								print($filehandleoutput1 "<li>PCR internal standard concentration: <ul>\n");
								foreach my $pcr_standard_conc (split(/;/, $experimenttable{'pcr_standard_conc'})) {
									print($filehandleoutput1 "<li>$pcr_standard_conc</li>\n");
								}
								print($filehandleoutput1 "</ul></li>\n");
								print($filehandleoutput1 "</ul>\n");
								if (exists($sampletable{'lat_lon'}) &&exists($sampletable{'worldmesh'}) && $sampletable{'lat_lon'} && $sampletable{'worldmesh'}) {
									$sampletable{'lat_lon'} =~ /^(\S+) +(\S+)$/;
									my ($latitude, $longitude) = ($1, $2);
									my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode2latlon($sampletable{'worldmesh'});
									if ($latitude && $longitude && $latN && $latS && $lonW && $lonE) {
										my $meshcode1 = substr($sampletable{'worldmesh'}, 0, 6);
										print($filehandleoutput1 "<h2>Geographic location</h2>[leaflet-map lat=\"$latitude\" lng=\"$longitude\" zoom=\"10\"][zoomhomemap][fullscreen]");
										open($filehandleinput1, "< $locus$meshcode1.html") or die(__LINE__ . ": Cannot open \"$locus$meshcode1.html\".\n");
										while (<$filehandleinput1>) {
											if (/(\[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" [^\[\]]+\][^\[\]]+\[\/leaflet-polygon\])/) {
												print($filehandleoutput1 $1);
											}
										}
										close($filehandleinput1);
									}
								}
								if ($wordcloudspecies) {
									print($filehandleoutput1 "<h2>Species composition</h2><p style=\"text-align: center;\"><a href=\"$wordcloudspecies\"><img src=\"$wordcloudspecies\" alt=\"Wordcloud image of species composition of $sample\" /></a></p>");
								}
								if ($wordcloudgenus) {
									print($filehandleoutput1 "<h2>Genus composition</h2><p style=\"text-align: center;\"><a href=\"$wordcloudgenus\"><img src=\"$wordcloudgenus\" alt=\"Wordcloud image of genus composition of $sample\" /></a></p>");
								}
								if ($wordcloudfamily) {
									print($filehandleoutput1 "<h2>Family composition</h2><p style=\"text-align: center;\"><a href=\"$wordcloudfamily\"><img src=\"$wordcloudfamily\" alt=\"Wordcloud image of family composition of $sample\" /></a></p>");
								}
								if ($wordcloudorder) {
									print($filehandleoutput1 "<h2>Order composition</h2><p style=\"text-align: center;\"><a href=\"$wordcloudorder\"><img src=\"$wordcloudorder\" alt=\"Wordcloud image of order composition of $sample\" /></a></p>");
								}
								if ($wordcloudclass) {
									print($filehandleoutput1 "<h2>Class composition</h2><p style=\"text-align: center;\"><a href=\"$wordcloudclass\"><img src=\"$wordcloudclass\" alt=\"Wordcloud image of class composition of $sample\" /></a></p>");
								}
								if ($wordcloudphylum) {
									print($filehandleoutput1 "<h2>Phylum composition</h2><p style=\"text-align: center;\"><a href=\"$wordcloudphylum\"><img src=\"$wordcloudphylum\" alt=\"Wordcloud image of phylum composition of $sample\" /></a></p>");
								}
								print($filehandleoutput1 "<h2>Complete community composition</h2><table>");
								$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/community_qc3nn_target.tsv.xz");
								{
									my $lineno = 1;
									my %rankname;
									my %switch;
									while (<$filehandleinput1>) {
										s/\r?\n?$//;
										my @row = split(/\t/, $_);
										if ($lineno == 1) {
											print($filehandleoutput1 "<thead>\n<tr>\n");
											for (my $i = 0; $i < scalar(@row); $i ++) {
												if ($row[$i] eq 'nreads' || $row[$i] eq 'ncopies' || $row[$i] eq 'ncopiesperml') {
													print($filehandleoutput1 "<th style=\"text-align: right;\">$row[$i]</th>\n");
													$switch{$i} = 1;
												}
												elsif ($row[$i] eq 'superkingdom' || $row[$i] eq 'kingdom' || $row[$i] eq 'phylum' || $row[$i] eq 'class' || $row[$i] eq 'order' || $row[$i] eq 'family' || $row[$i] eq 'genus' || $row[$i] eq 'species') {
													print($filehandleoutput1 "<th style=\"text-align: left;\">$row[$i]</th>\n");
													$rankname{$i} = $row[$i];
												}
												elsif ($row[$i] eq 'sequence') {
													print($filehandleoutput1 "<th style=\"text-align: left;\">$row[$i]</th>\n");
													$switch{$i} = 1;
												}
											}
											print($filehandleoutput1 "</tr>\n</thead>\n<tbody>\n");
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
											print($filehandleoutput1 "<tr>\n");
											my $parent;
											for (my $i = 0; $i < scalar(@row); $i ++) {
												if (exists($switch{$i}) && ($row[$i] =~ /^\d+$/ || $row[$i] =~ /^\d+\.\d+$/)) {
													print($filehandleoutput1 "<td style=\"white-space: nowrap; text-align: right;\">$row[$i]</td>\n");
												}
												elsif (exists($switch{$i})) {
													print($filehandleoutput1 "<td style=\"white-space: nowrap; text-align: left;\">$row[$i]</td>\n");
												}
												elsif (exists($rankname{$i}) && $row[$i] =~ /^unidentified /) {
													print($filehandleoutput1 "<td style=\"white-space: nowrap; text-align: left;\">$row[$i]</td>\n");
												}
												elsif (exists($rankname{$i})) {
													my $num = 1;
													while (!exists($combinedtaxonnames{$combinedtaxonname{$i}}) && exists($taxonnames{$row[$i] . '-' . $num})) {
														$num ++;
													}
													while (exists($combinedtaxonnames{$combinedtaxonname{$i}}) && exists($taxonnames{$row[$i] . '-' . $num}) && exists($parents{$row[$i] . '-' . $num}) && $parents{$row[$i] . '-' . $num} ne $parent) {
														$num ++;
													}
													my $taxonname = $row[$i] . '-' . $num;
													my $slug = lc($taxonname);
													$slug =~ s/[^a-z0-9\-]/-/g;
													$slug =~ s/-{2,}/-/g;
													print($filehandleoutput1 "<td style=\"white-space: nowrap; text-align: left;\"><a href=\"/map/$slug/\">$row[$i]</a></td>\n");
													$parent = $taxonname;
												}
											}
											print($filehandleoutput1 "</tr>\n");
										}
										$lineno ++;
									}
								}
								print($filehandleoutput1 "</tbody>\n");
								close($filehandleinput1);
								print($filehandleoutput1 "</table>");
								close($filehandleoutput1);
								# make or update sample page
								system("wp post update --post_author=1 --post_type='sample' --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$sample' --post_date_gmt='$sampletable{'collection_date_utc'}' --post_date='$sampletable{'collection_date_local'}' --post_excerpt='$locus metabarcoding sample by $sampletable{'correspondence'}' $sampleid $inputfolder/$locus/$team/$project/$run/$sample/wordpress.html");
								unlink("$inputfolder/$locus/$team/$project/$run/$sample/wordpress.html");
								exit;
							}
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
