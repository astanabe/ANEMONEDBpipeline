#use strict;
use File::Spec;
my $devnull = File::Spec->devnull();

my $filehandleinput1;
my $filehandleoutput1;

my $inputfolder = $ARGV[0];
if (!-d $inputfolder) {
	die(__LINE__);
}
my @loci;
foreach my $temp (glob("$inputfolder/*")) {
	if (-d $temp) {
		$temp =~ /([^\/]+$)/;
		push(@loci, $1);
	}
}
foreach my $locus (@loci) {
	my %meshcode2;
	my %meshcode1;
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
					# correspondence, lat_lon, worldmesh, collection_date_local, collection_date_utc, temp, salinity, samp_taxon_id, filter_prod_name, size_frac_low, samp_size, num_filter
					my %sampletable;
					$filehandleinput1 = &readFile("$inputfolder/$locus/$team/$project/$run/$sample/sample.tsv.xz");
					{
						my $lineno = 1;
						my @label;
						while (<$filehandleinput1>) {
							s/\r?\n?$//;
							my @row = split(/\t/, $_);
							if ($lineno > 1) {
								$sampletable{$row[1]} = $row[2];
							}
							$lineno ++;
						}
					}
					close($filehandleinput1);
					if (exists($sampletable{'worldmesh'}) && $sampletable{'worldmesh'} > 0) {
						$meshcode2{$sampletable{'worldmesh'}} += 1;
					}
				}
			}
		}
	}
	my $locusslug = lc($locus);
	$locusslug =~ s/[^a-z0-9\-]/-/g;
	$locusslug =~ s/-{2,}/-/g;
	{
		foreach my $meshcode2 (keys(%meshcode2)) {
			my $meshcode1 = substr($meshcode2, 0, 6);
			my $latitude;
			my $longitude;
			{
				my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode1latlon($meshcode1);
				($latitude, $longitude) = ($latC, $lonC);
			}
			if (exists($meshcode2{$meshcode2})) {
				$meshcode1{$meshcode1} += $meshcode2{$meshcode2};
			}
		}
		foreach my $meshcode1 (keys(%meshcode1)) {
			my $latitude;
			my $longitude;
			{
				my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode1latlon($meshcode1);
				($latitude, $longitude) = ($latC, $lonC);
			}
			open($filehandleoutput1, "> $locus$meshcode1.html") or die(__LINE__ . ": Cannot open \"$meshcode1.html\"\n");
			print($filehandleoutput1 "<h2>Browse all $locus metabarcoding samples of worldmesh $meshcode1</h2>[leaflet-map lat=\"$latitude\" lng=\"$longitude\" zoom=\"10\"][zoomhomemap][fullscreen]");
			close($filehandleoutput1);
		}
		foreach my $meshcode2 (keys(%meshcode2)) {
			my $meshcode1 = substr($meshcode2, 0, 6);
			open($filehandleoutput1, ">> $locus$meshcode1.html") or die(__LINE__ . ": Cannot open \"$meshcode1.html\"\n");
			if (exists($meshcode2{$meshcode2})) {
				my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode2latlon($meshcode2);
				if ($meshcode2{$meshcode2} > 160) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FC8E50\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode2{$meshcode2} > 80) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF7D93\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode2{$meshcode2} > 40) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF72C7\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode2{$meshcode2} > 20) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF74F1\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode2{$meshcode2} > 10) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#DF86FF\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode2{$meshcode2} > 5) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#9B9FFF\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode2{$meshcode2} > 1) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} samples are available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode2{$meshcode2} == 1) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.5\"]<a href=\"/project/$locusslug/?meshcode2=$meshcode2\">$meshcode2{$meshcode2} sample is available on this worldmesh. Link to sample list of this worldmesh.</a>[/leaflet-polygon]");
				}
			}
			close($filehandleoutput1);
		}
	}
	{
		my $nsamples = 0;
		open($filehandleoutput1, "> $locus.html") or die(__LINE__ . ": Cannot open \"$locus.html\"\n");
		print($filehandleoutput1 "<h2>Browse all $locus metabarcoding samples on the earth</h2>[leaflet-map fitbounds][zoomhomemap][fullscreen]");
		foreach my $meshcode1 (keys(%meshcode1)) {
			if (exists($meshcode1{$meshcode1})) {
				my ($latN, $latC, $latS, $lonW, $lonC, $lonE) = &meshcode1latlon($meshcode1);
				if ($meshcode1{$meshcode1} > 160) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FC8E50\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} samples are available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode1{$meshcode1} > 80) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF7D93\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} samples are available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode1{$meshcode1} > 40) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF72C7\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} samples are available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode1{$meshcode1} > 20) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#FF74F1\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} samples are available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode1{$meshcode1} > 10) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#DF86FF\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} samples are available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode1{$meshcode1} > 5) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#9B9FFF\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} samples are available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode1{$meshcode1} > 1) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} samples are available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				elsif ($meshcode1{$meshcode1} == 1) {
					print($filehandleoutput1 "[leaflet-polygon latlngs=\"$latN, $lonW; $latN, $lonE; $latS, $lonE; $latS, $lonW\" stroke=\"true\" color=\"#000000\" weight=\"1\" fill=\"true\" fillColor=\"#00B7FF\" fillOpacity=\"0.5\"]<a href=\"/map/$locusslug/$locusslug$meshcode1/\">$meshcode1{$meshcode1} sample is available on this worldmesh. Link to zoomed map of this worldmesh.</a>[/leaflet-polygon]");
				}
				$nsamples += $meshcode1{$meshcode1};
			}
		}
		print($filehandleoutput1 "<h2>List of submaps</h2>[subpages post_type=\"map\" depth=\"1\" post_status=\"publish\" sort_column=\"post_title\" sort_order=\"asc\"]");
		close($filehandleoutput1);
		my $maps = `wp post list --post_type='map' --format=csv --fields=ID,post_name,post_parent`;
		my $locusid;
		if ($maps !~ /\d+,$locusslug,/) {
			$locusid = `wp post create --porcelain --post_author=1 --post_type='map' --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$locus ($nsamples)' --post_name='$locusslug' --post_excerpt='' --menu_order=1 $locus.html`;
			$locusid =~ s/\r?\n?$//;
		}
		elsif ($maps =~ /(\d+),$locusslug,/) {
			$locusid = $1;
			system("wp post update --post_type='map' --post_author=1 --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$locus ($nsamples)' --post_name='$locusslug' --post_excerpt='' --menu_order=1 $locusid $locus.html \&");
		}
		foreach my $meshcode1 (keys(%meshcode1)) {
			if ($maps !~ /\d+,$locusslug$meshcode1,$locusid\n/) {
				system("wp post create --post_author=1 --post_type='map' --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$locus - $meshcode1 ($meshcode1{$meshcode1})' --post_name='$locusslug$meshcode1' --post_parent='$locusid' --post_excerpt='' $locus$meshcode1.html \&");
			}
			elsif ($maps =~ /(\d+),$locusslug$meshcode1,$locusid\n/) {
				my $locusmeshcode1id = $1;
				system("wp post update --post_type='map' --post_author=1 --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='$locus - $meshcode1 ($meshcode1{$meshcode1})' --post_name='$locusslug$meshcode1' --post_parent='$locusid' --post_excerpt='' $locusmeshcode1id $locus$meshcode1.html \&");
			}
		}
	}
}
{
	open($filehandleoutput1, "> Home.html") or die(__LINE__ . ": Cannot open \"Home.html\"\n");
	print($filehandleoutput1 "[if-logout]<p>ANEMONE DB is a database for biodiversity surveys based on environmental DNA metabarcoding. It is managed by Dr. Michio Kondoh, Ecological Integration Lab. at Tohoku University.</p><p><a href=\"/login/\">Link to login page.</a></p><p><a href=\"/register/\">Link to registration page.</a></p>");
	foreach my $locus (@loci) {
		open($filehandleinput1, "< $locus.html") or die(__LINE__ . ": Cannot open \"$locus.html\".");
		while (<$filehandleinput1>) {
			s/<h2>List of submaps<\/h2>.+//;
			s/<a href=[^>]+>//g;
			s/<\/a>//g;
			s/ Link to zoomed map of this worldmesh\.//g;
			s/<h2>.+?<\/h2>/<h2>Current distribution of all $locus metabarcoding samples<\/h2>/;
			print($filehandleoutput1 $_);
		}
		close($filehandleinput1);
	}
	print($filehandleoutput1 "[/if-logout][if-login]");
	foreach my $locus (@loci) {
		open($filehandleinput1, "< $locus.html") or die(__LINE__ . ": Cannot open \"$locus.html\".");
		while (<$filehandleinput1>) {
			s/<h2>List of submaps<\/h2>.+//;
			print($filehandleoutput1 $_);
		}
		close($filehandleinput1);
	}
	print($filehandleoutput1 "<h2>List of top-level metabarcoding maps</h2>[pagelist post_type=\"map\" depth=\"1\" post_status=\"publish\" sort_column=\"menu_order,post_title\" sort_order=\"asc\"][/if-login]");
	close($filehandleoutput1);
	my $pages = `wp post list --post_type='page' --format=csv --fields=ID,post_name`;
	if ($pages !~ /\d+,home\n/) {
		system("wp post create --post_author=1 --post_type='page' --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='Home' --post_excerpt='' Home.html \&");
	}
	elsif ($pages =~ /(\d+),home\n/) {
		my $homeid = $1;
		system("wp post update --post_type='page' --post_author=1 --comment_status='closed' --ping_status='closed' --post_status='publish' --post_title='Home' --post_excerpt='' $homeid Home.html \&");
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

sub meshcode1latlon {
	my $meshcode1 = shift(@_);
	my $code0 = substr($meshcode1, 0, 1);
	$code0 --;
	my $code12 = substr($meshcode1, 1, 3);
	if ($code12 =~ /^00/) {
		$code12 = substr($meshcode1, 3, 1);
	}
	elsif ($code12 =~ /^0/) {
		$code12 = substr($meshcode1, 2, 2);
	}
	else {
		$code12 = substr($meshcode1, 1, 3);
	}
	my $code34;
	if (substr($meshcode1, 4, 1) eq '0') {
		$code34 = substr($meshcode1, 5, 1);
	}
	else {
		$code34 = substr($meshcode1, 4, 2);
	}
	my $latwidth = 2 / 3;
	my $lonwidth = 1;
	my $z = $code0 % 2;
	my $y = (($code0 - $z) / 2) % 2;
	my $x = ($code0 - (2 * $y) - $z) / 4;
	my $latN = (($code12 - $x + 1) * 2) / 3;
	$latN *= (1 - (2 * $x));
	my $lonW = ($code34 + $y) + (100 * $z);
	$lonW *= (1 - (2 * $y));
	my $dlat = 2 / 3;
	my $dlon = 1;
	my $latS = sprintf("%.8f", ($latN - $dlat));
	my $latC = sprintf("%.8f", ($latN - ($dlat / 2)));
	my $lonE = sprintf("%.8f", ($lonW + $dlon));
	my $lonC = sprintf("%.8f", ($lonW + ($dlon / 2)));
	$latN = sprintf("%.8f", $latN);
	$lonW = sprintf("%.8f", $lonW);
	return($latN, $latC, $latS, $lonW, $lonC, $lonE);
}
