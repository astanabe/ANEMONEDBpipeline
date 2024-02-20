team=`pwd | perl -ne 'chomp;@path=split("/");print($path[-3])'`
project=`pwd | perl -ne 'chomp;@path=split("/");print($path[-2])'`
run=`pwd | perl -ne 'chomp;@path=split("/");print($path[-1])'`
runfolder=`pwd`
tempfolder="/work/$USER/temp/$team/$project/$run"
mkdir -p "$tempfolder"

cd "$runfolder" || exit $LINENO

# Set number of processor cores used for computation
export THREADS=`grep -c processor /proc/cpuinfo`

# Copy FASTQ files to temporary firectory
if ! test -d "$tempfolder/fastq_undemultiplexed"; then
cp -R \
fastq_undemultiplexed \
"$tempfolder/" || exit $LINENO
fi

# Demultiplex Type A (If you have undemultiplexed FASTQ files)
if ! test -d "$tempfolder/demultiplexed"; then
clsplitseq \
--runname="$project$run" \
--forwardprimerfile="$runfolder/forwardprimer.fasta" \
--reverseprimerfile="$runfolder/reverseprimer.fasta" \
--truncateN=enable \
--outputmultihit=enable \
--index1file="$runfolder/index1.fasta" \
--index2file="$runfolder/index2.fasta" \
--minqualtag=30 \
--seqnamestyle=illumina \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/fastq_undemultiplexed/"*_R1_*.fastq.?? \
"$tempfolder/fastq_undemultiplexed/"*_I1_*.fastq.?? \
"$tempfolder/fastq_undemultiplexed/"*_I2_*.fastq.?? \
"$tempfolder/fastq_undemultiplexed/"*_R2_*.fastq.?? \
"$tempfolder/demultiplexed" || exit $LINENO
fi

# Delete temporary FASTQ files
rm -rf "$tempfolder/fastq_undemultiplexed" || exit $LINENO

# Repeat analysis on each locus
for locus in `grep -P -o '^>\S+' "$runfolder/forwardprimer.fasta" | perl -npe 's/^>//'`
do \
# Set variables for each locus
if test "$locus" = 'MiFish' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='overall_species_wsp'
includetaxa='\t(Hyperoartia|Myxini|Chondrichthyes|Actinopterygii|Coelacanthiformes|Dipnomorpha)\t'
elif test "$locus" = 'MiMammal' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='overall_species_wsp'
includetaxa='\t(Mammalia)\t'
elif test "$locus" = 'MiBird' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='overall_species_wsp'
includetaxa='\t(Aves)\t'
else
exit $LINENO
fi
standardfasta="$runfolder/$locus/standard.fasta"
stdconctable="$runfolder/$locus/stdconctable.tsv"
solutionvoltable="$runfolder/$locus/solutionvoltable.tsv"
watervoltable="$runfolder/$locus/watervoltable.tsv"

if ! test -d "$tempfolder/$locus/demultiplexed" && ! test -d "$runfolder/$locus/demultiplexed"; then
# Make directory for each locus
mkdir -p \
"$tempfolder/$locus/demultiplexed" || exit $LINENO
mkdir -p \
"$runfolder/$locus" || exit $LINENO
# Move demultiplexed files
ls "$tempfolder/demultiplexed/$project$run"__*__"$locus".*.fastq.xz \
| xargs -I {} sh -c "mv -f {} \"$tempfolder/$locus/demultiplexed/\" || exit $?"
fi

# Concatenate pairs
if ! test -d "$tempfolder/$locus/concatenated" && ! test -d "$runfolder/$locus/concatenated"; then
clconcatpairv \
--mode="$ovlmode" \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/$locus/demultiplexed" \
"$tempfolder/$locus/concatenated" || exit $LINENO
fi

# Apply filtering out low quality sequences
if ! test -d "$tempfolder/$locus/filtered" && ! test -d "$runfolder/$locus/filtered"; then
clfilterseqv \
--maxqual=41 \
--minlen="$minlen" \
--maxlen="$maxlen" \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/$locus/concatenated" \
"$tempfolder/$locus/filtered" || exit $LINENO
fi

# Denoise using DADA2
if ! test -d "$tempfolder/$locus/denoised" && ! test -d "$runfolder/$locus/denoised"; then
cldenoiseseqd \
--pool=pseudo \
--numthreads="$THREADS" \
"$tempfolder/$locus/filtered" \
"$tempfolder/$locus/denoised" || exit $LINENO
fi

# Remove chimeras using UCHIME3
if ! test -d "$tempfolder/$locus/chimeraremoved1" && ! test -d "$runfolder/$locus/chimeraremoved1"; then
clremovechimev \
--mode=denovo \
--uchimedenovo=3 \
"$tempfolder/$locus/denoised" \
"$tempfolder/$locus/chimeraremoved1" || exit $LINENO
fi

if test -s "$standardfasta"; then
# Cluster internal standard sequences to otus
if ! test -d "$tempfolder/$locus/stdclustered" && ! test -d "$runfolder/$locus/stdclustered"; then
clclusterstdv \
--standardseq="$standardfasta" \
--numthreads="$THREADS" \
"$tempfolder/$locus/chimeraremoved1" \
"$tempfolder/$locus/stdclustered" || exit $LINENO
fi
# Remove chimeras using UCHIME3
if test -n "$referencedb"; then
if ! test -d "$tempfolder/$locus/chimeraremoved2" && ! test -d "$runfolder/$locus/chimeraremoved2"; then
clremovechimev \
--mode=ref \
--referencedb="$referencedb" \
--addtoref="$tempfolder/$locus/stdclustered/stdvariations.fasta" \
--numthreads="$THREADS" \
"$tempfolder/$locus/stdclustered" \
"$tempfolder/$locus/chimeraremoved2" || exit $LINENO
fi
fi
else
# Remove chimeras using UCHIME3
if test -n "$referencedb"; then
if ! test -d "$tempfolder/$locus/chimeraremoved2" && ! test -d "$runfolder/$locus/chimeraremoved2"; then
clremovechimev \
--mode=ref \
--referencedb="$referencedb" \
--numthreads="$THREADS" \
"$tempfolder/$locus/chimeraremoved1" \
"$tempfolder/$locus/chimeraremoved2" || exit $LINENO
fi
fi
fi

# Set variables
if test -d "$tempfolder/$locus/chimeraremoved2"; then
alltsv="$tempfolder/$locus/chimeraremoved2/nonchimeras.tsv"
fastafile="$tempfolder/$locus/chimeraremoved2/nonchimeras.fasta"
elif test -d "$tempfolder/$locus/stdclustered"; then
alltsv="$tempfolder/$locus/stdclustered/stdclustered.tsv"
fastafile="$tempfolder/$locus/stdclustered/stdclustered.fasta"
elif test -d "$tempfolder/$locus/chimeraremoved1"; then
alltsv="$tempfolder/$locus/chimeraremoved1/nonchimeras.tsv"
fastafile="$tempfolder/$locus/chimeraremoved1/nonchimeras.fasta"
elif test -d "$runfolder/$locus/chimeraremoved2"; then
alltsv="$runfolder/$locus/chimeraremoved2/nonchimeras.tsv"
fastafile="$runfolder/$locus/chimeraremoved2/nonchimeras.fasta"
elif test -d "$runfolder/$locus/stdclustered"; then
alltsv="$runfolder/$locus/stdclustered/stdclustered.tsv"
fastafile="$runfolder/$locus/stdclustered/stdclustered.fasta"
elif test -d "$runfolder/$locus/chimeraremoved1"; then
alltsv="$runfolder/$locus/chimeraremoved1/nonchimeras.tsv"
fastafile="$runfolder/$locus/chimeraremoved1/nonchimeras.fasta"
else
exit $LINENO
fi

if ! test -d "$tempfolder/$locus/taxonomy" && ! test -d "$runfolder/$locus/taxonomy"; then

echo "Taxonomy"

# Make directory for output
mkdir -p \
"$tempfolder/$locus/taxonomy" || exit $LINENO

# Set variables
if test -e "$runfolder/../../../$locus/qc_$blastdb.identdb"; then
qcidentdb="$runfolder/../../../$locus/qc_$blastdb.identdb"
fi
if test -e "$runfolder/../../../$locus/3nn_$blastdb.identdb"; then
threennidentdb="$runfolder/../../../$locus/3nn_$blastdb.identdb"
fi
if ! test -e "$runfolder/../../../$locus/qc_$blastdb.identdb" && ! test -e "$runfolder/../../../$locus/3nn_$blastdb.identdb"; then
mkdir -p \
"$runfolder/../../../$locus" || exit $LINENO
fi

# Make cachedb
if ! test -d "$tempfolder/$locus/taxonomy/cachedb_$blastdb"; then
if test -n "$qcidentdb"; then
if test -s "$standardfasta"; then
clmakecachedb \
--identdb="$qcidentdb" \
--blastdb="$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $LINENO
else
clmakecachedb \
--identdb="$qcidentdb" \
--blastdb="$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $LINENO
fi
else
if test -s "$standardfasta"; then
clmakecachedb \
--blastdb="$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $LINENO
else
clmakecachedb \
--blastdb="$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $LINENO
fi
fi
fi

# Assign taxonomy based on QCauto method using $blastdb
if ! test -e "$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt"; then
if test -n "$qcidentdb"; then
if test -s "$standardfasta"; then
clidentseq \
--method=QC \
--identdb="$qcidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $LINENO
else
clidentseq \
--method=QC \
--identdb="$qcidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $LINENO
fi
else
if test -s "$standardfasta"; then
clidentseq \
--method=QC \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $LINENO
else
clidentseq \
--method=QC \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $LINENO
fi
fi
fi

if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv"; then
classigntax \
--taxdb="$blastdb" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" || exit $LINENO
fi

# Assign taxonomy based on (95%-)3NN method using $blastdb
if ! test -e "$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt"; then
if test -n "$threennidentdb"; then
if test -s "$standardfasta"; then
clidentseq \
--method=3,95% \
--identdb="$threennidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $LINENO
else
clidentseq \
--method=3,95% \
--identdb="$threennidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $LINENO
fi
else
if test -s "$standardfasta"; then
clidentseq \
--method=3,95% \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $LINENO
else
clidentseq \
--method=3,95% \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $LINENO
fi
fi
fi

if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv"; then
classigntax \
--taxdb="$blastdb" \
--minnsupporter=3 \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" \
"$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv" || exit $LINENO
fi

# Merge taxonomic assignment results
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv"; then
clmergeassign \
--preferlower \
--priority=descend \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv" || exit $LINENO
fi

# Fill blanks in taxonomy
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv; then
clfillassign \
--fullfill=enable \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv; then
clfillassign \
--fullfill=enable \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv || exit $LINENO
fi

# Update identdb
clmakeidentdb \
--append \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" \
"$runfolder/../../../$locus/qc_$blastdb.identdb" &
clmakeidentdb \
--append \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" \
"$runfolder/../../../$locus/3nn_$blastdb.identdb" &
wait || exit $LINENO

fi

# Set variables
qctaxonomy="$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv
qc3nntaxonomy="$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv

if ! test -d "$tempfolder/$locus/community" && ! test -d "$runfolder/$locus/community"; then

echo "Community"

# Make directory for output
mkdir -p \
"$tempfolder/$locus/community" || exit $LINENO

# Combine replicates
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_all.tsv"; then
if test -e "$runfolder/$locus/replicatelist.txt"; then
clfiltersum \
--replicatelist="$runfolder/$locus/replicatelist.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_all.tsv" || exit $LINENO
else
cp \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_all.tsv" || exit $LINENO
fi
fi
alltsv="$tempfolder/$locus/community/sample_otu_matrix_all.tsv"

# Extract internal standard OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_standard.tsv"; then
if test -s "$standardfasta"; then
clfiltersum \
--otuseq="$standardfasta" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_standard.tsv" || exit $LINENO
fi
fi

# Pick jawless vertebrates, cartilaginous fishes, ray-finned fishes, coelacanths and lungfishes
grep -P "$includetaxa" \
"$qctaxonomy" \
| cut -f 1 \
> "$tempfolder/$locus/targetotu.txt"

# Pick nontarget taxa
grep -P -v "$includetaxa" \
"$qctaxonomy" \
| cut -f 1 \
> "$tempfolder/$locus/nontargetotu.txt"

# Delete nontarget OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_target.tsv"; then
if test -s "$tempfolder/$locus/targetotu.txt"; then
clfiltersum \
--otulist="$tempfolder/$locus/targetotu.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_target.tsv" || exit $LINENO
fi
fi

# Extract nontarget OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv"; then
if test -s "$tempfolder/$locus/nontargetotu.txt"; then
clfiltersum \
--otulist="$tempfolder/$locus/nontargetotu.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv" || exit $LINENO
fi
fi

# Set variables
if test -s "$standardfasta"; then
standardtsv="$tempfolder/$locus/community/sample_otu_matrix_standard.tsv"
fi
if test -s "$tempfolder/$locus/targetotu.txt"; then
targettsv="$tempfolder/$locus/community/sample_otu_matrix_target.tsv"
fi
if test -s "$tempfolder/$locus/nontargetotu.txt"; then
nontargettsv="$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv"
fi

# Make target community files
if test -n "$targettsv"; then
if ! test -e "$tempfolder/$locus/community/community_qc_target.tsv"; then
clsumtaxa \
--taxfile="$qctaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$targettsv" \
"$tempfolder/$locus/community/community_qc_target.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community/community_qc3nn_target.tsv"; then
clsumtaxa \
--taxfile="$qc3nntaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$targettsv" \
"$tempfolder/$locus/community/community_qc3nn_target.tsv" || exit $LINENO
fi
fi

# Make nontarget community files
if test -n "$nontargettsv"; then
if ! test -e "$tempfolder/$locus/community/community_qc_nontarget.tsv"; then
clsumtaxa \
--taxfile="$qctaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$nontargettsv" \
"$tempfolder/$locus/community/community_qc_nontarget.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community/community_qc3nn_nontarget.tsv"; then
clsumtaxa \
--taxfile="$qc3nntaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$nontargettsv" \
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" || exit $LINENO
fi
fi

# Estimate DNA concentration and make community files of target taxa
if test -n "$targettsv" && test -n "$standardtsv" && test -s "$stdconctable" && test -s "$solutionvoltable" && test -s "$watervoltable"; then
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv"; then
clestimateconc \
--stdconctable="$stdconctable" \
--stdtable="$standardtsv" \
--solutionvoltable="$solutionvoltable" \
--watervoltable="$watervoltable" \
--numthreads="$THREADS" \
"$targettsv" \
"$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community/community_qc_target_estimated.tsv"; then
clsumtaxa \
--taxfile="$qctaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv" \
"$tempfolder/$locus/community/community_qc_target_estimated.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv"; then
clsumtaxa \
--taxfile="$qc3nntaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv" \
"$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community_qc_target.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc_target.tsv" \
"$tempfolder/$locus/community/community_qc_target_estimated.tsv" \
"$tempfolder/$locus/community_qc_target.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_target.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc3nn_target.tsv" \
"$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv" \
"$tempfolder/$locus/community_qc3nn_target.tsv" || exit $LINENO
fi
else
if ! test -e "$tempfolder/$locus/community_qc_target.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc_target.tsv" \
"$tempfolder/$locus/community_qc_target.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_target.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc3nn_target.tsv" \
"$tempfolder/$locus/community_qc3nn_target.tsv" || exit $LINENO
fi
fi

# Estimate DNA concentration and make community files of nontarget taxa
if test -n "$nontargettsv" && test -n "$standardtsv" && test -s "$stdconctable" && test -s "$solutionvoltable" && test -s "$watervoltable"; then
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv"; then
clestimateconc \
--stdconctable="$stdconctable" \
--stdtable="$standardtsv" \
--solutionvoltable="$solutionvoltable" \
--watervoltable="$watervoltable" \
--numthreads="$THREADS" \
"$nontargettsv" \
"$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv"; then
clsumtaxa \
--taxfile="$qctaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv" \
"$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv"; then
clsumtaxa \
--taxfile="$qc3nntaxonomy" \
--targetrank=sequence \
--otuseq="$fastafile" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv" \
"$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community_qc_nontarget.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc_nontarget.tsv" \
"$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv" \
"$tempfolder/$locus/community_qc_nontarget.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_nontarget.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" \
"$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv" \
"$tempfolder/$locus/community_qc3nn_nontarget.tsv" || exit $LINENO
fi
else
if ! test -e "$tempfolder/$locus/community_qc_nontarget.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc_nontarget.tsv" \
"$tempfolder/$locus/community_qc_nontarget.tsv" || exit $LINENO
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_nontarget.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" \
"$tempfolder/$locus/community_qc3nn_nontarget.tsv" || exit $LINENO
fi
fi

# Make standard community files
if test -n "$standardtsv"; then
if ! test -e "$tempfolder/$locus/community_standard.tsv"; then
echo -e 'query\tspecies' > "$tempfolder/$locus/taxonomy_standard.tsv"
grep -P -o '^>\S+' "$standardfasta" | perl -ne 'chomp;s/^>//;print("$_\t$_\n")' >> "$tempfolder/$locus/taxonomy_standard.tsv"
clsumtaxa \
--taxfile="$tempfolder/$locus/taxonomy_standard.tsv" \
--targetrank=sequence \
--otuseq="$standardfasta" \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
"$standardtsv" \
"$tempfolder/$locus/community_standard.tsv" || exit $LINENO
perl -i -npe 's/\tspecies\t/\tstandard\t/' "$tempfolder/$locus/community_standard.tsv" || exit $LINENO
rm -f "$tempfolder/$locus/taxonomy_standard.tsv" || exit $LINENO
fi
fi

# Generate word cloud
if test -n "$targettsv"; then
cd "$tempfolder"
clplotwordcloud \
--taxfile="$qc3nntaxonomy" \
--append \
--numthreads="$THREADS" \
--anemone \
"$targettsv" \
samplename || exit $LINENO
ls -d "$project$run"__*__"$locus" \
| xargs -I {} sh -c "mv -f {} \"$runfolder/$locus/\" || exit $?"
cd "$runfolder"
fi

fi

# Move results from work to home and delete temporary files
rm \
-rf \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $LINENO
cp \
-R \
"$tempfolder/$locus" \
"$runfolder/" || exit $LINENO
rm \
-rf \
"$tempfolder/$locus" || exit $LINENO

# Unset variables
unset ovlmode
unset minlen
unset maxlen
unset referencedb
unset blastdb
unset includetaxa
unset standardfasta
unset stdconctable
unset solutionvoltable
unset watervoltable
unset alltsv
unset fastafile
unset qcidentdb
unset threennidentdb
unset qctaxonomy
unset qc3nntaxonomy

done

rm \
-rf \
"$tempfolder" || exit $LINENO
