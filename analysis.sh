team=`pwd | perl -ne 'chomp;@path=split("/");print($path[-3])'`
project=`pwd | perl -ne 'chomp;@path=split("/");print($path[-2])'`
run=`pwd | perl -ne 'chomp;@path=split("/");print($path[-1])'`
runfolder=`pwd`
tempfolder="/work/$USER/temp/$team/$project/$run"
mkdir -p "$tempfolder"

cd "$runfolder" || exit "$?"

# Set number of processor cores used for computation
export THREADS=`grep -c processor /proc/cpuinfo`

# Copy FASTQ files to temporary firectory
if ! test -d "$tempfolder/fastq_undemultiplexed"; then
cp -R \
fastq_undemultiplexed \
"$tempfolder/" || exit "$?"
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
"$tempfolder/demultiplexed" || exit "$?"
fi

# Delete temporary FASTQ files
rm -rf "$tempfolder/fastq_undemultiplexed" || exit "$?"

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
exit 1
fi
standardfasta="$runfolder/$locus/standard.fasta"
stdconctable="$runfolder/$locus/stdconctable.tsv"
solutionvoltable="$runfolder/$locus/solutionvoltable.tsv"
watervoltable="$runfolder/$locus/watervoltable.tsv"

if ! test -d "$tempfolder/$locus/demultiplexed"; then
# Make directory for each locus
mkdir -p \
"$tempfolder/$locus/demultiplexed" || exit "$?"
mkdir -p \
"$runfolder/$locus" || exit "$?"
# Move demultiplexed files
ls "$tempfolder/demultiplexed/$project$run"__*__"$locus".*.fastq.xz \
| xargs -I {} sh -c "mv -f {} \"$tempfolder/$locus/demultiplexed/\" || exit \"$?\""
fi

# Concatenate pairs
if ! test -d "$tempfolder/$locus/concatenated"; then
clconcatpairv \
--mode="$ovlmode" \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/$locus/demultiplexed" \
"$tempfolder/$locus/concatenated" || exit "$?"
fi

# Apply filtering out low quality sequences
if ! test -d "$tempfolder/$locus/filtered"; then
clfilterseqv \
--maxqual=41 \
--minlen="$minlen" \
--maxlen="$maxlen" \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/$locus/concatenated" \
"$tempfolder/$locus/filtered" || exit "$?"
fi

# Denoise using DADA2
if ! test -d "$tempfolder/$locus/denoised"; then
cldenoiseseqd \
--pool=pseudo \
--numthreads="$THREADS" \
"$tempfolder/$locus/filtered" \
"$tempfolder/$locus/denoised" || exit "$?"
fi

# Remove chimeras using UCHIME3
if ! test -d "$tempfolder/$locus/chimeraremoved1"; then
clremovechimev \
--mode=denovo \
--uchimedenovo=3 \
"$tempfolder/$locus/denoised" \
"$tempfolder/$locus/chimeraremoved1" || exit "$?"
fi

if test -s "$standardfasta"; then
# Cluster internal standard sequences to otus
if ! test -d "$tempfolder/$locus/stdclustered"; then
clclusterstdv \
--standardseq="$standardfasta" \
--numthreads="$THREADS" \
"$tempfolder/$locus/chimeraremoved1" \
"$tempfolder/$locus/stdclustered" || exit "$?"
fi
# Remove chimeras using UCHIME3
if test -n "$referencedb"; then
if ! test -d "$tempfolder/$locus/chimeraremoved2"; then
clremovechimev \
--mode=ref \
--referencedb="$referencedb" \
--addtoref="$tempfolder/$locus/stdclustered/stdvariations.fasta" \
--numthreads="$THREADS" \
"$tempfolder/$locus/stdclustered" \
"$tempfolder/$locus/chimeraremoved2" || exit "$?"
fi
fi
else
# Remove chimeras using UCHIME3
if test -n "$referencedb"; then
if ! test -d "$tempfolder/$locus/chimeraremoved2"; then
clremovechimev \
--mode=ref \
--referencedb="$referencedb" \
--numthreads="$THREADS" \
"$tempfolder/$locus/chimeraremoved1" \
"$tempfolder/$locus/chimeraremoved2" || exit "$?"
fi
fi
fi

# Set variables
if test -e "$tempfolder/$locus/chimeraremoved2"; then
alltsv="$tempfolder/$locus/chimeraremoved2/nonchimeras.tsv"
fastafile="$tempfolder/$locus/chimeraremoved2/nonchimeras.fasta"
elif test -e "$tempfolder/$locus/stdclustered"; then
alltsv="$tempfolder/$locus/stdclustered/stdclustered.tsv"
fastafile="$tempfolder/$locus/stdclustered/stdclustered.fasta"
elif test -e "$tempfolder/$locus/chimeraremoved1"; then
alltsv="$tempfolder/$locus/chimeraremoved1/nonchimeras.tsv"
fastafile="$tempfolder/$locus/chimeraremoved1/nonchimeras.fasta"
else
exit 1
fi

# Make directory for output
mkdir -p \
"$tempfolder/$locus/taxonomy" || exit "$?"

# Set variables
if test -e "$runfolder/../../../$locus/qc_$blastdb.identdb"; then
qcidentdb="$runfolder/../../../$locus/qc_$blastdb.identdb"
fi
if test -e "$runfolder/../../../$locus/3nn_$blastdb.identdb"; then
threennidentdb="$runfolder/../../../$locus/3nn_$blastdb.identdb"
fi
if ! test -e "$runfolder/../../../$locus/qc_$blastdb.identdb" && ! test -e "$runfolder/../../../$locus/3nn_$blastdb.identdb"; then
mkdir -p \
"$runfolder/../../../$locus" || exit "$?"
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
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit "$?"
else
clmakecachedb \
--identdb="$qcidentdb" \
--blastdb="$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit "$?"
fi
else
if test -s "$standardfasta"; then
clmakecachedb \
--blastdb="$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit "$?"
else
clmakecachedb \
--blastdb="$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit "$?"
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
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit "$?"
else
clidentseq \
--method=QC \
--identdb="$qcidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit "$?"
fi
else
if test -s "$standardfasta"; then
clidentseq \
--method=QC \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit "$?"
else
clidentseq \
--method=QC \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit "$?"
fi
fi
fi

if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv"; then
classigntax \
--taxdb="$blastdb" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" || exit "$?"
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
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit "$?"
else
clidentseq \
--method=3,95% \
--identdb="$threennidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit "$?"
fi
else
if test -s "$standardfasta"; then
clidentseq \
--method=3,95% \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit "$?"
else
clidentseq \
--method=3,95% \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit "$?"
fi
fi
fi

if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv"; then
classigntax \
--taxdb="$blastdb" \
--minnsupporter=3 \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" \
"$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv" || exit "$?"
fi

# Merge taxonomic assignment results
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv"; then
clmergeassign \
--preferlower \
--priority=descend \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv" || exit "$?"
fi

# Fill blanks in taxonomy
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv; then
clfillassign \
--fullfill=enable \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv || exit "$?"
fi
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv; then
clfillassign \
--fullfill=enable \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv || exit "$?"
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
wait || exit "$?"

# Set variables
qctaxonomy="$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv
qc3nntaxonomy="$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv

# Make directory for output
mkdir -p \
"$tempfolder/$locus/community" || exit "$?"

# Combine replicates
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_all.tsv"; then
if test -e "$runfolder/$locus/replicatelist.txt"; then
clfiltersum \
--replicatelist="$runfolder/$locus/replicatelist.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_all.tsv" || exit "$?"
else
cp \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_all.tsv" || exit "$?"
fi
fi
alltsv="$tempfolder/$locus/community/sample_otu_matrix_all.tsv"

# Extract internal standard OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_standard.tsv"; then
if test -s "$standardfasta"; then
clfiltersum \
--otuseq="$standardfasta" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_standard.tsv" || exit "$?"
fi
fi

# Pick jawless vertebrates, cartilaginous fishes, ray-finned fishes, coelacanths and lungfishes
grep -P "$includetaxa" \
"$qctaxonomy" \
| cut -f 1 \
> "$runfolder/$locus/targetotu.txt"

# Pick nontarget taxa
grep -P -v "$includetaxa" \
"$qctaxonomy" \
| cut -f 1 \
> "$runfolder/$locus/nontargetotu.txt"

# Delete nontarget OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_target.tsv"; then
if test -s "$runfolder/$locus/targetotu.txt"; then
clfiltersum \
--otulist="$runfolder/$locus/targetotu.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_target.tsv" || exit "$?"
fi
fi

# Extract nontarget OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv"; then
if test -s "$runfolder/$locus/nontargetotu.txt"; then
clfiltersum \
--otulist="$runfolder/$locus/nontargetotu.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv" || exit "$?"
fi
fi

# Set variables
if test -s "$standardfasta"; then
standardtsv="$tempfolder/$locus/community/sample_otu_matrix_standard.tsv"
fi
if test -s "$runfolder/$locus/targetotu.txt"; then
targettsv="$tempfolder/$locus/community/sample_otu_matrix_target.tsv"
fi
if test -s "$runfolder/$locus/nontargetotu.txt"; then
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
"$tempfolder/$locus/community/community_qc_target.tsv" || exit "$?"
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
"$tempfolder/$locus/community/community_qc3nn_target.tsv" || exit "$?"
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
"$tempfolder/$locus/community/community_qc_nontarget.tsv" || exit "$?"
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
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" || exit "$?"
fi
fi

# Generate word cloud
if test -n "$targettsv"; then
clplotwordcloud \
--taxfile="$qc3nntaxonomy" \
--append \
--numthreads="$THREADS" \
--anemone \
"$targettsv" \
samplename || exit "$?"
ls -d "$project$run"__*__"$locus" \
| xargs -I {} sh -c "mv -f {} \"$runfolder/$locus/\" || exit \"$?\""
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
"$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv" || exit "$?"
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
"$tempfolder/$locus/community/community_qc_target_estimated.tsv" || exit "$?"
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
"$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv" || exit "$?"
fi
if ! test -e "$runfolder/$locus/community_qc_target.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc_target.tsv" \
"$tempfolder/$locus/community/community_qc_target_estimated.tsv" \
"$runfolder/$locus/community_qc_target.tsv" || exit "$?"
fi
if ! test -e "$runfolder/$locus/community_qc3nn_target.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc3nn_target.tsv" \
"$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv" \
"$runfolder/$locus/community_qc3nn_target.tsv" || exit "$?"
fi
else
if ! test -e "$runfolder/$locus/community_qc_target.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc_target.tsv" \
"$runfolder/$locus/community_qc_target.tsv" || exit "$?"
fi
if ! test -e "$runfolder/$locus/community_qc3nn_target.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc3nn_target.tsv" \
"$runfolder/$locus/community_qc3nn_target.tsv" || exit "$?"
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
"$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv" || exit "$?"
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
"$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv" || exit "$?"
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
"$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv" || exit "$?"
fi
if ! test -e "$runfolder/$locus/community_qc_nontarget.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc_nontarget.tsv" \
"$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv" \
"$runfolder/$locus/community_qc_nontarget.tsv" || exit "$?"
fi
if ! test -e "$runfolder/$locus/community_qc3nn_nontarget.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" \
"$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv" \
"$runfolder/$locus/community_qc3nn_nontarget.tsv" || exit "$?"
fi
else
if ! test -e "$runfolder/$locus/community_qc_nontarget.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc_nontarget.tsv" \
"$runfolder/$locus/community_qc_nontarget.tsv" || exit "$?"
fi
if ! test -e "$runfolder/$locus/community_qc3nn_nontarget.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" \
"$runfolder/$locus/community_qc3nn_nontarget.tsv" || exit "$?"
fi
fi

# Make standard community files
if test -n "$standardtsv"; then
if ! test -e "$runfolder/$locus/community_standard.tsv"; then
echo -e 'query\tspecies' > "$runfolder/$locus/taxonomy_standard.tsv"
grep -P -o '^>\S+' "$standardfasta" | perl -ne 'chomp;s/^>//;print("$_\t$_\n")' >> "$runfolder/$locus/taxonomy_standard.tsv"
clsumtaxa \
--taxfile="$runfolder/$locus/taxonomy_standard.tsv" \
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
"$runfolder/$locus/community_standard.tsv" || exit "$?"
perl -i -npe 's/\tspecies\t/\tstandard\t/' "$runfolder/$locus/community_standard.tsv" || exit "$?"
rm -f "$runfolder/$locus/taxonomy_standard.tsv" || exit "$?"
fi
fi

# Move results from work to home and delete temporary files
rm \
-rf \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb"
mv \
-f \
"$tempfolder/$locus/demultiplexed" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/concatenated" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/filtered" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/denoised" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/chimeraremoved1" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/stdclustered" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/chimeraremoved2" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/taxonomy" \
"$runfolder/$locus/"
mv \
-f \
"$tempfolder/$locus/community" \
"$runfolder/$locus/"
rm \
-rf \
"$tempfolder/$locus"

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
"$tempfolder"
