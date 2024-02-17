team=ANEMONE
project=2024ANEMONE
projectfolder=/home/$USER/project/$team/$project
run=RUN01
runfolder=$projectfolder/$run
mkdir -p $runfolder
tempfolder=/work/$USER/temp/$team/$project/$run
mkdir -p $tempfolder

cd $runfolder || exit $?

# Set number of processor cores used for computation
export THREADS=128

# Copy FASTQ files to temporary firectory
cp -R \
fastq_undemultiplexed \
$tempfolder/ || exit $?

# Demultiplex Type A (If you have undemultiplexed FASTQ files)
clsplitseq \
--runname=$project$run \
--forwardprimerfile=$runfolder/forwardprimer.fasta \
--reverseprimerfile=$runfolder/reverseprimer.fasta \
--truncateN=enable \
--index1file=$runfolder/index1.fasta \
--index2file=$runfolder/index2.fasta \
--minqualtag=30 \
--compress=xz \
--numthreads=$THREADS \
--seqnamestyle=illumina \
$tempfolder/fastq_undemultiplexed/*_R1_*.fastq.?? \
$tempfolder/fastq_undemultiplexed/*_I1_*.fastq.?? \
$tempfolder/fastq_undemultiplexed/*_I2_*.fastq.?? \
$tempfolder/fastq_undemultiplexed/*_R2_*.fastq.?? \
$tempfolder/demultiplexed || exit $?

# Delete temporary FASTQ files
rm -rf $tempfolder/fastq_undemultiplexed || exit $?

# Repeat analysis on each locus
for locus in `grep -P -o '^>\S+' $runfolder/forwardprimer.fasta | perl -npe 's/^>//'`
do \
# Set variables for each locus
if test $locus = 'MiFish' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='overall_species_wsp'
includetaxa='\t(Hyperoartia|Myxini|Chondrichthyes|Actinopterygii|Coelacanthiformes|Dipnomorpha)\t'
elif test $locus = 'MiMammal' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='overall_species_wsp'
includetaxa='\t(Mammalia)\t'
elif test $locus = 'MiBird' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='overall_species_wsp'
includetaxa='\t(Aves)\t'
else
exit(1)
fi
standardfasta=$runfolder/$locus/standard.fasta
stdconctable=$runfolder/$locus/stdconctable.tsv
solutionvoltable=$runfolder/$locus/solutionvoltable.tsv
watervoltable=$runfolder/$locus/watervoltable.tsv

# Make directory for each locus
mkdir -p \
$tempfolder/$locus/demultiplexed || exit $?
mkdir -p \
$runfolder/$locus || exit $?

# Move demultiplexed files
ls $tempfolder/demultiplexed/$project$run\__*__$locus.fastq.xz \
| xargs -I {} sh -c "mv -f {} $tempfolder/$locus/demultiplexed/ || exit $?"

# Concatenate pairs
clconcatpairv \
--mode=$ovlmode \
--compress=xz \
--numthreads=$THREADS \
$tempfolder/$locus/demultiplexed \
$tempfolder/$locus/concatenated || exit $?

# Apply filtering out low quality sequences
clfilterseqv \
--maxqual=41 \
--minlen=$minlen \
--maxlen=$maxlen \
--maxnee=2.0 \
--maxnNs=0 \
--compress=xz \
--numthreads=$THREADS \
$tempfolder/$locus/concatenated \
$tempfolder/$locus/filtered || exit $?

# Denoise using DADA2
cldenoiseseqd \
--pool=pseudo \
--numthreads=$THREADS \
$tempfolder/$locus/filtered \
$tempfolder/$locus/denoised || exit $?

# Remove chimeras using UCHIME3
clremovechimev \
--mode=denovo \
--uchimedenovo=3 \
$tempfolder/$locus/denoised \
$tempfolder/$locus/chimeraremoved1 || exit $?

if test -s $standardfasta; then
# Cluster internal standard sequences to otus
clclusterstdv \
--standardseq=$standardfasta \
--numthreads=$THREADS \
$tempfolder/$locus/chimeraremoved1 \
$tempfolder/$locus/stdclustered || exit $?
# Remove chimeras using UCHIME3
if test -n $referencedb; then
clremovechimev \
--mode=ref \
--referencedb=$referencedb \
--addtoref=$tempfolder/$locus/stdclustered/stdvariations.fasta \
--numthreads=$THREADS \
$tempfolder/$locus/stdclustered \
$tempfolder/$locus/chimeraremoved2 || exit $?
fi
else
# Remove chimeras using UCHIME3
if test -n $referencedb; then
clremovechimev \
--mode=ref \
--referencedb=$referencedb \
--numthreads=$THREADS \
$tempfolder/$locus/chimeraremoved1 \
$tempfolder/$locus/chimeraremoved2 || exit $?
fi
fi

# Set variables
if test -e $tempfolder/$locus/chimeraremoved2; then
alltsv=$tempfolder/$locus/chimeraremoved2/nonchimeras.tsv
fastafile=$tempfolder/$locus/chimeraremoved2/nonchimeras.fasta
elif test -e $tempfolder/$locus/stdclustered; then
alltsv=$tempfolder/$locus/stdclustered/stdclustered.tsv
fastafile=$tempfolder/$locus/stdclustered/stdclustered.fasta
elif test -e $tempfolder/$locus/chimeraremoved1; then
alltsv=$tempfolder/$locus/chimeraremoved1/nonchimeras.tsv
fastafile=$tempfolder/$locus/chimeraremoved1/nonchimeras.fasta
else
exit(1)
fi

# Make directory for output
mkdir -p \
$tempfolder/$locus/taxonomy || exit $?

# Set variables
if test -e $projectfolder/../../$locus/qc_$blastdb.identdb; then
qcidentdb=$projectfolder/../../$locus/qc_$blastdb.identdb
fi
if test -e $projectfolder/../../$locus/3nn_$blastdb.identdb; then
3nnidentdb=$projectfolder/../../$locus/3nn_$blastdb.identdb
fi
if ! test -e $projectfolder/../../$locus/qc_$blastdb.identdb && ! test -e $projectfolder/../../$locus/3nn_$blastdb.identdb; then
mkdir -p \
$projectfolder/../../$locus || exit $?
fi

# Assign taxonomy based on QCauto method using $blastdb
if test -n $qcidentdb; then
if test -s $standardfasta; then
clmakecachedb \
--identdb=$qcidentdb \
--blastdb=$blastdb \
--ignoreotuseq=$standardfasta \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/cachedb_$blastdb || exit $?
clidentseq \
--method=QC \
--identdb=$qcidentdb \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--ignoreotuseq=$standardfasta \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt || exit $?
else
clmakecachedb \
--identdb=$qcidentdb \
--blastdb=$blastdb \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/cachedb_$blastdb || exit $?
clidentseq \
--method=QC \
--identdb=$qcidentdb \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt || exit $?
fi
else
if test -s $standardfasta; then
clmakecachedb \
--blastdb=$blastdb \
--ignoreotuseq=$standardfasta \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/cachedb_$blastdb || exit $?
clidentseq \
--method=QC \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--ignoreotuseq=$standardfasta \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt || exit $?
else
clmakecachedb \
--blastdb=$blastdb \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/cachedb_$blastdb || exit $?
clidentseq \
--method=QC \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt || exit $?
fi
fi

classigntax \
--taxdb=$blastdb \
$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt \
$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv || exit $?

# Assign taxonomy based on (95%-)3NN method using $blastdb
if test -n $3nnidentdb; then
if test -s $standardfasta; then
clidentseq \
--method=3,95% \
--identdb=$3nnidentdb \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--ignoreotuseq=$standardfasta \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt || exit $?
else
clidentseq \
--method=3,95% \
--identdb=$3nnidentdb \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt || exit $?
fi
else
if test -s $standardfasta; then
clidentseq \
--method=3,95% \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--ignoreotuseq=$standardfasta \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_3nn_species_wsp.txt || exit $?
else
clidentseq \
--method=3,95% \
--blastdb=$tempfolder/$locus/taxonomy/cachedb_$blastdb \
--numthreads=$THREADS \
$fastafile \
$tempfolder/$locus/taxonomy/neighborhoods_3nn_species_wsp.txt || exit $?
fi
fi

classigntax \
--taxdb=$blastdb \
--minnsupporter=3 \
$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt \
$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv || exit $?

# Merge taxonomic assignment results
clmergeassign \
--preferlower \
--priority=descend \
$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv \
$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv \
$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv || exit $?

# Fill blanks in taxonomy
clfillassign \
--fullfill=ENABLE \
$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv \
$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb_filled.tsv || exit $?
clfillassign \
--fullfill=ENABLE \
$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv \
$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb_filled.tsv || exit $?

# Update identdb
clmakeidentdb \
--append \
$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt \
$qcidentdb &
clmakeidentdb \
--append \
$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt \
$3nnidentdb &
wait || exit $?

# Set variables
qctaxonomy=$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb_filled.tsv
qc3nntaxonomy=$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb_filled.tsv

# Make directory for output
mkdir -p \
$tempfolder/$locus/community || exit $?

# Combine replicates
if test -e $runfolder/$locus/replicatelist.txt; then
clfiltersum \
--replicatelist=$runfolder/$locus/replicatelist.txt \
$alltsv \
$tempfolder/$locus/community/sample_otu_matrix_all.tsv || exit $?
else
cp \
$alltsv \
$tempfolder/$locus/community/sample_otu_matrix_all.tsv || exit $?
fi
alltsv=$tempfolder/$locus/community/sample_otu_matrix_all.tsv

# Extract internal standard OTUs
if test -s $standardfasta; then
clfiltersum \
--otuseq=$standardfasta \
$alltsv \
$tempfolder/$locus/community/sample_otu_matrix_standard.tsv || exit $?
fi

# Pick jawless vertebrates, cartilaginous fishes, ray-finned fishes, coelacanths and lungfishes
grep -P $includetaxa \
$qctaxonomy \
| cut -f 1 \
> $runfolder/$locus/targetotu.txt

# Pick nontarget taxa
grep -P -v $includetaxa \
$qctaxonomy \
| cut -f 1 \
> $runfolder/$locus/nontargetotu.txt

# Delete nontarget OTUs
if test -s $runfolder/$locus/targetotu.txt; then
clfiltersum \
--otulist=$runfolder/$locus/targetotu.txt \
$alltsv \
$tempfolder/$locus/community/sample_otu_matrix_target.tsv || exit $?
fi

# Extract nontarget OTUs
if test -s $runfolder/$locus/nontargetotu.txt; then
clfiltersum \
--otulist=$runfolder/$locus/nontargetotu.txt \
$alltsv \
$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv || exit $?
fi

# Set variables
if test -s $standardfasta; then
standardtsv=$tempfolder/$locus/community/sample_otu_matrix_standard.tsv
fi
if test -s $runfolder/$locus/targetotu.txt; then
targettsv=$tempfolder/$locus/community/sample_otu_matrix_target.tsv
fi
if test -s $runfolder/$locus/nontargetotu.txt; then
nontargettsv=$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv
fi

# Make target community files
if test -n $targettsv; then
clsumtaxa \
--taxfile=$qctaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$targettsv \
$tempfolder/$locus/community/community_qc_target.tsv || exit $?
clsumtaxa \
--taxfile=$qc3nntaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$targettsv \
$tempfolder/$locus/community/community_qc3nn_target.tsv || exit $?
fi

# Make nontarget community files
if test -n $nontargettsv; then
clsumtaxa \
--taxfile=$qctaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$nontargettsv \
$tempfolder/$locus/community/community_qc_nontarget.tsv || exit $?
clsumtaxa \
--taxfile=$qc3nntaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$nontargettsv \
$tempfolder/$locus/community/community_qc3nn_nontarget.tsv || exit $?
fi

# Generate word cloud
if test -n $targettsv; then
clplotwordcloud \
--taxfile=$qc3nntaxonomy \
--append \
--numthreads=$THREADS \
--anemone \
$targettsv \
samplename || exit $?
ls $project$run\__*__$locus \
| xargs -I {} sh -c "mv -f {} $runfolder/$locus/ || exit $?"
fi

# Estimate DNA concentration and make community files of target taxa
if test -n $targettsv && test -n $standardtsv && test -s $stdconctable && test -s $solutionvoltable && test -s $watervoltable; then
clestimateconc \
--stdconctable=$stdconctable \
--stdtable=$standardtsv \
--solutionvoltable=$solutionvoltable \
--watervoltable=$watervoltable \
--numthreads=$THREADS \
$targettsv \
$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv || exit $?
clsumtaxa \
--taxfile=$qctaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv \
$tempfolder/$locus/community/community_qc_target_estimated.tsv || exit $?
clsumtaxa \
--taxfile=$qc3nntaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv \
$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv || exit $?
perl \
$projectfolder/../../combinecommunity.pl \
$tempfolder/$locus/community/community_qc_target.tsv \
$tempfolder/$locus/community/community_qc_target_estimated.tsv \
$runfolder/$locus/community_qc_target.tsv || exit $?
perl \
$projectfolder/../../combinecommunity.pl \
$tempfolder/$locus/community/community_qc3nn_target.tsv \
$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv \
$runfolder/$locus/community_qc3nn_target.tsv || exit $?
else
cp \
$tempfolder/$locus/community/community_qc_target.tsv \
$runfolder/$locus/community_qc_target.tsv || exit $?
cp \
$tempfolder/$locus/community/community_qc3nn_target.tsv \
$runfolder/$locus/community_qc3nn_target.tsv || exit $?
fi

# Estimate DNA concentration and make community files of nontarget taxa
if test -n $nontargettsv && test -n $standardtsv && test -s $stdconctable && test -s $solutionvoltable && test -s $watervoltable; then
clestimateconc \
--stdconctable=$stdconctable \
--stdtable=$standardtsv \
--solutionvoltable=$solutionvoltable \
--watervoltable=$watervoltable \
--numthreads=$THREADS \
$nontargettsv \
$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv || exit $?
clsumtaxa \
--taxfile=$qctaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv \
$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv || exit $?
clsumtaxa \
--taxfile=$qc3nntaxonomy \
--targetrank=sequence \
--otuseq=$fastafile \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv \
$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv || exit $?
perl \
$projectfolder/../../combinecommunity.pl \
$tempfolder/$locus/community/community_qc_nontarget.tsv \
$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv \
$runfolder/$locus/community_qc_nontarget.tsv || exit $?
perl \
$projectfolder/../../combinecommunity.pl \
$tempfolder/$locus/community/community_qc3nn_nontarget.tsv \
$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv \
$runfolder/$locus/community_qc3nn_nontarget.tsv || exit $?
else
cp \
$tempfolder/$locus/community/community_qc_nontarget.tsv \
$runfolder/$locus/community_qc_nontarget.tsv || exit $?
cp \
$tempfolder/$locus/community/community_qc3nn_nontarget.tsv \
$runfolder/$locus/community_qc3nn_nontarget.tsv || exit $?
fi

# Make standard community files
if test -n $standardtsv; then
echo -e 'query\tspecies' > $runfolder/$locus/taxonomy_standard.tsv
grep -P -o '^>\S+' $standardfasta | perl -ne 'chomp;s/^>//;print("$_\t$_\n")' >> $runfolder/$locus/taxonomy_standard.tsv
clsumtaxa \
--taxfile=$runfolder/$locus/taxonomy_standard.tsv \
--targetrank=sequence \
--otuseq=$standardfasta \
--taxnamereplace=disable \
--taxranknamereplace=disable \
--fuseotu=enable \
--numbering=disable \
--base62encoding=disable \
--sortkey=abundance \
--tableformat=column \
$standardtsv \
$runfolder/$locus/community_standard.tsv || exit $?
perl -i -npe 's/\tspecies\t/\tstandard\t/' $runfolder/$locus/community_standard.tsv || exit $?
rm -f $runfolder/$locus/taxonomy_standard.tsv || exit $?
fi

# Move results from work to home and delete temporary files
rm \
-rf \
$tempfolder/$locus/taxonomy/cachedb_$blastdb
mv \
-f \
$tempfolder/$locus/demultiplexed \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/concatenated \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/filtered \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/denoised \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/chimeraremoved1 \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/stdclustered \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/chimeraremoved2 \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/taxonomy \
$runfolder/$locus/
mv \
-f \
$tempfolder/$locus/community \
$runfolder/$locus/
rm \
-rf \
$tempfolder/$locus

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
unset 3nnidentdb
unset qctaxonomy
unset qc3nntaxonomy

done

rm \
-rf \
$tempfolder
