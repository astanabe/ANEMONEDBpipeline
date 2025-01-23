#!/bin/bash
#SBATCH --export=ALL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -o runlog_%j.txt

date --rfc-3339=seconds

echo HOSTNAME=$HOSTNAME
echo USER=$USER
echo HOME=$HOME
echo SHELL=$SHELL
echo LANG=$LANG
echo PATH=$PATH

team=`pwd | perl -ne 'chomp;@path=split("/");print($path[-3])'`
project=`pwd | perl -ne 'chomp;@path=split("/");print($path[-2])'`
run=`pwd | perl -ne 'chomp;@path=split("/");print($path[-1])'`
runfolder=`pwd`
tempfolder="/work/$USER/temp/$team/$project/$run"
mkdir -p "$tempfolder" || exit $?

echo $LINENO

cd "$runfolder" || exit $?

echo $LINENO

# Set number of processor cores used for computation
export THREADS=`grep -c processor /proc/cpuinfo`

echo $LINENO

if test -d "$runfolder/fastq_undemultiplexed"; then

echo $LINENO

# Copy FASTQ files to temporary firectory
if ! test -d "$tempfolder/fastq_undemultiplexed"; then
cp \
-RLf \
"$runfolder/fastq_undemultiplexed" \
"$tempfolder/" || exit $?
fi
# Demultiplex Type A (If you have undemultiplexed FASTQ files)
if ! test -d "$tempfolder/demultiplexed"; then
for p in `ls "$tempfolder/fastq_undemultiplexed/"*_I1_*.fastq "$tempfolder/fastq_undemultiplexed/"*_I1_*.fastq.gz "$tempfolder/fastq_undemultiplexed/"*_I1_*.fastq.bz2 "$tempfolder/fastq_undemultiplexed/"*_I1_*.fastq.xz | grep -P -o '.+(?=_I1_.+\.fastq)'`
do \
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
--append \
"$p"_R1_*.fastq \
"$p"_I1_*.fastq \
"$p"_I2_*.fastq \
"$p"_R2_*.fastq \
"$tempfolder/demultiplexed" || \
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
--append \
"$p"_R1_*.fastq.gz \
"$p"_I1_*.fastq.gz \
"$p"_I2_*.fastq.gz \
"$p"_R2_*.fastq.gz \
"$tempfolder/demultiplexed" || \
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
--append \
"$p"_R1_*.fastq.bz2 \
"$p"_I1_*.fastq.bz2 \
"$p"_I2_*.fastq.bz2 \
"$p"_R2_*.fastq.bz2 \
"$tempfolder/demultiplexed" || \
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
--append \
"$p"_R1_*.fastq.xz \
"$p"_I1_*.fastq.xz \
"$p"_I2_*.fastq.xz \
"$p"_R2_*.fastq.xz \
"$tempfolder/demultiplexed" || exit $?
done
fi
# Delete temporary FASTQ files
rm -rf "$tempfolder/fastq_undemultiplexed" || exit $?

echo $LINENO

elif test -d "$runfolder/fastq_demultiplexed"; then

echo $LINENO

# Copy FASTQ files to temporary firectory
if ! test -d "$tempfolder/fastq_demultiplexed"; then
cp \
-RLf \
"$runfolder/fastq_demultiplexed" \
"$tempfolder/" || exit $?
fi
# Demultiplex Type B (If FASTQ files have been already demultiplexed)
if ! test -d "$tempfolder/demultiplexed"; then
cltruncprimer \
--runname="$project$run" \
--forwardprimerfile="$runfolder/forwardprimer.fasta" \
--reverseprimerfile="$runfolder/reverseprimer.fasta" \
--truncateN=enable \
--outputmultihit=enable \
--index1file="$runfolder/index1.fasta" \
--index2file="$runfolder/index2.fasta" \
--seqnamestyle=illumina \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/fastq_demultiplexed" \
"$tempfolder/demultiplexed" || exit $?
fi
# Delete temporary FASTQ files
rm -rf "$tempfolder/fastq_demultiplexed" || exit $?

echo $LINENO

else
echo $LINENO
exit 1
fi

echo $LINENO

# Repeat analysis on each locus
for locus in `grep -P -o '^>\S+' "$runfolder/forwardprimer.fasta" | perl -npe 's/^>//'`
do \
# Set variables for each locus
if test "$locus" = 'MiFish' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='animals_12S_species'
includetaxa='\t(Hyperoartia|Myxini|Chondrichthyes|Actinopterygii|Coelacanthiformes|Dipnomorpha)\t'
elif test "$locus" = 'MiMammal' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='animals_12S_species'
includetaxa='\t(Mammalia)\t'
elif test "$locus" = 'MiBird' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='animals_12S_species'
includetaxa='\t(Aves)\t'
elif test "$locus" = 'MiDeca' ; then
ovlmode='OVL'
minlen=100
maxlen=250
referencedb='cdu12s'
blastdb='animals_12S_species'
includetaxa='\t(Decapoda)\t'
elif test "$locus" = 'Amph16S' ; then
ovlmode='OVL'
minlen=100
maxlen=400
referencedb='cdu16s'
blastdb='animals_16S_species'
includetaxa='\t(Amphibia)\t'
elif test "$locus" = 'MtInsects-16S' ; then
ovlmode='OVL'
minlen=100
maxlen=300
referencedb='cdu16s'
blastdb='animals_16S_species'
includetaxa='\t(Hexapoda)\t'
else
echo $LINENO
exit 1
fi
standardfasta="$runfolder/$locus/standard.fasta"
stdconctable="$runfolder/$locus/stdconctable.tsv"
solutionvoltable="$runfolder/$locus/solutionvoltable.tsv"
watervoltable="$runfolder/$locus/watervoltable.tsv"

echo $LINENO

if ! test -d "$tempfolder/$locus/demultiplexed" && ! test -d "$runfolder/$locus/demultiplexed"; then
# Make directory for each locus
mkdir -p \
"$tempfolder/$locus/demultiplexed" || exit $?
mkdir -p \
"$runfolder/$locus" || exit $?
# Move demultiplexed files
ls "$tempfolder/demultiplexed/$project$run"__*__"$locus".*.fastq.xz \
| xargs -I {} sh -c "mv -f {} \"$tempfolder/$locus/demultiplexed/\" || exit $?"
fi

echo $LINENO

# Concatenate pairs
if ! test -d "$tempfolder/$locus/concatenated" && ! test -d "$runfolder/$locus/concatenated"; then
clconcatpairv \
--mode="$ovlmode" \
--output=directory \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/$locus/demultiplexed" \
"$tempfolder/$locus/concatenated" || exit $?
fi

echo $LINENO

# Apply filtering out low quality sequences
if ! test -d "$tempfolder/$locus/filtered" && ! test -d "$runfolder/$locus/filtered"; then
clfilterseqv \
--maxqual=41 \
--minlen="$minlen" \
--maxlen="$maxlen" \
--maxnee=2.0 \
--maxnNs=0 \
--output=directory \
--compress=xz \
--numthreads="$THREADS" \
"$tempfolder/$locus/concatenated" \
"$tempfolder/$locus/filtered" || exit $?
fi

echo $LINENO

# Denoise using DADA2
if ! test -d "$tempfolder/$locus/denoised" && ! test -d "$runfolder/$locus/denoised"; then
cldenoiseseqd \
--pool=pseudo \
--numthreads="$THREADS" \
"$tempfolder/$locus/filtered" \
"$tempfolder/$locus/denoised" || exit $?
fi

echo $LINENO

# Remove chimeras using UCHIME3
if ! test -d "$tempfolder/$locus/chimeraremoved1" && ! test -d "$runfolder/$locus/chimeraremoved1"; then
clremovechimev \
--mode=denovo \
--uchimedenovo=3 \
"$tempfolder/$locus/denoised" \
"$tempfolder/$locus/chimeraremoved1" || exit $?
fi

echo $LINENO

if test -s "$standardfasta"; then
# Cluster internal standard sequences to otus
if ! test -d "$tempfolder/$locus/stdclustered" && ! test -d "$runfolder/$locus/stdclustered"; then
clclusterstdv \
--standardseq="$standardfasta" \
--numthreads="$THREADS" \
"$tempfolder/$locus/chimeraremoved1" \
"$tempfolder/$locus/stdclustered" || exit $?
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
"$tempfolder/$locus/chimeraremoved2" || exit $?
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
"$tempfolder/$locus/chimeraremoved2" || exit $?
fi
fi
fi

echo $LINENO

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
echo $LINENO
exit 1
fi

echo $LINENO

if ! test -d "$tempfolder/$locus/taxonomy" && ! test -d "$runfolder/$locus/taxonomy"; then

echo $LINENO

# Make directory for output
mkdir -p \
"$tempfolder/$locus/taxonomy" || exit $?

echo $LINENO

# Set variables
if test -e "$runfolder/../../../$locus/qc_$blastdb.identdb"; then
qcidentdb="$runfolder/../../../$locus/qc_$blastdb.identdb"
fi
if test -e "$runfolder/../../../$locus/3nn_$blastdb.identdb"; then
threennidentdb="$runfolder/../../../$locus/3nn_$blastdb.identdb"
fi
if ! test -e "$runfolder/../../../$locus/qc_$blastdb.identdb" && ! test -e "$runfolder/../../../$locus/3nn_$blastdb.identdb"; then
mkdir -p \
"$runfolder/../../../$locus" || exit $?
fi

echo $LINENO

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
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $?
else
clmakecachedb \
--identdb="$qcidentdb" \
--blastdb="$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $?
fi
else
if test -s "$standardfasta"; then
clmakecachedb \
--blastdb="$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $?
else
clmakecachedb \
--blastdb="$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $?
fi
fi
fi

echo $LINENO

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
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $?
else
clidentseq \
--method=QC \
--identdb="$qcidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $?
fi
else
if test -s "$standardfasta"; then
clidentseq \
--method=QC \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $?
else
clidentseq \
--method=QC \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" || exit $?
fi
fi
fi

echo $LINENO

if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv"; then
classigntax \
--taxdb="$blastdb" \
--maxpopposer=0.05 \
--minsoratio=19 \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" || exit $?
fi

echo $LINENO

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
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $?
else
clidentseq \
--method=3,95% \
--identdb="$threennidentdb" \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $?
fi
else
if test -s "$standardfasta"; then
clidentseq \
--method=3,95% \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--ignoreotuseq="$standardfasta" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $?
else
clidentseq \
--method=3,95% \
--blastdb="$tempfolder/$locus/taxonomy/cachedb_$blastdb" \
--numthreads="$THREADS" \
"$fastafile" \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" || exit $?
fi
fi
fi

echo $LINENO

if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv"; then
classigntax \
--taxdb="$blastdb" \
--minnsupporter=3 \
--maxpopposer=0.05 \
--minsoratio=19 \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" \
"$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv" || exit $?
fi

echo $LINENO

# Merge taxonomic assignment results
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv"; then
clmergeassign \
--preferlower \
--priority=descend \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_3nn_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv" || exit $?
fi

echo $LINENO

# Fill blanks in taxonomy
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv; then
clfillassign \
--fullfill=enable \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv || exit $?
fi
if ! test -e "$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv; then
clfillassign \
--fullfill=enable \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb.tsv" \
"$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv || exit $?
fi

echo $LINENO

# Update identdb
clmakeidentdb \
--append \
"$tempfolder/$locus/taxonomy/neighborhoods_qc_$blastdb.txt" \
"$runfolder/../../../$locus/qc_$blastdb.identdb" &
clmakeidentdb \
--append \
"$tempfolder/$locus/taxonomy/neighborhoods_3nn_$blastdb.txt" \
"$runfolder/../../../$locus/3nn_$blastdb.identdb" &
wait || exit $?

echo $LINENO

fi

echo $LINENO

# Set variables
qctaxonomy="$tempfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv
qc3nntaxonomy="$tempfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv

echo $LINENO

if ! test -d "$tempfolder/$locus/community" && ! test -d "$runfolder/$locus/community"; then

echo $LINENO

# Make directory for output
mkdir -p \
"$tempfolder/$locus/community" || exit $?

echo $LINENO

# Combine replicates
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_all.tsv"; then
if test -e "$runfolder/$locus/replicatelist.txt"; then
clfiltersum \
--replicatelist="$runfolder/$locus/replicatelist.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_all.tsv" || exit $?
else
cp \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_all.tsv" || exit $?
fi
fi
alltsv="$tempfolder/$locus/community/sample_otu_matrix_all.tsv"

echo $LINENO

# Extract internal standard OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_standard.tsv"; then
if test -s "$standardfasta"; then
clfiltersum \
--otuseq="$standardfasta" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_standard.tsv" || exit $?
fi
fi

echo $LINENO

# Pick target taxa
grep -P "$includetaxa" \
"$qctaxonomy" \
| cut -f 1 \
> "$tempfolder/$locus/targetotu.txt"

echo $LINENO

# Pick nontarget taxa
grep -P -v "$includetaxa" \
"$qctaxonomy" \
| cut -f 1 \
> "$tempfolder/$locus/nontargetotu.txt"

echo $LINENO

# Delete nontarget OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_target.tsv"; then
if test -s "$tempfolder/$locus/targetotu.txt"; then
clfiltersum \
--otulist="$tempfolder/$locus/targetotu.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_target.tsv" || exit $?
fi
fi

echo $LINENO

# Extract nontarget OTUs
if ! test -e "$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv"; then
if test -s "$tempfolder/$locus/nontargetotu.txt"; then
clfiltersum \
--otulist="$tempfolder/$locus/nontargetotu.txt" \
"$alltsv" \
"$tempfolder/$locus/community/sample_otu_matrix_nontarget.tsv" || exit $?
fi
fi

echo $LINENO

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

echo $LINENO

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
"$tempfolder/$locus/community/community_qc_target.tsv" || exit $?
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
"$tempfolder/$locus/community/community_qc3nn_target.tsv" || exit $?
fi
fi

echo $LINENO

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
"$tempfolder/$locus/community/community_qc_nontarget.tsv" || exit $?
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
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" || exit $?
fi
fi

echo $LINENO

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
"$tempfolder/$locus/community/sample_otu_matrix_target_estimated.tsv" || exit $?
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
"$tempfolder/$locus/community/community_qc_target_estimated.tsv" || exit $?
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
"$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv" || exit $?
fi
if ! test -e "$tempfolder/$locus/community_qc_target.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc_target.tsv" \
"$tempfolder/$locus/community/community_qc_target_estimated.tsv" \
"$tempfolder/$locus/community_qc_target.tsv" || exit $?
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_target.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc3nn_target.tsv" \
"$tempfolder/$locus/community/community_qc3nn_target_estimated.tsv" \
"$tempfolder/$locus/community_qc3nn_target.tsv" || exit $?
fi
else
if ! test -e "$tempfolder/$locus/community_qc_target.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc_target.tsv" \
"$tempfolder/$locus/community_qc_target.tsv" || exit $?
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_target.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc3nn_target.tsv" \
"$tempfolder/$locus/community_qc3nn_target.tsv" || exit $?
fi
fi

echo $LINENO

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
"$tempfolder/$locus/community/sample_otu_matrix_nontarget_estimated.tsv" || exit $?
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
"$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv" || exit $?
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
"$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv" || exit $?
fi
if ! test -e "$tempfolder/$locus/community_qc_nontarget.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc_nontarget.tsv" \
"$tempfolder/$locus/community/community_qc_nontarget_estimated.tsv" \
"$tempfolder/$locus/community_qc_nontarget.tsv" || exit $?
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_nontarget.tsv"; then
perl \
"$runfolder/../../../combinecommunity.pl" \
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" \
"$tempfolder/$locus/community/community_qc3nn_nontarget_estimated.tsv" \
"$tempfolder/$locus/community_qc3nn_nontarget.tsv" || exit $?
fi
else
if ! test -e "$tempfolder/$locus/community_qc_nontarget.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc_nontarget.tsv" \
"$tempfolder/$locus/community_qc_nontarget.tsv" || exit $?
fi
if ! test -e "$tempfolder/$locus/community_qc3nn_nontarget.tsv"; then
cp \
"$tempfolder/$locus/community/community_qc3nn_nontarget.tsv" \
"$tempfolder/$locus/community_qc3nn_nontarget.tsv" || exit $?
fi
fi

echo $LINENO

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
"$tempfolder/$locus/community_standard.tsv" || exit $?
perl -i -npe 's/\tspecies\t/\tstandard\t/' "$tempfolder/$locus/community_standard.tsv" || exit $?
rm -f "$tempfolder/$locus/taxonomy_standard.tsv" || exit $?
fi
fi

echo $LINENO

# Generate word cloud
if test -n "$targettsv"; then
cd "$tempfolder/$locus"
clplotwordcloud \
--taxfile="$qc3nntaxonomy" \
--append \
--numthreads="$THREADS" \
--anemone \
"$targettsv" \
samplename || exit $?
cd "$runfolder"
fi

echo $LINENO

fi

echo $LINENO

# Move results from work to home and delete temporary files
rm \
-rf \
"$tempfolder/$locus/taxonomy/cachedb_$blastdb" || exit $?
cp \
-RLf \
"$tempfolder/$locus" \
"$runfolder/" || exit $?
rm \
-rf \
"$tempfolder/$locus" || exit $?

echo $LINENO

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
unset standardtsv
unset targettsv
unset nontargettsv

echo $LINENO

done

echo $LINENO

rm \
-rf \
"$tempfolder" || exit $?

date --rfc-3339=seconds
