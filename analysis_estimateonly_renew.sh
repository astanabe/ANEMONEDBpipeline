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

rm -rf \
"$runfolder/$locus/community/"*_estimated.tsv \
"$tempfolder/$locus/community/"*_estimated.tsv \
"$runfolder/$locus/community_"*.tsv \
"$tempfolder/$locus/community_"*.tsv || exit $?

echo $LINENO

if ! test -d "$tempfolder/$locus/community"; then
# Make directory for each locus
mkdir -p \
"$tempfolder/$locus" || exit $?
# Copy community files
if test -d "$runfolder/$locus/community"; then
cp \
-RLf \
"$runfolder/$locus/community" \
"$tempfolder/$locus/" || exit $?
else
echo $LINENO
exit 1
fi
fi

echo $LINENO

# Set variables
if test -d "$runfolder/$locus/chimeraremoved2"; then
fastafile="$runfolder/$locus/chimeraremoved2/nonchimeras.fasta"
elif test -d "$runfolder/$locus/stdclustered"; then
fastafile="$runfolder/$locus/stdclustered/stdclustered.fasta"
elif test -d "$runfolder/$locus/chimeraremoved1"; then
fastafile="$runfolder/$locus/chimeraremoved1/nonchimeras.fasta"
else
echo $LINENO
exit 1
fi

echo $LINENO

# Set variables
if test -s "$runfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv; then
qctaxonomy="$runfolder/$locus/taxonomy/taxonomy_qc_$blastdb"_filled.tsv
else
echo $LINENO
exit 1
fi
if test -s "$runfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv; then
qc3nntaxonomy="$runfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv
else
echo $LINENO
exit 1
fi

echo $LINENO

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

# Move results from work to home and delete temporary files
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
