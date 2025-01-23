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

ls -d "$runfolder/$locus/$project$run"__*__"$locus" \
| xargs -L 1 -P $THREADS rm -rf
ls -d "$tempfolder/$locus/$project$run"__*__"$locus" \
| xargs -L 1 -P $THREADS rm -rf

echo $LINENO

if ! test -d "$tempfolder/$locus"; then
# Make directory for each locus
mkdir -p \
"$tempfolder/$locus" || exit $?
fi

echo $LINENO

# Set variables
if test -s "$runfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv; then
qc3nntaxonomy="$runfolder/$locus/taxonomy/taxonomy_qc3nn_$blastdb"_filled.tsv
else
echo $LINENO
exit 1
fi

echo $LINENO

# Set variables
if test -s "$runfolder/$locus/community/sample_otu_matrix_target.tsv"; then
targettsv="$runfolder/$locus/community/sample_otu_matrix_target.tsv"
else
echo $LINENO
exit 1
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
