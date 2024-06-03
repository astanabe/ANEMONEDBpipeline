team=`pwd | perl -ne 'chomp;@path=split("/");print($path[-3])'`
project=`pwd | perl -ne 'chomp;@path=split("/");print($path[-2])'`
run=`pwd | perl -ne 'chomp;@path=split("/");print($path[-1])'`
runfolder=`pwd`

# make samplenames.txt
tail \
-n +2 \
community/sample_otu_matrix_target.tsv \
| grep -P -o '^\S+' \
> samplenames.txt

# compose dist files
perl \
$runfolder/../../../composedistfiles.pl \
Metadata.tsv \
samplenames.txt \
/home/$USER/dist
