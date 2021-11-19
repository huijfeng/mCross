##get transfac
if [ $# -ne 5 ]
then
	echo "Usage: sh motifcluster.sh motifDir resultPrefix resultDir mCrossDir StampscoreFile"
	exit 1
fi
files="$1/*.mat"
f=`ls $files`
b=$2
if [ ! -d "$3" ];then
	mkdir $3
fi
echo $b
echo ${f[@]}
for i in ${f[@]}
do
echo $i
perl $4/summarizeMotif.pl $i >>  $3/$b.consensus.txt
linenum=$(wc -l < $i)
# remove file that has trunctaed transfac file format to avoid errors in the following perl script
if [ "$linenum" -gt 11 ]; then
echo $linenum
perl $4/formatMotifs.pl $i $3/$b.transfac 0.2 >> $3/$b.log
fi
done

stamp -tf $3/$b.transfac  -sd $5   -cc PCC -align SWU  -forwardonly -printpairwise -chp  -out $3/$b > $3/$b.summary.txt
awk '/Pairwise alignment scores:/,/Alignments Finished/' $3/$b.summary.txt | egrep -v "Pairwise alignment scores:|Alignments Finished" | grep -v '^$' > $3/${b}.pairwisealignmentdist.txt
awk '/Calinski & Harabasz:/,/Tree Built/' $3/$b.summary.txt | egrep -v "Calinski & Harabasz:|Tree Built" |  head -n -1 | sed -e 's/^[ \t]//' > $3/$b.CHidx.txt
sort $3/$b.consensus.txt | uniq | awk '
{
        if(!match($2,",") && $1!="name") 
        {
                print ">" $1; print $2
        } 
        if(match($2,",") && $1!="name") 
        {
                n=split($2, a,","); 
                print ">" $1; print a[1]
        }
}' >  $3/$b.consensus.seq.txt


