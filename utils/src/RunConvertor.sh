lndir=$LNDIR

#source root
source $lndir/fitter/lh-soft.sh
#set code directory
fitdir=$FITDIR/utils/src/
echo $fitdir

#set in/out file names
infile=$1
outfile=$2

echo $outfile
#.root in file has to be in current directory
mydir=$(pwd)
infile=$mydir/$infile

#make a directory for your output
if [ ! -d fitLLdata ]
    then mkdir fitLLdata
fi

#set out filename
outfile=$mydir/fitLLdata/$outfile
#echo $outfile

echo $infile
echo $outfile

#link code to current directory
ln -s $fitdir/treeConvert.* $mydir

root -l <<EOF
.L treeConvert.cc
treeConvert("$infile","default","$outfile","dpTT","truedpTT","selmu_mom","selmu_truemom")
.q
EOF
echo "Complete"

if [ $mydir != $fitdir ]; then
rm treeConvert*
fi