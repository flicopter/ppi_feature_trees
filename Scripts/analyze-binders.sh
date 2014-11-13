#!/bin/bash

### usage
# bash analyze-binders.sh RECEPTOR_PROTEIN_DIR PARTNER_PDB [CONFIG_FILE]
GITHOME=$(git rev-parse --show-toplevel)
if [ $# -lt 3 ]; then CONFIG_FILE=$GITHOME/Configs/$USER.cfg; else CONFIG_FILE=$3; fi
source $CONFIG_FILE 

ANCHOR=$PWD
RECDIR=$1
RECNAME=$( echo $1 | sed 's/\/$//' | awk -F '/' '{print $NF}')
PARTNER_PDB=$( readlink -f $2 )

cd $RECDIR/
rm -f *.pdb
cp $PDBDIR/*.pdb .
grep "^$RECNAME" $REFFILE | awk '{print $2}' > binders.txt
for i in $RECNAME-*.out
  do python $SCRIPTDIR/contactmap.py $i $PARTNER_PDB
done
bash $SCRIPTDIR/movefile.sh $RECNAME

Rscript $SCRIPTDIR/analyze-binders.R Binders No_Binders > classify_results.txt

rm -f *.pdb
cd $ANCHOR 
