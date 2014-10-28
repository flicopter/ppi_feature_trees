#!/bin/bash

### usage
# bash analyze-binders.sh CONFIG_FILE RECEPTOR_PROTEIN_DIR PARTNER_PDB
source $1 

ANCHOR=$PWD
RECDIR=$2
RECNAME=$( echo $2 | sed 's/\/$//' | awk -F '/' '{print $NF}')
PARTNER_PDB=$( readlink -f $3 )

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
