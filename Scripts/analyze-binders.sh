#!/bin/bash

### usage
# bash analyze-binders.R RECEPTOR_PROTEIN_DIR PARTNER_PDB

PDBDIR=/home/yuri/Depot/PDB/benchmark4_structures/unbound/
REFFILE=/home/yuri/Depot/rel/bench4_mono_bu_correct.tsv
SCRIPTDIR=/home/yuri/src/ppi_feature_trees/Scripts/
ANCHOR=$PWD
RECNAME=$( echo $1 | sed 's/\///' )
PDBID=$( echo $1 | sed 's/\///' ) 
PARTNER_PDB=$2

cd $RECNAME/
rm -f *.pdb
cp $PDBDIR/*.pdb .
grep "^$PDBID" $REFFILE | awk '{print $2}' > binders.txt
for i in $RECNAME-*.out
  do python $SCRIPTDIR/contactmap.py $i $PARTNER_PDB
done
bash $SCRIPTDIR/movefile.sh $PDBID

Rscript $SCRIPTDIR/analyze-binders.R Binders No_Binders > classify_results.txt

rm -f *.pdb
cd $ANCHOR 
