#!/bin/bash

### usage
# bash analyze-binders.sh RECEPTOR_PROTEIN_DIR PARTNER_PDB [CONFIG_FILE]
GITHOME=$(git rev-parse --show-toplevel)
if [ $# -lt 3 ]
then
  if [ -e "$HOME/.ppi_feature_trees.cfg" ]
  then
    CONFIG_FILE="$HOME/.ppi_feature_trees.cfg"
  else
    CONFIG_FILE=$GITHOME/Configs/$USER.cfg
  fi
else 
  CONFIG_FILE=$3
fi
source $CONFIG_FILE 

if [ ! -d "$SCRIPTDIR" ]; then
  echo "Directory SCRIPTDIR ($SCRIPTDIR) does not exist."
  exit 1
fi

if [ ! -d "$PDBDIR" ]; then
  echo "Directory PDBDIR ($PDBDIR) does not exist."
  exit 1
fi

if [ ! -e "$REFFILE" ]; then
  echo "File REFFILE ($REFFILE) does not exist."
  exit 1
fi

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
