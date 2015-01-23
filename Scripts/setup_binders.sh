#!/bin/bash

### usage
# bash setup_binders.sh RECEPTOR_PROTEIN_DIR PARTNER_PDB [CONFIG_FILE]
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
if [ $# -lt 3 ]
then
  source $SCRIPTDIR/find_config.sh
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
RECNAME=$( basename $1 )
PARTNER_PDB=$( readlink -f $2 )

cd $RECDIR/
rm -f *.pdb
cp $PDBDIR/*.pdb .
grep "^$RECNAME" $REFFILE | awk '{print $2}' > binders.txt
for i in $RECNAME-*.out
  do python $SCRIPTDIR/contactmap.py $i $PARTNER_PDB
done
bash $SCRIPTDIR/movefile.sh $RECNAME

rm -f *.pdb
cd $ANCHOR 
