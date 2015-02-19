#!/bin/bash

if [ "$#" -lt 1 ]
    then echo "Usage: analyze_binders.sh <DIR_HAVING_BINDERS_NOBINDERS>"
    exit 1
fi

SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
RECDIR="$1"

if [ ! -d "$RECDIR/Binders" ]; then 
    echo "The directory '$RECDIR/Binders' does not exist."
    exit 1
fi
if [ ! -d "$RECDIR/No_Binders" ]; then
    echo "The directory '$RECDIR/No_Binders' does not exist."
    exit 1
fi

pushd $RECDIR
Rscript $SCRIPTDIR/analyze_binders.R Binders No_Binders > classify_results.txt
echo "Classification results written to '$RECDIR/classify_results.txt'."
popd
