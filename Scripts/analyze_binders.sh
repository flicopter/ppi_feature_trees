#!/bin/bash

if [ "$#" -lt 1 ]
    then echo "Usage: analyze_binders.sh <DIR_HAVING_BINDERS_NOBINDERS>"
    exit 1
fi

SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
RECDIR="$1"

pushd $RECDIR
Rscript $SCRIPTDIR/analyze_binders.R Binders No_Binders > classify_results.txt
popd
