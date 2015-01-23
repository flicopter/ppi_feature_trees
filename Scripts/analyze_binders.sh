#!/bin/bash

SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

Rscript $SCRIPTDIR/analyze_binders.R Binders No_Binders > classify_results.txt
