if [ $# -lt 1 ]; then echo $#; echo "enter receptor"; exit 0; fi
if [ ! -d Binders ]; then mkdir Binders; fi
if [ ! -d No_Binders ]; then mkdir No_Binders; fi
cat binders.txt | while read line; do filename=${1}-${line}.recint.score.csv; if [ -e $filename ] ; then mv $filename Binders/.; fi; done
mv *.recint.score.csv No_Binders/.
