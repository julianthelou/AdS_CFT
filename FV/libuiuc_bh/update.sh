#!/bin/bash
#
# This is a small script to copy the UIUC_BH source code to the local application.
# We actually manage the UIUC_BH code at https://bitbucket.org/svek/uiucinitialdata
# and want to keep the code there. ExaHyPE lacks a proper way of managing external
# repos...
# SvenK, 2018-07-13

UIUC_BlackHole="../../../ExternalLibraries/UIUC_BlackHole"

[[ -e $UIUC_BlackHole ]] || { echo "Cannot find $UIUC_BlackHole"; exit -1; }

uiuc="$UIUC_BlackHole/uiuc"

[[ -e $uiuc ]] || $UIUC_BlackHole/download-uiuc.sh || { echo "Cannot download UIUC"; exit -1; }

for remotef in $uiuc/standalone/*; do
   localf=$(basename $remotef)
   if [[ -e $localf ]] && ! diff $localf $remotef ; then
       echo "$localf has local changes, compared to $remotef , please merge back!"
       exit -1
   else
       cp -v $remotef $localf
   fi
done

echo "Done"
