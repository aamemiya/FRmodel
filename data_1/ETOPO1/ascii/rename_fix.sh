#!/bin/sh

files=`ls -x ??N-??N`
for file in $files ; do
 ista=`echo $file | cut -c 1-2`
 iend=`echo $file | cut -c 5-6`
 newfile=${iend}S-${ista}S
 mv $file $newfile
done


files=`ls -x ??S-??S`
for file in $files ; do
 ista=`echo $file | cut -c 1-2`
 iend=`echo $file | cut -c 5-6`
 newfile=${iend}N-${ista}N
 mv $file $newfile
done

