#!/bin/sh

for ilat in `seq 1 180`;do
 clat=`printf %03d $ilat`
 ivlatl=`expr 90 - $ilat`
 ivlatr=`expr 91 - $ilat`
 cvlatl=`printf %03d $ivlatl | cut -c 2-3`
 cvlatr=`printf %03d $ivlatr | cut -c 2-3`
 [ $ivlatl -ge 0 ] && cvlatl=${cvlatl}N || cvlatl=${cvlatl}S 
 [ $ivlatr -ge 0 ] && cvlatr=${cvlatr}N || cvlatr=${cvlatr}S
 mv dir$clat $cvlatl-$cvlatr
done