#!/bin/sh

londirs=`ls -x -d  ???-???`

for dir in $londirs ; do
for ilon in `seq 1 360`;do
 clon=`printf %03d $ilon` 
 ilonl=`expr $ilon - 181`
 ilonr=`expr $ilon - 180`
 clonl=`printf %04d $ilonl | cut -c 2-4`  ### abs value
 clonr=`printf %04d $ilonr | cut -c 2-4`  
 [ $ilonl -ge 0 ] && clonl=${clonl}E || clonl=${clonl}W
 [ $ilonr -ge 0 ] && clonr=${clonr}E || clonr=${clonr}W
 mv $dir/$clon $dir/$clonl-$clonr  
#echo $clon $clonl-$clonr 
done
done