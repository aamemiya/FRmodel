#!/bin/sh

split -l 1296000 ETOPO1_Ice_c_int.xyz

echo 'first split done.'

files=`ls -x ???`
  ilat=1
  for file in $files ; do
   clat=`printf %03d $ilat`
   mv $file $clat
   mkdir dirm$clat
   ilat=`expr $ilat + 1`
  done 

echo 'first rename done.'

for ilat in `seq 1 180` ; do
echo 'ilat = '$ilat' /  180'
 clat=`printf %03d $ilat`
 mkdir dir$clat
 cd dir$clat
 split -l 21600 ../$clat
  files=`ls -x ???`
  imin=1
  for file in $files ; do
   cmin=`printf %02d $imin`
   mv $file m$cmin
   mkdir dirm$cmin
   imin=`expr $imin + 1`
  done 
   for imin in `seq 1 60`;do
    cmin=`printf %02d $imin`
    cd dirm$cmin
    split -l 60 ../m$cmin
    files=`ls -x ???`
    ilon=1
    for file in $files ; do
     clon=`printf %03d $ilon`
     mv $file m$cmin$clon
     cat m$cmin$clon >> ../$clon
     ilon=`expr $ilon + 1`
    done 
    cd ..
   done
   cd ..
done

#SRC=ETOPO1_Ice_c_int.xyz
#for ilat in `seq 1 180`;do
# mkdir `printf %03d $ilat`
# head $SRC 
#for ilon in `seq 1 360`;do
 
 
#done
#done