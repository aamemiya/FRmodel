#!/bin/sh

invdir=./inv

lat_dirs=`ls -x -d ???-???`
for lat_dir in $lat_dirs ;do

mkdir $invdir/$lat_dir

lon_files=`ls -x $lat_dir/????-????`

for lon_file in $lon_files ;do

file_inv=$invdir/$lon_file

cp $lon_file ./work/temp
file=work/temp
file_new=work/temp_new

split -l 60 $file
temp_files=`ls -x -r x??`
for temp_file in $temp_files ; do
 cat $temp_file >> $file_new
done
cp $file_new $file_inv
rm $file_new
rm  $temp_files

done
done