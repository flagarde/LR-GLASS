#!/bin/bash 
shopt -s nullglob
for f in $1*.txt
do
  filename=$(basename "$f")
  extension="${filename##*.}"
  filename="${filename%.*}"
	echo "Launch lrGlass for $f - Creating $filename.root "
        ./lrGlass $filename.root $f 
done
