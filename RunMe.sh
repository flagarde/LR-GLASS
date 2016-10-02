#!/bin/bash 
shopt -s nullglob
for f in $1*.txt
do
  filename=$(basename "$f")
  extension="${filename##*.}"
  filename="${filename%.*}"
  sep=_
  name=$filename$sep$2
  a="$(echo $1|sed 's#/#\_#g')"
  b=Results/
  if [ ! -d "$a$b" ]; then
  	mkdir $a$b
  fi
  if [ "$#" -ne 1 ]; 
	then
		echo "Launch lrGlass for $f - Creating $name.root in $a$b "
  		./LRGRPC $name.root $f
		mv $name.root $b$a$name.root
		mv $filename  $b$a$filename
	else 
		echo "Launch lrGlass for $f - Creating $filename.root in $a$b"
		./LRGRPC $filename.root $f
                mv $filename.root $b$a$filename.root
		mv $filename $b$a$filename
  fi 
done
