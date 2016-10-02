#!/bin/bash 
VerifyResultsExists()
{
  if [ -d Results ];
    then 
      echo "Copying into Results"
    else 
      mkdir ./Results
  fi
}

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
    echo "Please give me the folder where you have your datacards"
    echo "Syntax : ./RunMe.sh path_to_datacards {optional string}"
  else
  shopt -s nullglob
  for f in $1*.txt
  do
    filename=$(basename "$f")
    filename="${filename%.*}"
    sep=_
    if [ $# -ne 1 ]; 
      then name=$filename$sep$2
      else name=$filename
    fi
    c=./Results/
    if [ "$#" -ne 1 ]; 
	  then
		  echo "Launch lrGlass for $f - Creating $name.root in $name "
  	  ./LRGRPC $name.root $f
  	  VerifyResultsExists
  	  #mkdir $c$a
  	  mv $filename  $c$name
		  mv $name.root $c$name.root
		else 
		  echo "Launch lrGlass for $f - Creating $filename.root in $name"
		  ./LRGRPC $filename.root $f
		  VerifyResultsExists
		  #mkdir $c$a
		  mv $filename $c$name
      mv $filename.root $c$name.root
    fi 
  done
fi
