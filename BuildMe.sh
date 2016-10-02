#!/bin/bash
CreateBuildDir()
{
  if [ -d build ]; 
  then 
    rm -r build
  fi
  mkdir build
}

CreateMakeFile()
{
  echo "========================================================================="
  echo "=                                                                       ="
  echo "=                    Build Makefile for compilation                     ="
  echo "=                                                                       ="
  echo "========================================================================="
  echo "=->Creating the new makefile                                            ="
  echo "========================================================================="
  cmake ..
}

Compile()
{
  echo "========================================================================="
  echo "=                                                                       ="
  echo "=                              Compiling!                               ="
  echo "=                                                                       ="
  echo "========================================================================="
  make 
  make install
  echo "========================================================================="
  echo "=->Done!!!                                                              ="
  echo "========================================================================="
}

CreateBuildDir
cd build
CreateMakeFile
Compile

