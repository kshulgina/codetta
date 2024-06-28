#!/bin/bash
# Setup script for Codetta
# Yekaterina Shulgina - Nov 2021

# check that Codetta required Python packages are present
python3 check_requirements.py
if [ $? -ne 0 ]
then
    exit
fi

echo Checking wget
if ! command -v wget &> /dev/null
then
    echo "ERROR: wget could not be found"
    exit
fi

echo Checking gzip
if ! command -v gzip &> /dev/null
then
    echo "ERROR: gzip could not be found"
    exit
fi

# untar HMMER distribution
tar xf hmmer.tar.gz

# compile HMMER
echo Installing local HMMER
cd hmmer-3.3.2
pwd | xargs -I {} ./configure --prefix={}
if [ $? -ne 0 ]
then
    'ERROR with HMMER installation'
    exit
fi

make
if [ $? -ne 0 ]
then
    'ERROR with HMMER installation'
    exit
fi

make install
if [ $? -ne 0 ]
then
    'ERROR with HMMER installation'
    exit
fi

# compile Easel
echo Installing local Easel
cd easel; make install 
if [ $? -ne 0 ]
then
    'ERROR with Easel installation'
    exit
fi

cd ../..

# clean up
# rm hmmer-3.3.2.tar.gz 

echo Codetta setup complete!
