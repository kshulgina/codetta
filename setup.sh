#!/bin/bash
# Setup script for Codetta
# Yekaterina Shulgina - Nov 2021

# check that Codetta required Python packages are present
python check_requirements.py

# untar HMMER distribution
tar xf hmmer.tar.gz

# compile HMMER
echo Installing local HMMER
cd hmmer-3.3.2
pwd | xargs -I {} ./configure --prefix={}
make --quiet
make install --quiet

# compile Easel
echo Installing local Easel
cd easel; make install --quiet
cd ../..

# download Pfam database built with --enone option
# TBD

# clean up
# rm hmmer-3.3.2.tar.gz 