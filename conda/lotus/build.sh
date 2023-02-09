#!/bin/bash

echo "Start"
mkdir -p $PREFIX/bin
echo ${PREFIX}/bin
cp -r python_scripts $PREFIX/bin
cp lotus.py $PREFIX/bin/lotus
chmod +x $PREFIX/bin/*
echo "End"
