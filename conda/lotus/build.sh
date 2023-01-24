#!/bin/bash

mkdir -p $PREFIX/bin
cp -r python_scripts $PREFIX/bin
cp lotus.py $PREFIX/bin
chmod +x $PREFIX/bin/*
