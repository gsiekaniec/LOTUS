#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import uuid
from python_scripts.path_modification import true_stem


#############
# true stem #
#############

file1 = '/a/b/file1.vcf'
file2 = '/a/b/file2.txt.svg'
file3 = '/a/b/file3'
no_file = ''

################################################################################################################
################################################# Tests ########################################################
################################################################################################################


#############
# true stem #
#############

def test_true_stem():
	assert true_stem(file1) == 'file1'
	assert true_stem(file2) == 'file2'
	assert true_stem(file3) == 'file3'
	assert true_stem(no_file) == ''
