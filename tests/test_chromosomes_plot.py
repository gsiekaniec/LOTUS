#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import re
import uuid
import logging
import pytest
from pathlib import Path
from python_scripts.path_modification import true_stem
from python_scripts.chromosomes_plot import extract_num, count_impacted_genes, get_value_for_chr, get_plot, scientific, cross_product, sorted_nicely, create_plot, create_chromosomes_plot, COLOR_G_STAINING

logger = logging.getLogger('This is not a logger')


###############
# extract num #
###############

pattern = re.compile(r'chr(\d+)')
pattern2 = re.compile(r'bob(\d+)')


########################
# count impacted genes #
########################

positions = 10
span = 6
genes_positions = [(2,3), (6,8), (12,16), (17,18)]
genes_positions2 = [(2,3), (15,16), (17,18)]


#####################
# get value for chr #
#####################

chr = 'chrP'
length = 100
span = 5
step = 5
positions_for_chr = {'chrP' : [(2,3), (15,16), (17,18)]}
x_value = [0, 2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97, 98, 103]
y_value = [1, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
x_value_0 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95]
y_value_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


############
# get plot #
############

chr_size = {'chrP' : 100, 'chrG' : 100}


##############
# scientific #
##############

number1 = 126531123671
number2 = 123
number3 = 0


#################
# cross product #
#################

list_values = [1,4,18,25]
max1 = 50
max2 = 100
list_values2 = [2.0,8.0,36.0,50.0]


#################
# sorted nicely #
#################

list_not_ok = ['chr1', 'chr11', 'chr2', 'chrY', 'chrX']
list_ok = ['chr1', 'chr2', 'chr11', 'chrX', 'chrY']


###########################
# create chromosomes plot #
###########################

gene_positions = {'chr1' : [(125600,136700),(145600,150500),(20000000,20003000),(20003020,20003025),(20003100,20003150)], 'chr2' : [(125600000,125600050)], 'chr17' :  [(43044295, 43170245), (43044295, 43170245)]}
gene_positions2 = {'chr1' : [(125600,136700),(20000000,20003000),(20003020,20003025),(20003100,20003150)], 'chr2' : [(125600000,125600050)], 'chr19' :  [(43044295, 43170245), (43044295, 43170245)]}
gene_positions3 = {'chr1' : [(125600,136700),(145600,150500),(20000000,20003000),(20003020,20003025),(20003100,20003150)], 'chr2' : [(125600000,125600050)], 'chr17' :  [(43044295, 43170245), (43044295, 43170245)], 'chr19' :  [(43044295, 43170245), (43044295, 43170245)]}
cytoband_file = 'LOTUS_external_files/hg38_cytoband.tsv'
chr_plot_name =  str(uuid.uuid4())+'.svg'


################################################################################################################
################################################# Tests ########################################################
################################################################################################################


###############
# extract num #
###############

def test_extract_num():
	assert extract_num('chr34', pattern, float('inf')) == 34
	assert extract_num('chr23and2', pattern, float('inf')) == 23
	assert extract_num('chr', pattern, float('inf')) == float('inf')
	assert extract_num('bob1', pattern, None) == None
	assert extract_num('bob5', pattern2, float('inf')) == 5


########################
# count impacted genes #
########################

def test_count_impacted_genes():
	assert count_impacted_genes(positions, span, genes_positions) == 2
	assert count_impacted_genes(positions, span, genes_positions2) == 0


#####################
# get value for chr #
#####################

def test_get_value_for_chr():
	assert get_value_for_chr(positions_for_chr, chr, length, span, step) == (x_value, y_value)


############
# get plot #
############

def test_get_plot():
	assert get_plot(positions_for_chr, chr, chr_size, span, step) == (x_value, y_value)
	assert get_plot(positions_for_chr, 'chrG', chr_size, span, step) == (x_value_0, y_value_0)


##############
# scientific #
##############

def test_scientific():
	assert scientific(number1, None) == '1.27e+11'
	assert scientific(number2, None) == '1.23e+02'
	assert scientific(number3, None) == '0.00e+00'


#################
# cross product #
#################

def test_cross_product():
	assert cross_product(max1, max2, list_values) == list_values2


#################
# sorted nicely #
#################

def test_sorted_nicely():
	assert sorted_nicely(list_not_ok) == list_ok


###########################
# create chromosomes plot #
###########################

def test_create_chromosomes_plot():
	create_chromosomes_plot(gene_positions, gene_positions2, gene_positions3, cytoband_file, chr_plot_name)	
	assert Path(chr_plot_name).exists()
	os.remove(chr_plot_name)
