#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import uuid
import logging
import pytest
from pathlib import Path
from python_scripts.toppgene_api import ToppGene_GEOA
from python_scripts.path_modification import true_stem

logger = logging.getLogger('This is not a logger')

#################
# ToppGene GEOA #
#################

list_genes = ['BRCA1', 'BRCA2']
toppgene_name = str(uuid.uuid4())+'.tsv'
pvalue = 0.05
maxres = 100
ids = ['ID', 'Name', 'PValue', 'QValueFDRBH', 'QValueFDRBY', 'QValueBonferroni', 'TotalGenes', 'GenesInTerm', 'GenesInQuery', 'GenesInTermInQuery', 'URL', 'Genes']


################################################################################################################
################################################# Tests ########################################################
################################################################################################################

#################
# ToppGene GEOA #
#################

def test_ToppGene_GEOA():
	ToppGene_GEOA(list_genes, toppgene_name, logger, pvalue, maxres)
	assert Path(true_stem(toppgene_name)+'.tsv').exists()
	assert Path(true_stem(toppgene_name)+'.xlsx').exists()
	with open(Path(true_stem(toppgene_name)+'.tsv'), 'r') as f:
		first_line_list = f.readlines()[0].strip().split('\t') 
	assert first_line_list == ids
	os.remove(Path(true_stem(toppgene_name)+'.tsv'))
	os.remove(Path(true_stem(toppgene_name)+'.xlsx'))
