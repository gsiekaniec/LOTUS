#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import uuid
import logging
import pytest
from pathlib import Path
from python_scripts.panther_api import Panther_GEOA
from python_scripts.path_modification import true_stem

logger = logging.getLogger('This is not a logger')

################
# Panther GEOA #
################

list_genes = ['BRCA1', 'BRCA2']
panther_name = str(uuid.uuid4())+'.tsv'
pvalue = 0.05
maxres = 100
ids = ["Id", "name", "PValue", "FDR", "Genes number"]
#ids = 


################################################################################################################
################################################# Tests ########################################################
################################################################################################################

################
# Panther GEOA #
################

def test_Panther_GEOA():
	Panther_GEOA(list_genes, panther_name, logger, pvalue, maxres)
	assert Path(true_stem(panther_name)+'.tsv').exists()
	assert Path(true_stem(panther_name)+'.xlsx').exists()
	with open(Path(true_stem(panther_name)+'.tsv'), 'r') as f:
		first_line_list = f.readlines()[0].strip().split('\t')
	assert first_line_list == ids
	os.remove(Path(true_stem(panther_name)+'.tsv'))
	os.remove(Path(true_stem(panther_name)+'.xlsx'))
