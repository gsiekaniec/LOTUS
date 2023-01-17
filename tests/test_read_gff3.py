#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import uuid
import logging
import pytest
from pathlib import Path
from python_scripts.read_gff3 import read_gff3
from python_scripts.path_modification import true_stem

logger = logging.getLogger('This is not a logger')

#############
# read gff3 #
#############

gff3 = 'LOTUS_external_files/Homo_sapiens.GRCh38.108.chr.gff3.gz'
not_gff3 = str(uuid.uuid4())+'.gff3'
with open(not_gff3, 'w') as o:
	o.write("I'm not a gff3 file !")

################################################################################################################
################################################# Tests ########################################################
################################################################################################################

#############
# read gff3 #
#############

def test_read_gff3():
	dico = read_gff3(gff3)
	gene_name_dico = dico[0]
	gene_id_dico = dico[1]
	transcript_dico = dico[2]
	assert gene_name_dico['BRCA1'] == [('ENSG00000012048', 'chr17', (43044295, 43170245))]
	assert gene_id_dico['ENSG00000012048'] == [('BRCA1', 'chr17', (43044295, 43170245))]
	assert transcript_dico['ENST00000357654'] == [('BRCA1', 'chr17', (43044295, 43170245))]	
	assert Path('LOTUS_external_files/Homo_sapiens.GRCh38.108.chr.gff3.gz').with_suffix('.pk').exists()
	with pytest.raises(ValueError, match=str(not_gff3)+' is not a gff3 file !'):
		read_gff3(not_gff3)
	os.remove(not_gff3)

