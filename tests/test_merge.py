#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import numpy as np
import pandas as pd
import os
import logging
import uuid
from pathlib import Path
import pytest
from python_scripts.merge import read_config, get_nb_files_and_file_names, create_upsetplot
from python_scripts.path_modification import true_stem


logger = logging.getLogger('This is not a logger')



###############
# read config #
###############

false_tsv1 = str(uuid.uuid4())+'.tsv'
with open(false_tsv1, 'w') as o:
	o.write('Gene symbol\tVariant weakness (in %)\tTumor burden (symmetrical difference)\tsample11\tsample12\nGene1\t0.0\t4\t2\t2\nGene2\t100.0\t3\t1\t2')
false_tsv2 = str(uuid.uuid4())+'.tsv'
with open(false_tsv2, 'w') as o:
	o.write('Gene symbol\tVariant weakness (in %)\tTumor burden (symmetrical difference)\tsample21\tsample22\nGene2\t30.0\t7\t3\t4\nGene3\t70.0\t1\t1\t0')
dt_xlsx = np.dtype([('Gene symbol', np.unicode_, 16),('Variant weakness (in %)', float), ('Tumor burden (symmetrical difference', int),('sample31', int),('sample32', int)])
df_for_xlsx = pd.DataFrame(np.array([('Gene2', 100.0, 4, 2, 2)], dtype=dt_xlsx), columns=['Gene symbol', 'Variant weakness (in %)', 'Tumor burden (symmetrical difference)', 'sample31', 'sample32'])
df_for_xlsx = df_for_xlsx.set_index('Gene symbol')
false_xlsx = str(uuid.uuid4())+'.xlsx'
df_for_xlsx.to_excel(false_xlsx)
config_file = str(uuid.uuid4())+'.vcf'
with open(config_file, 'w') as o:
	o.write(f'{false_tsv1},sample1\n{false_tsv2},sample2\n{false_xlsx},sample3\n')

dt_tsv1 = np.dtype([('Gene symbol', np.unicode_, 16),('Variant weakness (in %)', float), ('Tumor burden (symmetrical difference)', int),('sample11', int),('sample12', int)])
df_tsv1 =  pd.DataFrame(np.array([('Gene1', 0.0, 4, 2, 2),('Gene2', 100.0, 3, 1, 2)], dtype=dt_tsv1), columns=['Gene symbol', 'Variant weakness (in %)', 'Tumor burden (symmetrical difference)', 'sample11', 'sample12'])
df_tsv1 = df_tsv1.set_index('Gene symbol')
dt_tsv2 = np.dtype([('Gene symbol', np.unicode_, 16),('Variant weakness (in %)', float), ('Tumor burden (symmetrical difference)', int),('sample21', int),('sample22', int)])
df_tsv2 =  pd.DataFrame(np.array([('Gene2', 30.0, 7, 3, 4),('Gene3', 70.0, 1, 1, 0)], dtype=dt_tsv2), columns=['Gene symbol', 'Variant weakness (in %)', 'Tumor burden (symmetrical difference)', 'sample21', 'sample22'])
df_tsv2 = df_tsv2.set_index('Gene symbol')
res = [('sample1',df_tsv1),('sample2',df_tsv2),('sample3',df_for_xlsx)]


####################
# create upsetplot #
####################

upset_name = str(uuid.uuid4())+'.svg'
data_for_upsetplot = {'sample1' : {'BRCA1', 'BRCA2', 'R2D2', 'HDFA1'}, 'sample2' : {'BRCA1', 'R2D2'}, 'sample3' : {'R2D2', 'HDFA1', 'MUC3A'}}
category_for_upsetplot =  [['sample1'], ['sample2'], ['sample3'], ['sample1', 'sample2'], ['sample1', 'sample3'], ['sample2', 'sample3'], ['sample1', 'sample2', 'sample3']]
names_for_upsetplot = {'sample1', 'sample2', 'sample3'}



################################################################################################################
################################################# Tests ########################################################
################################################################################################################


###############
# read config #
###############

def test_read_config():
	for i,data in enumerate(list(read_config(config_file))):
		assert data[0] == res[i][0]
		assert all(data[1] == res[i][1])
	os.remove(false_tsv1)
	os.remove(false_tsv2)	
	os.remove(false_xlsx)


###############################
# get nb files and file names #
###############################

def test_get_nb_files_and_file_names():
	assert get_nb_files_and_file_names(config_file) == (3,['sample1','sample2','sample3'])
	os.remove(config_file)	


####################
# create upsetplot #
####################

def test_create_upsetplot():
	create_upsetplot(data_for_upsetplot, category_for_upsetplot, upset_name, names_for_upsetplot, 1, 0, 1, 0, 101, logger)
	assert Path(upset_name).exists()
	os.remove(upset_name)





