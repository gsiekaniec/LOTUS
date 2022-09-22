#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import logging
import uuid
import pandas as pd
import numpy as np
import pytest
from pathlib import Path
from python_scripts.summary import create_snp_dict, create_dataframe_from_gene, create_ordered_dataframe, graph_snp



###################
# create snp dict #
###################

empty_count_dictionary = {'T>G': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0}, 'C>G': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0}, 'C>T': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0}, 'T>C': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0}, 'C>A': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0}, 'T>A': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0}}

##############################
# create dataframe from gene #
##############################

gene_dictionary = {'LINC01410': [1, [1, 0, 0, 0, 0, 0], 'chr9', ['G'], [['T']], ['62801195']], 'PLXNB3': [1, [1, 0, 0, 0, 0, 0], 'chr6', ['C'], [['A']], ['143507595']], 'CADM3': [1, [0, 0, 0, 0, 1, 0], 'chr1', ['G'], [['GC']], ['159192649']]}
dt = np.dtype([('tumor_burden', int),('details_(snp,dnp,tnp,np,insertion,deletion)', np.unicode_, 16),('chromosome', np.unicode_, 16),('ref', np.unicode_, 16),('alt_variant(s)', np.unicode_, 16),('position(s)', np.unicode_, 16)])

df_gene = pd.DataFrame(np.array([(1, '0,0,0,0,1,0', 'chr1', 'G', 'GC', 159192649),(1,'1,0,0,0,0,0','chr9','G','T',62801195),(1,'1,0,0,0,0,0','chr6','C','A',143507595)], dtype=dt),
		columns=['tumor_burden', 'details_(snp,dnp,tnp,np,insertion,deletion)', 'chromosome', 'ref', 'alt_variant(s)', 'position(s)'], 
		index=['CADM3', 'LINC01410', 'PLXNB3'])
df_gene.index.name = 'gene_name'

############################
# create ordered dataframe #
############################

profile_dictionary = {'C>A' : {'ACA' : 0.01}, 'T>G': {'ATA': 0.8, 'ATT': 0.15, 'ATC': 0.04}}
tuples_profile = [('C>A', 'ACA'),('T>G', 'ATA'),('T>G', 'ATC'),('T>G', 'ATT')]
index_profile = pd.MultiIndex.from_tuples(tuples_profile)
array_profile = np.array([0.01, 0.8, 0.15, 0.04])
df_profile = pd.DataFrame(array_profile, index=index_profile)

#############
# graph snp #
#############

complete_profile_dictionary = {'T>G': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'GTA': 0.0, 'GTT': 0.0, 'GTC': 0.0, 'GTG': 0.0}, 'C>A': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.5, 'CCG': 0.3, 'GCA': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}, 'T>C': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.01, 'CTG': 0.04, 'GTA': 0.0, 'GTT': 0.05, 'GTC': 0.0, 'GTG': 0.0}, 
		'C>G': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.0, 'CCG': 0.0, 'GCA': 0.05, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}, 
		'T>A': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'GTA': 0.0, 'GTT': 0.0, 'GTC': 0.0, 'GTG': 0.0}, 
		'C>T': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.05, 'CCG': 0.0, 'GCA': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}
		}
snp_plot_name = str(uuid.uuid4())+'.svg'
logger = logging.getLogger('Logger de la muerte')


################################################################################################################
################################################# Tests ########################################################
################################################################################################################


###################
# create snp dict #
###################

def test_create_snp_dict():
	assert create_snp_dict() == empty_count_dictionary
	
##############################
# create dataframe from gene #
##############################

def test_create_dataframe_from_gene():
	assert (create_dataframe_from_gene(gene_dictionary) == df_gene).all

############################
# create ordered dataframe #
############################

def test_create_ordered_dataframe():
	assert (create_ordered_dataframe(profile_dictionary) == df_profile).all

#############
# graph snp #
#############


def test_graph_snp():
	graph_snp(complete_profile_dictionary, snp_plot_name , logger)	
	assert Path(snp_plot_name).exists()
	os.remove(snp_plot_name)	




