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
from python_scripts.compare import get_variants, add_to_genes, modify_variants_pass_and_get_genes, create_graph_snp, create_graph_indel, compare_vcf
from python_scripts.path_modification import true_stem


logger = logging.getLogger('This is not a logger')


################
# get variants #
################

false_vcf = str(uuid.uuid4())+'.vcf'
with open(false_vcf, 'w') as out_vcf:
	out_vcf.write('##fileformat=VCFv4.2\n##FILTER=<ID=FAIL,Description=\"Fail the site if all alleles fail but for different reasons.\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXAMPLE\nchr1\t5\t.\tA\tT\t.\tPASS\tAS_FilterStatus=SITE;AS_SB_TABLE=17,31|57,57;DP=182;ECNT=1;FUNCOTATION=[CADM3|hg38|chr9|39609857|39609857|MISSENSE||INSERTION|G|G|GC|false|false];GERMQ=49;MBQ=35,33;MFRL=168,162;MMQ=60,60;MPOS=27;POPAF=7.30;ROQ=93;RPA=3,4;RU=C;STR;STRQ=93;TLOD=234.51;LOTUS_FILTER=PASS\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/1:48,114:0.669:162:12,35:13,27:42,86:17,31,57,57')


################
# add to genes #
################

genes = {}
one_gene = 'BRCA1'
type_ok = 'weak'
type_not_ok = 'i\'m not a variant type'


######################################
# modify variants pass and get genes #
######################################

variants_ok =  {('chr1', '5', 'A', 'T')}
variants_2 = {('chr1', '5', 'A', 'T'), ('chr2', '7', 'G', 'C')}
genes_from_vcf = {'CADM3': {'weak': 1, 'strong': 0}}


####################
# create graph snp #
####################

out_snp_graph = str(uuid.uuid4())+'.svg'

dt = np.dtype([('Unnamed: 0', np.unicode_, 16),('Unnamed: 1', np.unicode_, 16),('sample', float)])
df_snp = pd.DataFrame(np.array([('C>A', 'ACA', 1.0)], dtype=dt), columns=['Unnamed: 0', 'Unnamed: 1', 'sample'])
df_snp2 = pd.DataFrame(np.array([('C>T', 'ACT', 0.7),('C>A', 'ACA', 0.3)], dtype=dt), columns=['Unnamed: 0', 'Unnamed: 1', 'sample'])


######################
# create graph indel #
######################

out_indel_graph = str(uuid.uuid4())+'.svg'
out_insertion_only = str(uuid.uuid4())+'.svg'
out_deletion_only = str(uuid.uuid4())+'.svg'
out_no_indel_graph = str(uuid.uuid4())+'.svg'

ins1 = str(uuid.uuid4())+'.tsv'
with open(ins1 ,'w') as o:
	o.write('\tsample1\n1\t6\n2\t3\n72\t1\n')
ins2 = str(uuid.uuid4())+'.tsv'
with open(ins2 ,'w') as o:
	o.write('\tsample2\n1\t2\n2\t8\n45\t1\n')
del1 = str(uuid.uuid4())+'.tsv'
with open(del1 ,'w') as o:
	o.write('\tsample1\n1\t3\n3\t9\n5\t10\n')
del2 = str(uuid.uuid4())+'.tsv'
with open(del2 ,'w') as o:
	o.write('\tsample2\n1\t2\n2\t3\n5\t4\n')


################################################################################################################
################################################# Tests ########################################################
################################################################################################################



################
# get variants #
################

def test_get_variants():	
	assert get_variants(false_vcf, set()) == variants_ok
	assert not get_variants(false_vcf, set()) == {}


################
# add to genes #
################

def test_add_to_genes():
	add_to_genes(genes, one_gene, type_ok)
	assert genes == {'BRCA1' : {'strong': 0, 'weak': 1}}
	add_to_genes(genes, one_gene, type_ok)
	assert genes == {'BRCA1' : {'strong': 0, 'weak': 2}}
	with pytest.raises(KeyError, match=type_not_ok):
		add_to_genes(genes, one_gene, type_not_ok)


######################################
# modify variants pass and get genes #
######################################

def test_modify_variants_pass_and_get_genes():
	outfile = true_stem(false_vcf)+'_compare_to_'+true_stem('file2')+'.pass.vcf'
	assert modify_variants_pass_and_get_genes(false_vcf, 'file2', variants_2, variants_ok, {}, logger) == genes_from_vcf
	assert Path(outfile).exists()
	os.remove(outfile)
	os.remove(false_vcf)	

	
####################
# create graph snp #
####################

def test_create_graph_snp():
	print(df_snp,'\n',df_snp2,'\n')
	create_graph_snp(df_snp, df_snp2, out_snp_graph, logger)
	assert Path(out_snp_graph).exists()
	os.remove(out_snp_graph)


######################
# create_graph_indel #
######################

def test_create_graph_indel():
	create_graph_indel(del1, del2, ins1, ins2, out_indel_graph, logger)
	assert Path(out_indel_graph).exists()
	os.remove(out_indel_graph)
	
	create_graph_indel(del1, del2, ins1, None, out_insertion_only, logger)
	assert Path(true_stem(out_insertion_only)+'_sample1.svg').exists()
	os.remove(true_stem(out_insertion_only)+'_sample1.svg')

	create_graph_indel(del1, None, ins1, ins2, out_deletion_only, logger)
	assert Path(true_stem(out_deletion_only)+'_sample1.svg').exists()
	os.remove(true_stem(out_deletion_only)+'_sample1.svg')

	create_graph_indel(del1, None, ins1, None, out_no_indel_graph, logger)
	assert not Path(true_stem(out_no_indel_graph)+'_sample1.svg').exists()

	os.remove(ins1)
	os.remove(ins2)
	os.remove(del1)
	os.remove(del2)














