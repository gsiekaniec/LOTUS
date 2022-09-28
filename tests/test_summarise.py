#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import pickle
import pyfastx
import os
import logging
import uuid
import pandas as pd
import numpy as np
import pytest
from collections import Counter
from pathlib import Path
from python_scripts.summarise import create_snp_dict, create_dataframe_from_gene, create_ordered_dataframe, graph_snp, graph_indel, is_fasta, is_pickle, get_genome_dict, read_vcf, add_snp, size_indel, variants_count, is_snp, is_dnp, is_tnp, is_np, is_deletion, is_insertion, write_stats, write_impacted_genes


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

complete_profile_dictionary = {'T>G': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'GTA': 0.0, 'GTT': 0.0, 'GTC': 0.0, 'GTG': 0.0}, 
		'C>A': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.5, 'CCG': 0.3, 'GCA': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}, 
		'T>C': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.01, 'CTG': 0.04, 'GTA': 0.0, 'GTT': 0.05, 'GTC': 0.0, 'GTG': 0.0}, 
		'C>G': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.0, 'CCG': 0.0, 'GCA': 0.05, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}, 
		'T>A': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'GTA': 0.0, 'GTT': 0.0, 'GTC': 0.0, 'GTG': 0.0}, 
		'C>T': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.05, 'CCG': 0.0, 'GCA': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}
		}
snp_plot_name = str(uuid.uuid4())+'.svg'
logger = logging.getLogger('Logger de la muerte')

###############
# graph indel #
###############

empty_counter = Counter()
deletion_counter = Counter({1:7,2:3,8:2,65:2})
insertion_counter = Counter({1:14,2:1,5:2,68:3})
indel_plot_name = str(uuid.uuid4())+'.svg'
insertion_plot_name = str(uuid.uuid4())+'.svg'
deletion_plot_name = str(uuid.uuid4())+'.svg'

##############################################
# is fasta and is pickle and get genome dict #
##############################################

fasta_ok = str(uuid.uuid4())+'.fasta'
with open(fasta_ok , 'w') as fa:
        fa.write('>chr1\nATTTAGAG\n>chr2\nATTTAGGCAGCAGC\n')
fasta_not_fasta = str(uuid.uuid4())+'.fasta'
with open(fasta_not_fasta , 'w') as fa:
	fa.write('May the force be with you !')
empty_fasta = str(uuid.uuid4())+'.fasta'
Path(empty_fasta).touch()

file_not_exists = str(uuid.uuid4())+'.txt'

pickle_ok = str(uuid.uuid4())+'.pk'
dict_pickle = {'chr1':'ATTTAGAG','chr2':'ATTTAGGCAGCAGC'}
with open(pickle_ok, 'wb') as out_pk:
	pickle.dump(dict_pickle, out_pk)
pickle_empty = str(uuid.uuid4())+'.pk'
Path(pickle_empty).touch()

############
# read vcf #
############

false_vcf = str(uuid.uuid4())+'.vcf'
with open(false_vcf, 'w') as out_vcf:
	out_vcf.write('##fileformat=VCFv4.2\n##FILTER=<ID=FAIL,Description=\"Fail the site if all alleles fail but for different reasons.\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXAMPLE\nchr1\t5\t.\tA\tT\t.\tPASS\tAS_FilterStatus=SITE;AS_SB_TABLE=17,31|57,57;DP=182;ECNT=1;FUNCOTATION=[CADM3|hg38|chr9|39609857|39609857|MISSENSE||INSERTION|G|G|GC|false|false];GERMQ=49;MBQ=35,33;MFRL=168,162;MMQ=60,60;MPOS=27;POPAF=7.30;ROQ=93;RPA=3,4;RU=C;STR;STRQ=93;TLOD=234.51;LOTUS_FILTER=PASS\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/1:48,114:0.669:162:12,35:13,27:42,86:17,31,57,57')
vcf_generator = (i for i in [{'idx_chr':0,'idx_pos':1,'idx_ref':3,'idx_alts':4,'idx_filt':6,'idx_info':7,'idx_format':8,'idx_values':9},['chr1','5','.','A','T','.','PASS','AS_FilterStatus=SITE;AS_SB_TABLE=17,31|57,57;DP=182;ECNT=1;FUNCOTATION=[CADM3|hg38|chr9|39609857|39609857|MISSENSE||INSERTION|G|G|GC|false|false];GERMQ=49;MBQ=35,33;MFRL=168,162;MMQ=60,60;MPOS=27;POPAF=7.30;ROQ=93;RPA=3,4;RU=C;STR;STRQ=93;TLOD=234.51;LOTUS_FILTER=PASS','GT:AD:AF:DP:F1R2:F2R1:FAD:SB','0/1:48,114:0.669:162:12,35:13,27:42,86:17,31,57,57']])

###########
# add snp #
###########

empty_profile_dictionary = {'T>G': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0}, 
		'C>A': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0},
                'T>C': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0},
                'C>G': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0},
                'T>A': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0},
                'C>T': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0}
                }
new_profile_dictionary = {'T>G': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0},
                'C>A': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0},
                'T>C': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 1, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0},
                'C>G': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0},
                'T>A': {'ATA': 0, 'ATT': 0, 'ATC': 0, 'ATG': 0, 'TTA': 0, 'TTT': 0, 'TTC': 0, 'TTG': 0, 'CTA': 0, 'CTT': 0, 'CTC': 0, 'CTG': 0, 'GTA': 0, 'GTT': 0, 'GTC': 0, 'GTG': 0},
                'C>T': {'ACA': 0, 'ACT': 0, 'ACC': 0, 'ACG': 0, 'TCA': 0, 'TCT': 0, 'TCC': 0, 'TCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'CCG': 0, 'GCA': 0, 'GCT': 0, 'GCC': 0, 'GCG': 0}
                }
reference = 'A'
alternative = 'G'
triplet = 'AAG'

##############
# size indel #
##############

ref = 'ATTTG'
deletion = 'ATTG'
insertion = 'ATTTGG'
deletion3 = 'AG'
insertion3 = 'ATTTGCTG'

##################
# variants count #
##################

expected_result = (Counter(), Counter(), {'T>C': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'GTA': 0.0, 'GTT': 0.0, 'GTC': 0.0, 'GTG': 0.0}, 'C>A': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.0, 'CCG': 0.0, 'GCA': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}, 'C>G': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.0, 'CCG': 0.0, 'GCA': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}, 'C>T': {'ACA': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACG': 0.0, 'TCA': 0.0, 'TCT': 0.0, 'TCC': 0.0, 'TCG': 0.0, 'CCA': 0.0, 'CCT': 0.0, 'CCC': 0.0, 'CCG': 0.0, 'GCA': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCG': 0.0}, 'T>G': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 0.0, 'CTT': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'GTA': 0.0, 'GTT': 0.0, 'GTC': 0.0, 'GTG': 0.0}, 'T>A': {'ATA': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'TTC': 0.0, 'TTG': 0.0, 'CTA': 1.0, 'CTT': 0.0, 'CTC': 0.0, 'CTG': 0.0, 'GTA': 0.0, 'GTT': 0.0, 'GTC': 0.0, 'GTG': 0.0}}, {'CADM3': [1, [1, 0, 0, 0, 0, 0], 'chr1', ['A'], [['T']], [5]]}, {'PASS': 1, 'SNP': [1, {'CADM3'}], 'DNP': [0, set()], 'TNP': [0, set()], 'NP': [0, set()], 'INSERTION': [0, set()], 'DELETION': [0, set()]})

########
# is * #
########

ref_snp = 'A'
ref_dnp = 'AT'
ref_tnp = 'ATC'
ref_np = 'ATCGTA'
ref_indel = 'ATG'
snp = 'T'
dnp = 'TC'
tnp = 'TCG'
np = 'TCGAAT'
alt_insertion = 'ATCG'
alt_deletion = 'AG'

###############
# write stats #
###############

count_stats = Counter({'Total': 5, 'germline': 0, 'PON': 0, 'non functional': 0, 'germline+PON': 0, 'germline+non functional': 0, 'PON+non functional': 0, 'germline+PON+non functional': 0, 'PASS': 3, 'SNP': [2, {'LINC01410', 'PLXNB3'}], 'DNP': [0, set()], 'TNP': [0, set()], 'NP': [0, set()], 'INSERTION': [1, {'CADM3'}], 'DELETION': [0, set()]})
stats_file = str(uuid.uuid4())+'.txt'

########################
# write impacted genes #
########################

out_genes_file = str(uuid.uuid4())+'.xlsx'
list_genes = ['BRCA1', 'BRCA2']

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

###############
# graph indel #
###############

def test_graph_indel():
	graph_indel(deletion_counter, insertion_counter, indel_plot_name, logger)
	assert Path(indel_plot_name).exists()
	os.remove(indel_plot_name)

def test_graph_insertion_only():
	graph_indel(empty_counter, insertion_counter, insertion_plot_name, logger)
	assert Path(insertion_plot_name).exists()
	os.remove(insertion_plot_name)

def test_graph_deletion_only():
	graph_indel(deletion_counter, empty_counter, deletion_plot_name, logger)
	assert Path(deletion_plot_name).exists()
	os.remove(deletion_plot_name)

def test_graph_no_count():
        assert graph_indel(empty_counter, empty_counter, indel_plot_name, logger) == None

############
# is fasta #
############

def test_is_fasta_ok():
	assert is_fasta(fasta_ok)

def test_is_fasta_pickle():
	assert not is_fasta(pickle_ok)

def test_is_fasta_not_exists():
	with pytest.raises(FileExistsError, match='File '+file_not_exists+' doesn\'t exists'):
		is_fasta(file_not_exists)

def test_is_fasta_not_ok():
	assert not is_fasta(fasta_not_fasta)

def test_is_fasta_empty():
	assert is_fasta(empty_fasta) == False
	os.remove(empty_fasta)

#############
# is pickle #
#############


def test_is_pickle_ok():
	assert is_pickle(pickle_ok)

def test_is_pickle_empty():
	assert not is_pickle(pickle_empty)	
	os.remove(pickle_empty)

def test_is_pickle_not_ok():
	assert not is_pickle(fasta_ok)

def test_is_pickle_not_exists():
	with pytest.raises(FileNotFoundError, match='File '+file_not_exists+' doesn\'t exists'):
		is_pickle(file_not_exists)


###################
# get genome dict #
###################


def test_get_genome_dict_from_fasta():
	assert get_genome_dict(fasta_ok, logger) == dict_pickle
	os.remove(fasta_ok)
	os.remove(Path(fasta_ok).with_suffix('.pk'))

def test_get_genome_dict_from_pickle():
	assert get_genome_dict(pickle_ok, logger) == dict_pickle
	os.remove(pickle_ok)

def test_get_genome_dict_not_fasta_or_pickle():
	with pytest.raises(ValueError, match=str(fasta_not_fasta)+' is not a fasta or a pickle file !'):
		get_genome_dict(fasta_not_fasta, logger)
	os.remove(fasta_not_fasta)	

############
# read vcf #
############

def test_read_vcf():
	assert list(read_vcf(false_vcf)) == list(vcf_generator)

###########
# add snp #
###########

def test_add_vcf():
	assert add_snp(empty_profile_dictionary, reference, alternative, triplet) == new_profile_dictionary


##############
# size indel #
##############

def test_size_insertion ():
	assert size_indel(ref, insertion) == 1

def test_size_deletion_3 ():
	assert size_indel(ref, deletion3) == 3

def test_size_insertion_3 ():
	assert size_indel(ref, insertion3) == 3


def test_size_deletion ():
	assert size_indel(ref, deletion) == 1


##################
# variants count #
##################

def test_variants_count():
	assert variants_count(False, false_vcf, dict_pickle, logger) == expected_result


########
# is * #
########


def test_is_snp():
	assert is_snp(len(ref_snp), len(snp))
	assert not is_snp(len(ref_dnp), len(dnp))
	assert not is_snp(len(ref_indel), len(alt_deletion))

def test_is_dnp():
	assert is_dnp(len(ref_dnp), len(dnp))
	assert not is_dnp(len(ref_snp), len(snp))
	assert not is_dnp(len(ref_indel), len(alt_deletion))

def test_is_tnp():
	assert is_tnp(len(ref_tnp), len(tnp))
	assert not is_tnp(len(ref_snp), len(snp))
	assert not is_tnp(len(ref_indel), len(alt_deletion))

def test_is_np():
	assert is_np(len(ref_np), len(np))
	assert not is_np(len(ref_snp), len(snp))
	assert not is_np(len(ref_indel), len(alt_deletion))

def test_is_insertion():
	assert is_insertion(len(ref_indel), len(alt_insertion))
	assert not is_insertion(len(ref_indel), len(alt_deletion))
	assert not is_insertion(len(ref_snp), len(snp))

def test_is_deletion():
	assert is_deletion(len(ref_indel), len(alt_deletion))
	assert not is_deletion(len(ref_indel), len(alt_insertion))
	assert not is_deletion(len(ref_snp), len(snp))


###############
# write stats #
###############


def test_write_stats():
	write_stats(False, false_vcf, stats_file, count_stats, out_genes_file, logger)
	assert Path(stats_file).exists()			
	os.remove(stats_file)
	os.remove(false_vcf)	


########################
# write impacted genes #
########################

def test_write_impacted_genes():
	write_impacted_genes(out_genes_file, df_gene, logger)
	assert Path(out_genes_file).exists()
	assert Path(out_genes_file).with_suffix('.tsv').exists()
	os.remove(out_genes_file)
	os.remove(Path(out_genes_file).with_suffix('.tsv'))



