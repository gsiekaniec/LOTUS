#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import uuid
from python_scripts.read_vcf import read_vcf, get_vcf_header


############
# read vcf #
############

false_vcf = str(uuid.uuid4())+'.vcf'
with open(false_vcf, 'w') as out_vcf:
	out_vcf.write('##fileformat=VCFv4.2\n##FILTER=<ID=FAIL,Description=\"Fail the site if all alleles fail but for different reasons.\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXAMPLE\nchr1\t5\t.\tA\tT\t.\tPASS\tAS_FilterStatus=SITE;AS_SB_TABLE=17,31|57,57;DP=182;ECNT=1;FUNCOTATION=[CADM3|hg38|chr9|39609857|39609857|MISSENSE||INSERTION|G|G|GC|false|false];GERMQ=49;MBQ=35,33;MFRL=168,162;MMQ=60,60;MPOS=27;POPAF=7.30;ROQ=93;RPA=3,4;RU=C;STR;STRQ=93;TLOD=234.51;LOTUS_FILTER=PASS\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/1:48,114:0.669:162:12,35:13,27:42,86:17,31,57,57')
vcf_generator = (i for i in [{'idx_chr':0,'idx_pos':1,'idx_ref':3,'idx_alts':4,'idx_filt':6,'idx_info':7,'idx_format':8,'idx_values':9},['chr1','5','.','A','T','.','PASS','AS_FilterStatus=SITE;AS_SB_TABLE=17,31|57,57;DP=182;ECNT=1;FUNCOTATION=[CADM3|hg38|chr9|39609857|39609857|MISSENSE||INSERTION|G|G|GC|false|false];GERMQ=49;MBQ=35,33;MFRL=168,162;MMQ=60,60;MPOS=27;POPAF=7.30;ROQ=93;RPA=3,4;RU=C;STR;STRQ=93;TLOD=234.51;LOTUS_FILTER=PASS','GT:AD:AF:DP:F1R2:F2R1:FAD:SB','0/1:48,114:0.669:162:12,35:13,27:42,86:17,31,57,57']])

##################
# get vcf header #
##################

header_list = ['##fileformat=VCFv4.2','##FILTER=<ID=FAIL,Description=\"Fail the site if all alleles fail but for different reasons.\">','#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXAMPLE']

################################################################################################################
################################################# Tests ########################################################
################################################################################################################


############
# read vcf #
############

def test_read_vcf():
	assert list(read_vcf(false_vcf)) == list(vcf_generator)


##################
# get vcf header #
##################

def test_get_vcf_header():
	assert list(get_vcf_header(false_vcf)) == header_list
	os.remove(false_vcf)
