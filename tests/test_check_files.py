#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import uuid
from pathlib import Path
import pandas as pd
import numpy as np
import pytest
from python_scripts.check_files import verif_input_vcf, verif_output, verif_input, verif_input_xlsx, verif_input_tsv, verif_input_config_merge, verif_input_config, verif_supplementary_information_file


###################
# Verif Input VCF #
###################

test1 = str(uuid.uuid4())+'.vcf'
with open(test1 , 'w') as t1:
        t1.write('##fileformat=VCFv.001\n')
test2 = str(uuid.uuid4())+'.vcf'
Path(test2).touch()
test3 = str(uuid.uuid4())+'.vcf'
if Path(test3).exists():
        os.chmod(test3, 711)
with open(test3 , 'w') as t3:
        t3.write('##fileformat=VCFv.001\n')
os.chmod(test3, 000)
test4 = str(uuid.uuid4())
with open(test4 , 'w') as t4:
        t4.write('##fileformat=VCFv.001\n')
test5 = str(uuid.uuid4())+'.vcf'
with open(test5 , 'w') as t5:
        t5.write('This is not a vcf file !\n')


###############
# Verif Ouput #
###############

output = str(uuid.uuid4())+'.vcf'
Path(output).touch()
output2 = str(uuid.uuid4())+'.vcf'
Path(output2).touch()
output3 = str(uuid.uuid4())+'.vcf'
if Path(output3).exists():
        os.chmod(output3, 711)
else:
	Path(output3).touch()
os.chmod(output3, 000)


###############
# Verif Input #
###############

file1 = str(uuid.uuid4())+'.txt'
with open(file1 , 'w') as f1:
        f1.write('I\'m not empty\n')
file2 = str(uuid.uuid4())+'.txt'
Path(file2).touch()
file3 = str(uuid.uuid4())+'.vcf'
if Path(file3).exists():
        os.chmod(file3, 711)
with open(file3 , 'w') as f3:
        f3.write('Con-fidential\n')
os.chmod(file3, 000)
file4 = str(uuid.uuid4())+'.txt'


####################
# Verif Input XLSX #
####################

dt = np.dtype([('Name', np.unicode_, 16),('Column 1', int)])
xlsx_df = pd.DataFrame(np.array([('M', 1), ('B', 0)], dtype=dt), columns=['Name', 'Column 1'])
xlsx_df = xlsx_df.set_index('Name')
xlsx_file = str(uuid.uuid4())+'.xlsx'
xlsx_df.to_excel(xlsx_file)
false_xlsx_file = str(uuid.uuid4())+'.xlsx'
with open(false_xlsx_file , 'w') as f:
	f.write('Not an xlsx or tsv file\n')


###################
# Verif Input TSV #
###################

tsv_file = str(uuid.uuid4())+'.tsv'
with open(tsv_file , 'w') as f:
	f.write('Name\tColumn 1\nM\t1\nB\t0\n')


############################
# Verif Input Config Merge #
############################

tsv1 = str(uuid.uuid4())+'.tsv'
with open(tsv1, 'w') as o:
	o.write('Gene symbol\tVariant weakness (in %)\tTumor burden (symmetrical difference)\tsample11\tsample12\nGene1\t0.0\t4\t2\t2\nGene2\t100.0\t3\t1\t2')
tsv2 = str(uuid.uuid4())+'.tsv'
with open(tsv2, 'w') as o:
	o.write('Gene symbol\tVariant weakness (in %)\tTumor burden (symmetrical difference)\tsample21\tsample22\nGene2\t30.0\t7\t3\t4\nGene3\t70.0\t1\t1\t0')
dt_xlsx = np.dtype([('Gene symbol', np.unicode_, 16),('Variant weakness (in %)', float), ('Tumor burden (symmetrical difference', int),('sample31', int),('sample32', int)])
df_for_xlsx = pd.DataFrame(np.array([('Gene2', 100.0, 4, 2, 2)], dtype=dt_xlsx), columns=['Gene symbol', 'Variant weakness (in %)', 'Tumor burden (symmetrical difference)', 'sample31', 'sample32'])
df_for_xlsx = df_for_xlsx.set_index('Gene symbol')
xlsx = str(uuid.uuid4())+'.xlsx'
df_for_xlsx.to_excel(xlsx)
config_file_all1 = str(uuid.uuid4())+'.txt'
with open(config_file_all1, 'w') as o:
	o.write(f'{tsv1},sample1\n{tsv2},sample2\n{xlsx},sample3\n')
config_file_all2 = str(uuid.uuid4())+'.txt'
with open(config_file_all2, 'w') as o:
	o.write(f'{tsv1},sample1\n{tsv2},sample2\n{file4},sample3\n')
false_config = str(uuid.uuid4())+'.txt'
with open(false_config, 'w') as o:
	o.write(f'Not a good file\n')


######################
# Verif Input Config #
######################

# Empty file
empty_file = str(uuid.uuid4())+'.vcf'
Path(empty_file).touch()

# Vcf
vcf1 = str(uuid.uuid4())+'.vcf'
with open(vcf1, 'w') as o:
	o.write('##fileformat=VCFv4.2\n##FILTER=<ID=FAIL,Description=\"Fail the site if all alleles fail but for different reasons.\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXAMPLE\nchr1\t5\t.\tA\tT\t.\tPASS\tAS_FilterStatus=SITE;AS_SB_TABLE=17,31|57,57;DP=182;ECNT=1;FUNCOTATION=[CADM3|hg38|chr9|39609857|39609857|MISSENSE||INSERTION|G|G|GC|false|false];GERMQ=49;MBQ=35,33;MFRL=168,162;MMQ=60,60;MPOS=27;POPAF=7.30;ROQ=93;RPA=3,4;RU=C;STR;STRQ=93;TLOD=234.51;LOTUS_FILTER=PASS\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/1:48,114:0.669:162:12,35:13,27:42,86:17,31,57,57')
vcf2 = str(uuid.uuid4())+'.vcf'
with open(vcf2, 'w') as o:
	o.write('##fileformat=VCFv4.2\n##FILTER=<ID=FAIL,Description=\"Fail the site if all alleles fail but for different reasons.\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXAMPLE\nchr1\t5\t.\tA\tT\t.\tPASS\tAS_FilterStatus=SITE;AS_SB_TABLE=17,31|57,57;DP=182;ECNT=1;FUNCOTATION=[CADM3|hg38|chr9|39609857|39609857|MISSENSE||INSERTION|G|G|GC|false|false];GERMQ=49;MBQ=35,33;MFRL=168,162;MMQ=60,60;MPOS=27;POPAF=7.30;ROQ=93;RPA=3,4;RU=C;STR;STRQ=93;TLOD=234.51;LOTUS_FILTER=PASS\tGT:AD:AF:DP:F1R2:F2R1:FAD:SB\t0/1:48,114:0.669:162:12,35:13,27:42,86:17,31,57,57')	

# Profile file
profile1 = str(uuid.uuid4())+'.tsv'
with open(profile1, 'w') as o:
	o.write('\t\t01A-041-2-20ESM00382A-2_S5_R1_001.pass\nC>A\tACA\t0.0\nC>A\tACC\t0.01470588')

profile2 = str(uuid.uuid4())+'.tsv'
with open(profile2, 'w') as o:
	o.write('\t\t01A-041-2-20ESM00382A-2_S5_R1_001.pass\nC>A\tACA\t0.0\nC>A\tACC\t0.01470588')

# Deletion file
del1 = str(uuid.uuid4())+'.tsv'
with open(del1, 'w') as o:
	o.write('\t\t01A-041-2-20ESM00382A-2_S5_R1_001.pass\n1\t2\n3\t2\n4\t1\n9\t1\n12\t1\n72\t1')

del2 = str(uuid.uuid4())+'.tsv'
with open(del2, 'w') as o:
	o.write('\t\t01A-041-2-20ESM00382A-2_S5_R1_001.pass\n1\t2\n3\t2\n4\t1\n9\t1\n12\t1\n72\t1')

# Insertion file
ins1 = str(uuid.uuid4())+'.tsv'
with open(ins1, 'w') as o:
	o.write('\t\t01A-041-2-20ESM00382A-2_S5_R1_001.pass\n1\t2\n2\t2\n5\t1\n6\t2\n7\t1\n8\t1\n13\t1')

ins2 = str(uuid.uuid4())+'.tsv'
with open(ins2, 'w') as o:
	o.write('\t\t01A-041-2-20ESM00382A-2_S5_R1_001.pass\n1\t2\n2\t2\n5\t1\n6\t2\n7\t1\n8\t1\n13\t1')

config1 = str(uuid.uuid4())+'.txt'
with open(config1, 'w') as o:
	o.write(f'2\n\n### VCF files ###\n# Filtered.vcf files (from g-LOTUS filter)\n\n- filtered:\n{vcf1}\n{vcf2}\n\n# Corresponding Pass.vcf files in the same order (from g-LOTUS filter)\n\n- pass:\n{vcf1}\n{vcf2}\n\n### TSV files ###\n# snp profile tsv (same order than vcf files)\n\n- profile:\n{profile1}\n{profile2}\n\n# insertion count tsv (None if no file)\n\n- insertion:\n{ins1}\n{ins2}\n\n# deletion count tsv (None if no file)\n\n- deletion:\n{del1}\n{del2}')
config2 = str(uuid.uuid4())+'.txt'
with open(config2, 'w') as o:
	o.write(f'2\n\n### VCF files ###\n# Filtered.vcf files (from g-LOTUS filter)\n\n- filtered:\n{empty_file}\n{vcf2}\n\n# Corresponding Pass.vcf files in the same order (from g-LOTUS filter)\n\n- pass:\n{vcf1}\n{vcf2}\n\n### TSV files ###\n# snp profile tsv (same order than vcf files)\n\n- profile:\n{profile1}\n{profile2}\n\n# insertion count tsv (None if no file)\n\n- insertion:\n{ins1}\n{ins2}\n\n# deletion count tsv (None if no file)\n\n- deletion:\n{del1}\n{del2}')	
config3 = str(uuid.uuid4())+'.txt'
with open(config3, 'w') as o:
	o.write(f'2\n\n### VCF files ###\n# Filtered.vcf files (from g-LOTUS filter)\n\n- filtered:\n{vcf1}\n{vcf2}\n\n# Corresponding Pass.vcf files in the same order (from g-LOTUS filter)\n\n- pass:\n{vcf1}\n{vcf2}\n\n### TSV files ###\n# snp profile tsv (same order than vcf files)\n\n- profile:\n{empty_file}\n{profile2}\n\n# insertion count tsv (None if no file)\n\n- insertion:\n{ins1}\n{ins2}\n\n# deletion count tsv (None if no file)\n\n- deletion:\n{del1}\n{del2}')
config4 = str(uuid.uuid4())+'.txt'
with open(config4, 'w') as o:
	o.write(f'2\n\n### VCF files ###\n# Filtered.vcf files (from g-LOTUS filter)\n\n- filtered:\n{vcf1}\n{vcf2}\n\n# Corresponding Pass.vcf files in the same order (from g-LOTUS filter)\n\n- pass:\n{vcf1}\n{vcf2}\n\n### TSV files ###\n# snp profile tsv (same order than vcf files)\n\n- profile:\n{profile1}\n{profile2}\n\n# insertion count tsv (None if no file)\n\n- insertion:\n{empty_file}\n{ins2}\n\n# deletion count tsv (None if no file)\n\n- deletion:\n{del1}\n{del2}')
config5 = str(uuid.uuid4())+'.txt'
with open(config5, 'w') as o:
	o.write(f'2\n\n### VCF files ###\n# Filtered.vcf files (from g-LOTUS filter)\n\n- filtered:\n{vcf1}\n{vcf2}\n\n# Corresponding Pass.vcf files in the same order (from g-LOTUS filter)\n\n- pass:\n{vcf1}\n{vcf2}\n\n### TSV files ###\n# snp profile tsv (same order than vcf files)\n\n- profile:\n{profile1}\n{profile2}\n\n# insertion count tsv (None if no file)\n\n- insertion:\n{ins1}\n{ins2}\n\n# deletion count tsv (None if no file)\n\n- deletion:\n{empty_file}\n{del2}')
config6 = str(uuid.uuid4())+'.txt'
Path(config6).touch()


########################################
# Verif Supplementary Information File #
########################################

infos_file = 'LOTUS_external_files/Lotus_ExternalBases_202301.xlsx'
wrong1 = str(uuid.uuid4())+'.xlsx'
Path(wrong1).touch()
wrong2 = str(uuid.uuid4())+'.xlsx'


################################################################################################################
################################################# Tests ########################################################
################################################################################################################


###################
# Verif Input VCF #
###################

def test_verif_input_ok() :
        assert verif_input_vcf(test1) == None
        os.remove(test1)

def test_verif_input_empty() :
        with pytest.raises(ValueError, match='File '+str(test2)+' is empty !'):
                verif_input_vcf(test2)
        os.remove(test2)

def test_verif_input_no_permission() :
        with pytest.raises(ValueError, match='File '+str(test3)+' is not available for reading !'):
                verif_input_vcf(test3)
        os.remove(test3)

def test_verif_input_extension() :
        with pytest.raises(ValueError, match='Extension of file '+str(test4)+' is not .vcf !'):
                verif_input_vcf(test4)
        os.remove(test4)

def test_verif_input_format() :
        with pytest.raises(ValueError, match='File '+str(test5)+' is not a vcf file or the header ##fileformat=VCF doesn\'t exist !'):
                verif_input_vcf(test5)
        os.remove(test5)


###############
# Verif Ouput #
###############

def test_verif_ouput_yes(monkeypatch) :
	monkeypatch.setattr('builtins.input', lambda _: "y")
	assert verif_output(output) == None	
	os.remove(output)

def test_verif_output_no(monkeypatch, capsys) :
	monkeypatch.setattr('builtins.input', lambda _: "n")
	with pytest.raises(SystemExit) as e:
		verif_output(output2)
	captured_stdout, captured_stderr = capsys.readouterr()
	assert captured_stdout == "Warning output file "+str(output2)+" already exist !\n\nRetry with another output file !\n\n"
	assert e.type == SystemExit
	assert e.value.code == 1
	os.remove(output2)

def test_verif_output_no_permission(monkeypatch) :
	monkeypatch.setattr('builtins.input', lambda _: "y")
	with pytest.raises(ValueError, match='File '+str(output3)+' is not available in writing !'):
		verif_output(output3)
	os.remove(output3)


###############
# Verif Input #
###############

def test_verif_input():
	assert verif_input(file1) == None	
	with pytest.raises(ValueError, match='File '+file2+' is empty !'):
		verif_input(file2)
	with pytest.raises(ValueError, match='File '+file3+' is not available for reading !'):
		verif_input(file3)
	with pytest.raises(ValueError, match='File '+file4+' doesn\'t exists !'):
		verif_input(file4)
	os.remove(file1)
	os.remove(file2)
	os.remove(file3)


####################
# Verif Input XLSX # 
####################

def test_verif_input_xlsx():
	assert verif_input_xlsx(xlsx_file)	
	assert not verif_input_xlsx(false_xlsx_file)
	os.remove(xlsx_file)
	os.remove(false_xlsx_file)


###################
# Verif Input TSV #
###################

def test_verif_input_tsv():
	assert verif_input_tsv(tsv_file)
	os.remove(tsv_file)


############################
# Verif Input Config Merge #
############################

def test_verif_input_config_merge():
	assert verif_input_config_merge(config_file_all1) == None
	with pytest.raises(ValueError, match='File '+file4+' doesn\'t exists !'):
		assert verif_input_config_merge(config_file_all2)
	with pytest.raises(ValueError, match='File Not a good file doesn\'t exists !'):
		assert verif_input_config_merge(false_config)
	os.remove(false_config)
	os.remove(config_file_all1)
	os.remove(config_file_all2)
	os.remove(tsv1)
	os.remove(tsv2)
	os.remove(xlsx)


######################
# Verif Input Config #
######################

def test_verif_input_config():
	assert verif_input_config(config1) == ([vcf1,vcf2],[vcf1,vcf2],[profile1,profile2],[ins1,ins2],[del1,del2])
	with pytest.raises(ValueError, match='File '+empty_file+' is empty !'):
		verif_input_config(config2)
	with pytest.raises(ValueError, match='File '+empty_file+' is empty !'):
		verif_input_config(config3)
	with pytest.raises(ValueError, match='File '+empty_file+' is empty !'):
		verif_input_config(config4)
	with pytest.raises(ValueError, match='File '+empty_file+' is empty !'):
		verif_input_config(config5)
	with pytest.raises(ValueError, match='File '+config6+' is empty !'):
		verif_input_config(config6)
	os.remove(vcf1)
	os.remove(vcf2)
	os.remove(profile1)
	os.remove(profile2)
	os.remove(ins1)
	os.remove(ins2)
	os.remove(del1)
	os.remove(del2)
	os.remove(config1)
	os.remove(config2)
	os.remove(config3)
	os.remove(config4)
	os.remove(config5)
	os.remove(config6)
	os.remove(empty_file)	


########################################
# Verif Supplementary Information File #
########################################

def test_verif_supplementary_information_file():
	#assert verif_input(wrong_infos_file) == None
	with pytest.raises(ValueError, match='File '+wrong1+' is empty !'):
		verif_supplementary_information_file(wrong1, '.')
	verif_supplementary_information_file(wrong2, '.')
	os.remove(wrong1)

