#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../')
import os
import uuid
from pathlib import Path
import pytest
from python_scripts.check_files import verif_input_vcf, verif_output


###############
# Verif Input #
###############

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




################################################################################################################
################################################# Tests ########################################################
################################################################################################################


###############
# Verif Input #
###############

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


