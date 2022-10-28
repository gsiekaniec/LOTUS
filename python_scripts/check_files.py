#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path
from os import access, R_OK, W_OK


def verif_input_vcf (vcf : str):
	'''
	Test if vcf file exist and have a correct format else raise error
	'''
	verif_input(vcf)
	path = Path(vcf) 
	if vcf.upper().endswith('.VCF'):
		with open(vcf ,'r') as f:
			if f.readline().startswith('##fileformat=VCF'):
				pass
			else:
				raise ValueError(f'File {vcf} is not a vcf file or the header ##fileformat=VCF doesn\'t exist !')
	else:
		raise ValueError(f'Extension of file {vcf} is not .vcf !')


def verif_output (output : str):
        '''
        Test if output file exist and if it exist ask user to overwrite it
        '''
        path = Path(output)
        if path.exists():
                authorized_answers = {'Y', 'YES', 'N', 'NO'}
                print(f'Warning output file {output} already exist !\n')
                answer = input('Do you want to overwrite it ? (y/n)\n').upper()
                while answer not in authorized_answers:
                        print(f'{answer} is not a correct answer !\n')
                        answer = input(f'Do you want to overwrite the {output} file ? (y/n)\n').upper()
                if answer == 'N' or  answer == 'NO':
                        print('Retry with another output file !\n')
                        exit(1)
                else:
                        if not access(output, W_OK) :
                                raise ValueError(f'File {output} is not available in writing !')

def verif_input (file : str):	
	'''
	Basic test on input files
	'''
	path = Path(file)
	if path.exists():
		if path.stat().st_size > 0:
			if access(file, R_OK) :
				pass
			else :
				raise ValueError(f'File {file} is not available for reading !')
		else:
			raise ValueError(f'File {file} is empty !')
	else:
		raise ValueError(f'File {file} doesn\'t exists !')

def verif_input_xlsx (file : str):
	try:
		pd.read_excel(file, index_col=0)
		return True
	except ValueError:
		return False

def verif_input_tsv (file : str):
	try:
		pd.read_csv(file, sep='\t')
		return True
	except ValueError:
		return False

def verif_input_config_merge (config : str):
	with open(config, 'r') as f:
		for line in f:
			line = line.strip()
			line = line.split(',')[0]
			if line != '':
				verif_input(line)
				result_xlsx = verif_input_xlsx(line)
				result_tsv = verif_input_tsv(line)
				if result_xlsx:
					pass
				elif result_tsv:
					pass
				else:
					raise ValueError(f'File {line} is not a tsv or excel file !')


def verif_input_config (config : str):
	'''
	Test if the config file and file it contains exists and are ok 
	'''
	verif_input(config)
	config_path = Path(config)
	list_to_fill = None
	vcf_filtered = []
	vcf_pass = []
	snp = []
	insertion = []
	deletion = []
	with open(config ,'r') as f:
		try:
			nb_sample = int(f.readline().strip())
		except ValueError:
			raise ValueError(f'Config file {config} does\'t start with the number of samples!')
		for line in f:
			line = line.strip()
			if line != '':
				if line[0] != '#':
					if line[0] == '-':
						if '- filtered:' in line:
							list_to_fill = 'filtered'
						elif '- pass:' in line:
							list_to_fill = 'pass'
						elif '- profile:' in line:
							list_to_fill = 'profile'
						elif '- insertion:' in line:
							list_to_fill = 'insertion'
						elif '- deletion:' in line:
							list_to_fill = 'deletion'
					if line[0] != '-':
						if list_to_fill == 'filtered':
							verif_input_vcf(line)
							vcf_filtered.append(line)
						elif list_to_fill == 'pass':
							verif_input_vcf(line)
							vcf_pass.append(line)
						elif list_to_fill == 'profile':
							verif_input(line)
							snp.append(line)
						elif list_to_fill == 'insertion':
							verif_input(line)
							insertion.append(line)
						elif list_to_fill == 'deletion':			
							verif_input(line)							
							deletion.append(line)
	if len(vcf_filtered) != nb_sample: 
		raise ValueError(f'Wrong number of filtered vcf files in config file {config} !')
	elif len(vcf_pass) != nb_sample:
		raise ValueError(f'Wrong number of pass vcf files in config file {config} !')
	elif len(snp) != nb_sample:
		raise ValueError(f'Wrong number of snp profile tsv files in config file {config} !')
	elif len(insertion) != nb_sample:
		raise ValueError(f'Wrong number of insertion count tsv files in config file {config} !')
	elif len(deletion) != nb_sample:
		raise ValueError(f'Wrong number of deletion count tsv files in config file {config} !')
	return vcf_filtered, vcf_pass, snp, insertion, deletion










