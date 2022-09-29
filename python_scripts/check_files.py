#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from os import access, R_OK, W_OK


def verif_input_vcf (vcf : str):
        '''
        Test if vcf file exist and have a correct format else raise error
        '''
        path = Path(vcf)
        if path.exists():
                if path.stat().st_size > 0:
                        if access(vcf, R_OK) :
                                if vcf.upper().endswith('.VCF'):
                                        with open(vcf ,'r') as f:
                                                if f.readline().startswith('##fileformat=VCF'):
                                                        pass
                                                else:
                                                        raise ValueError(f'File {vcf} is not a vcf file or the header ##fileformat=VCF doesn\'t exist !')
                                else:
                                        raise ValueError(f'Extension of file {vcf} is not .vcf !')
                        else :
                                raise ValueError(f'File {vcf} is not available for reading !')
                else:
                        raise ValueError(f'File {vcf} is empty !')
        else:
                raise ValueError(f'File {vcf} doesn\'t exists !')


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


def verif_input_config (config : str):
	'''
	Test if the config file and file it contains exists and are ok 
	'''
	config_path = Path(config)
	if config_path.exists():
		if config_path.stat().st_size > 0:
			if access(config, R_OK) :
				with open(config ,'r') as f:
					try:
						line1 = int(f.readline().strip())
					except ValueError:
						raise ValueError(f'Config file {config} does\'t start with the number of samples!')
					#if line1.startswith('###VCF files') or line1.startswith('### VCF files') :
					#	pass
					#else:
					#	raise ValueError(f'File {vcf} is not a vcf file or the header ##fileformat=VCF doesn\'t exist !')
			else :
				raise ValueError(f'Config file {config} is not available for reading !')
		else:
			raise ValueError(f'Config file {config} is empty !')
	else:
		raise ValueError(f'Config file {config} doesn\'t exists !')










