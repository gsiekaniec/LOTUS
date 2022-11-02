#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle as pk
from pathlib import Path
import sys
import os
import numpy as np
import logging
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
from more_itertools import powerset, unique_everseen
from upsetplot import from_memberships
from python_scripts.check_files import verif_output, verif_input_config_merge, verif_input_xlsx, verif_input_tsv
from python_scripts.ToppGene_api import ToppGene_GEOA
from python_scripts.Panther_api import Panther_GEOA
from python_scripts.read_vcf import read_vcf
from python_scripts.path_modification import true_stem
from upsetplot import from_memberships
from upsetplot import UpSet
import matplotlib.pyplot as plt


def read_config(config_file):
	with open(config_file, 'r') as f:
		for line in f:
			line = line.strip()
			if line != '':
				line = line.split(',')
				if verif_input_xlsx(line[0]):
					data = pd.read_excel(line[0], index_col=0)					 
				elif verif_input_tsv(line[0]):
					data = pd.read_csv(line[0], sep='\t', index_col=0)
				if len(line) > 1:
					yield line[1], data	
				else:
					yield line[0], data


def get_nb_files_and_file_names(config_file):
	names = [i.strip().split(',')[1] if i.strip() != '' and len(i.strip().split(',')) > 1 else true_stem(i.split(',')[0]) for i in open(config_file, 'r').readlines()]
	nb = len(names)
	return nb, names		


def merge_results(it_files, category, output, upset_output, weakness_threshold, min_subset_size, max_subset_size, min_degree, max_degree, logger):
	
	names = set()
	d = {}
	data = {}
	dico = {'weakness':[], 'burden':[], 'B1':[], 'B2':[], 'samples':[]}
	for name, df in it_files:
		name = true_stem(name.strip())
		if not name in names:
			names.add(name)
		b1 = df.columns.values.tolist()[2]
		b2 = df.columns.values.tolist()[3]
		print(f'Processing of samples B1: {b1} and B2: {b2}')
		logger.info(f'Processing of samples B1: {b1} and B2: {b2}')
		for index, row in df.iterrows():
			if index not in d.keys():
				d[index]=deepcopy(dico)
			d[index]['B1'].append(int(row.iloc[2]))
			d[index]['B2'].append(int(row.iloc[3]))
			d[index]['burden'].append(int(row['Tumor burden (symmetrical difference)']))
			d[index]['weakness'].append(float(row['Variant weakness (in %)']))
			d[index]['samples'].append(name)		

	column_names = ['Gene', 'Nb samples', 'Mean weakness', 'Nb B1 variants (mean)', 'Nb b2 variants (mean)', 'Sum B1+B2']

	d2 = {}
	data = {}
	for k, v in d.items():
		if not k in d2.keys():
			d2[k]={}
		if np.mean(v['weakness']) <= weakness_threshold:
			for sample in v['samples']:
				if not sample in data.keys():
					data[sample]=set()
				data[sample].add(k)
		d2[k]['Nb samples'] = int(len(v['samples']))
		d2[k]['Mean weakness'] = np.mean(v['weakness'])
		d2[k]['Nb B1 variants (mean)'] = np.mean(v['B1'])
		d2[k]['Nb B2 variants (mean)'] = np.mean(v['B2'])
		d2[k]['Sum B1+B2'] = sum([np.mean(v['B1']),np.mean(v['B2'])])

	##### Create the UpSetPlot
	create_upsetplot(data, category, upset_output, names, min_subset_size, max_subset_size, min_degree, max_degree, weakness_threshold, logger)
	#####

	del d
	df = pd.DataFrame(d2).T														# get pandas Dataframe from d2
	del d2
	df.index.name = 'Gene'														# add a name to the index column

	df = df.sort_values(by=['Nb samples', 'Mean weakness', 'Sum B1+B2', 'Gene'], ascending = [False, True, False, True]) 		# sort columns
	df['Nb samples'] = df['Nb samples'].astype('int')												 
	df[['Mean weakness', 'Nb B1 variants (mean)', 'Nb B2 variants (mean)', 'Sum B1+B2']] = np.round(df[['Mean weakness', 'Nb B1 variants (mean)', 'Nb B2 variants (mean)', 'Sum B1+B2']], 2)

	# Save gene list files

	if not str(output).endswith('.xlsx'):
		output = Path(output).with_suffix('.xlsx')
	df.to_excel(output)
	if not str(output).endswith('.tsv'):
		output = Path(output).with_suffix('.tsv')
	df.to_csv(output, sep='\t')


def create_upsetplot(data, category, upset_output, names, min_subset_size, max_subset_size, min_degree, max_degree, weakness_threshold, logger):

	if max_subset_size == 0:
		max_subset_size = max([len(v) for v in data.values()])

	no_gene = names-set([k for k in data.keys()])
	if no_gene != set():
		for sample in no_gene:
			data[sample]=set()
			print(f'Warning : sample {sample} don\'t have gene !')
			logger.warning(f'Sample {sample} don\'t have gene !')

	category.sort(key=len)
	if max_degree == 0:
		max_degree = len(category[-1])
	category = [c for c in category if len(c) >= min_degree and len(c) <= max_degree]

	u = set()
	for v in data.values():
	        if u == set():
	                u = v
	        else:
	                u = u.union(v)
	union_data = dict([(k, u) for k in data.keys()])

	data2 = []
	cat = []
	to_suppr = set()
	for c in category[::-1]:
		save_set = None
		for v in c:
			if save_set == None:
				save_set = data[v]
			else:
				save_set = save_set.intersection(data[v])
		if save_set-to_suppr != set() and len(save_set-to_suppr) >= min_subset_size:
			if max_subset_size != 0:
				if len(save_set-to_suppr) <= max_subset_size:
					cat.append(c)
					data2.append(len(save_set-to_suppr))
			else :
				cat.append(c)
				data2.append(len(save_set-to_suppr))
		to_suppr = to_suppr.union(save_set)

	df = from_memberships(cat,data=data2)

	UpSet(df, show_counts=True).plot()

	plt.title('Number of genes specific to the samples sets')
	plt.ylabel('Gene number')

	info = f'minsb : {min_subset_size}, maxsb : {max_subset_size}, mind : {min_degree}, maxd : {max_degree}, wt : {weakness_threshold}'
	plt.figtext(0.5, 0.01, info, ha="center", fontsize=7, bbox={"facecolor":"lightgray", "alpha":0.5, "pad":2})

	if not str(upset_output).endswith('.svg'):
		upset_output = Path(upset_output).with_suffix('.svg')
	print(f'Create UpSetPlot in {upset_output} !')
	logger.info(f'Create UpSetPlot in {upset_output} !')
	
	plt.savefig(upset_output,format='svg')


def main(args):

	config = args.conf
	output = args.output
	upset_output = args.upset
	min_subset_size_upset_plot = args.min_subset_size
	max_subset_size_upset_plot = args.max_subset_size
	min_degree = args.min_degree
	max_degree = args.max_degree
	weakness_threshold = args.threshold
	working_directory = Path(config).parent.absolute()

	# Logger configuration

	logger = logging.getLogger('g-LOTUS merge')
	logger.setLevel(logging.DEBUG)

	fh = logging.FileHandler(args.log)
	fh.setLevel(logging.DEBUG)

	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)

	logger.addHandler(fh)

	# Verification of given arguments

	verif_input_config_merge(config)
	try:
		logger.info('Verification of {output}')
		verif_output(output)
		logger.info('- Output file ok -')
	except ValueError:
		print (f'Problem with {output}: ', sys.exc_info()[0])
		logger.error('- Problem with {output} -')
		exit(1)
	try:
		logger.info('Verification of {upset_output}')
		verif_output(upset_output)
		logger.info('- Output file ok -')
	except ValueError:
		print (f'Problem with {upset_output}: ', sys.exc_info()[0])
		logger.error('- Problem with {upset_output} -')
		exit(1)

	# Start

	logger.info('**************************************************************************************************************')
	logger.info('*** g-LOTUS merging module ***')
	logger.info(f'** cmd line : python lotus.py merge -c {str(config)} -o {str(output)} -u {str(upset_output)} -wt {str(weakness_threshold)} -minsb {str(min_subset_size_upset_plot)} -maxsb {str(max_subset_size_upset_plot)} -mind {str(min_degree)} -maxd {str(max_degree)} **')
	logger.info('* Start merging *')
	logger.info(f'Working directory (config file folder) : {working_directory}')
	logger.info(f'Current directory : {Path().absolute()}')

	nb_files, names = get_nb_files_and_file_names(config)
	category = [list(i) for i in list(powerset(unique_everseen(names))) if i != ()]
	print(f'Merging {nb_files} files...')
	logger.info(f'Merging {nb_files} files !')
	it_files = read_config(config)
	
	print('Start merging...')
	merge_results(it_files, category, output, upset_output, weakness_threshold, min_subset_size_upset_plot, max_subset_size_upset_plot, min_degree, max_degree, logger)

	logger.info('* End merging *')
	logger.info('**************************************************************************************************************')

	# End






