#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle as pk
from pathlib import Path
import sys
import os
import numpy as np
import logging
import warnings
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
from more_itertools import powerset, unique_everseen
from upsetplot import from_memberships
from upsetplot import UpSet
from python_scripts.check_files import verif_output, verif_input_config_merge, verif_input_xlsx, verif_input_tsv, verif_input, verif_supplementary_information_file
from python_scripts.toppgene_api import ToppGene_GEOA
from python_scripts.panther_api import Panther_GEOA
from python_scripts.read_vcf import read_vcf
from python_scripts.chromosomes_plot import create_chromosomes_plot
from python_scripts.path_modification import true_stem
import matplotlib.pyplot as plt


def get_informations_for_genes(info_file, logger):
	'''
	Extract informations from the external cancer databases file.
	'''
	df = pd.read_excel(info_file, index_col=1)
	df = df.drop(['Ordre'], axis=1)
	df.set_axis([source.split(' Info')[0] for source in df.columns], axis="columns", inplace=True)
	print(f'Extract informations from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	logger.info(f'Extract informations from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	return (df)	
		

def read_config(config_file):
	'''
	Read the merge config file (composed of one gene.tsv compare file per line).
	'''
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
	'''
	Extract files number and name from the config file.
	'''
	names = [i.strip().split(',')[1] if i.strip() != '' and len(i.strip().split(',')) > 1 else true_stem(i.split(',')[0]) for i in open(config_file, 'r').readlines()]
	nb = len(names)
	return nb, names		


def add_variants_to_dictionnary (d : {}, gene_name : str, pos : int, type : str, row):
	'''
	Add a variant in count dictionnary.
	'''
	for i in row.iloc[pos].split('|'):
		if not i in d[gene_name][type].keys():
			d[gene_name][type][i]=1
		else:
			d[gene_name][type][i]+=1


def merge_results(it_files, category, output, upset_output, infos, cytoband_file, chromosomes_output, step, weakness_threshold, min_subset_size, max_subset_size, min_degree, max_degree, nb_files, enrichment, logger):
	'''
	
	'''

	# Get genes infos if not None
	if infos:
		infos_df = get_informations_for_genes(infos, logger)

	names = set()
	d = {}										# gene dictionnary to create tsv union file (list format)
	dico = {'chr':'', 'start':'', 'end':'', 'weakness':[], 'gb1':{}, 'cb1':{}, 'pb1':{}, 'gb2':{}, 'cb2':{}, 'pb2':{}, 'samples':[]}
	print('Read files...')
	pbar_files = tqdm(total=nb_files)
	for name, df in it_files:
		name = true_stem(name.strip())
		if not name in names:
			names.add(name)
		b1 = df.columns.values.tolist()[5]
		b2 = df.columns.values.tolist()[9]
		logger.info(f'Processing of samples B1: {b1} and B2: {b2}')
		for index, row in df.iterrows():
			
			if index not in d.keys():
				d[index]=deepcopy(dico)
			if row.iloc[5] != 0:
				add_variants_to_dictionnary (d, index, 6, 'gb1', row)
				if not pd.isna(row.iloc[7]) and ''.join(row.iloc[7].split('|')) != '':
					add_variants_to_dictionnary (d, index, 7, 'cb1', row)
				if not pd.isna(row.iloc[8]) and ''.join(row.iloc[8].split('|')) != '':
					add_variants_to_dictionnary (d, index, 8, 'pb1', row)
			if row.iloc[9] != 0:
				add_variants_to_dictionnary (d, index, 10, 'gb2', row)
				if not pd.isna(row.iloc[11]) and ''.join(row.iloc[11].split('|')) != '':
					add_variants_to_dictionnary (d, index, 11, 'cb2', row)
				if not pd.isna(row.iloc[12]) and ''.join(row.iloc[12].split('|')) != '':
					add_variants_to_dictionnary (d, index, 12, 'pb2', row)
			d[index]['weakness'].append(float(row['Gene weakness (in %)']))
			d[index]['samples'].append(name)
			if d[index]['chr'] == '':
				d[index]['chr']=row['Chromosome']
			if d[index]['start'] == '':
				d[index]['start'] = int(row['Gene position start'])
			if d[index]['end'] == '':	
				d[index]['end'] = int(row['Gene position end'])
		pbar_files.update(1)
	pbar_files.close()


	######################################
	# Columns name for the union.tsv file #
	column_names = ['Gene Symbol', 'Chromosome', 'Gene position start', 'Gene position end', 'Nb samples', 'Sample names',  'Tumour burden union', 'Mean gene weakness', 'Nb Mutation b1', 'g.TPn union', 'c.TPn union', 'p.TPn union', 'Nb Mutation TPn+1', 'g.TPn+1 union', 'c.TPn+1 union', 'p.TPn+1 union']

	d2 = {}										# gene dictionnary to create tsv union file (mean computation for each list)
	data = {}									# dictionnary containing genes list for each sample - use to create the upset plot
	genes_pos_for_chromosomes = {}
	genes_pos_for_chromosomes_b1 = {}
	genes_pos_for_chromosomes_b2 = {}
	for k, v in d.items():
		if np.mean(v['weakness']) <= weakness_threshold:
			if not v['chr'] in genes_pos_for_chromosomes.keys():
				genes_pos_for_chromosomes[v['chr']]=[]
			if not v['chr'] in genes_pos_for_chromosomes_b1.keys():
				genes_pos_for_chromosomes_b1[v['chr']]=[]
			if not v['chr'] in genes_pos_for_chromosomes_b2.keys():
				genes_pos_for_chromosomes_b2[v['chr']]=[]

			genes_pos_for_chromosomes[v['chr']] += [(v['start'],v['end'])]*int(len(v['samples']))
			if v['gb1'] != {} and v['gb2'] != {}:
				genes_pos_for_chromosomes_b1[v['chr']] += [(v['start'],v['end'])]*int(sum(list(v['gb1'].values()))) 
				genes_pos_for_chromosomes_b2[v['chr']] += [(v['start'],v['end'])]*int(sum(list(v['gb2'].values())))
			elif v['gb1'] != {}:
				genes_pos_for_chromosomes_b1[v['chr']] += [(v['start'],v['end'])]*int(sum(list(v['gb1'].values())))
			elif v['gb2'] != {}:
				genes_pos_for_chromosomes_b2[v['chr']] += [(v['start'],v['end'])]*int(sum(list(v['gb2'].values())))
		

		if not k in d2.keys():
			d2[k]={}
		if np.mean(v['weakness']) <= weakness_threshold:
			for sample in v['samples']:
				if not sample in data.keys():
					data[sample]=set()
				data[sample].add(k)
		d2[k]['Chromosome'] = v['chr']
		d2[k]['Gene position start'] = v['start']
		d2[k]['Gene position end'] = v['end']
		d2[k]['Nb samples'] = int(len(v['samples']))
		d2[k]['Sample names'] = ';'.join(v['samples'])
		d2[k]['Tumour burden union'] = len(set(v['gb1'].keys()).union(set(v['gb2'].keys())))
		d2[k]['Mean gene weakness'] = np.mean(v['weakness'])
		d2[k]['Nb Mutation TPn'] = len(set(v['gb1'].keys()))
		d2[k]['g.TPn union'] = '|'.join([str(k)+'('+str(v)+')' for k,v in v['gb1'].items()])
		d2[k]['c.TPn union'] = '|'.join([str(k)+'('+str(v)+')' for k,v in v['cb1'].items()])
		d2[k]['p.TPn union'] = '|'.join([str(k)+'('+str(v)+')' for k,v in v['pb1'].items()])
		d2[k]['Nb Mutation TPn+1'] = len(set(v['gb2'].keys()))
		d2[k]['g.TPn+1 union'] = '|'.join([str(k)+'('+str(v)+')' for k,v in v['gb2'].items()])
		d2[k]['c.TPn+1 union'] = '|'.join([str(k)+'('+str(v)+')' for k,v in v['cb2'].items()])
		d2[k]['p.TPn+1 union'] = '|'.join([str(k)+'('+str(v)+')' for k,v in v['pb2'].items()])
	

	if cytoband_file:
		##### Create the Chromosome plot
		create_chromosomes_plot(genes_pos_for_chromosomes_b1, genes_pos_for_chromosomes_b2, genes_pos_for_chromosomes, cytoband_file, chromosomes_output, step, logger)
		#####

	if upset_output:
		if len(names) < 16:					# Actually upset plot can not be calculated for more than 15 samples
			##### Create the UpSetPlot
			create_upsetplot(data, category, upset_output, names, min_subset_size, max_subset_size, min_degree, max_degree, weakness_threshold, logger)
			#####

	del d
	df = pd.DataFrame(d2).T														# get pandas Dataframe from d2
	del d2
	df.index.name = 'Gene'														# add a name to the index column


	df = df.sort_values(by=['Nb samples', 'Mean gene weakness', 'Tumour burden union', 'Gene'], ascending = [False, True, False, True]) 		# sort columns
	df['Nb samples'] = df['Nb samples'].astype('int')												 
	df[['Tumour burden union', 'Mean gene weakness', 'g.TPn union', 'c.Tpn union', 'p.TPn union', 'Nb Mutation TPn', 'g.TPn+1 union', 'c.TPn+1 union', 'p.TPn+1 union', 'Nb Mutation TPn+1']] = np.round(df[['Tumour burden union','Mean gene weakness', 'g.TPn union', 'c.TPn union', 'p.TPn union', 'Nb Mutation TPn', 'g.TPn+1 union', 'c.TPn+1 union', 'p.TPn+1 union', 'Nb Mutation TPn+1']], 2)

	if infos:
		# Add additional cancer centric genes informations
		df = df.join(infos_df)
		# Drop empty informational columns
		empty_cols = [col for col in df.columns if df[col].isnull().all()]
		df.drop(empty_cols, axis=1, inplace=True)

	genes_list = list(df.index)

	##########################################
	# Save genes list files (.xlsx and .tsv) #

	print(f'Save genes list files in {Path(output).with_suffix(".MutatedGenes.xlsx")} and {Path(output).with_suffix(".MutatedGenes.tsv")}')
	logger.info(f'Save genes list files in {Path(output).with_suffix(".MutatedGenes.xlsx")} and {Path(output).with_suffix(".MutatedGenes.tsv")}')

	if not str(output).endswith('.xlsx'):
		output = Path(output).with_suffix('.xlsx')
	df.to_excel(output)
	if not str(output).endswith('.tsv'):
		output = Path(output).with_suffix('.tsv')
	df.to_csv(output, sep='\t')

	################################################
	# Biological process enrichment using the genes list with the ToppGene and Panther API

	if enrichment:
		toppgene_name = str(Path(output).resolve().parent)+'/'+str(true_stem(output))+'_ToppGene_enrichment'
		panther_name =  str(Path(output).resolve().parent)+'/'+str(true_stem(output))+'_Panther_enrichment'
		ToppGene_GEOA(genes_list, toppgene_name, logger)
		Panther_GEOA(genes_list, panther_name, logger)


def create_upsetplot(data, category, upset_output, names, min_subset_size, max_subset_size, min_degree, max_degree, weakness_threshold, logger):
	'''
	Create an Upsetplot showing the number of common genes to the different sample set
	'''

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

	info = f'minsb: {min_subset_size}; maxsb: {max_subset_size}; mind: {min_degree}; maxd: {max_degree}; w: {weakness_threshold}'
	plt.figtext(0.5, 0.01, info, ha="center", fontsize=7, bbox={"facecolor":"lightgray", "alpha":0.5, "pad":2})

	if not str(upset_output).endswith('.svg'):
		upset_output = Path(upset_output).with_suffix('.svg')
	print(f'Create UpSetPlot in {upset_output} !')
	logger.info(f'Create UpSetPlot in {upset_output} !')
	plt.savefig(upset_output,format='svg')	

	if not str(upset_output).endswith('.png'):
		upset_output = Path(upset_output).with_suffix('.png')
	print(f'Create UpSetPlot in {upset_output} !')
	logger.info(f'Create UpSetPlot in {upset_output} !')
	plt.savefig(upset_output,format='png')


def main(args):

	config = args.conf
	output = args.output
	upset_output = args.upset
	min_subset_size_upset_plot = args.min_subset_size
	max_subset_size_upset_plot = args.max_subset_size
	min_degree = args.min_degree
	max_degree = args.max_degree
	weakness_threshold = args.threshold
	enrichment = args.enrichment

	sup_infos = args.agi

	# chromosomes plot variables
	cytoband_file = args.cyto
	chromosomes_output = args.chromosomes_output
	step = args.step

	current_directory = os.getcwd()
	working_directory = Path(config).parent.absolute()


	# Logger configuration

	logger = logging.getLogger('LOTUS merge')
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
		if not output.endswith('MutatedGenes.tsv'):
			output = Path(output).with_suffix('.MutatedGenes.tsv')
		verif_output(Path(output).with_suffix('.xlsx'))
		verif_output(Path(output).with_suffix('.tsv'))
		logger.info('- Output file ok -')
	except ValueError:
		print (f'Problem with {output}: ', sys.exc_info()[0])
		logger.error('- Problem with {output} -')
		exit(1)
	if upset_output:
		try:
			logger.info('Verification of {upset_output}')
			verif_output(upset_output)
			logger.info('- Output file ok -')
		except ValueError:
			print (f'Problem with {upset_output}: ', sys.exc_info()[0])
			logger.error('- Problem with {upset_output} -')
			exit(1)
	if cytoband_file:
		try:
			logger.info('Verification of {cytoband_file}')
			verif_input(cytoband_file)
			try:
				logger.info('Verification of {chromosomes_output}')
				verif_output(chromosomes_output)
				logger.info('- Output file ok -')
			except ValueError:
				print (f'Problem with {chromosomes_output}: ', sys.exc_info()[0])
				logger.error('- Problem with {chromosomes_output} -')
				exit(1)
			logger.info('- Input file ok -')
		except ValueError:
			print (f'Problem with {cytoband_file}: ', sys.exc_info()[0])
			logger.error('- Problem with {cytoband_file} -')
			exit(1)

	if step < 500000:
		warnings.warn("A step value below 500000 may cause a long calculation time !", RuntimeWarning)
		logger.warning('A step value below 500000 may cause a long calculation time !')

	infos = None
	if sup_infos:
		print('\nSearch for file Lotus_ExternalBases_202301.xlsx ...\n')
		infos = 'Lotus_ExternalBases_202301.xlsx'
		logger.info('Verification of {infos}')
		try:
			infos = verif_supplementary_information_file(infos, current_directory)
		except ValueError:
			print (f'Problem with {infos}: ', sys.exc_info()[0])
			logger.error('- Problem with {infos} -')			
			exit(1)

	# Start
	print('--Start--')

	logger.info('**************************************************************************************************************')
	logger.info('*** LOTUS merging module ***')
	no_argument = ''
	if sup_infos:
		no_argument+=' --additional_gene_information'
	if enrichment:
		no_argument+=' --enrichment'
	logger.info(f'** cmd line : python lotus.py merge -c {str(config)} -o {str(output)} -cyto {cytoband_file} -step {step} -co {chromosomes_output} -u {str(upset_output)} -w {str(weakness_threshold)} -minsb {str(min_subset_size_upset_plot)} -maxsb {str(max_subset_size_upset_plot)} -mind {str(min_degree)} -maxd {str(max_degree)}'+str(no_argument)+' **')
	logger.info('* Start merging *')
	logger.info(f'Working directory (config file folder) : {working_directory}')
	logger.info(f'Current directory : {Path().absolute()}')

	nb_files, names = get_nb_files_and_file_names(config)

	category = None
	if nb_files < 16:
		if upset_output:
			category = [list(i) for i in list(powerset(unique_everseen(names))) if i != ()]
	else:
		warnings.warn("Upset plot not computed because of the combination explosion !", RuntimeWarning)
		logger.warning('Upset plot not computed because of the combination explosion !')

	print(f'Merging {nb_files} files...')
	logger.info(f'Merging {nb_files} files !')
	it_files = read_config(config)
	
	print('Start merging...')
	merge_results(it_files, category, output, upset_output, infos, cytoband_file, chromosomes_output, step, weakness_threshold, min_subset_size_upset_plot, max_subset_size_upset_plot, min_degree, max_degree, nb_files, enrichment, logger)

	logger.info('* End merging *')
	logger.info('**************************************************************************************************************')

	# End








