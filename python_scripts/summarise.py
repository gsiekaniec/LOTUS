#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pyfastx
import pickle as pk
from pathlib import Path
import sys
import os
import numpy as np
import logging
from tqdm import tqdm
from itertools import dropwhile
from collections import OrderedDict, Counter
from copy import deepcopy
import warnings
import pandas as pd
from python_scripts.check_files import verif_input_vcf, verif_output 
from python_scripts.toppgene_api import ToppGene_GEOA
from python_scripts.panther_api import Panther_GEOA
from python_scripts.read_vcf import read_vcf
import matplotlib.pyplot as plt

TRANSLATE = {'A>C': 'T>G', 'A>G': 'T>C', 'A>T': 'T>A', 'C>A': 'C>A', 'C>G': 'C>G', 'C>T': 'C>T',
             'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A', 'T>A': 'T>A', 'T>C': 'T>C', 'T>G': 'T>G'}
TAB = str.maketrans("ACGT", "TGCA")
	
def create_snp_dict():
	'''
	Create the dictionnary to count the SNP variant and create the vcf SNP profile. Only the C>* and T>* value are used because the variant sens isn't known.
	Input : None:
	Output : dictionnary with empty counts
	'''
	c_words = {}
	t_words = {}
	for i in 'ATCG':
		word = i
		for j in 'CT':
			word2 = word + j
			for k in 'ATCG':
				word3 = word2 + k
				if word3[1] == 'C':
					c_words[word3] = 0
				else:
					t_words[word3] = 0
	save = {}
	for snp in set(TRANSLATE.values()):
		if snp[0] == 'C':
			save[snp] = deepcopy(c_words)
		else:
			save[snp] = deepcopy(t_words)
	return save


def create_dataframe_from_gene(d):
	'''
	Take a dictionary containing the genes list with information for each genes and transform it in pandas dataframe
	Input : dictionary containing genes
	Output : dataframe corresponding to the dictionary
	'''
	
	col = ['tumor_burden','details_(snp,dnp,tnp,np,insertion,deletion)','chromosome','ref','alt_variant(s)','position(s)'] 				# Columns name
	id = 'gene_name'																# Index name
	d = OrderedDict(sorted(d.items()))
	df = pd.DataFrame.from_dict(d, orient='index', columns=col)
	df.index.name = id
	
	#Modification of the dataframe using a copy of it
	
	df2 = pd.DataFrame()
	df2['details_(snp,dnp,tnp,np,insertion,deletion)'] = [','.join(map(str, l)) for l in df['details_(snp,dnp,tnp,np,insertion,deletion)']]
	df['details_(snp,dnp,tnp,np,insertion,deletion)'] = df2['details_(snp,dnp,tnp,np,insertion,deletion)'].values
	del df2
	df2 = pd.DataFrame()
	df2['ref'] = [','.join(map(str, l)) for l in df['ref']]
	df['ref'] = df2['ref'].values
	del df2
	df2 = pd.DataFrame()
	df2['alt_variant(s)'] = [ ','.join([str(l2) if len(l2) > 1 else str(l2[0]) for l2 in l]) for l in df['alt_variant(s)']]
	df['alt_variant(s)'] = df2['alt_variant(s)'].values
	del df2
	df2 = pd.DataFrame()
	df2['position(s)'] = [','.join(map(str, l)) for l in df['position(s)']]
	df['position(s)'] = df2['position(s)'].values
	del df2
	
	return df


def create_ordered_dataframe(d):
	'''
	Take a dictionary contaning snp counts, ordered it using it and return the pandas dataframe corresponding
	Input : dictionary contaning snp counts
	Output : dataframe corresponding to the dictionary
	'''
	for k, v in d.items():
		d[k] = OrderedDict(sorted(v.items(), key=lambda x: x[0]))
	d = OrderedDict(sorted(d.items(), key=lambda t: t[0]))
	df = pd.DataFrame.from_dict(d, orient="index").stack().to_frame()
	return df


def graph_snp(d, name, vcf_name, logger):
	'''
        Creation of the profile graph 
        Input : a dictionnary containing the count for all kind of SNP (d), the output file name and the logger
        Output : write a .svg file 
        '''

	print(f'Draw profile in {name}...')
	logger.info(f'Draw profile in {name}')

	df = create_ordered_dataframe(d)
	df.columns = [str(Path(vcf_name).stem)]
	
	#color for graph

	white_96 = ['white']*96
	color_6 = ['darkblue', 'blue', 'lightblue', 'darkgreen', 'green', 'lightgreen']
	color_96 = []
	for i in color_6:
		color_96 += [i]*16
	bars = [index[1] for index, row in df.iterrows()]
	height = [float(row) for index, row in df.iterrows()]
	group = []
	for index, row in df.iterrows():
		if index[0] not in group:
			group.append(index[0])
	x_pos = np.arange(len(bars))

	# Create bars
	plt.figure(figsize=(15, 10))
	ax1 = plt.subplot(1,1,1)
	plt.bar(x_pos, height, color=color_96)
	
	# Create names on the x-axis
	plt.xticks(x_pos, bars,rotation=90, fontsize=8)
	ax1.set_ylabel('Fraction of each mutation-type')

	# Set scond x-axis
	ax2 = ax1.twiny()
	newlabel = group
	newpos = [8, 24, 40, 56, 72, 88]
	ax2.set_xticks(newpos)
	ax2.set_xticklabels(newlabel)
	# set the position of the second x-axis to bottom
	ax2.xaxis.set_ticks_position('bottom')
	ax2.xaxis.set_label_position('bottom')
	ax2.spines['bottom'].set_position(('outward', 36))
	ax2.set_xlabel('Mutation types')
	ax2.set_xlim(ax1.get_xlim())
	plt.tick_params(
		axis='x',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		bottom=False,      # ticks along the bottom edge are off
		top=False,         # ticks along the top edge are off
		labelbottom=True)  # label along the bottom edge are on
	for xtick, color in zip(ax2.get_xticklabels(), color_6):
		xtick.set_color(color)
	ax2.spines['bottom'].set_visible(False)
	plt.savefig(name)

	name = Path(name).with_suffix('.png')
	plt.savefig(name)

	plt.close()

	print(f'Save profile values in {Path(name).with_suffix(".tsv")}...')
	logger.info(f'Save profile values in {Path(name).with_suffix(".tsv")}')

	df.to_csv(Path(name).with_suffix(".tsv"), sep='\t')


def graph_indel(deletion, insertion, name, vcf_name, logger):
	'''
	Creation of the indel size graph
	Input : two dictionnaries containing the count for insertion/deletion of different size, the output file name and the logger
	Output : write a .svg file
	'''

	insert = True
	delet = True
	if deletion == Counter() and insertion == Counter():						# If no Indel
		print('Warning ! No indel in vcf !')
		logger.warning(f'No indel in vcf !')
		return None 
	elif insertion == Counter():									# If no insertion
		print('Warning ! No insertion in vcf !')
		logger.warning(f'No insertion in vcf !')
		insert = False
	elif deletion == Counter():									# If no deletion
		print('Warning ! No deletion in vcf !')
		logger.warning(f'No deletion in vcf !')
		delet = False
	
	print(f'Draw indel size barplot in {name}...')
	logger.info(f'Draw indel size barplot in {name}')
	
	width = 0.25

	sum_del = 0
	sum_ins = 0
	if delet:
		sum_del = sum(deletion.values())
		df_del = pd.DataFrame.from_dict(deletion, orient='index').sort_index()
		height_del = list([float(i)/sum_del for i in df_del[0]])
		r_del = list(df_del.index)
		max_del = max([float(i) for i in r_del])
	if insert:
		sum_ins = sum(insertion.values())
		df_ins = pd.DataFrame.from_dict(insertion, orient='index').sort_index()
		height_ins = list([float(i)/sum_ins for i in df_ins[0]])	
		r_ins = list(df_ins.index)
		max_ins = max([float(i) for i in r_ins])

	if delet and insert:
		maximum = max([int(max_del)+1, int(max_ins)+1])
	elif delet:
		maximum = int(max_del)+1
	elif not delet:
		maximum = int(max_ins)+1
	x = [0]+[i+1 for i in range(maximum)]
	
	if delet and insert:
		plt.bar([float(i)-(width*0.65) for i in r_del], height_del, color = 'k', width = width, edgecolor = 'k', label='Deletion')
		plt.bar([float(i)+(width*0.65) for i in r_ins], height_ins, color = 'r', width = width, edgecolor = 'r', label='Insertion')
	elif not insert:
		plt.bar([float(i) for i in r_del], height_del, color = 'k', width = width, edgecolor = 'k', label='Deletion')
	elif not delet:
		plt.bar([float(i) for i in r_ins], height_ins, color = 'r', width = width, edgecolor = 'r', label='Insertion')
			

	plt.xticks(x, fontsize=9)
	plt.yticks()

	if  maximum > 10:
	
		##############################################
		# Adaptation of the figure size to the x range

		plt.gca().margins(x=0)
		plt.gcf().canvas.draw()
		tl = plt.gca().get_xticklabels()
		maxsize = max([t.get_window_extent().width for t in tl])+1
		m = 0.5 # inch margin
		s = maxsize/plt.gcf().dpi*maximum+2*m
		margin = m/plt.gcf().get_size_inches()[0]
		plt.gcf().subplots_adjust(left=margin, right=1.-margin)
		plt.gcf().set_size_inches(s*1.5, plt.gcf().get_size_inches()[1])

		plot_margin = 1
		x0, x1, y0, y1 = plt.axis()
		plt.axis((x0 - plot_margin, x1 + plot_margin*2, y0, y1))

	plt.xlabel("Indel size")
	plt.ylabel("Percentage")

	plt.legend()
	plt.savefig(name)

	name = Path(name).with_suffix('.png')
	plt.savefig(name)

	plt.close()	


	print(f'Save indel counts in {Path(name).with_suffix(".tsv")}...')
	logger.info(f'Save indel counts in {Path(name).with_suffix(".tsv")}')

	if delet and insert:
		df_del.columns = [str(Path(vcf_name).stem)]
		df_ins.columns = [str(Path(vcf_name).stem)]
		df_del.to_csv(Path(name).with_suffix(".deletion.tsv"), sep='\t')
		df_ins.to_csv(Path(name).with_suffix(".insertion.tsv"), sep='\t')
	elif delet:
		df_del.columns = [str(Path(vcf_name).stem)]
		df_del.to_csv(Path(name).with_suffix(".deletion.tsv"), sep='\t')
	elif insert:
		df_ins.columns = [str(Path(vcf_name).stem)]
		df_ins.to_csv(Path(name).with_suffix(".insertion.tsv"), sep='\t')


def is_fasta(filename : str) -> bool:
	'''
	Is the file a fasta file?
	Input : path to the file
	Output : True, False (if the file is a pickle file for example) or raise an error if file doesn't exist
	'''
	try:
		fa = pyfastx.Fastx(filename)
		fasta = [name for name,seq,comment in fa]
		return any(fasta)
	except RuntimeError as runerr:
		return False
	except UnicodeDecodeError as uderr:
		return False
	except FileNotFoundError as fnferr:
		print(f'\nFile {filename} doesn\'t exists: {fnferr}\n')
		raise FileNotFoundError(f'\nFile {filename} doesn\'t exists\n')
	except FileExistsError as feerr:
		print(f'\nFile {filename} doesn\'t exists: {feerr}\n')
		raise FileExistsError(f'\nFile {filename} doesn\'t exists\n')
	except : 
		print(f'\nUnexpected error: {sys.exc_info()[0]}\n')
		raise


def is_pickle(filename : str) -> bool:
	'''
	Is the file a pickle file?
	Input : path to the file
	Output : True, False (if the file is a fatsa file for example) or raise an error if file doesn't exist
	'''
	try:
		with open(filename, 'rb') as f:
			genome = pk.load(f)
			return any(genome)
	except EOFError as eoferr:
		return False
	except UnicodeDecodeError as uderr:
		return False
	except FileNotFoundError as fnferr:
		print(f'\nFile {filename} doesn\'t exist: {fnferr}\n')
		raise FileNotFoundError(f'\nFile {filename} doesn\'t exists\n')
	except pk.UnpicklingError:
		return False
	except :
		print(f'\nUnexpected error: {sys.exc_info()[0]}\n')
		raise


def get_genome_dict(genome_file, logger):
	'''
	Check the genome file and (1) read the fatsa file and create the pickle dictionary or (2) load the pickle dictionary and return this dictionary
	Input : path to the genome file and logger
	Output : dictionnary containing the genome
	'''
	print(f'Check input genome file {genome_file}...')
	logger.info(f'Check input genome file {genome_file}')
	if is_fasta(genome_file):
		genome = {}
		path = Path(genome_file).with_suffix('.pk')
		fa = pyfastx.Fastx(genome_file)
		print("chr", "length")
		for name,seq,comment in fa:
			genome[name] = seq
			print(name, len(seq))

		if not path.exists():
			print(f'Create the pickle genome dictionary in {path}')
			logger.info(f'Create the pickle genome dictionary in {path}')
			with open(path, 'wb') as out_pickle:
				pk.dump(genome, out_pickle)
		else:
			print(f'File {path} already exists !')
			logger.warning(f'File {path} already exists !')
	elif is_pickle(genome_file):
		print('Charge the pickle genome dictionary...')
		with open(genome_file, 'rb') as f:
			genome = pk.load(f)
	else: 
		raise ValueError(f'{genome_file} is not a fasta or a pickle file !')

	return genome



def add_snp (snp_count, ref, alt, triplet):
	'''
	Add 1 to the SNP in counter 
	As in VCF the variant is be given for the leading (forward, 5'-3') strand, if we get G>* or A>* the reference, the variant and the triplet are reversed complement 
	Input : SNP counter, the reference allele, the variant and the triplet containing the reference 
	Output : SNP counter
	'''
	if TRANSLATE[str(ref + ">" + alt)] == str(ref + ">" + alt):
		snp_count[str(ref + ">" + alt)][triplet] += 1
	else:
		snp_count[TRANSLATE[str(ref + ">" + alt)]][triplet.translate(TAB)[::-1]] += 1
	return snp_count


def size_indel(ref, alt):
	'''
	Get the indel size
	Input : reference allele and the alternative variant 
	Output : size of the indel variant
	'''
	size = 0	
	if len(ref) > len(alt):
		size = len(ref)-len(alt)
	elif len(alt) > len(ref):
		size = len(alt)-len(ref)
	else:
		print('Not an indel !')
		raise ValueError('Not an indel !')
	return size


def variants_count(vcf_file, vcf_file_pass, genome, logger):
	'''
	Read the vcf file(s) and create the different counters to summary the files
	Input : path to the vcf file from filter module (optional), path to the vcf file containing only pass variants from filter module, genome in a dictionnary and the output path file and the logger
	Output : counters of (1) the deletion size, (2) the insertion size, (3) the snp (in tripet), the impacted genes list and a dictiionnary contaning general stats from the vcf
	'''
	print('Counting variants and get genes...')	
	logger.info(f'Compute stats on variants and genes')

	##################
	#Counter/list creation

	stats = {}

	if  vcf_file:
		stats['Total'] = 0
		stats['germline'] = 0
		stats['PON'] = 0
		stats['non functional'] = 0
		stats['germline+PON'] = 0
		stats['germline+non functional'] = 0
		stats['PON+non functional'] = 0
		stats['germline+PON+non functional'] = 0
	stats['PASS'] = 0
	stats['SNP'] = [0,set()]
	stats['DNP'] = [0,set()]
	stats['TNP'] = [0,set()]
	stats['NP'] = [0,set()]
	stats['INSERTION'] = [0,set()]
	stats['DELETION'] = [0,set()]
	genes_list = {}
	gene_info = [0,[0,0,0,0,0,0],'',[], [], []] # [gene, [variants number (tumor burden),[snp, dnp, tnp, np, insertion, deletion]], chromosome, reference, alternative, position]
	idx = {}
	snp_count = create_snp_dict()
	counter_deletion_size = Counter()
	counter_insertion_size = Counter()

	f = open(vcf_file_pass, 'r', encoding='latin1')
	nb_lines_pass = len(list(dropwhile(lambda x: x[0]=='#', (line for line in f))))
	f.close()
	if vcf_file:
		f = open(vcf_file, 'r', encoding='latin1')
		nb_lines = len(list(dropwhile(lambda x: x[0]=='#', (line for line in f))))
	f.close()

	####################
	# Read pass vcf file	

	print(f'Read {vcf_file_pass}...')
	logger.info(f'Read {vcf_file_pass}')
	with tqdm(total=nb_lines_pass) as pbar:						# progress bar
		vcf_reader_pass = read_vcf(vcf_file_pass)				# create the read file generator
		for line in vcf_reader_pass:
			if type(line) == type({}):					# get index position
				idx = line
			else:

				# Get values
				chr = line[idx['idx_chr']]
				ref = line[idx['idx_ref']]
				alts = line[idx['idx_alts']].split(',')
				nb_variants = len(alts)
				pos = int(line[idx['idx_pos']])			
			
				# Get triplet by using the genome and variant position			
				triplet = str(genome[chr][pos - 2:pos + 1])
				

				infos = [tuple(infos.split('=')) if len(infos.split('='))>1 else (infos,'') for infos in line[idx['idx_info']].split(';')]
				if 'FUNCOTATION' in list(zip(*infos))[0]:
					idx_funcotation = list(zip(*infos))[0].index('FUNCOTATION')
					gene = list(zip(*infos))[1][idx_funcotation].split('|')[0].lstrip('[')
					if not gene in genes_list.keys():
						genes_list[gene]=deepcopy(gene_info)
					if genes_list[gene][2] == '':
						genes_list[gene][2] = chr
					genes_list[gene][3].append(ref)
					genes_list[gene][4].append(alts)
					genes_list[gene][5].append(pos)
				else:
					warnings.warn("Warning! Some variants aren\'t annotated by Funcotation", UserWarning)					
					logger.warning(f'Some variants aren\'t annotated by Funcotation !')
				
				for i, alt in enumerate(alts):

					stats['PASS']+=1
					if 'FUNCOTATION' in list(zip(*infos))[0]:
						genes_list[gene][0]+=1
					if is_snp(len(ref), len(alt)):
						snp_count = add_snp(snp_count, ref, alt, triplet)
						stats['SNP'][0]+=1
						if 'FUNCOTATION' in list(zip(*infos))[0]:
							genes_list[gene][1][0]+=1
							stats['SNP'][1].add(gene)
					elif is_dnp(len(ref), len(alt)):
						stats['DNP'][0]+=1
						if 'FUNCOTATION' in list(zip(*infos))[0]:
							genes_list[gene][1][1]+=1
							stats['DNP'][1].add(gene)
					elif is_tnp(len(ref), len(alt)):
						stats['TNP'][0]+=1
						if 'FUNCOTATION' in list(zip(*infos))[0]:
							genes_list[gene][1][2]+=1
							stats['TNP'][1].add(gene)
					elif is_np(len(ref), len(alt)):
						stats['NP'][0]+=1
						if 'FUNCOTATION' in list(zip(*infos))[0]:
							genes_list[gene][1][3]+=1
							stats['NP'][1].add(gene)
					elif is_insertion(len(ref), len(alt)):
						counter_insertion_size.update([size_indel(ref, alt)])
						stats['INSERTION'][0]+=1
						if 'FUNCOTATION' in list(zip(*infos))[0]:
							genes_list[gene][1][4]+=1
							stats['INSERTION'][1].add(gene)
					elif is_deletion(len(ref), len(alt)):
						counter_deletion_size.update([size_indel(ref, alt)])
						stats['DELETION'][0]+=1
						if 'FUNCOTATION' in list(zip(*infos))[0]:
							genes_list[gene][1][5]+=1
							stats['DELETION'][1].add(gene)	
				pbar.update(1)

	#####################################################
	# If whole vcf (not just pass variants) file is given

	if vcf_file:
		print(f'Read {vcf_file}...')
		logger.info(f'Read {vcf_file}')
		with tqdm(total=nb_lines) as pbar: 					# progress bar
			vcf_reader = read_vcf(vcf_file)					# create the read file generator				
			for line in vcf_reader:
				if type(line) != type({}):
					non_functional = True
					alts = line[idx['idx_alts']].split(',')
					nb_variants = len(alts)
					filters = line[idx['idx_filt']].split(',')
	
					for i, alt in enumerate(alts):
						stats['Total']+=1	
	
						non_functional = any([True for nf in [information.split('=')[1].split(',')[-nb_variants:] for information in line[idx['idx_info']].split(';') if information.split('=')[0] == 'OTHER_FILTER'] if 'NOT_FUNCTIONAL' in nf])
	
						if 'germline' in filters and 'panel_of_normals' in filters and non_functional:	
							stats['germline+PON+non functional']+=1
						elif 'germline' in filters and 'panel_of_normals' in filters:
							stats['germline+PON']+=1
						elif 'germline' in filters and non_functional:
							stats['germline+non functional']+=1
						elif 'panel_of_normals' in filters and non_functional:
							stats['PON+non functional']+=1
						elif 'germline' in filters:
							stats['germline']+=1
						elif 'panel_of_normals' in filters:
							stats['PON']+=1
						elif non_functional:
							stats['non functional']+=1
					pbar.update(1)

	# Get percentage from count for pass SNP
	snp_count_pct = deepcopy(snp_count)
	tot = sum([sum(val.values()) for key, val in snp_count.items()])
	for key, val in snp_count.items():
		for k, v in val.items():
			snp_count_pct[key][k] = round((v / tot), 8)

	return counter_deletion_size, counter_insertion_size, snp_count_pct, genes_list, stats


def is_snp(ref_length, alt_length):
	'''ex: A > T '''
	if ref_length == 1:
		return ref_length == alt_length
	else:
		return False


def is_dnp(ref_length, alt_length):
	'''ex: AT > GA '''
	if ref_length == 2:
		return ref_length == alt_length
	else:
		return False


def is_tnp(ref_length, alt_length):
	'''ex: AAG > TGC '''
	if ref_length == 3:
		return ref_length == alt_length
	else:
		return False


def is_np(ref_length, alt_length):
	'''ex: AATGG > GTCAA '''
	if ref_length > 3:
		return ref_length == alt_length
	else:
		return False


def is_deletion(ref_length, alt_length):
	'''ex: TAAG > T '''
	return ref_length > alt_length


def is_insertion(ref_length, alt_length):
	'''ex: T > TAAG '''
	return ref_length < alt_length


def write_stats(vcf_file : str, vcf_file_pass : str, out_stats : str, stats : Counter, out_genes : str, logger):
	'''
	Writes a simple statistics file on the variants of the VCF file
	Input : path to the vcf file from filter module (optional), path to the vcf file containing only pass variants from filter module, counter containing vcf stats and the output path file 
	Output : Write the vcf stats file
	'''
	print('Write stats file...')
	logger.info('Write stats file')
	with open(out_stats ,'w') as o:
		if vcf_file:
			o.write(f'########################### {vcf_file} ###########################\n')		
			o.write(f'\n###########################\nTotal variants : {stats["Total"]}\n###########################\n')
			o.write(f'germline: {stats["germline"]}\t|\tPON: {stats["PON"]}\t|\tnot functional: {stats["non functional"]}\n')
			o.write(f'germline & PON: {stats["germline+PON"]}\t|\tgermline &  not functional: {stats["germline+non functional"]}\t|\tPON &  not functional: {stats["PON+non functional"]}\n')
			o.write(f'germline & PON & not functional: {stats["germline+PON+non functional"]}\n')
		else:
			o.write(f'########################### {vcf_file_pass} ###########################\n')
		o.write(f'\n###########################\nPASS : {stats["PASS"]}\n###########################\n')
		o.write(f'---\nvariants : {stats["SNP"][0]+stats["DNP"][0]+stats["TNP"][0]+stats["NP"][0]+stats["INSERTION"][0]+stats["DELETION"][0]}\n---\n')
		o.write(f'SNP : {stats["SNP"][0]}\tDNP : {stats["DNP"][0]}\tTNP : {stats["TNP"][0]}\tNP : {stats["NP"][0]}\n')
		o.write(f'INDEL : {stats["INSERTION"][0]+stats["DELETION"][0]}\tINSERTION : {stats["INSERTION"][0]}, DELETION : {stats["DELETION"][0]}\n')
		o.write(f'---\nImpacted genes : { len(stats["SNP"][1]) + len(stats["DNP"][1]) + len(stats["TNP"][1]) + len(stats["NP"][1]) + len(stats["DELETION"][1]) + len(stats["INSERTION"][1]) } (list in {out_genes})\n---\n')
		o.write(f'SNP : {len(stats["SNP"][1])}\tDNP : {len(stats["DNP"][1])}\tTNP : {len(stats["TNP"][1])}\tNP : {len(stats["NP"][1])}\n')
		o.write(f'INDEL : {len(stats["INSERTION"][1])+len(stats["DELETION"][1])} \tINSERTION : {len(stats["INSERTION"][1])}, DELETION : {len(stats["DELETION"][1])}\n')

def write_impacted_genes(out_genes : str, genes_list : pd.DataFrame, logger):
	'''
	Writes an xlsx file containing the list of genes impacted by the variants from the VCF file
	Input : name use for the output files (.xlsx), the genes list and the logger
	Output : Write the genes list in two different format
	'''
	print('Write genes list file... (.xlsx and .tsv)')
	logger.info('Write genes list file (.xlsx and .tsv)')
	genes_list.to_excel(out_genes)
	genes_list.to_csv(Path(out_genes).with_suffix('.tsv'), sep='\t')


def summary(vcf_file : str, vcf_file_pass : str, genome_file : str, out_stats : str, out_genes : str, out_profile : str, out_indel : str, logger, enrichment : bool):
	'''
	Summarizes the vcf file
	Input : path to the vcf file from filter module (optional), path to the vcf file containing only pass variants from filter module, path to genome (fasta or pickle),
	outputs path, logger and boolean for writing gene list enrichments analysis 
	Output : stat summary file, file containing the impacted genes, mutation profile graph, indel size graph 
	'''

	genome = get_genome_dict(genome_file, logger)

	####################
	# Get stats and data

	counter_deletion_size, counter_insertion_size, snp_count_pct, genes_list, stats = variants_count(vcf_file, vcf_file_pass, genome, logger)

	##################
	# Write stats file

	write_stats(vcf_file, vcf_file_pass, out_stats, stats, out_genes, logger)

	##########################################
	# Write genes file in tsv and wlxs formats

	genes_df = create_dataframe_from_gene(genes_list)
	write_impacted_genes(out_genes, genes_df, logger)
	################################################
	# Biological process enrichment using the genes list with the ToppGene and Panther API
	if enrichment:
		toppgene_name = str(Path(out_genes).resolve().parent)+'/'+str(Path(vcf_file_pass).stem.replace('.','_'))+'_ToppGene_enrichment'
		panther_name =  str(Path(out_genes).resolve().parent)+'/'+str(Path(vcf_file_pass).stem.replace('.','_'))+'_Panther_enrichment'
		ToppGene_GEOA(genes_list, toppgene_name, logger)
		Panther_GEOA(genes_list, panther_name, logger)

	#################
	# Create snp mutation types plot and indel size plot

	graph_snp(snp_count_pct, out_profile, vcf_file_pass, logger)
	if not graph_indel(counter_deletion_size, counter_insertion_size, out_indel, vcf_file_pass, logger):
		logger.warning(f'No figures will be plotted !')


def main(args):

	vcf_file = args.vcf
	vcf_file_pass = args.vcf_pass
	working_directory = Path(vcf_file_pass).parent.absolute()
	genome_file = args.genome
	out_stats = args.stats
	out_genes = args.genes
	out_profile = args.profile
	out_indel = args.indel
	enrichment = args.enrichment

	# Logger configuration

	logger = logging.getLogger('LOTUS summarise')
	logger.setLevel(logging.DEBUG)

	fh = logging.FileHandler(args.log)
	fh.setLevel(logging.DEBUG)

	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)

	logger.addHandler(fh)	

	# Verification of given arguments

	try:
		if vcf_file:
			vcf_file_stats = os.stat(vcf_file)
			vcf_file_pass_stats = os.stat(vcf_file_pass)
			if vcf_file_stats.st_size < vcf_file_pass_stats.st_size:
				raise ValueError(f'Classique vcf file cannot be smaller than pass vcf file !')
			
	except ValueError:
		print (f'{vcf_file} cannot be smaller than {vcf_file_pass} !:', sys.exc_info()[0])
		logger.error('- Problem with a vcf file -')
		exit(1)

	try:
		logger.info('Verification of vcf files')
		if vcf_file:
			verif_input_vcf(vcf_file)
		verif_input_vcf(vcf_file_pass)
		logger.info('- Vcf file ok -')
	except ValueError:
		print (f'Problem with a vcf file {vcf_file}:', sys.exc_info()[0])
		logger.error('- Problem with a vcf file -')
		exit(1)

	try:
		logger.info('Verification of outputs file')
		verif_output(out_stats)
		verif_output(out_genes)
		if not out_profile.endswith('.svg'):
                	out_profile = Path(out_profile).with_suffix('.svg')
		verif_output(out_profile)
		if not out_indel.endswith('.svg'):
			out_indel = Path(out_indel).with_suffix('.svg')
		verif_output(out_indel)
		if enrichment:
			verif_output(str(Path(out_genes).resolve().parent)+'/'+str(Path(vcf_file_pass).stem.replace('.','_'))+'_ToppGene_enrichment.xlsx')
			verif_output(str(Path(out_genes).resolve().parent)+'/'+str(Path(vcf_file_pass).stem.replace('.','_'))+'_Panther_enrichment.xlsx')
		logger.info('- Outputs file ok -')
	except ValueError:
		print ('Problem with one or more output files: ', sys.exc_info()[0])
		logger.error('- Problem with output files -')
		exit(1)

	# Start

	logger.info('**************************************************************************************************************')
	logger.info('*** LOTUS summarise module ***')
	no_argument = ''
	if enrichment:
		no_argument+=' --enrichment'
	logger.info(f'** cmd line : python lotus.py summarise -v {str(vcf_file)} -vp {str(vcf_file_pass)} -g {str(genome_file)} -s {str(out_stats)} -genes {str(out_genes)} -p {str(out_profile)} -i {str(out_indel)}'+str(no_argument)+' **')

	logger.info('* Start summarizing *')
	logger.info(f'Working directory (vcf files folder) : {working_directory}')
	logger.info(f'Current directory : {Path().absolute()}')

	summary(vcf_file, vcf_file_pass, genome_file, out_stats, Path(out_genes).with_suffix('.xlsx'), out_profile, out_indel, logger, enrichment)

	logger.info('* End summarizing *')
	logger.info(f'File created :\t{out_stats}\t{out_genes}\t{out_profile}')
	logger.info('**************************************************************************************************************')

	# End







