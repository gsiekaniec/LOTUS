#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import logging
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from pathlib import Path
from python_scripts.check_files import verif_input_vcf, verif_output, verif_input_config
from python_scripts.enrichment_api import extract_biological_process
from python_scripts.read_vcf import read_vcf
from python_scripts.path_modification import true_stem


def get_variants(file, variants_save):
	'''

	'''
	for line in read_vcf(file):
		if type(line) == type({}):
			intfield = line
		else:
			chr = line[intfield['idx_chr']]
			pos = line[intfield['idx_pos']]
			ref = line[intfield['idx_ref']]
			alts = line[intfield['idx_alts']]
			for alt in alts.split(','):
				variants_save.add((chr,pos,ref,alt))
	return variants_save


def add_to_genes(genes, gene, type):
	'''
	
	'''
	if not gene in genes.keys():
		genes[gene]={}
		genes[gene]['weak']=0
		genes[gene]['strong']=0
	genes[gene][type]+=1


def modify_variants_pass_and_get_genes(file1, file2, variants, weak, strong, logger):
	'''
	
	'''
	outfile = true_stem(file1)+'_compare_to_'+true_stem(file2)+'.pass.vcf'
	genes = {}
	with open(outfile, 'w') as o:
		for line in read_vcf(file1):
			if type(line) == type({}):
				intfield = line
			else:
				infos = [tuple(infos.split('=')) if len(infos.split('='))>1 else (infos,'') for infos in line[intfield['idx_info']].split(';')]
				if 'FUNCOTATION' in list(zip(*infos))[0]:
					idx_funcotation = list(zip(*infos))[0].index('FUNCOTATION')
					gene = list(zip(*infos))[1][idx_funcotation].split('|')[0].lstrip('[')
				variant_type = []
				for i, alt in enumerate(line[intfield['idx_alts']]):
					variant = (line[intfield['idx_chr']],line[intfield['idx_pos']],line[intfield['idx_ref']],alt)
					if variant in weak:
						variant_type.append('weak')
						add_to_genes(genes, gene, 'weak')
					elif variant in strong:
						variant_type.append('strong')
						add_to_genes(genes, gene, 'strong')
					else:
						variant_type.append('common')
				line[intfield['idx_info']] = line[intfield['idx_info']]+';VARIANT_TYPE='+','.join(variant_type)
				o.write("\t".join(line)+'\n')
	return genes


def create_graph_snp(df_snp, df_snp2, outname, logger):
	'''
	
	'''
	#Colors
	white_96 = ['white']*96
	color_6 = ['darkblue', 'blue', 'lightblue', 'darkgreen', 'green', 'lightgreen']
	color_96 = []
	for i in color_6:
		color_96 += [i]*16

	name1 = true_stem(df_snp.columns[2])
	name2 = true_stem(df_snp2.columns[2])
	outname = true_stem(outname)+'.svg'

	logger.info(f'Save profile comparison graph in {outname}')

	bars = [row[1] for index, row in df_snp.iterrows()]
	height = [float(i) for i in df_snp.iloc[:,2]-df_snp2.iloc[:,2]]
	height2 = [float(i) for i in df_snp.iloc[:,2]]
	height3 = [float(-i) for i in df_snp2.iloc[:,2]]
	group = []
	for index, row in df_snp.iterrows():
		if row[0] not in group:
			group.append(row[0])
	x_pos = np.arange(len(bars))

	#Create bars
	
	plt.figure(figsize=(15, 10))
	ax1 = plt.subplot(1,1,1)
	plt.bar(x_pos, height2, color=white_96, edgecolor=color_96, linestyle="--")	
	plt.bar(x_pos, height3, color=white_96, edgecolor=color_96, linestyle="--")
	plt.bar(x_pos, height, color=color_96)
	
	#Create names on the x-axis 

	plt.xticks(x_pos, bars, rotation=90, fontsize=8)
	ax1.set_ylabel(f'Difference between the percentages of each mutation-type')	
	ax2 = ax1.twiny()
	newpos = [8, 24, 40, 56, 72, 88]
	ax2.set_xticks(newpos)
	ax2.set_xticklabels(group)
	ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
	ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
	ax2.spines['bottom'].set_position(('outward', 36))
	ax2.set_xlabel('Mutation types')
	ax2.set_xlim(ax1.get_xlim())
	plt.tick_params(
	axis='x',          # changes apply to the x-axis
	which='both',      # both major and minor ticks are affected
	bottom=False,      # ticks along the bottom edge are off
	top=False,         # ticks along the top edge are off
	labelbottom=True)  # labels along the bottom edge are off
	for xtick, color in zip(ax2.get_xticklabels(), color_6):
		xtick.set_color(color)
	ax2.spines['bottom'].set_visible(False)
	plt.annotate(name1, xy=(0, ax1.get_ylim()[1]-0.08*ax1.get_ylim()[1]))
	plt.annotate(name2, xy=(0, ax1.get_ylim()[0] - 0.05 * ax1.get_ylim()[0]))

	plt.savefig(outname)
	plt.close()	


def graph_snp(snp_profile_files, out_profile, logger):
	'''
	
	'''
	print('Save profile comparison graph...')
	pbar_snp = tqdm(total=len(snp_profile_files)-1)
	for num in range(1,len(snp_profile_files)):
		df1 = pd.read_csv(snp_profile_files[num-1], sep='\t')
		df2 = pd.read_csv(snp_profile_files[num], sep='\t')
		create_graph_snp(df1,df2,true_stem(out_profile), logger)
		pbar_snp.update(1)
	pbar_snp.close()


def create_graph_indel(deletion1, deletion2, insertion1, insertion2, outname, logger):
	'''

	'''
	insert = True
	delet = True
	if (deletion1 == None or deletion2 == None) and (insertion1 == None or insertion2 == None):
		print('Warning ! No indel files ! No graphs will be created !')
		logger.warning(f'No indel files ! No graphs will be created !')
		return None
	elif insertion1 == None or insertion2 == None:
		print('Warning ! No insertion file !')
		logger.warning(f'No insertion file !')
		insert = False
	elif deletion1 == None or deletion2 == None:
		print('Warning ! No deletion file !')
		logger.warning(f'No deletion file !')
		delet = False

	# Get dataframe from tsv file
		
	if insert:
		ins1 = pd.read_csv(insertion1, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_ins1 = sum(ins1)
		ins2 = pd.read_csv(insertion2, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_ins2 = sum(ins2)
		height_ins1 = list([float(i)/sum_ins1 for i in ins1])
		height_ins2 = list([-(float(i)/sum_ins2) for i in ins2])
		bars_ins1 = list(ins1.index)
		bars_ins2 = list(ins2.index)
		name_ins1 = ins1.name
		name_ins2 = ins2.name
	if delet:
		del1 = pd.read_csv(deletion1, sep='\t', header=0, index_col=0).iloc[:, 0]
		del2 = pd.read_csv(deletion2, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_del1 = sum(del1)
		sum_del2 = sum(del2)
		height_del1 = list([float(i)/sum_del1 for i in del1])
		height_del2 = list([-(float(i)/sum_del2) for i in del2])
		bars_del1 = list(del1.index)
		bars_del2 = list(del2.index)
		name_del1 = del1.name
		name_del2 = del2.name	

	if delet and insert:
		maximum = max([max(bars_ins1)+1, max(bars_ins2)+1, max(bars_del1)+1, max(bars_del2)+1])
	elif delet:
		maximum = max([max(bars_ins1)+1, max(bars_ins2)+1])
	elif not delet:
		maximum = max([max(bars_del1)+1, max(bars_del2)+1])
	x1 = [0]+[i+1 for i in range(maximum)]	
	x2 = [0]+[i+1 for i in range(maximum)]	

	width = 0.25

	plt.figure(figsize=(15, 10))
	ax1 = plt.subplot(1,1,1)
	ax2 = ax1.twiny()

	if delet and insert:
		plt.bar([float(i)+(width*0.65) for i in bars_ins1], height_ins1, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+true_stem(name_ins1))
		plt.bar([float(i)-(width*0.65) for i in bars_del1], height_del1, color = 'k', width = width, edgecolor = 'k', label='Deletion_'+true_stem(name_del1))
		plt.bar([float(i)+(width*0.65) for i in bars_ins2], height_ins2, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+true_stem(name_ins2))
		plt.bar([float(i)-(width*0.65) for i in bars_del2], height_del2, color = 'k', width = width, edgecolor = 'k', label='Deletion_'+true_stem(name_del2))
		outname = true_stem(outname)+'.svg'
		x1 = sorted(list(set(bars_del1).union(set(bars_ins1))))
		x2 = sorted(list(set(bars_ins2).union(set(bars_del2))))
	elif not insert:
		plt.bar([float(i) for i in bars_del1], height_del, color = 'k', width = width, edgecolor = 'k', label='Deletion_'+name_del1)
		plt.bar([float(i) for i in bars_del2], height_del2, color = 'k', width = width, edgecolor = 'k', label='Deletion_'+name_del2)
		outname = true_stem(outname)+'_'+true_stem(name_del1)+'.svg'
		x1 = bars_del1
		x2 = bars_del2
	elif not delet:
		plt.bar([float(i) for i in bars_ins1], height_ins1, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+name_ins1)
		plt.bar([float(i) for i in bars_ins2], height_ins2, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+name_ins2)
		outname = true_stem(outname)+'_'+true_stem(name_ins1)+'.svg'
		x1 = bars_ins1
		x2 = bars_ins2

	logger.info(f'Draw indel size barplot in {outname}')

	ax1.set_xticks(x2)
	ax2.set_xticks(x1)
	ax1.set_xlim(0,max([max(x2), max(x1)]))
	ax2.set_xlim(0,max([max(x2), max(x1)]))
	plt.yticks(fontsize=8)

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
		plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])
	
		plot_margin = 0.25
		x0, x1, y0, y1 = plt.axis()
		plt.axis((x0 - plot_margin, x1 + plot_margin*2, y0, y1))

	plt.hlines(0, 0, maximum+1, color='black', linewidth=1)

	ax1.set_title("Indel size")
	plt.ylabel("Percentage")
		
	plt.legend()
	plt.savefig(outname)
	plt.close()


def graph_indel(insertions_count_files, deletions_count_files, out_indel, logger):
	'''

	'''
	if len(insertions_count_files) != len(deletions_count_files):
		logger.error(f'Not the same number of insertions and deletions count files. No graph will be produced !')
		raise ValueError(f'Not the same number of insertions and deletions count files. No graph will be produced !')
		exit(1)
	print('Draw indel size barplot...')
	pbar_indel = tqdm(total=len(insertions_count_files)-1)
	for num in range(1,len(insertions_count_files)):
		create_graph_indel(insertions_count_files[num-1], insertions_count_files[num], deletions_count_files[num-1], deletions_count_files[num], out_indel, logger)
		pbar_indel.update(1)
	pbar_indel.close()

def compare_vcf(vcf_pass_files, vcf_filtered_files, enrichment, logger):
	'''
	
	'''
	print('Save genes lists from comparisons...')
	pbar_compare = tqdm(total=len(vcf_pass_files)-1)
	for num in range(1,len(vcf_pass_files)):
		logger.info(f'Start comparing {true_stem(vcf_pass_files[num-1])} and {true_stem(vcf_pass_files[num])} !')

		variants_pass_save = set()
		variants_pass_save2 = set()
		variants_filtered_save = set()
		variants_filtered_save2 = set()
			
		variants_pass_save = get_variants(vcf_pass_files[num-1], variants_pass_save)		
		variants_pass_save2 = get_variants(vcf_pass_files[num], variants_pass_save2)
		
		variants_filtered_save = get_variants(vcf_filtered_files[num-1], variants_filtered_save)
		variants_filtered_save2 = get_variants(vcf_filtered_files[num], variants_filtered_save2)
		
		croissant1 = variants_pass_save-variants_pass_save2
		croissant2 = variants_pass_save2-variants_pass_save
		strong1 = croissant1-variants_filtered_save2
		strong2 = croissant2-variants_filtered_save
		weak1 = croissant1.intersection(variants_filtered_save2)
		weak2 = croissant2.intersection(variants_filtered_save)

		dfGenes = pd.DataFrame(columns=['gene', '%weakness', 'charge', str(true_stem(vcf_pass_files[num-1])), str(true_stem(vcf_pass_files[num]))])		
		weakness_list = []
		charge_list = []
		tumor1_list = []
		tumor2_list = []		
		genes1 = modify_variants_pass_and_get_genes(vcf_pass_files[num-1], vcf_pass_files[num], croissant1, weak1, strong1, logger)
		genes2 = modify_variants_pass_and_get_genes(vcf_pass_files[num], vcf_pass_files[num-1], croissant2, weak2, strong2, logger)
		genes_list = list(set(genes1.keys()).union(set(genes2.keys())))
		
		for gene in genes_list:
			
			if gene in genes1.keys() and gene in genes2.keys():
				weakness_list.append(((genes1[gene]['weak']+genes2[gene]['weak'])/(genes1[gene]['weak']+genes2[gene]['weak']+genes1[gene]['strong']+genes2[gene]['strong']))*100)
				charge_list.append(genes1[gene]['weak']+genes2[gene]['weak']+genes1[gene]['strong']+genes2[gene]['strong'])
				tumor1_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
				tumor2_list.append(genes2[gene]['weak']+genes2[gene]['strong'])	
			elif gene in genes1.keys():
				weakness_list.append((genes1[gene]['weak']/(genes1[gene]['weak']+genes1[gene]['strong']))*100)
				charge_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
				tumor1_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
				tumor2_list.append(0)
			elif gene in genes2.keys():
				weakness_list.append((genes2[gene]['weak']/(genes2[gene]['weak']+genes2[gene]['strong']))*100)
				charge_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
				tumor1_list.append(0)
				tumor2_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
		dfGenes['gene']=genes_list
		dfGenes['%weakness']=weakness_list
		dfGenes['charge']=charge_list
		dfGenes[str(true_stem(vcf_pass_files[num-1]))]=tumor1_list
		dfGenes[str(true_stem(vcf_pass_files[num]))]=tumor2_list
		
		dfGenes = dfGenes.sort_values(by=['charge', '%weakness', 'gene'], ascending = [False, True, True])		

		dfGenes = dfGenes.set_index('gene')		

		path_gene_list_out = str(true_stem(vcf_pass_files[num-1]))+'_'+str(true_stem(vcf_pass_files[num]))+'_genes_list.tsv'
		
		logger.info(f'Save genes list in {path_gene_list_out}')
		dfGenes.to_csv(path_gene_list_out, sep='\t')
		logger.info(f'Save genes list in {Path(path_gene_list_out).with_suffix(".xlsx")}')
		dfGenes.to_excel(Path(path_gene_list_out).with_suffix('.xlsx'))

		################################################
		# Biological process enrichment using the genes list with the ToppGene and Panther API
		if enrichment:
			toppgene_name = str(true_stem(vcf_pass_files[num-1]))+'_compare_to_'+str(true_stem(vcf_pass_files[num]))+'_ToppGene_enrichment'
			panther_name =  str(true_stem(vcf_pass_files[num-1]))+'_compare_to_'+str(true_stem(vcf_pass_files[num]))+'_Panther_enrichment'
			extract_biological_process(genes_list, toppgene_name, panther_name, logger)
		pbar_compare.update(1)
	pbar_compare.close()

def main(args):

	out_stats = args.stats
	out_profile = args.profile
	out_indel = args.indel
	enrichment = args.enrichment
	out_gene = args.out

	# Logger configuration

	logger = logging.getLogger('g-LOTUS compare')
	logger.setLevel(logging.DEBUG)

	fh = logging.FileHandler(args.log)
	fh.setLevel(logging.DEBUG)

	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	
	logger.addHandler(fh)

	# Verification of config file

	try:
		logger.info('Verification of input in config file')
		vcf_filtered_files, vcf_pass_files, snp_profile_files, insertions_count_files, deletions_count_files = verif_input_config(args.config)
		logger.info('- Inputs file ok -')
	except ValueError:
		print ('Problem with config file: ', sys.exc_info()[0])
		logger.error('- Problem with config file -')
		exit(1)

	# Verification of given arguments

	try:
		logger.info('Verification of outputs file')
		for i in range(1,len(vcf_filtered_files)):
			file_one = true_stem(vcf_filtered_files[i-1])
			file_two = true_stem(vcf_filtered_files[i])
				
			out_stats = Path(file_one+'_'+file_two+'_'+out_stats).with_suffix('.txt')
			#verif_output(out_stats)
		
			out_profile = Path(true_stem(out_profile)+'_'+file_one+'_'+file_two).with_suffix('.svg')
			#verif_output(out_profile)

			out_indel = Path(true_stem(out_indel)+'_'+file_one+'_'+file_two).with_suffix('.svg')
			#verif_output(out_indel)	
		
			out_gene = file_one+'_'+file_two+'_genes_list.xlsx'
			#verif_output(out_gene)
			out_gene = file_one+'_'+file_two+'_genes_list.tsv'
			#verif_output(out_gene)

			vcf_strong_weak1 = Path(file_one+'_compare_to_'+file_two).with_suffix('.pass.vcf')
			#verif_output(vcf_strong_weak1)
			vcf_strong_weak2 = Path(file_two+'_compare_to_'+file_one).with_suffix('.pass.vcf')
			#verif_output(vcf_strong_weak2)

			if enrichment:
				verif_output(Path(file_one+'_compare_to_'+file_two+'_ToppGene_enrichment.xlsx'))
				verif_output(Path(file_one+'_compare_to_'+file_two+'_ToppGene_enrichment.tsv'))
				verif_output(Path(file_one+'_compare_to_'+file_two+'_Panther_enrichment.xlsx'))
				verif_output(Path(file_one+'_compare_to_'+file_two+'_Panther_enrichment.tsv'))

		logger.info('- Outputs file ok -')

	except ValueError:
		print ('Problem with one or more output files: ', sys.exc_info()[0])
		logger.error('- Problem with output files -')
		exit(1)

	# Start

	logger.info('********************************************************************************************************************')
	logger.info('*** g-LOTUS compare module ***')
	logger.info('* Start comparing *')
	logger.info(f'Current directory : {Path().absolute()}')
	
	##############################
	# Comparison of both vcf files

	compare_vcf(vcf_pass_files, vcf_filtered_files, enrichment, logger)

	####################################################
	# Create snp mutation types plot and indel size plot

	graph_snp(snp_profile_files, out_profile, logger)
	graph_indel(insertions_count_files, deletions_count_files, out_indel, logger)

	logger.info('* End comparing *')
	logger.info('********************************************************************************************************************')

	# End








