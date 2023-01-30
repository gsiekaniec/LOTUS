#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import logging
import sys
import os
import pickle as pk
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from pathlib import Path
from python_scripts.check_files import verif_input_vcf, verif_output, verif_input_config, verif_input, verif_supplementary_information_file
from python_scripts.toppgene_api import ToppGene_GEOA
from python_scripts.panther_api import Panther_GEOA
from python_scripts.read_vcf import read_vcf
from python_scripts.path_modification import true_stem
from python_scripts.read_gff3 import read_gff3


def get_informations_for_genes(info_file, logger):
	df = pd.read_excel(info_file, index_col=1)
	df = df.drop(['Ordre'], axis=1)
	df.set_axis([source.split(' Info')[0] for source in df.columns], axis="columns", inplace=True)
	print(f'Extract informations from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	logger.info(f'Extract informations from {len(list(df.columns))} sources: {", ".join(list(df.columns))}')
	#print('DF',df)
	return (df)


def get_variants(file, variants_save):
	'''
	Take a vcf file and return all its variants in a set (or a dictionnary with the index fields)
	Input : vcf file and empty set to save variants
	Output : set containing the variants
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


def add_to_genes(genes : dict, gene : str, type : str, chr : str, gene_position : (), mg : str, mc : str, mp : str):
	'''
	Add 1 to a gene in the genes dictionary according to its type
        Input : genes dictionary, gene name, type of the gene ('weak' or 'strong') + add chromosome, gene position, mutation in genomic sequence, mutation in coding sequence, mutation in proteic sequence
	'''
	if not gene in genes.keys():
		genes[gene]={}
		genes[gene]['weak']=0
		genes[gene]['strong']=0
		genes[gene]['mg']=[]
		genes[gene]['mc']=[]
		genes[gene]['mp']=[]
	genes[gene][type]+=1
	genes[gene]['chr']=chr
	genes[gene]['gene_pos']=gene_position
	genes[gene]['mg'].append(mg)
	genes[gene]['mc'].append(mc)
	genes[gene]['mp'].append(mp)


def modify_variants_pass_and_get_genes(file1, file2, variants, weak, strong, gene_name_dico : {}, gene_id_dico : {}, transcript_dico : {}, parent_name : str, logger):
	'''
	Get the pass variants lines from a vcf and add the variant type ('weak'/'strong'/'common') + get the count of types weak and strong for every gene 
	Input : vcf file 1 to modify, vcf file 2 (just used for the name), set of specific variants from file 1, set of weak variants from file 1, set of strong variants from file 1 and the logger
	Output : dictionary containing genes and their number of weak and strong variants (specific to vcf file 1)
	'''

	outfile = parent_name+'/'+str(true_stem(file1)+'_'+true_stem(file2)+'.passed.vcf')
	genes = {}
	with open(outfile, 'w') as o:
		for line in read_vcf(file1):
			if type(line) == type({}):
				intfield = line
			else:
				chromosome = line[intfield['idx_chr']]
				ID_ensemble = []
				mg = []
				mc = []
				mp = []
				infos = [tuple(infos.split('=')) if len(infos.split('='))>1 else (infos,'') for infos in line[intfield['idx_info']].split(';')]
				nb_alt = len(line[intfield['idx_alts']].split(','))
				if 'FUNCOTATION' in list(zip(*infos))[0]:
					idx_funcotation = list(zip(*infos))[0].index('FUNCOTATION')
					info_funcotation = list(zip(*infos))[1][idx_funcotation].split(',')
					gene = info_funcotation[0].split('|')[0].lstrip('[')
					id_ensembl = info_funcotation[0].split('|')[12]
					try:
						if len(gene_name_dico[gene]) > 1:
							gene_position = transcript_dico[id_ensembl.split('.')[0]][0]
						else:
							gene_position = gene_name_dico[gene][0]
					except KeyError:
						try:
							gene_position = transcript_dico[id_ensembl.split('.')[0]][0]
						except KeyError:
							gene_position = ('',chromosome,'')
					for i in range(nb_alt):
						mg.append(info_funcotation[i].split('|')[11])
						mc.append(info_funcotation[i].split('|')[16])
						mp.append(info_funcotation[i].split('|')[18])
						ID_ensemble.append(info_funcotation[i].split('|')[12])
				variant_type = []
				for i, alt in enumerate(line[intfield['idx_alts']].split(',')):							#For each variants
					variant = (line[intfield['idx_chr']], line[intfield['idx_pos']], line[intfield['idx_ref']], alt)
					# Add variant type to the info field
					if variant in weak:
						variant_type.append('weak')
						add_to_genes(genes, gene, 'weak', chromosome, gene_position, mg[i], mc[i], mp[i])
					elif variant in strong:
						variant_type.append('strong')
						add_to_genes(genes, gene, 'strong', chromosome, gene_position, mg[i], mc[i], mp[i])
					else:
						variant_type.append('common')
				# Recreate and write the vcf line
				line[intfield['idx_info']] = line[intfield['idx_info']]+';VARIANT_TYPE='+','.join(variant_type)		#Add the variant type to the vcf Info field
				o.write("\t".join(line)+'\n')
	return genes


def create_graph_snp(df_snp, df_snp2, outname, logger):
	'''
	Snp count plot creation
        Input : snp profile of file 1, snp profile of file 2, output name for the plot and the logger
	'''

	#Colors
	white_96 = ['white']*96
	color_6 = ['darkblue', 'blue', 'lightblue', 'darkgreen', 'green', 'lightgreen']
	color_96 = []
	for i in color_6:
		color_96 += [i]*16

	name1 = true_stem(df_snp.columns[2])
	name2 = true_stem(df_snp2.columns[2])
	outname = str(Path(outname).resolve().parent)+'/'+str(true_stem(outname)+'.svg')

	logger.info(f'Save profile comparison graph in {outname}')

	bars = [row[1] for index, row in df_snp.iterrows()]
	height = [float(i) for i in df_snp.iloc[:,2]-df_snp2.iloc[:,2]]
	height2 = [float(i) for i in df_snp.iloc[:,2]]
	height3 = [float(-i) for i in df_snp2.iloc[:,2]]
	group = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
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

	outname = Path(outname).with_suffix('.png')
	plt.savefig(outname)

	plt.close()	


def graph_snp(snp_profile_files, out_profile, logger):
	'''
	Create the snp profile plot for each snp profile file (if exist)
	Input : snp profile files, snp profile files, output name for the plot and the logger	
	'''
	print('Save profile comparison graph...')
	pbar_snp = tqdm(total=len(snp_profile_files)-1)
	for num in range(1,len(snp_profile_files)):
		df1 = pd.read_csv(snp_profile_files[num-1], sep='\t')
		df2 = pd.read_csv(snp_profile_files[num], sep='\t')
		create_graph_snp(df1,df2,out_profile, logger)
		pbar_snp.update(1)
	pbar_snp.close()


def create_graph_indel(deletion1, deletion2, insertion1, insertion2, outname, logger):
	'''
	Create the indel count comparison plot
        Input : Insertion count file 1, Insertion count file 2, Deletion count file 1, Deletion count file 2, output name for the plot and the logger
	'''
	
	insert = True
	delet = True
	if (deletion1 == None or deletion2 == None) and (insertion1 == None or insertion2 == None):  # no insertion nor deltion counts file
		print('Warning ! No indel files ! No graphs will be created !')
		logger.warning(f'No indel files ! No graphs will be created !')
		return None
	elif insertion1 == None or insertion2 == None:  # insertion counts file does not exist
		print('Warning ! No insertion file !')
		logger.warning(f'No insertion file !')
		insert = False
	elif deletion1 == None or deletion2 == None:  # deltion counts file does not exist
		print('Warning ! No deletion file !')
		logger.warning(f'No deletion file !')
		delet = False

	# Get dataframe from tsv file

	if insert:
		ins1 = pd.read_csv(insertion1, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_ins1 = sum(ins1)
		ins2 = pd.read_csv(insertion2, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_ins2 = sum(ins2)
		height_ins1 = list([float(i) / sum_ins1 for i in ins1])
		height_ins2 = list([-(float(i) / sum_ins2) for i in ins2])
		bars_ins1 = list(ins1.index)
		bars_ins2 = list(ins2.index)
		name_ins1 = ins1.name
		name_ins2 = ins2.name
	if delet:
		del1 = pd.read_csv(deletion1, sep='\t', header=0, index_col=0).iloc[:, 0]
		del2 = pd.read_csv(deletion2, sep='\t', header=0, index_col=0).iloc[:, 0]
		sum_del1 = sum(del1)
		sum_del2 = sum(del2)
		height_del1 = list([float(i) / sum_del1 for i in del1])
		height_del2 = list([-(float(i) / sum_del2) for i in del2])
		bars_del1 = list(del1.index)
		bars_del2 = list(del2.index)
		name_del1 = del1.name
		name_del2 = del2.name

	if delet and insert:
		maximum = max([max(bars_ins1)+1, max(bars_ins2)+1, max(bars_del1)+1, max(bars_del2)+1])
	elif delet:
		maximum = max([max(bars_del1)+1, max(bars_del2)+1])
	elif not delet:
		maximum = max([max(bars_ins1)+1, max(bars_ins2)+1])
	
	x1 = [0]+[i+1 for i in range(maximum)]
	x2 = [0]+[i+1 for i in range(maximum)]

	width = 0.25

	fig = plt.figure(figsize=(15, 10))
	ax1 = plt.subplot(1,1,1)

	# create plot according to deletion and insertion counts

	if delet and insert:
		plt.bar([float(i)+(width*0.65) for i in bars_ins1], height_ins1, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+true_stem(name_ins1))
		plt.bar([float(i)-(width*0.65) for i in bars_del1], height_del1, color = 'darkred', width = width, edgecolor = 'darkred', label='Deletion_'+true_stem(name_del1))
		plt.bar([float(i)+(width*0.65) for i in bars_ins2], height_ins2, color = 'b', width = width, edgecolor = 'b', label='Insertion_'+true_stem(name_ins2))
		plt.bar([float(i)-(width*0.65) for i in bars_del2], height_del2, color = 'darkblue', width = width, edgecolor = 'darkblue', label='Deletion_'+true_stem(name_del2))
		x1 = sorted(list(set(bars_del1).union(set(bars_ins1))))
		x2 = sorted(list(set(bars_ins2).union(set(bars_del2))))
	elif not insert:
		plt.bar([float(i) for i in bars_del1], height_del1, color = 'darkred', width = width, edgecolor = 'darkred', label='Deletion_'+true_stem(name_del1))
		plt.bar([float(i) for i in bars_del2], height_del2, color = 'darkblue', width = width, edgecolor = 'darkblue', label='Deletion_'+true_stem(name_del2))
		x1 = bars_del1
		x2 = bars_del2
	elif not delet:
		plt.bar([float(i) for i in bars_ins1], height_ins1, color = 'r', width = width, edgecolor = 'r', label='Insertion_'+true_stem(name_ins1))
		plt.bar([float(i) for i in bars_ins2], height_ins2, color = 'b', width = width, edgecolor = 'b', label='Insertion_'+true_stem(name_ins2))
		x1 = bars_ins1
		x2 = bars_ins2

	plt.legend()
	ax1.grid(axis='x', color='lightgrey', linestyle='-', linewidth=0.5)
	max_x = max(x1+x2)

	ax1.set_xlim(0, max_x+1)
	ax1.set_xlabel("Indel size (bp)")
	ax1.set_ylabel("Indel percentage")
	ax1.set_xticks(x2)
	ax1.set_xticklabels(x2, fontsize=10)

	ax2 = ax1.twiny()
	ax2.set_xlim(0, max_x+1)
	ax2.set_xticks(x1)
	ax2.set_xticklabels(x1, fontsize=10)
	ax2.grid(axis='x', color='lightgrey', linestyle='-', linewidth=0.5)

	fig.set_size_inches(max_x/3, 10)

	y_pos = list(np.arange(-1, 1.1, 0.1))
	y_value = [round(i,1) for i in (list(np.arange(1, 0, -0.1))+list(np.arange(0,1.1, 0.1)))]
	ax1.set_yticks(y_pos)
	ax1.set_yticklabels(y_value)

	plt.hlines(0, 0, maximum+1, color='black', linewidth=1)

	logger.info(f'Draw indel size barplot in {outname}')

	plt.savefig(outname)

	outname = Path(outname).with_suffix('.png')
	plt.savefig(outname)

	plt.close()


def graph_indel(insertions_count_files, deletions_count_files, out_indel, logger):
	'''
	Create the indel count comparison plot for each indel count file (if exist) 
	Input : Insertion count files, Deletion count files, output name for the plot and the logger 
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


def compare_vcf(vcf_pass_files : list, vcf_filtered_files : list, gene_name_dico : {}, gene_id_dico : {}, transcript_dico : {}, out_gene : str, infos : str, enrichment : bool, logger):
	'''
	Compare each vcf files with the n-1 vcf file (starting from the second file)
	Input : vcf files (pass and filtered from the filter module), boolean to know if GOEA need to be done and the logger
	'''

	# Get genes infos if not None
	if infos:
		infos_df = get_informations_for_genes(infos, logger)

	print('Save genes lists from comparisons...')
	pbar_compare = tqdm(total=len(vcf_pass_files)-1)
	for num in range(1,len(vcf_pass_files)):					# Comparison file by file starting from the second with is n-1

		logger.info(f'Start comparing {true_stem(vcf_pass_files[num-1])} and {true_stem(vcf_pass_files[num])} !')

		# Get variants for all vcf files

		variants_pass_save = set()
		variants_pass_save2 = set()
		variants_filtered_save = set()
		variants_filtered_save2 = set()
		variants_pass_save = get_variants(vcf_pass_files[num-1], variants_pass_save)		
		variants_pass_save2 = get_variants(vcf_pass_files[num], variants_pass_save2)
		variants_filtered_save = get_variants(vcf_filtered_files[num-1], variants_filtered_save)
		variants_filtered_save2 = get_variants(vcf_filtered_files[num], variants_filtered_save2)
		
		# Get variants sets specific to each file and separate this sets in weak and strong variants using the filtered file from the other vcf

		croissant1 = variants_pass_save-variants_pass_save2				# variants from pass 1 without pass 2
		croissant2 = variants_pass_save2-variants_pass_save				# variants from pass 2 without pass 1
		strong1 = croissant1-variants_filtered_save2					# variants from pass 1 without pass 2 and filter 2
		strong2 = croissant2-variants_filtered_save					# variants from pass 2 without pass 1 and filter 1
		weak1 = croissant1.intersection(variants_filtered_save2)			# variants from pass 1 without pass 2 but present in filter 2
		weak2 = croissant2.intersection(variants_filtered_save)				# variants from pass 2 without pass 1 but present in filter 1

		# Creation of the dataframe containing impacted gene and creation of the new vcf files containing the genes type ('strong'/'weak'/'common')  

		dfGenes = pd.DataFrame(columns=['Gene symbol', 'Chromosome', 'Gene position start', 'Gene position end', 'Tumour burden (symmetrical difference)', 'Gene weakness (in %)', str(true_stem(vcf_pass_files[num-1])), 'g.'+str(true_stem(vcf_pass_files[num-1])), 'c.'+str(true_stem(vcf_pass_files[num-1])), 'p.'+str(true_stem(vcf_pass_files[num-1])), str(true_stem(vcf_pass_files[num])), 'g.'+str(true_stem(vcf_pass_files[num])), 'c.'+str(true_stem(vcf_pass_files[num])), 'p.'+str(true_stem(vcf_pass_files[num]))])
		chromosomes_list = []
		gene_position_list = []
		weakness_list = []
		charge_list = []
		tumor1_list = []
		tumor2_list = []
		genomic_variant_annotation_list = []
		coding_variant_annotation_list = []
		proteomic_variant_annotation_list = []
		genomic_variant2_annotation_list = []
		coding_variant2_annotation_list = []
		proteomic_variant2_annotation_list = []
		genes1 = modify_variants_pass_and_get_genes(vcf_pass_files[num-1], vcf_pass_files[num], croissant1, weak1, strong1, gene_name_dico, gene_id_dico, transcript_dico, str(Path(out_gene).resolve().parent), logger)
		genes2 = modify_variants_pass_and_get_genes(vcf_pass_files[num], vcf_pass_files[num-1], croissant2, weak2, strong2, gene_name_dico, gene_id_dico, transcript_dico, str(Path(out_gene).resolve().parent), logger)
		genes_list = list(set(genes1.keys()).union(set(genes2.keys())))
		
		for gene in genes_list:
			
			if gene in genes1.keys() and gene in genes2.keys():
				chromosomes_list.append(genes1[gene]['chr'])
				gene_position_list.append(genes1[gene]['gene_pos'])
				weakness_list.append(round(((genes1[gene]['weak']+genes2[gene]['weak'])/(genes1[gene]['weak']+genes2[gene]['weak']+genes1[gene]['strong']+genes2[gene]['strong']))*100,2))
				charge_list.append(genes1[gene]['weak']+genes2[gene]['weak']+genes1[gene]['strong']+genes2[gene]['strong'])
				tumor1_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
				tumor2_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
				genomic_variant_annotation_list.append('|'.join(genes1[gene]['mg']))
				coding_variant_annotation_list.append('|'.join(genes1[gene]['mc']))
				proteomic_variant_annotation_list.append('|'.join(genes1[gene]['mp']))
				genomic_variant2_annotation_list.append('|'.join(genes2[gene]['mg']))
				coding_variant2_annotation_list.append('|'.join(genes2[gene]['mc']))
				proteomic_variant2_annotation_list.append('|'.join(genes2[gene]['mp']))
			elif gene in genes1.keys():
				chromosomes_list.append(genes1[gene]['chr'])
				gene_position_list.append(genes1[gene]['gene_pos'])
				weakness_list.append(round((genes1[gene]['weak']/(genes1[gene]['weak']+genes1[gene]['strong']))*100,2))
				charge_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
				tumor1_list.append(genes1[gene]['weak']+genes1[gene]['strong'])
				tumor2_list.append(0)
				genomic_variant_annotation_list.append('|'.join(genes1[gene]['mg']))
				coding_variant_annotation_list.append('|'.join(genes1[gene]['mc']))
				proteomic_variant_annotation_list.append('|'.join(genes1[gene]['mp']))
				genomic_variant2_annotation_list.append('')
				coding_variant2_annotation_list.append('')
				proteomic_variant2_annotation_list.append('')
			elif gene in genes2.keys():
				chromosomes_list.append(genes2[gene]['chr'])
				gene_position_list.append(genes2[gene]['gene_pos'])
				weakness_list.append(round((genes2[gene]['weak']/(genes2[gene]['weak']+genes2[gene]['strong']))*100,2))
				charge_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
				tumor1_list.append(0)
				genomic_variant_annotation_list.append('')
				coding_variant_annotation_list.append('')
				proteomic_variant_annotation_list.append('')
				tumor2_list.append(genes2[gene]['weak']+genes2[gene]['strong'])
				genomic_variant2_annotation_list.append('|'.join(genes2[gene]['mg']))
				coding_variant2_annotation_list.append('|'.join(genes2[gene]['mc']))
				proteomic_variant2_annotation_list.append('|'.join(genes2[gene]['mp']))
		dfGenes['Gene symbol']=genes_list
		dfGenes['Chromosome']=chromosomes_list
		dfGenes['Gene position start']=[i[2][0] for i in gene_position_list]
		dfGenes['Gene position end']=[i[2][1] for i in gene_position_list]
		dfGenes['Tumour burden (symmetrical difference)']=charge_list
		dfGenes['Gene weakness (in %)']=weakness_list
		b1 = str(true_stem(vcf_pass_files[num-1]))
		b2 = str(true_stem(vcf_pass_files[num]))
		dfGenes[b1]=tumor1_list
		dfGenes['g.'+b1]=genomic_variant_annotation_list
		dfGenes['c.'+b1]=coding_variant_annotation_list
		dfGenes['p.'+b1]=proteomic_variant_annotation_list
		dfGenes[b2]=tumor2_list
		dfGenes['g.'+b2]=genomic_variant2_annotation_list
		dfGenes['c.'+b2]=coding_variant2_annotation_list
		dfGenes['p.'+b2]=proteomic_variant2_annotation_list

		dfGenes = dfGenes.sort_values(by=['Tumour burden (symmetrical difference)', 'Gene weakness (in %)', 'Gene symbol'], ascending = [False, True, True])		
		dfGenes = dfGenes.set_index('Gene symbol')		
		
		# Save impacted genes informations

		if infos:
			dfGenes = dfGenes.join(infos_df)
			# Drop empty informational columns
			empty_cols = [col for col in dfGenes.columns if dfGenes[col].isnull().all()]
			dfGenes.drop(empty_cols, axis=1, inplace=True)

		logger.info(f'Save genes list in {out_gene}')
		print(f'Save genes list in {out_gene}')
		dfGenes.to_csv(out_gene, sep='\t')
		logger.info(f'Save genes list in {Path(out_gene).with_suffix(".xlsx")}')
		print(f'Save genes list in {Path(out_gene).with_suffix(".xlsx")}')
		dfGenes.to_excel(Path(out_gene).with_suffix('.xlsx'))

		################################################
		# Biological process enrichment using the genes list with the ToppGene and Panther API
		if enrichment:
			toppgene_name = str(Path(out_gene).resolve().parent)+'/'+str(true_stem(vcf_pass_files[num-1]))+'_'+str(true_stem(vcf_pass_files[num]))+'_ToppGene_enrichment'
			panther_name =  str(Path(out_gene).resolve().parent)+'/'+str(true_stem(vcf_pass_files[num-1]))+'_'+str(true_stem(vcf_pass_files[num]))+'_Panther_enrichment'
			ToppGene_GEOA(genes_list, toppgene_name, logger)
			Panther_GEOA(genes_list, panther_name, logger)
		pbar_compare.update(1)
	pbar_compare.close()


def main(args):

	out_profile = args.profile
	out_indel = args.indel
	enrichment = args.enrichment
	out_gene = args.out
	gff3 = args.gff3
	pk_gff3 = args.pk_gff3
	sup_infos = args.agi

	current_directory = os.getcwd()

	# Logger configuration

	logger = logging.getLogger('LOTUS compare')
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
		print ('Problem with config file : ', sys.exc_info()[1])
		logger.error('- Problem with config file -')
		exit(1)

	try:
		verif_input(gff3)
		if not pk_gff3:
			genes_positions = read_gff3(gff3)
		else:
			with open(gff3, 'rb') as f:
				genes_positions = pk.load(f)
		transcript_dico = genes_positions[2]
		gene_id_dico = genes_positions[1]
		gene_name_dico = genes_positions[0]
	except pk.UnpicklingError:
		print (f'{gff3} is not a pickle file : ', sys.exc_info()[1])
		logger.error(f'- {gff3} is not a pickle file ! -')
		exit(1)	
	except UnicodeDecodeError:
		print (f'{gff3} can not be read as a classic gff3 file : ', sys.exc_info()[1])
		logger.error(f'- {gff3} can not be read as a classic gff3 file ! -')
		exit(1)
	except ValueError:
		print (f'{gff3} is not a gff3 file (or the first # line missing) : ', sys.exc_info()[1])
		logger.error(f'- {gff3} is not a gff3 file ! -')
		exit(1)

	# Verification of given arguments

	try:
		logger.info('Verification of outputs file')
		for i in range(1,len(vcf_filtered_files)):
			file_one = true_stem(vcf_filtered_files[i-1])
			file_two = true_stem(vcf_filtered_files[i])

			out_profile = str(Path(out_profile).resolve().parent)+'/'+str(Path(true_stem(out_profile)+'_'+file_one+'_'+file_two).with_suffix('.svg'))
			verif_output(out_profile)

			out_indel = str(Path(out_indel).resolve().parent)+'/'+str(Path(true_stem(out_indel)+'_'+file_one+'_'+file_two).with_suffix('.svg'))
			verif_output(out_indel)	

			out_gene_test = str(Path(out_gene).resolve().parent)+'/'+str(file_one+'_'+file_two+'_'+true_stem(out_gene)+'.MutatedGenes.xlsx')
			verif_output(out_gene_test)
			out_gene = str(Path(out_gene).resolve().parent)+'/'+str(file_one+'_'+file_two+'_'+true_stem(out_gene)+'.MutatedGenes.tsv')
			verif_output(out_gene)

			vcf_strong_weak1 = str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two).with_suffix('.passed.vcf'))
			verif_output(vcf_strong_weak1)
			vcf_strong_weak2 = str(Path(out_gene).resolve().parent)+'/'+str(Path(file_two+'_'+file_one).with_suffix('.passed.vcf'))
			verif_output(vcf_strong_weak2)

			if enrichment:
				verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_ToppGene_enrichment.xlsx')))
				verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_ToppGene_enrichment.tsv')))
				verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_Panther_enrichment.xlsx')))
				verif_output(str(Path(out_gene).resolve().parent)+'/'+str(Path(file_one+'_'+file_two+'_Panther_enrichment.tsv')))

		logger.info('- Outputs file ok -')

	except ValueError:
		print ('Problem with one or more output files: ', sys.exc_info()[1])
		logger.error('- Problem with output files -')
		exit(1)


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

	logger.info('********************************************************************************************************************')
	logger.info('*** LOTUS compare module ***')
	no_argument = ''
	if pk_gff3:
		no_argument+=' --pickle_gff3'
	if enrichment:
		no_argument+=' --enrichment'
	if sup_infos:
		no_argument+=' --additional_gene_information'
	logger.info(f'** cmd line : python lotus.py compare -c {str(args.config)} -o {str(out_gene)} -p {str(out_profile)} -i {str(out_indel)}'+str(no_argument)+' **')
	logger.info('* Start comparing *')
	logger.info(f'Current directory : {Path().absolute()}')
	
	##############################
	# Comparison of both vcf files

	compare_vcf(vcf_pass_files, vcf_filtered_files, gene_name_dico, gene_id_dico, transcript_dico, out_gene, infos, enrichment, logger)

	####################################################
	# Create snp mutation types plot and indel size plot

	graph_snp(snp_profile_files, out_profile, logger)
	graph_indel(insertions_count_files, deletions_count_files, out_indel, logger)

	logger.info('* End comparing *')
	logger.info('********************************************************************************************************************')

	# End








