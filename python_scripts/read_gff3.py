#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle as pk
import gzip
from pathlib import Path

def read_gff3(file):
	'''
	Read a gff3 file and return a dictionary containing chromosome and positions for each gene
	input : gff3 file path
	output (return + pickle file) : {('gene name', 'ensembl id') : ('chr', (start_pos, end_pos)), ...}
	'''
	gz = False
	if file.endswith('gz') or file.endswith('gzip'):
		f = gzip.open(file, "rb")
		gz = True
	else:
		f = open(file, 'r')
	gene_id_dico = {}
	gene_name_dico = {}
	transcript_dico = {}

	first_line = True
	for line in f:
		if gz:
			line = str(line, 'utf-8')
		if first_line:
			if not line.startswith('##gff-version'):
				raise ValueError(f'{file} is not a gff3 file !')
			first_line = False
		if line[0] != '#':                              # Don't care about comments lines
			line = line.strip().split('\t')

			if 'gene' in line[2]:
				if not 'chr' in line[0]:
					chr = 'chr' + str(line[0])
				else:
					chr = str(line[0])
				strand = str(line[6])
				pos_start = int(line[3])
				pos_end = int(line[4])
				dico_infos = {}
				for i in line[8].split(';'):
					dico_infos[i.split('=')[0]] = i.split('=')[1]
				ensemble_id = dico_infos['ID'].split(':')[-1]
				if 'Name' in dico_infos.keys():
					name = dico_infos['Name']
				else:
					name = ''
				if not ensemble_id in gene_id_dico.keys():
					gene_id_dico[ensemble_id] = []
				gene_id_dico[ensemble_id].append((name, chr, (pos_start, pos_end)))
				if name != '':
					if not name in gene_name_dico.keys():
						gene_name_dico[name] = []
					gene_name_dico[name].append((ensemble_id, chr, (pos_start, pos_end)))
				
			# Transcript

			if 'Parent=gene' in line[8]:
				for i in line[8].split(';'):
					if 'Parent=gene:' in i:
						g = i.split('Parent=gene:')[1]
					elif 'ID=transcript:' in i:
						t = i.split('ID=transcript:')[1]
				transcript_dico[t] = gene_id_dico[g]

	f.close()
	list_gff3 = [gene_name_dico, gene_id_dico, transcript_dico]
	pickle_file = Path(file).with_suffix('.pk')
	with open(pickle_file, 'wb') as o:
		pk.dump(list_gff3, o)
	return list_gff3



