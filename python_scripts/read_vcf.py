#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def read_vcf(vcf_file, get_variant_idx = False):
	'''
	Create a generator to read the vcf file
	Input : path to the vcf file
	Output : generator returning a dictionnary containing the index position or the strip and split line
	'''
	idx = {}
	with open(vcf_file, 'r', encoding='latin-1') as vcf:
		if not get_variant_idx:
			for line in vcf:
				if line != '':
					if line [0] != '#':
						yield (line.strip().split('\t'))
					elif line.startswith('#') and not line.startswith('##'):
						line = line.lstrip('#').split('\t')
						idx['idx_chr'] = line.index('CHROM')
						idx['idx_pos'] = line.index('POS')
						idx['idx_ref'] = line.index('REF')
						idx['idx_alts'] = line.index('ALT')
						idx['idx_filt'] = line.index('FILTER')
						idx['idx_info'] = line.index('INFO')
						idx['idx_format'] = line.index('FORMAT')
						idx['idx_values'] = idx['idx_format']+1
						yield (idx)
		else:
			pos = vcf.tell()
			line = vcf.readline()
			if line.startswith('#') and not line.startswith('##'):
				line = line.lstrip('#').split('\t')
				idx['idx_chr'] = line.index('CHROM')
				idx['idx_pos'] = line.index('POS')
				idx['idx_ref'] = line.index('REF')
				idx['idx_alts'] = line.index('ALT')
				idx['idx_filt'] = line.index('FILTER')
				idx['idx_info'] = line.index('INFO')
				idx['idx_format'] = line.index('FORMAT')
				idx['idx_values'] = idx['idx_format']+1
				yield(pos, idx)
			while line != '':
				pos = vcf.tell()
				line = vcf.readline()
				if line != '':
					if line [0] != '#':
						yield (pos,line.strip().split('\t'))
					elif line.startswith('#') and not line.startswith('##'):
						line = line.lstrip('#').split('\t')
						idx['idx_chr'] = line.index('CHROM')
						idx['idx_pos'] = line.index('POS')
						idx['idx_ref'] = line.index('REF')
						idx['idx_alts'] = line.index('ALT')
						idx['idx_filt'] = line.index('FILTER')
						idx['idx_info'] = line.index('INFO')
						idx['idx_format'] = line.index('FORMAT')
						idx['idx_values'] = idx['idx_format']+1
						yield (pos, idx)
