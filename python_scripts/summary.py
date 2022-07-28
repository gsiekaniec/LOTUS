#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import pickle as pk
from pathlib import Path
import vcf
from collections import OrderedDict, Counter
from copy import deepcopy
import warnings
import pandas as pd
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig


TRANSLATE = {'A>C': 'T>G', 'A>G': 'T>C', 'A>T': 'T>A', 'C>A': 'C>A', 'C>G': 'C>G', 'C>T': 'C>T',
             'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A', 'T>A': 'T>A', 'T>C': 'T>C', 'T>G': 'T>G'}
TAB = str.maketrans("ACGT", "TGCA")
	

def create_dataframe_from_gene_list(d):
	col = ['nb variant(s)','nb details (snp,dnp,tnp,np,insertion,deletion)','chromosome','ref','alt variant(s)','position(s)']
	id = 'gene_name'
	d = OrderedDict(sorted(d.items()))
	df = pd.DataFrame.from_dict(d, orient='index', columns=col)
	df.index.name = id
	
	df2 = pd.DataFrame()
	df2['nb details (snp,dnp,tnp,np,insertion,deletion)'] = [','.join(map(str, l)) for l in df['nb details (snp,dnp,tnp,np,insertion,deletion)']]
	df['nb details (snp,dnp,tnp,np,insertion,deletion)'] = df2['nb details (snp,dnp,tnp,np,insertion,deletion)'].values
	del df2
	df2 = pd.DataFrame()
	df2['ref'] = [','.join(map(str, l)) for l in df['ref']]
	df['ref'] = df2['ref'].values
	del df2
	df2 = pd.DataFrame()
	df2['alt variant(s)'] = [ ','.join([str(l2) if len(l2) > 1 else str(l2[0]) for l2 in l]) for l in df['alt variant(s)']]
	df['alt variant(s)'] = df2['alt variant(s)'].values
	del df2
	df2 = pd.DataFrame()
	df2['position(s)'] = [','.join(map(str, l)) for l in df['position(s)']]
	df['position(s)'] = df2['position(s)'].values
	del df2
	return df


def is_fasta(filename : str) -> bool:
	try:
		with open(filename, "r") as handle:
			fasta = SeqIO.parse(handle, "fasta")
			return any(fasta)
	except UnicodeDecodeError as uderr:
		return False
	except FileNotFoundError as fnferr:
		print(f'\nFile {filename} doesn\'t exist: {fnferr}\n')
		raise
	except : 
		print(f'\nUnexpected error: {sys.exc_info()[0]}\n')
		raise


def is_pickle(filename : str) -> bool:
	try:
		with open(filename, 'rb') as f:
			genome = pk.load(f)
			return any(genome)
	except UnicodeDecodeError as uderr:
		return False
	except FileNotFoundError as fnferr:
		print(f'\nFile {filename} doesn\'t exist: {fnferr}\n')
		raise
	except :
		print(f'\nUnexpected error: {sys.exc_info()[0]}\n')
		raise


def get_genome_dict(genome_file):
	print(f'Check input genome file {genome_file}...')
	if is_fasta(genome_file):
		print('Ok, create the pickle genome dictionary...')
		genome = {}
		with open(genome_file, 'r') as handle:
			with open(Path(genome_file).with_suffix('.pk'), 'wb') as out_pickle:
				print("chr", "length")
				for record in SeqIO.parse(handle, "fasta"):
					genome[record.id]=record.seq
					print(record.id,len(record.seq))
				pk.dump(genome, out_pickle)
	else:
		if is_pickle(genome_file):
			print('Charge the pickle genome dictionary...')
			with open(genome_file, 'rb') as f:
				genome = pk.load(f)
	return genome


def check_vcf_file(vcf_file):
	print(f'Check input vcf file {vcf_file}...')
	vcf_reader = vcf.Reader(open(vcf_file, 'r', encoding='latin-1'))
	meta = vcf_reader.metadata
	info = vcf_reader.infos
	filter = vcf_reader.filters
	format = vcf_reader.formats
	samples = vcf_reader.samples
	if meta == OrderedDict() or info == OrderedDict() or filter == OrderedDict() or format == OrderedDict() or samples == OrderedDict():
		raise IOError(f'Invalid vcf file {vcf_file}!')


def variants_count(vcf_file, genome, pass_only):

	print('Counting variants...')	

	stats = Counter()

	if not pass_only:
		stats['Total'] = 0
		stats['germline'] = 0
		stats['PON'] = 0
		stats['germline+PON'] = 0
		stats['functional'] = 0
		stats['non functional'] = 0
	stats['PASS'] = 0
	stats['SNP'] = 0
	stats['DNP'] = 0
	stats['TNP'] = 0
	stats['NP'] = 0
	stats['INSERTION'] = 0
	stats['DELETION'] = 0
	genes_list = {}
	gene_info = [0,[0,0,0,0,0,0],'',[], [], []]

	vcf_reader = vcf.Reader(open(vcf_file, 'r', encoding='latin-1'))
	#sample = vcf_reader.samples[0] #Only one (the first) sample from vcf are processed
	for record in vcf_reader:
		chr = record.CHROM
		pos = record.POS
		ref = record.REF
		alts = record.ALT
		if pass_only:
			if 'FUNCOTATION' in record.INFO.keys():
				gene = record.INFO['FUNCOTATION'][0].split('|')[0].lstrip('[')
				if not gene in genes_list.keys():
					genes_list[gene]=deepcopy(gene_info)
				if genes_list[gene][2] == '':
					genes_list[gene][2] = chr
				genes_list[gene][3].append(ref)
				genes_list[gene][4].append(alts)
				genes_list[gene][5].append(pos)
			for i, alt in enumerate(alts):
				stats.update(['PASS'])
				if 'FUNCOTATION' in record.INFO.keys():
					genes_list[gene][0]+=1
				if is_snp(len(ref), len(alt)):
					stats.update(['SNP'])
					if 'FUNCOTATION' in record.INFO.keys():
						genes_list[gene][1][0]+=1
				elif is_dnp(len(ref), len(alt)):
					stats.update(['DNP'])
					if 'FUNCOTATION' in record.INFO.keys():
						genes_list[gene][1][1]+=1
				elif is_tnp(len(ref), len(alt)):
					stats.update(['TNP'])
					if 'FUNCOTATION' in record.INFO.keys():
						genes_list[gene][1][2]+=1
				elif is_np(len(ref), len(alt)):
					stats.update(['NP'])
					if 'FUNCOTATION' in record.INFO.keys():
						genes_list[gene][1][3]+=1
				elif is_insertion(len(ref), len(alt)):
					stats.update(['INSERTION'])
					genes_list[gene][1][4]+=1
				elif is_deletion(len(ref), len(alt)):
					stats.update(['DELETION'])
					if 'FUNCOTATION' in record.INFO.keys():
						genes_list[gene][1][5]+=1
	s = 0
	for k, v in genes_list.items():
		s+=v[0]
	if s != stats['PASS']:
		warnings.warn("Warning! Some variants aren\'t annotated by Funcotation", UserWarning)	
	return genes_list, stats


def is_snp(ref_length, alt_length):
        if ref_length == 1:
                return ref_length == alt_length
        else:
                return False


def is_dnp(ref_length, alt_length):
        if ref_length == 2:
                return ref_length == alt_length
        else:
                return False


def is_tnp(ref_length, alt_length):
        if ref_length == 3:
                return ref_length == alt_length
        else:
                return False


def is_np(ref_length, alt_length):
        if ref_length > 3:
                return ref_length == alt_length
        else:
                return False


def is_deletion(ref_length, alt_length):
        return ref_length > alt_length


def is_insertion(ref_length, alt_length):
        return ref_length < alt_length


def write_stats(vcf_file : str, out_stats : str, stats : Counter):
	'''
	Writes a simple statistics file on the variants of the VCF file
	'''
	print('Write stats file...')
	with open(out_stats ,'w') as o:
		o.write(f'########### {vcf_file} ##########\nPASS : {stats["PASS"]}\n----------\n')
		o.write(f'SNP : {stats["SNP"]}\nDNP : {stats["DNP"]}\nTNP : {stats["TNP"]}\nNP : {stats["NP"]}\n----------\n')
		o.write(f'INDEL : {stats["INSERTION"]+stats["DELETION"]} \tINSERTION : {stats["INSERTION"]}, DELETION : {stats["DELETION"]}\n')


def write_impacted_genes(out_genes : str, genes_list : pd.DataFrame):
	'''
	Writes an xlsx file containing the list of genes impacted by the variants from the VCF file
	'''
	print('Write genes list file... (.xlsx)')
	genes_list.to_excel(out_genes)

def summary(vcf_file : str , genome_file : str, out_stats : str, out_genes : str, out_profil : str, pass_only : bool):
	'''
	'''

	check_vcf_file(vcf_file)
	#genome = get_genome_dict(genome_file)
	genome = []
	genes_list, stats = variants_count(vcf_file, genome, pass_only)
	genes_list = create_dataframe_from_gene_list(genes_list)
	
	write_stats(vcf_file, out_stats, stats)
	write_impacted_genes(out_genes, genes_list)

	#test installation SigProfiler (to move in compare.py)
	#genInstall.install('GRCh38', offline_files_path='./vcf/')
	#matrices = matGen.SigProfilerMatrixGeneratorFunc("lotus_res", "GRCh38", "../testSig", plot=True)
	#print(matrices)
	#sig.sigProfilerExtractor("matrix", "Sig_output", matrices, reference_genome="GRCh38", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100)


def main(args):

	vcf_file = args.vcf
	genome_file = args.genome
	out_stats = args.stats
	out_genes = args.genes
	out_profil = args.profil
	pass_only = args.pass_only

	summary(vcf_file , genome_file, out_stats, Path(out_genes).with_suffix('.xlsx'), out_profil, pass_only)









