#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import vcf
from copy import deepcopy
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts, _Filter as VcfFilter
from pathlib import Path

def line_without_fail(record, to_suppr : str):
	'''
	Takes a record and deletes all variants that do not pass the filters, if no variants pass the filters returns False otherwise returns True
	Input : record and list of variants to suppressed
	Output : True (Record modified) or False (No record left)
	'''

	nb_alt = len(record.ALT)
	record.ALT=[alt for i, alt in enumerate(record.ALT) if i not in to_suppr]

	if record.ALT == []:
		return False

	for id, info in record.INFO.items():
		#print(id, info, type(info))
		if id == 'OTHER_FILTER':
			record.INFO['OTHER_FILTER']='PASS'
		elif id == 'AS_FilterStatus':
			record.INFO['AS_FilterStatus']=['|'.join([status for i, status in enumerate(info[0].split('|')) if i not in to_suppr])]
		elif id == 'AS_SB_TABLE':
			record.INFO['AS_SB_TABLE']= (','.join(info[:2]).split('|')[0]+'|'+'|'.join([fr for i, fr in enumerate(','.join(info).split('|')[1:]) if i not in to_suppr])).split(',')
		elif type(info) != type(1) and type(info) != type(1.1) and type(info) != type("1") and type(info) != type(True):
			if len(info) == nb_alt:
				record.INFO[id] = [elmt for i, elmt in enumerate(info) if i not in to_suppr]
			elif len(info) == nb_alt+1:
				record.INFO[id] = [info[0]]+[elmt for i, elmt in enumerate(info[1:]) if i not in to_suppr]
	
	names = record.FORMAT.split(':')
	f_vals = [record.samples[0].data[vx] for vx in range(len(names))]
	handy_dict = dict(zip(names, f_vals))
	for id, data in handy_dict.items():
		if id == 'GT':
			handy_dict[id] = '/'.join([str(i) for i in range(nb_alt-len(to_suppr)+1)])
		elif type(data) == type([]):
			if len(data) == nb_alt:
				handy_dict[id] = [d for i, d in enumerate(data) if i not in to_suppr]
			elif len(data) == nb_alt+1:
				handy_dict[id] = [data[0]]+[d for i, d in enumerate(data[1:]) if i not in to_suppr]
	new_vals = [handy_dict[x] for x in names]
	record.samples[0].data = record.samples[0].data._make(new_vals)

	return True


def fail_filters(info : {}, AD : [], AF : [], nb_alt : int):
	'''
	Do variant(s) at a position passes the filters ?
	Input : info, AD and AF fields and the number of variants
	Output : False (don't fail) or a dictionnary containing the failed filters
	'''
	
	mutation_to_save = set({'MISSENSE', 'NONSENSE', 'NONSTOP', 'RNA', 'LINCRNA', 'START_CODON_SNP', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'FRAME_SHIFT_INS', 'FRAME_SHIFT_DEL', 'START_CODON_INS', 'START_CODON_DEL', 'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME'})
	splice_site = 'SPLICE_SITE'
	fail = {'all' : set()}
	#print(info, AD, AF)
	for i in range(nb_alt):
		fail[i]=set()
	if type(AF) != type([]):
		AF = [AF]
	
	# Minimum allele sequencing depth
	for i, ad in enumerate(AD[1:]):
		if int(ad) < 5:
			fail[i].add('AD')
	# Minimum allele frequency
	for i, af in enumerate(AF):
		if float(af) < 0.1:
			fail[i].add('AF')
	
	for id, elmt in info.items():
		
		# global sequencing depth
		if id == 'DP':
			if type(elmt) == type([]):
				if int(elmt[0]) < 10:
					fail['all'].add('DP')
			else:
				if int(elmt) < 10:
					fail['all'].add('DP')
		# median fragment length by allele
		if id == 'MFRL':
			for i, mfrl in enumerate(elmt[1:]):
				if float(mfrl) == 0:
					fail[i].add('MFRL')
		# median base quality by allele
		if id == 'MBQ':
			for i, mbq in enumerate(elmt[1:]):
				if float(mbq) < 20:
					fail[i].add('MBQ')
		# population allele frequencies
		if id == 'POPAF':
			for i, popaf in enumerate(elmt):
				popaf = 10**(-float(popaf))
				if float(popaf) >= 0.00005:
					fail[i].add('POPAF')
		# only allele with a potential functional effect (not silent)
		if id == 'FUNCOTATION':
			for i, funco in enumerate(elmt):
				t = funco.split('|')[5]
				t2 = funco.split('|')[6]
				if not (t in mutation_to_save or (t == 'SPLICE_SITE' and t2 in mutation_to_save)):
					fail[i].add('NOT_FUNCTIONAL')			
	
	if all([True if v == set() else False for v in fail.values()]) :
		return(False)
	else :
		return(fail)


def filter(vcf_file : str, output : str):
	'''
	Filters out variants:
	Input : vcf file 
	Output : vcf file containing variant with new filter column and vcf without variants that don't pass filters
	'''

	print(f'{vcf_file}')
	output2 = Path(output).with_suffix('.pass.vcf')
	print(f'Read file : {vcf_file}\nWrite in {output} and {output2}')

	vcf_reader = vcf.Reader(open(vcf_file, 'r', encoding='latin-1')) #read in latin-1 because some char from funcotator anotation are not in utf-8
	vcf_reader.infos['OTHER_FILTER'] = VcfInfo('OTHER_FILTER', vcf_field_counts['A'], 'String', 'Other filter that dont pass.', None, None) #New info field
	vcf_reader.infos['AS_SB_TABLE'] = VcfInfo('AS_SB_TABLE', vcf_field_counts['R'], 'String', vcf_reader.infos['AS_SB_TABLE'][4], None, None) #Modify the info.AS_SB_TABLE field
	vcf_reader.filters['LOTUS_filter'] = VcfFilter('LOTUS_filter', 'Mutation does not pass g-LOTUS filters') #Add the g-LOTUS filter to other filters

	sample = vcf_reader.samples[0] #Only one (the first) sample from vcf are processed. g-LOTUS compute only single sample vcf.

	#Output writer
	vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader)
	vcf_writer_pass = vcf.Writer(open(output2, 'w'), vcf_reader)


	#Reading vcf
	for record in vcf_reader:
		info = record.INFO
		
		call = record.genotype(sample)
		AD = call['AD']
		AF = call['AF']
		old_filter = record.FILTER
		fail = fail_filters(info, AD, AF, len(record.ALT)) #False if the variant does not pass the filters and the name of the filter that does not pass otherwise 
		print(fail)
		if fail: # suppress variant that don't pass filters
			to_suppr = set()
			filters_failed = ''
			if 'all' in fail.keys() and fail['all'] != set(): 
				filters_failed += 'DP,'
			
			for i in range(len(fail)-1):
				if fail[i] != set():
					filters_failed += str(i)+':'+':'.join([f for f in fail[i]])+','
					to_suppr.add(i)
				else:
					filters_failed += str(i)+':PASS'+','
			filters_failed = filters_failed.rstrip(',')
			record.add_info('OTHER_FILTER', str(filters_failed))
			record.add_filter('LOTUS_filter')
			vcf_writer.write_record(record)
			if not old_filter and not 'DP' in filters_failed:
				if line_without_fail(record, to_suppr): #Modification of the current record to keep only variants that pass filters
					vcf_writer_pass.write_record(record)
		else:
			record.add_info('OTHER_FILTER', 'PASS')
			vcf_writer.write_record(record)
			if not old_filter :
				vcf_writer_pass.write_record(record)
					
			
def main(args):
	
	vcf_file = args.vcf    
	output = args.out
	
	filter(vcf_file , output)


























