#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import time
import logging
import sys
from copy import deepcopy
from itertools import takewhile
from tqdm import tqdm
from python_scripts.check_files import verif_input_vcf, verif_output
from python_scripts.path_modification import true_stem

def line_without_fail(record, to_suppr : str, intIndex : dict, strLotusFilterCode : str, strFilterPassKey : str, strFilterPassOk : str):
	'''
	Takes a record and deletes all variants that do not pass the filters, if no variants pass the filters returns False otherwise returns True
	Input : record and list of variants to suppressed
	Output : True (Record modified) or False (No record left)
	'''

	# Get the number of alternatif variants
	nb_alt = len(record[intIndex['Alt']].split(','))

	# Suppression of variants that don't pass LOTUS filters
	alternative=[alt for i, alt in enumerate(record[intIndex['Alt']].split(',')) if i not in to_suppr]
	if alternative == []:
		return False
	else:
		record[intIndex['Alt']] = ','.join(alternative)

	# Filter modification

	record[intIndex['Filter']] = 'PASS'

	# Suppresion of the useless modification

	informations = [[info.split('=')[0], info.split('=')[1]] if (len(info.split('=')) == 2) else tuple([info.split('=')[0],'']) for info in [info for info in record[intIndex['Info']].split(';')]]

	for j, information in enumerate(informations):
		id, info = information
		if id == 'OTHER_FILTER':
			informations[j][1]=strFilterPassOk
		
		elif id == 'AS_FilterStatus':
			informations[j][1]='|'.join([status for i, status in enumerate(info.split('|')) if i not in to_suppr])

		elif id == 'AS_SB_TABLE':
			informations[j][1]='|'.join([fr for i, fr in enumerate(info.split('|')) if i-1 not in to_suppr])
	
		else:
			if len(info) > 1 and not info.isnumeric():
				if len(info.split(',')) == nb_alt:
					informations[j][1] = ','.join([elmt for i, elmt in enumerate(info.split(',')) if i not in to_suppr])
				elif len(info.split(',')) == nb_alt+1:
					informations[j][1] = ','.join([elmt for i, elmt in enumerate(info.split(',')) if i-1 not in to_suppr])
	record[intIndex['Info']] = ';'.join(['='.join([id, info]) for id, info in informations])

	names = record[intIndex['Format']].split(':')
	
	values = [value for value in record[intIndex['Values']].split(':')]
	for i, v in enumerate(values):
		if names[i] == 'GT':
			values[i] = '/'.join([str(i) for i in range(nb_alt-len(to_suppr)+1)]) 
		elif len(v.split(',')) == nb_alt:
			values[i] = '/'.join( [char for i, char in enumerate(v.split(',')) if i not in to_suppr] )
		elif len(v.split(',')) == nb_alt+1:
			values[i] = '/'.join( [char for i, char in enumerate(v.split(',')) if i-1 not in to_suppr] )
	record[intIndex['Values']] = ':'.join(values)

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


def tumor_sample_header_logging(x, logger):
	'''
	Is the tumor sample ID exists in the header
        Input : a record (x) and logger 
	Output : False if tumor sample ID if does not exist or True otherwise
	'''
	if x.startswith('##tumor_sample='):
		strSampleCode = x.split('=')[1].rstrip()
		logger.info(f'Sample id : {strSampleCode}')
		return True
	return False


def get_values (record : str, strFieldSplit : dict, intIndex : dict, InfoFieldsToProcess : set):
	'''
        Extract useful informations too filter a line of the VCF
        Input : the record and dictionnary needed to find and split line fields
        Output : information needed to filter the line (AD, AF, Info fields such as Funcotation...)
        '''

	# Extract data from the record to feed further functions
	# 1. the Info field -> from string or list to Dictionary
	
	recInfoOrig = record[intIndex['Info']].split(strFieldSplit['Info1'])
	recInfo = {}
	for r in recInfoOrig:
		# To avoid empty values:
		if strFieldSplit['Info2'] in r:
			info = r.split(strFieldSplit['Info2'])
			# Special case in Info field that needs treatment
			if info[0] in InfoFieldsToProcess:
				info[1] = info[1].split(strFieldSplit['Funcota'])
				if (info[0] == 'FUNCOTATION'):
					for k, v in enumerate(info[1]):
						info[1][k] = v.strip("[]")
			recInfo[info[0]]=info[1]
		else:
			# Key without a value (eg: PON)
			recInfo[r]=''

	# 2. the AD and AF fields -> into List of float

	recFormat = record[intIndex['Format']].split(strFieldSplit['Format1'])
	recValues = record[intIndex['Values']].split(strFieldSplit['Format1'])
	recAD = recValues[recFormat.index('AD')].split(strFieldSplit['Format2'])
	recAF = recValues[recFormat.index('AF')].split(strFieldSplit['Format2'])


	# 3. get the filter field in case record dont pass LOTUS filter

	recOriginalFilter = False if record[intIndex['Filter']] == 'PASS' else record[intIndex['Filter']]

	return (recInfo, recFormat, recValues, recAD, recAF, recOriginalFilter)


def create_new_records(record, strFieldSplit, intIndex, InfoFieldsToProcess, strCRLF):
	'''
	Creation of new records containing new g-Lotus filter or without variants that dont pass filters  
        Input : record and dictionnary needed to find and split line fields
        Output : boolean (is the record pass filters) and the two new records
	'''

	recInfo, recFormat, recValues, recAD, recAF, recOriginalFilter = get_values (record, strFieldSplit, intIndex, InfoFieldsToProcess)

	# #####################
	# Treatment of failure
	fail = fail_filters(recInfo, recAD, recAF, len(record[intIndex['Alt']].split(','))) #False if the variant does not pass the filters and the name of the filter that does not pass otherwise

	# ####################
	# Rebuilt and storage of filtered data

	# Variable
	strLotusFilterCode = 'LOTUS_filter'
	strFilterPassKey = 'OTHER_FILTER='
	strFilterPassOk = 'PASS'
	blnPass = False                                                 # set as default : the filter was not successful

	# Treatment
	if fail:                                                        # suppress variant that don't pass filters
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

		# Add a new item (new filters) to Info field
		record[intIndex['Info']] = record[intIndex['Info']]+str(';')+strFilterPassKey+str(filters_failed)

		if filters_failed != '':
			if not recOriginalFilter:
				# Add the 'Lotus' filter item in Filter field
				record[intIndex['Filter']] = strLotusFilterCode
			else:	
				record[intIndex['Filter']] = record[intIndex['Filter']]+';'+strLotusFilterCode

		if not recOriginalFilter and not 'DP' in filters_failed:
			cleanRecord = deepcopy(record)
			if line_without_fail(cleanRecord, to_suppr, intIndex, strLotusFilterCode, strFilterPassKey, strFilterPassOk): # Modification of the current record to keep only variants that pass filter
				blnPass = True
		else:
			cleanRecord = record
	else:
		#Add a new item to Info field
		record[intIndex['Info']] = record[intIndex['Info']]+str(';')+strFilterPassKey+strFilterPassOk
		cleanRecord = record
		if not recOriginalFilter :
			blnPass = True

	# Rebuilt the new full record
	newRecord = strFieldSplit[1].join(record)+strCRLF
	newRecord2 = strFieldSplit[1].join(cleanRecord)+strCRLF
	return blnPass, newRecord, newRecord2


def filter(vcf_file : str, logger : str, output : str, working_method : str):
	'''
	Filters out variants:
	Input : vcf file
	Output : 2 vcf files : one containing variant with new filter column and another one without variants that don't pass filters
	logfile : the full path with file name for the log file to write in
	working_method : 2 possibilities : 'InMemory' (more speed but higher memory consumption) or 'Direct' (slow speed but low memory consumption)
	'''

	output = Path(true_stem(output)).with_suffix('.filtered.vcf')
	output2 = Path(true_stem(output)).with_suffix('.pass.vcf')
	print(f'Read file :\t{vcf_file}\nWrite in {output} and {output2}\n')
	

	# log data
	logger.info(f'Working method chosen by user : {str(working_method)}')
	logger.info(f'Read input file : {str(vcf_file)}')
	logger.info(f'Write to output files : {str(output)} and {str(output2)}')

	# ######################################################################################
	# Opening the file, and read it line by line, store it in a list for further treatments

	try:
		with open(vcf_file, mode='r', encoding='latin-1', errors='backslashreplace') as obFileInput:		 #read in latin-1 because some char from funcotator anotation are not in utf-8

			n = 0				# total number of lines in the VCF
			n_header = 0			# number of header lines from the VCF

			if working_method == 'InMemory':
				lstInput1=[]	# The original data
				lstInput1 = obFileInput.readlines()
				n = len(lstInput1)
				n_header = len(list(takewhile(lambda x: x[0]=='#', lstInput1)))
	
			elif working_method == 'Direct':
				# read the open file in r mode for counting lines, then close and reopen it
				n = sum(1 for _ in obFileInput)
				obFileInput.seek(0,0)
				n_header = len(list(takewhile(lambda x: x[0]=='#', (line for line in obFileInput))))
				obFileInput.seek(0,0)

			print(f'Number of lines to read : {n}\nNumber of header lines : {n_header}')
			logger.info(f'Input file contains : {str(n)} lines ({n_header} header and {n-n_header} data)')
			

			# ####################################################
			# Opening for writing outputs
			obFileOuput1 = open(output,'w', encoding='latin-1')		# the modified data -> with all variants and a new complementary filter
			obFileOuput2 = open(output2,'w', encoding='latin-1')		# the modified data -> without variants that don't pass filters


			# Criteria for header (only 6 (intlenghtCrit) char long for searching frame)
			intlenghtCrit = 6
			strHeaderFilter = '##FILT'				# To locate filter fields
			strHeaderInfo = '##INFO'				# To locate info fields

			# New lines to add in header
			strNewFilter = '##FILTER=<ID=LOTUS_filter,Description="Mutation does not pass g-LOTUS filters">'
			strNewInfo = '##INFO=<ID=OTHER_FILTER,Number=A,Type=String,Description="Other filter that dont pass.">'
			
			# Line ending
			strCRLF = '\n'

			# Counters
			intCountHeaderIns = 0					# Inserted header lines counter
			intCountHeaderMod = 0					# Modified header lines counter
			count = 0						# Total lines counter

			id_column = ''						# Saved the id of colunm for later

			# Presetting variable for reading
			sample_id = []
			if working_method=='InMemory':
				strLinePrevious = lstInput1[0]
			elif working_method == 'Direct':
				strLinePrevious = '#'

			# searching through the header

			if working_method=='InMemory':
				print('Header treatment')
				id_column = lstInput1[n_header-1].lstrip('#')
				lstSave = deepcopy(lstInput1)
				for i in tqdm(range(n_header)):
					x = lstSave[i]
					# 1. Search for a specific information without affecting data integrity
					sample_id.append(tumor_sample_header_logging(x, logger))
					# 2. Search for the end of a theme : to insert a line at the end of the theme
					strLinePreviousStart = strLinePrevious[0:intlenghtCrit]			
					if x[0:intlenghtCrit] != strLinePreviousStart:
						if strLinePreviousStart == strHeaderFilter:
							lstInput1 = lstInput1[0:i+intCountHeaderIns]+[strNewFilter+strCRLF]+lstInput1[i+intCountHeaderIns:]
							intCountHeaderIns +=1
						elif strLinePreviousStart == strHeaderInfo:
							lstInput1 = lstInput1[0:i+intCountHeaderIns]+[strNewInfo+strCRLF]+lstInput1[i+intCountHeaderIns:]
							intCountHeaderIns +=1
					strLinePrevious = x
				del lstSave

			elif working_method=='Direct':
				x = obFileInput.readline().strip()
				print ('Header treatment')
				pbar = tqdm(total=n_header)
				while x[0] == '#':
					# 1 Write unmodified lines	
					if strLinePrevious != '#':
						obFileOuput1.write(str(strLinePrevious+strCRLF))                # Simple copy of the actual value
						obFileOuput2.write(str(strLinePrevious+strCRLF))
					# 2. Search for a specific information without affecting data integrity
					sample_id.append(tumor_sample_header_logging(x, logger))
					# 3. Search for the end of a theme : to insert a line at the end of the theme
					strLinePreviousStart = strLinePrevious[0:intlenghtCrit]
					if x[0:intlenghtCrit] != strLinePreviousStart:
						if strLinePreviousStart == strHeaderFilter:
							obFileOuput1.write(str(strNewFilter+strCRLF))
							obFileOuput2.write(str(strNewFilter+strCRLF))
							intCountHeaderIns +=1
						elif strLinePreviousStart == strHeaderInfo :						
							obFileOuput1.write(str(strNewInfo+strCRLF))
							obFileOuput2.write(str(strNewInfo+strCRLF))
							intCountHeaderIns +=1
						if x[0] == '#' and x[1] != '#':
							intDataStartPos=obFileInput.tell()
							id_column = x.lstrip('#')
							obFileOuput1.write(str(x+strCRLF))
							obFileOuput2.write(str(x+strCRLF))						
					strLinePrevious = x
					count += 1
					pbar.update(1)
					x = obFileInput.readline().strip()
				pbar.close()
			if not any(sample_id):
				logger.warning('No sample id !')
			logger.info(f'{intCountHeaderIns} line(s) added to header')
			logger.info(f'Actual number of header line is {n_header+intCountHeaderIns} ({n_header} before)')
			n_header = n_header+intCountHeaderIns
			
					
			# Criteria for data
			# To know how to split each sub-fields of datas
			strFieldSplit = {1 : '\t', 'Info1' : ';', 'Info2' : '=', 'Funcota' : ',', 'Format1' : ':', 'Format2' : ','}

			# index for fields
			intIndex = {'Chr' : id_column.split(strFieldSplit[1]).index('CHROM'), 'Alt' : id_column.split(strFieldSplit[1]).index('ALT'), 'Filter' : id_column.split(strFieldSplit[1]).index('FILTER'), 'Info' : id_column.split(strFieldSplit[1]).index('INFO'), 'Format' : id_column.split(strFieldSplit[1]).index('FORMAT')}
			intIndex['Values'] = intIndex['Format']+1

			InfoFieldsToProcess = {'DP', 'MFRL', 'MBQ', 'POPAF', 'FUNCOTATION'}

			# Counter for new files
			intCountData = 0
			intCountPass = 0

			# ###################################################
                        # Analyze and extract genotype results, line by line

			if working_method=='InMemory':
				lstInput2 = deepcopy(lstInput1)
				for id_line, x in enumerate(lstInput1[n_header:]):
					record = x.strip().split(strFieldSplit[1])

					blnPass, newRecord, newRecord2 = create_new_records(record, strFieldSplit, intIndex, InfoFieldsToProcess, strCRLF)					

					lstInput1[n_header+id_line]=str(newRecord)
					intCountData +=1
					if (blnPass):
						lstInput2[n_header+id_line]=str(newRecord2)
						blnPass = False		# reset to default
						intCountPass +=1
					else:
						lstInput2[n_header+id_line]=None

			elif working_method=='Direct':
				obFileInput.seek(intDataStartPos)
				x = obFileInput.readline().strip()
				count = 1
				print(f'Saving {str(output)} and {str(output2)} on disk', flush=True)
				pbar = tqdm(total=n-(n_header-intCountHeaderIns))
				while x != '':
					record = x.split(strFieldSplit[1])
				
					blnPass, newRecord, newRecord2 = create_new_records(record, strFieldSplit, intIndex, InfoFieldsToProcess, strCRLF)
				
					obFileOuput1.write(newRecord)
					intCountData +=1
					if (blnPass):
						obFileOuput2.write(newRecord2)
						blnPass = False         # reset to default
						intCountPass +=1
					
					x = obFileInput.readline().strip()
					count += 1
					pbar.update(1)
				pbar.close()
			
			if working_method == 'InMemory':				
				print(f'Saving {str(output)} and {str(output2)} on disk', flush=True)
				for count, record in enumerate(tqdm(lstInput1)):
					obFileOuput1.write(record)
					if lstInput2[count]:
						obFileOuput2.write(lstInput2[count])
					count += 1

			print(f'Number of lines containing variants : {str(intCountData)}')
			print(f'Number of lines flagged as PASS : {str(intCountPass)}')
			
			logger.info(f'Number of lines for input header :\t{str(n_header)}')
			logger.info(f'Number of lines for output header :\t{str(n_header)} of which {str(n_header-intCountHeaderIns)} lines added and {str(intCountHeaderIns)} line modified')
			logger.info(f'Number of variants flagged as PASS :\t{str(intCountPass)}')
			

			# ##############################################
			# Close the 2 destination files
			obFileOuput1.close()
			obFileOuput2.close()
			del obFileOuput1
			del obFileOuput2
			if working_method == 'InMemory':
				del lstInput1
				del lstInput2


			logger.info(f' -> file {str(output)} :\t{str(n_header+intCountData)} total lines of which {str(intCountData)} genotype data lines')
			logger.info(f' -> file {str(output2)} :\t{str(n_header+intCountPass)} total lines of which {str(intCountPass)} genotype data lines')
			logger.info(f'* End of filtering *')
			logger.info('**************************************************************************************************************')



	except IOError:
		print ('Error raised : could not read file: '+str(vcf_file)+'. Please check your parameters or filename.', sep='')

	except:
		print ('Unexpected error:', sys.exc_info()[0])
		raise 


def main(args):

	vcf_file = args.vcf
	output = args.out

	working_method = args.working_method	# working_method = 'InMemory' (more speed but higher memory consumption) or 'Direct' (slow speed but low memory consumption)
	working_directory = Path(vcf_file).parent.absolute()

	# Logger configuration
	
	logger = logging.getLogger('g-LOTUS filter')
	logger.setLevel(logging.DEBUG)

	fh = logging.FileHandler(args.log)
	fh.setLevel(logging.DEBUG)
	
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)

	logger.addHandler(fh)

	# Verification of given arguments

	try:
		logger.info('Verification of vcf file')
		verif_input_vcf(vcf_file)
		logger.info('- Vcf file ok -')
	except ValueError:
		print ('Problem with vcf file {vcf_file}:', sys.exc_info()[0])
		logger.error('- Problem with vcf file -')
		exit(1)
	try:
		logger.info('Verification of output file')
		verif_output(output)
		logger.info('- Output file ok -')
	except ValueError:
		print ('Problem with output file {output}:', sys.exc_info()[0])
		logger.error('- Problem with output file -')
		exit(1)

	# Start

	logger.info('**************************************************************************************************************')
	logger.info('*** g-LOTUS filtering module ***')
	logger.info('* Start filtering *')
	logger.info(f'Working directory (vcf files folder) : {working_directory}') 
	logger.info(f'Current directory : {Path().absolute()}')

	filter(vcf_file , logger, output, working_method)

        # End
